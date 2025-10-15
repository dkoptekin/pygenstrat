import numpy as np
import pandas as pd
import logging
import sys
import os
import shutil
from tqdm import tqdm
from .config import EigenConfig
from .io_handlers import (
    EigenstratReader,
    read_ind_file,
    read_snp_file,
    OutputHandler)
from .processors import (
    filter_individuals,
    filter_snps,
    filter_by_geno,
    filter_by_maf,
    calculate_ind_missing,
    calculate_snp_stats,
    update_individual_populations,
    update_snp_ids,
    update_genetic_distances,
    random_haploidise,
    load_regions,
    load_chromosomes,
    load_genetic_map,
    polarise,
    apply_sex_chr_missing,
    _get_sex_chr_ids,
    _normalize_chr_ids,
    complement_snp_alleles,
    apply_update,
    log_sample_sex_counts,
    log_snp_counts_with_sex,
    detect_sex_chromosomes,
    log_missingness,
    log_chr_missingness,
    log_freq,
    log_female_y,
    log_polarisation)

class EigenDataset:
    def __init__(self, config: EigenConfig):
        self.config = config
        self.config.validate_parameters()
        self.config.validate_files()
        self.ind_data = read_ind_file(config.ind_path)
        self.snp_data = read_snp_file(config.snp_path)
        log_sample_sex_counts(self.ind_data, logger=logging)
        self.genetic_distances_updated = False
        self.ind_data, self.snp_data, _ = apply_update(self.ind_data, self.snp_data, self.config, in_place=False)
        if hasattr(config, 'interpolators') and config.interpolators is not None:
            self.interpolators = config.interpolators
            if config.out is not None:
                logging.info(f"Updating genetic distances using genetic map (unit: {config.map_unit})")
                self.snp_data = update_genetic_distances(
                    self.snp_data,
                    interpolators=self.interpolators,
                    output_unit=config.map_unit,
                    force_zero=getattr(config, 'genetic_distance', None) == 'zero')
                self.genetic_distances_updated = True
        elif hasattr(config, 'genetic_distance') and config.genetic_distance:
            if config.genetic_distance == 'zero':
                self.interpolators = None
            else:
                logging.info(f"Loading genetic distance map: {config.genetic_distance}")
                _, self.interpolators = load_genetic_map(config.genetic_distance)
            if config.out is not None:
                logging.info(f"Updating genetic distances using genetic distance map (unit: {config.map_unit})")
                self.snp_data = update_genetic_distances(
                    self.snp_data,
                    interpolators=getattr(self.config, 'interpolators', None),
                    output_unit=config.map_unit,
                    force_zero=getattr(self.config, 'genetic_distance', None) == 'zero')
                self.genetic_distances_updated = True
        self.geno_reader = EigenstratReader(config.geno_path, config.chunk_size)
        has_sex_chromosomes = detect_sex_chromosomes(self.snp_data, self.config)
        if has_sex_chromosomes:
            log_snp_counts_with_sex(self.snp_data, self.config, logger=logging)
            filtering_applied = any([
                self.config.geno is not None,
                self.config.min_maf is not None,
                self.config.max_maf is not None,
                self.config.missing,
                self.config.freq,
                self.config.sex_chr_missing])
            if filtering_applied:
                logging.info(
                    f"Filtering will be applied appropriately for each chromosome. If you want to ignore sex chromosomes use --ignore-sex-chr flag")

    def process_chunks(self) -> tuple:
        chunk_count = 0
        total_processed_snps = 0
        final_n_snps = 0
        indv_mask = np.zeros(0, dtype=np.bool_)
        all_stats = []
        total_n_snps = 0
        total_n_missing = 0
        chr_missing_data = {}
        first_chunk = True
        snp_filter_mask = None
        pol_df = None
        snp_annot_full = None
        is_polarise_file = False
        random_haploidise_diploid_mask = None
        try:
            n_chunks = (len(self.snp_data) + self.config.chunk_size - 1) // self.config.chunk_size
            pbar = tqdm(total=n_chunks,
                        desc="Processing chunks",
                        unit="chunk",
                        file=sys.stdout)

            needs_indv_filter = any([
                self.config.keep_indv, self.config.remove_indv,
                self.config.keep_pop, self.config.remove_pop])
            if needs_indv_filter:
                ind_data_filtered, indv_mask = filter_individuals(self.ind_data, self.config)
                if indv_mask is not None and indv_mask.any():
                    self.ind_data = ind_data_filtered
            else:
                indv_mask = np.zeros(0, dtype=np.bool_)

            needs_snp_static_mask = any([
                self.config.keep_snps, self.config.remove_snps
            ])
            if needs_snp_static_mask:
                snp_mask_file = self.config.remove_snps if self.config.remove_snps else self.config.keep_snps
                snp_filter_set = set(np.loadtxt(snp_mask_file, dtype=str))
                mask_values = self.snp_data['snpID'].isin(snp_filter_set).values
                if self.config.keep_snps:
                    snp_mask = mask_values
                else:
                    snp_mask = ~mask_values

                if snp_mask.sum() > 0:
                    self.snp_data = self.snp_data.loc[snp_mask].reset_index(drop=True)
                else:
                    logging.warning(f"All SNPs filtered out by SNP ID static filtering - output will be empty.")

            output = OutputHandler(self.config, self.ind_data)
            snp_filter_mask = None
            regions = self._load_regions()
            chromosomes = self._load_chromosomes()
            if self.config.polarise and os.path.isfile(str(self.config.polarise)):
                pol_df = pd.read_csv(self.config.polarise, sep=None, engine='python', names=['snpID', 'anc'])
                self.snp_data['snpID'] = self.snp_data['snpID'].astype(pol_df['snpID'].dtype)
                pol_df['snpID'] = pol_df['snpID'].astype(self.snp_data['snpID'].dtype)
                snp_annot_full = self.snp_data.merge(pol_df, on='snpID', how='left')
                is_polarise_file = True
            if self.config.random_haploidise:
                random_haploidise_diploid_mask = np.zeros(self.ind_data.shape[0], dtype=bool)
            for geno_chunk in self.geno_reader:
                chunk_count += 1
                chunk_size = len(geno_chunk)
                start_idx = (chunk_count - 1) * self.config.chunk_size
                end_idx = start_idx + chunk_size
                snp_chunk = self.snp_data.iloc[start_idx:end_idx]
                total_processed_snps += chunk_size
                geno_chunk, snp_chunk = self._apply_filters_to_chunk(
                    geno_chunk, snp_chunk, start_idx, end_idx, indv_mask, snp_filter_mask, pol_df, snp_annot_full, is_polarise_file, regions, chromosomes, random_haploidise_diploid_mask)
                final_n_snps, total_n_snps, total_n_missing, chr_missing_data, all_stats = self._update_statistics(
                    geno_chunk, snp_chunk, final_n_snps, total_n_snps, total_n_missing, chr_missing_data, all_stats, indv_mask)
                self._write_outputs(
                    output, first_chunk, geno_chunk, snp_chunk, ref_alt_alleles=getattr(self.geno_reader, 'ref_alt_alleles', None))
                first_chunk = False
                pbar.update(1)
                if n_chunks > 0:
                    progress = (chunk_count / n_chunks) * 100
                    pbar.set_postfix({"Progress": f"{progress:.2f}%"})
            self._finalize_all_outputs(output, final_n_snps, total_processed_snps,
                total_n_snps, total_n_missing, chr_missing_data, all_stats, indv_mask)
            return final_n_snps, total_processed_snps
        except Exception as e:
            logging.error(f"Error processing files: {e}")
            raise

    def _filter_individuals(self) -> np.ndarray:
        if any([self.config.keep_indv, self.config.remove_indv,
                self.config.keep_pop, self.config.remove_pop]):
            _, indv_mask = filter_individuals(self.ind_data, self.config)
            return indv_mask
        return np.zeros(0, dtype=np.bool_)

    def _create_snp_filter_mask(self) -> np.ndarray:
        if any([self.config.keep_snps, self.config.remove_snps]):
            snp_mask_file = self.config.remove_snps if self.config.remove_snps else self.config.keep_snps
            snp_filter_set = set(np.loadtxt(snp_mask_file, dtype=str))
            mask_values = self.snp_data['snpID'].isin(snp_filter_set).values
            return mask_values if self.config.keep_snps else ~mask_values
        return None

    def _load_regions(self):
        if self.config.keep_region:
            return load_regions(self.config.keep_region)
        elif self.config.remove_region:
            return load_regions(self.config.remove_region)
        return None

    def _load_chromosomes(self):
        if self.config.keep_chr:
            return load_chromosomes(self.config.keep_chr)
        elif self.config.remove_chr:
            return load_chromosomes(self.config.remove_chr)
        return None

    def _apply_filters_to_chunk(self, geno_chunk, snp_chunk, start_idx, end_idx, indv_mask, snp_filter_mask, pol_df, snp_annot_full, is_polarise_file, regions, chromosomes, random_haploidise_diploid_mask):
        config = self.config
        if config.polarise:
            if is_polarise_file:
                snp_annot_chunk = snp_annot_full.iloc[start_idx:end_idx]
                geno_chunk, snp_chunk = polarise(geno_chunk, snp_annot_chunk, polarise_is_file=True, config=config)
            else:
                geno_chunk, snp_chunk = polarise(geno_chunk, snp_chunk, polarise_is_file=False, ind_data=self.ind_data, config=config, sample_id=config.polarise)
        if config.random_haploidise:
            if random_haploidise_diploid_mask is not None:
                random_haploidise_diploid_mask |= (geno_chunk == 1).any(axis=0)
            geno_chunk = random_haploidise(geno_chunk)
        if config.sex_chr_missing:
            geno_chunk = apply_sex_chr_missing(geno_chunk, snp_chunk, self.ind_data, config)
        if indv_mask.size > 0 and indv_mask.sum() > 0:
            geno_chunk = geno_chunk[:, indv_mask]
        if snp_filter_mask is not None:
            chunk_mask = snp_filter_mask[start_idx:end_idx]
            geno_chunk = geno_chunk[chunk_mask]
            snp_chunk = snp_chunk[chunk_mask]
        if any([config.keep_region, config.remove_region, config.keep_chr, config.remove_chr]):
            geno_chunk, snp_chunk = filter_snps(
                geno_chunk, snp_chunk, config=config, regions=regions, chromosomes=chromosomes)
        if config.geno is not None:
            geno_chunk, snp_chunk = filter_by_geno(
                geno_chunk, snp_chunk, config.geno, config, self.ind_data)
        if config.min_maf is not None or config.max_maf is not None:
            geno_chunk, snp_chunk = filter_by_maf(
                geno_chunk, snp_chunk, config.min_maf, config.max_maf, config, self.ind_data)
        return geno_chunk, snp_chunk

    def _update_statistics(self, geno_chunk, snp_chunk, final_n_snps, total_n_snps, total_n_missing, chr_missing_data, all_stats, indv_mask):
        config = self.config
        final_n_snps += len(snp_chunk)
        # individual-level missing
        if config.missing:
            n_snps_chunk, n_missing_chunk = calculate_ind_missing(
                geno_chunk, self.ind_data, snp_chunk, config)
            total_n_snps += n_snps_chunk
            total_n_missing += n_missing_chunk
        # missing-by-chrom
        if config.missing_by_chr:
            current_ind_data = self.ind_data if indv_mask.size == 0 or not indv_mask.any() else self.ind_data[indv_mask]
            for chr_id in snp_chunk['chr'].unique():
                if chr_id not in chr_missing_data:
                    chr_missing_data[chr_id] = {'n_snps': np.zeros(len(current_ind_data), dtype=int),
                                                'n_missing': np.zeros(len(current_ind_data), dtype=int)}
                chr_mask = (snp_chunk['chr'] == chr_id).values
                chr_geno = geno_chunk[chr_mask]
                chr_snp = snp_chunk[chr_mask]
                if len(chr_snp) > 0:
                    chr_n_snps, chr_n_missing = calculate_ind_missing(
                        chr_geno, current_ind_data, chr_snp, config)
                    chr_missing_data[chr_id]['n_snps'] += chr_n_snps
                    chr_missing_data[chr_id]['n_missing'] += chr_n_missing
        # freq
        if config.freq:
            stats_df, _ = calculate_snp_stats(
                geno_chunk, snp_chunk, self.ind_data, config)
            all_stats.append(stats_df)
        return final_n_snps, total_n_snps, total_n_missing, chr_missing_data, all_stats

    def _write_outputs(self, output: OutputHandler, first_chunk: bool, geno_chunk, snp_chunk, ref_alt_alleles=None):
        if first_chunk:
            if self.config.out_ind is not None:
                output.write_ind(self.ind_data, self.config.out_ind)
        if self.config.out_snp is not None:
            output.write_snp(snp_chunk, self.config.out_snp, first_chunk, ref_alt_alleles)
        if self.config.out_geno is not None:
            output.write_geno(geno_chunk, first_chunk)

    def _finalize_all_outputs(self, output: OutputHandler, final_n_snps, total_processed_snps,
            total_n_snps, total_n_missing, chr_missing_data, all_stats, indv_mask):
        config = self.config
        if config.missing:
            if config.out_missing is not None:
                n_total = total_n_snps + total_n_missing
                missing_rate = np.where(n_total > 0, total_n_missing / n_total, np.nan)
                miss_df = pd.DataFrame({
                    'iid': self.ind_data['iid'],
                    'pop': self.ind_data['population'],
                    'n_snps': total_n_snps,
                    'n_missing': total_n_missing,
                    'n_total': n_total,
                    'missing_rate': np.round(missing_rate,4)})
                miss_df.to_csv(
                    config.out_missing,
                    sep='\t',
                    index=False)
                log_missingness(miss_df, logger=logging)
        if hasattr(config, 'missing_by_chr') and config.missing_by_chr:
            if config.out_missing_by_chr is not None:
                current_ind_data = self.ind_data if indv_mask.size == 0 or not indv_mask.any() else self.ind_data[indv_mask]
                for chr_id, chr_data in chr_missing_data.items():
                    n_total_chr = chr_data['n_snps'] + chr_data['n_missing']
                    missing_rate_chr = np.divide(chr_data['n_missing'], n_total_chr,
                                                out=np.full_like(chr_data['n_missing'], np.nan, dtype=float),
                                                where=(n_total_chr > 0))
                    _, y_chr_id = _get_sex_chr_ids(config)
                    if str(chr_id) == str(y_chr_id):
                        male_mask = (current_ind_data['sex'] == 'M').values
                        chr_miss_df = pd.DataFrame({
                            'iid': current_ind_data[male_mask]['iid'],
                            'population': current_ind_data[male_mask]['population'],
                            'chr': chr_id,
                            'n_snps': chr_data['n_snps'][male_mask],
                            'n_missing': chr_data['n_missing'][male_mask],
                            'n_total': n_total_chr[male_mask],
                            'missing_rate': np.round(missing_rate_chr[male_mask], 4)})
                    else:
                        chr_miss_df = pd.DataFrame({
                            'iid': current_ind_data['iid'],
                            'population': current_ind_data['population'],
                            'chr': chr_id,
                            'n_snps': chr_data['n_snps'],
                            'n_missing': chr_data['n_missing'],
                            'n_total': n_total_chr,
                            'missing_rate': np.round(missing_rate_chr, 4)})
                    chr_filename = f"{config.out_missing_by_chr}.chr{chr_id}.missing"
                    chr_miss_df.to_csv(chr_filename, sep='\t', index=False)
                    log_chr_missingness(chr_id, missing_rate_chr, y_chr_id, logger=logging)
        if config.freq:
            if config.out_freq is not None:
                combined_stats = pd.concat(all_stats, axis=0)
                combined_stats.to_csv(
                    config.out_freq,
                    sep='\t',
                    index=False)
                log_freq(combined_stats, logger=logging)
        # update header of anc geno
        if config.out_type == 'anc' and config.out_geno is not None:
            if os.path.exists(config.out_geno):
                output.finalise_packed_geno(config.out_geno)
            else:
                logging.error(f"Expected genotype file {config.out_geno} does not exist")
                raise FileNotFoundError(f"Genotype file missing: {config.out_geno}")
        log_female_y(config)
        if config.polarise:
            log_polarisation(config)
        if config.random_haploidise:
            logging.info("Random haploidisation: check log for details.")
        output_snps = final_n_snps
        output_individuals = len(self.ind_data)
        logging.info(f"{output_snps} SNPs and {output_individuals} individuals in output.")