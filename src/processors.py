import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
import logging
import intervaltree
import shutil
from .io_handlers import backup_file, safe_remove

#### update
def apply_update(ind_data, snp_data, config, in_place=False):
    """
    Apply updates to individual and SNP data.
    
    Handles: population updates, SNP ID updates, genetic distance updates, strand flipping.
    Can operate in-place (modifying files) or return updated data.
    """
    changes_made = False
    if config.update_ind:
        updated_ind_data = update_individual_populations(ind_data, config.update_ind)
        if in_place:
            backup_file(config.ind_path)
            updated_ind_data.to_csv(
                config.ind_path, sep='\t', header=False, index=False,
                columns=['iid', 'sex', 'population'])
            logging.info(f"Successfully updated {config.ind_path}")
        ind_data = updated_ind_data
        changes_made = True
    if config.update_snp:
        updated_snp_data = update_snp_ids(snp_data, config.update_snp)
        if in_place:
            backup_file(config.snp_path)
            updated_snp_data.to_csv(
                config.snp_path, sep='\t', header=False, index=False)
            logging.info(f"Successfully updated {config.snp_path}")
        snp_data = updated_snp_data
        changes_made = True
    if hasattr(config, 'genetic_distance') and config.genetic_distance:
        updated_snp_data = update_genetic_distances(
            snp_data,
            interpolators=getattr(config, 'interpolators', None),
            output_unit=config.map_unit,
            force_zero=getattr(config, 'genetic_distance', None) == 'zero')
        if in_place:
            backup_file(config.snp_path)
            updated_snp_data.to_csv(
                config.snp_path, sep='\t', header=False, index=False)
            logging.info(f"Successfully updated genetic distances in {config.snp_path}")
        snp_data = updated_snp_data
        changes_made = True
    if hasattr(config, 'flip_strand') and config.flip_strand:
        with open(config.flip_strand) as f:
            flip_ids = [line.strip() for line in f if line.strip()]
        updated_snp_data = complement_snp_alleles(snp_data, flip_ids)
        if in_place:
            backup_file(config.snp_path)
            updated_snp_data.to_csv(
                config.snp_path, sep='\t', header=False, index=False)
            logging.info(f"Successfully complemented alleles for {len(flip_ids)} SNPs in {config.snp_path}")
        snp_data = updated_snp_data
        changes_made = True
    if in_place:
        return None, None, changes_made
    else:
        return ind_data, snp_data, changes_made

def update_individual_populations(ind_data, update_file_path):
    """
    Update population labels for individuals.
    File format: two columns (individual_id, new_population_name).
    """
    try:
        update_df = pd.read_csv(
            update_file_path,
            sep=r'\s+',
            header=None,
            names=['iid', 'new_population'])
        update_map = dict(zip(update_df['iid'], update_df['new_population']))
        to_update = ind_data['iid'].isin(update_map.keys()).sum()
        # If updatePOP has more rows than are updated, warn with missing sampleIDs
        if to_update < len(update_df['iid']):
            not_found = [iid for iid in update_df['iid'] if iid not in set(ind_data['iid'])]
            logging.warning(f"{len(not_found)} sample IDs from updatePOP file not found in IND file (not updated): {', '.join(map(str, not_found))}")
        
        if to_update == 0:
            logging.warning("No individuals found to update. Check your IDs.")
            return ind_data
        ind_data_copy = ind_data.copy()
        ind_data_copy.loc[ind_data_copy['iid'].isin(update_map.keys()), 'population'] = \
            ind_data_copy.loc[ind_data_copy['iid'].isin(update_map.keys()), 'iid'].map(update_map)
        logging.info(f"Population information updated for {to_update} individuals")
        return ind_data_copy

    except Exception as e:
        logging.error(f"Error updating individual population information: {e}")
        raise

def update_snp_ids(snp_data, update_file_path):
    """
    Update SNP ID
    File format: two columns (old_snp_id, new_snp_id).
    """
    try:
        update_df = pd.read_csv(
            update_file_path,
            sep=r'\s+',
            header=None,
            names=['old_snpID', 'new_snpID'])
        update_map = dict(zip(update_df['old_snpID'], update_df['new_snpID']))
        to_update = snp_data['snpID'].isin(update_map.keys()).sum()

        if to_update == 0:
            logging.warning("No SNPs found to update. Check your SNP IDs.")
            return snp_data

        snp_data_copy = snp_data.copy()
        snp_data_copy.loc[snp_data_copy['snpID'].isin(update_map.keys()), 'snpID'] = \
            snp_data_copy.loc[snp_data_copy['snpID'].isin(update_map.keys()), 'snpID'].map(update_map)

        logging.info(f"SNP IDs updated for {to_update} SNPs")
        return snp_data_copy

    except Exception as e:
        logging.error(f"Error updating SNP IDs: {e}")
        raise

def load_genetic_map(genetic_map):
    """
    Load genetic map from file and create interpolation functions.
    
    Returns chromosome-specific genetic maps and interpolation functions
    for converting physical positions (bp) to genetic distances (cM).
    """
    try:
        with open(genetic_map, 'r') as f:
            first_line = f.readline().strip()
            header_present = "Position" in first_line or "position" in first_line

        map_df = pd.read_csv(
            genetic_map,
            sep=r'\s+',
            header=0 if header_present else None,
            names=None if header_present else ['chr', 'position', 'rate', 'map'],
            dtype=str)

        if header_present:
            col_map = {}
            for col in map_df.columns:
                if 'chr' in col.lower():
                    col_map[col] = 'chr'
                elif 'position' in col.lower() or 'bp' in col.lower():
                    col_map[col] = 'position'
                elif 'rate' in col.lower() or 'cm/mb' in col.lower():
                    col_map[col] = 'rate'
                elif 'map' in col.lower() or 'cm' in col.lower():
                    col_map[col] = 'map'
            map_df.rename(columns=col_map, inplace=True)

        map_df['position'] = pd.to_numeric(map_df['position'], errors='coerce')
        map_df['rate'] = pd.to_numeric(map_df['rate'], errors='coerce')
        map_df['map'] = pd.to_numeric(map_df['map'], errors='coerce')
        map_df['chr'] = map_df['chr'].str.replace('chr', '', regex=False)
        map_df = map_df.dropna(subset=['position', 'map'])
        chromosome_maps = {}
        interp_funcs = {}

        for chrom, group in map_df.groupby('chr'):
            group = group.sort_values('position')
            chromosome_maps[str(chrom)] = group

            bps = group['position'].values
            cMs = group['map'].values

            if len(bps) > 1:
                interp_funcs[str(chrom)] = interp1d(bps, cMs, kind='linear', fill_value="extrapolate")
            else:
                logging.warning(f"Not enough points for interpolation on chromosome {chrom}")

        logging.info(f"Created genetic distance map interpolation functions for {len(interp_funcs)} chromosomes")
        return chromosome_maps, interp_funcs

    except Exception as e:
        logging.error(f"Error loading genetic map: {e}")
        raise

def update_genetic_distances(snp_chunk, interpolators=None, output_unit='cM', force_zero=False):
    """
    Update genetic distances using interpolation from a genetic map.
    
    Interpolates centiMorgan (cM) positions from physical positions (bp)
    either in centiMorgans (cM) or Morgans (M).
    """
    if force_zero:
        snp_chunk['cM'] = 0
        logging.info(f"Set genetic distances (cM) to 0 for all {len(snp_chunk)} SNPs (force_zero mode)")
        return snp_chunk
    if interpolators is None:
        logging.error("No interpolation functions provided")
        return snp_chunk

    chrom_map = {}
    interpolated_chroms = set(interpolators.keys()) if interpolators is not None else set()
    snp_chr_list = snp_chunk['chr'].astype(str).str.lower().str.replace('chr', '').unique()
    missing_chr_in_map = [chr_val for chr_val in snp_chr_list if chr_val not in set([c.lower().replace('chr','') for c in interpolated_chroms])]
    if missing_chr_in_map:
        logging.warning(f"Genetic map does not contain chromosomes: {', '.join(missing_chr_in_map)}; genetic distances for these chromosomes will be set to zero.")

    for chrom in interpolators.keys():
        chrom_clean = chrom.lower().replace('chr', '')
        chrom_map[chrom_clean] = chrom
        chrom_map[chrom] = chrom
        chrom_map[f"chr{chrom_clean}"] = chrom

    updated_count = 0

    for i, row in snp_chunk.iterrows():
        chrom = str(row['chr']).lower().replace('chr', '')
        pos = int(row['pos'])
        if chrom not in chrom_map or chrom_map[chrom] not in interpolators:
            snp_chunk.at[i, 'cM'] = 0
            continue
        orig_chrom = chrom_map[chrom]
        interp_func = interpolators[orig_chrom]
        try:
            genetic_pos_cM = float(interp_func(pos))
            if output_unit and str(output_unit).lower() == 'm':
                out_val = genetic_pos_cM / 100.0
            else:
                out_val = genetic_pos_cM
            snp_chunk.at[i, 'cM'] = round(out_val, 6)
            updated_count += 1
        except Exception as e:
            logging.warning(f"Error calculating genetic distance for SNP at position {pos} on chromosome {row['chr']}: {e}")

    unit_description = "Morgans (M)" if output_unit and str(output_unit).lower() == 'm' else "centiMorgans (cM)"
    logging.info(f"Updated genetic distances for {updated_count} SNPs using {unit_description}")
    return snp_chunk

#### sex chr-related functions
def detect_sex_chromosomes(snp_data, config):
    if config is None or getattr(config, 'ignore_sex_chr', False):
        return False
    from .processors import _get_sex_chr_ids, _normalize_chr_ids
    x_chr_id, y_chr_id = _get_sex_chr_ids(config)
    unique_chrs = set(_normalize_chr_ids(snp_data['chr']))
    has_x_chr = x_chr_id in unique_chrs
    has_y_chr = y_chr_id in unique_chrs
    return has_x_chr or has_y_chr

def _get_sex_chr_ids(config):
    if not config or not hasattr(config, 'sex_chr') or not config.sex_chr:
        return '23', '24'

    x_chr_id = config.sex_chr[0] if len(config.sex_chr) > 0 else '23'
    y_chr_id = config.sex_chr[1] if len(config.sex_chr) > 1 else config.sex_chr[0]
    return x_chr_id, y_chr_id

def _normalize_chr_ids(chr_values):
    return [str(chr_val).replace('chr', '').replace('Chr', '') for chr_val in chr_values]

def _detect_y_chromosome_snps(snp_chunk, config):
    if not config or config.ignore_sex_chr:
        return False, np.array([])
    _, y_chr_id = _get_sex_chr_ids(config)
    current_chr_ids = _normalize_chr_ids(snp_chunk['chr'])
    y_chr_mask = np.array([chr_id == y_chr_id for chr_id in current_chr_ids])
    return np.any(y_chr_mask), y_chr_mask

def _get_sex_chr_info(snp_chunk, config):
    if not config or config.ignore_sex_chr:
        return None

    x_chr_id, y_chr_id = _get_sex_chr_ids(config)
    current_chr_ids = _normalize_chr_ids(snp_chunk['chr'])

    x_chr_mask = np.array([chr_id == x_chr_id for chr_id in current_chr_ids])
    y_chr_mask = np.array([chr_id == y_chr_id for chr_id in current_chr_ids])

    has_x_chr = np.any(x_chr_mask)
    has_y_chr = np.any(y_chr_mask)
    has_sex_chr = has_x_chr or has_y_chr

    if not has_sex_chr:
        return None

    x_indices = np.where(x_chr_mask)[0] if has_x_chr else np.array([])
    y_indices = np.where(y_chr_mask)[0] if has_y_chr else np.array([])

    sex_chr_mask = x_chr_mask | y_chr_mask
    non_sex_chr_indices = np.where(~sex_chr_mask)[0]

    return {
        'has_x_chr': has_x_chr,
        'has_y_chr': has_y_chr,
        'has_sex_chr': has_sex_chr,
        'x_chr_mask': x_chr_mask,
        'y_chr_mask': y_chr_mask,
        'x_indices': x_indices,
        'y_indices': y_indices,
        'non_sex_chr_indices': non_sex_chr_indices}

def check_and_log_female_y_genotypes(geno, snp_chunk, ind_data, config, sex_chr_info=None):
    if (config is not None and not config.ignore_sex_chr and
            ind_data is not None and snp_chunk is not None):

        if sex_chr_info is None:
            sex_chr_info = _get_sex_chr_info(snp_chunk, config)

        if sex_chr_info and sex_chr_info['has_y_chr']:
            if not hasattr(config, '_female_y_count'):
                female_mask = (ind_data['sex'] == 'F').values
                config._female_y_count = 0
                config._total_females = np.sum(female_mask)
                config._females_checked = set()

            female_mask = (ind_data['sex'] == 'F').values
            y_indices = sex_chr_info['y_indices']

            for idx in y_indices:
                female_y_genos = geno[idx, female_mask]
                non_missing_female_positions = np.where(female_y_genos != 9)[0]

                female_indices = np.where(female_mask)[0]
                non_missing_female_indices = female_indices[non_missing_female_positions]
                config._females_checked.update(non_missing_female_indices)

            config._female_y_count = len(config._females_checked)

def log_female_y(config):
    if hasattr(config, '_female_y_count') and config._female_y_count > 0:
        logging.info(
            f"WARNING: {config._female_y_count} of {config._total_females} female(s) and unknown sex individuals have non-missing Y chromosome genotypes.")

    if hasattr(config, '_male_x_het_count') and config._male_x_het_count > 0:
        logging.info(
            f"INFO: Set {config._male_x_het_count} non-female heterozygous X chromosome genotypes to missing.")

def apply_sex_chr_missing(geno, snp_data, ind_data, config, sex_chr_info=None):
    """
    - X chromosome: set male heterozygous genotypes to missing
    - Y chromosome: set non-male genotypes to missing 
    """
    if (config is None or not config.sex_chr_missing or config.ignore_sex_chr or
            ind_data is None or snp_data is None):
        return geno.copy()

    if sex_chr_info is None:
        sex_chr_info = _get_sex_chr_info(snp_data, config)

    if not (sex_chr_info and sex_chr_info['has_sex_chr']):
        return geno.copy()

    geno_processed = geno.copy()
    x_indices = sex_chr_info['x_indices']
    y_indices = sex_chr_info['y_indices']

    check_and_log_female_y_genotypes(geno_processed, snp_data, ind_data, config, sex_chr_info)

    male_mask = (ind_data['sex'] == 'M').values
    female_mask = (ind_data['sex'] == 'F').values
    non_female_mask = ~female_mask

    if len(y_indices) > 0:
        non_male_y_genos = geno_processed[np.ix_(y_indices, ~male_mask)]
        non_missing_mask = non_male_y_genos != 9
        if np.any(non_missing_mask):
            geno_processed[np.ix_(y_indices, ~male_mask)] = np.where(
                non_missing_mask, 9, non_male_y_genos
            )
            n_changed = np.sum(non_missing_mask)
            logging.info(f"Set {n_changed} non-male Y chromosome genotypes to missing")

    if np.any(non_female_mask) and len(x_indices) > 0:
        non_female_x_genos = geno_processed[np.ix_(x_indices, non_female_mask)]
        het_mask = (non_female_x_genos == 1) & (non_female_x_genos != 9)  # Heterozygous and not already missing
        if np.any(het_mask):
            geno_processed[np.ix_(x_indices, non_female_mask)] = np.where(
                het_mask, 9, non_female_x_genos
            )
            n_changed = np.sum(het_mask)
            if not hasattr(config, '_male_x_het_count'):
                config._male_x_het_count = 0
            config._male_x_het_count += n_changed
            logging.info(f"Set {n_changed} non-female heterozygous X chromosome genotypes to missing")

    return geno_processed

#### indv filtering
def filter_individuals(ind_data, config):
    """
    Filter individuals based on keep/remove lists for individuals or populations.
    """
    indv_mask = np.ones(len(ind_data), dtype=bool)
    # individual level
    if config.keep_indv:
        with open(config.keep_indv) as f:
            keep_ids = set(line.strip() for line in f)
        indv_mask &= ind_data['iid'].isin(keep_ids)
    elif config.remove_indv:
        with open(config.remove_indv) as f:
            remove_ids = set(line.strip() for line in f)
        indv_mask &= ~ind_data['iid'].isin(remove_ids)
    # population level
    if config.keep_pop:
        with open(config.keep_pop) as f:
            keep_pops = set(line.strip() for line in f)
        indv_mask &= ind_data['population'].isin(keep_pops)
    elif config.remove_pop:
        with open(config.remove_pop) as f:
            remove_pops = set(line.strip() for line in f)
        indv_mask &= ~ind_data['population'].isin(remove_pops)
    return ind_data[indv_mask], indv_mask

#### snp filtering
def filter_snps(geno_chunk, snp_chunk, config, snp_filter_list=None, regions=None, chromosomes=None):
    """
    Filter SNPs based on chromosome, SNP ID, or genomic region criteria.
    """
    need_filtering = False
    n_snps_original = len(snp_chunk)
    if config.keep_chr and chromosomes is not None:
        chr_mask = filter_by_chromosomes(snp_chunk, chromosomes, keep=True)
        snp_chunk = apply_mask_and_log(snp_chunk, chr_mask, "CHR", config)
        geno_chunk = apply_mask_and_log(geno_chunk, chr_mask, "CHR", config, snp_chunk)
        need_filtering = True
    elif config.remove_chr and chromosomes is not None:
        chr_mask = filter_by_chromosomes(snp_chunk, chromosomes, keep=False)
        snp_chunk = apply_mask_and_log(snp_chunk, chr_mask, "CHR", config)
        geno_chunk = apply_mask_and_log(geno_chunk, chr_mask, "CHR", config, snp_chunk)
        need_filtering = True

    if len(snp_chunk) == 0:
        logging.debug(f"All SNPs filtered out by chromosome filter")
        return geno_chunk, snp_chunk

    if config.keep_snps and snp_filter_list is not None:
        if not isinstance(snp_filter_list, set):
            snp_filter_list = set(snp_filter_list)
        snp_mask = snp_chunk['snpID'].isin(snp_filter_list).values
        snp_chunk = apply_mask_and_log(snp_chunk, snp_mask, "SNP_ID", config)
        geno_chunk = apply_mask_and_log(geno_chunk, snp_mask, "SNP_ID", config, snp_chunk)
        need_filtering = True
    elif config.remove_snps and snp_filter_list is not None:
        if not isinstance(snp_filter_list, set):
            snp_filter_list = set(snp_filter_list)

        snp_mask = ~snp_chunk['snpID'].isin(snp_filter_list).values
        snp_chunk = apply_mask_and_log(snp_chunk, snp_mask, "SNP_ID", config)
        geno_chunk = apply_mask_and_log(geno_chunk, snp_mask, "SNP_ID", config, snp_chunk)
        need_filtering = True

    if len(snp_chunk) == 0:
        logging.debug(f"All SNPs filtered out by SNP ID filter")
        return geno_chunk, snp_chunk

    if config.keep_region and regions is not None:
        region_mask = filter_by_regions(snp_chunk, regions, keep=True)
        snp_chunk = apply_mask_and_log(snp_chunk, region_mask, "REGION", config)
        geno_chunk = apply_mask_and_log(geno_chunk, region_mask, "REGION", config, snp_chunk)
        need_filtering = True
    elif config.remove_region and regions is not None:
        region_mask = filter_by_regions(snp_chunk, regions, keep=False)
        snp_chunk = apply_mask_and_log(snp_chunk, region_mask, "REGION", config)
        geno_chunk = apply_mask_and_log(geno_chunk, region_mask, "REGION", config, snp_chunk)
        need_filtering = True

    if need_filtering:
        n_snps_filtered = n_snps_original - len(snp_chunk)
        if n_snps_filtered > 0:
            logging.debug(f"Filtered out {n_snps_filtered} SNPs ({len(snp_chunk)} remaining)")

    return geno_chunk, snp_chunk

def load_regions(region_input):
    try:
        regions = []
        chr_to_tree = {}
        if os.path.exists(region_input):
            with open(region_input, 'r') as f:
                region_lines = [line.strip() for line in f if line.strip()]
            logging.info(f"Loaded region specifications from file: {region_input}")
        elif ':' in region_input and '-' in region_input and ',' not in region_input:
            region_lines = [region_input]
            logging.info(f"Processing direct region specification: {region_input}")
        else:
            region_lines = region_input.split(',')
            logging.info(f"Processing {len(region_lines)} comma-separated region specifications")

        for line in region_lines:
            try:
                if ':' in line and '-' in line:
                    chrom, pos_range = line.split(':')
                    start, end = map(int, pos_range.split('-'))
                    parts = [chrom, start, end]
                else:
                    parts = line.split()
                    if len(parts) < 3:
                        logging.warning(f"Skipping invalid region: {line}")
                        continue
                    parts[1] = int(parts[1])
                    parts[2] = int(parts[2])

                chrom = parts[0]
                start = parts[1]
                end = parts[2]
                chrom_clean = str(chrom).lower().replace('chr', '')

                if start == end:
                    regions.append((chrom_clean, start, end))
                    if chrom_clean not in chr_to_tree:
                        chr_to_tree[chrom_clean] = intervaltree.IntervalTree()
                    chr_to_tree[chrom_clean].add(intervaltree.Interval(start, start + 1))
                else:
                    regions.append((chrom_clean, start, end))
                    if chrom_clean not in chr_to_tree:
                        chr_to_tree[chrom_clean] = intervaltree.IntervalTree()
                    chr_to_tree[chrom_clean].add(intervaltree.Interval(start + 1, end + 1))

            except Exception as e:
                logging.warning(f"Skipping invalid region: {line} - {str(e)}")

        for tree in chr_to_tree.values():
            tree.merge_overlaps()

        logging.info(f"Loaded {len(regions)} regions")
        return (regions, chr_to_tree)

    except Exception as e:
        logging.error(f"Error loading regions from {region_input}: {e}")
        return ([], {})

def load_chromosomes(chr_input):
    try:
        chromosomes = set()
        if os.path.exists(chr_input):
            with open(chr_input, 'r') as f:
                chromosomes = set([line.strip() for line in f if line.strip()])
        else:
            chs = chr_input.split(',')
            for ch in chs:
                ch = ch.strip()
                if '-' in ch and ch.replace('-', '').isdigit():
                    start, end = map(int, ch.split('-'))
                    for i in range(start, end + 1):
                        chromosomes.add(str(i))
                else:
                    chromosomes.add(ch)

        if not chromosomes:
            logging.warning("No chromosomes loaded from the input")

        logging.info(f"Loaded {len(chromosomes)} chromosomes from {chr_input}")
        return chromosomes

    except Exception as e:
        logging.error(f"Error loading chromosomes from {chr_input}: {e}")
        raise

def filter_by_regions(snp_chunk, region_data, keep=True):
    if not region_data or len(region_data) != 2:
        return np.ones(len(snp_chunk), dtype=bool) if keep else np.zeros(len(snp_chunk), dtype=bool)

    regions, chr_to_tree = region_data

    if not regions or not chr_to_tree:
        return np.ones(len(snp_chunk), dtype=bool) if keep else np.zeros(len(snp_chunk), dtype=bool)

    mask = np.zeros(len(snp_chunk), dtype=bool)

    snp_chrs = snp_chunk['chr'].astype(str).str.lower().str.replace('chr', '')

    for i, (chr_val, pos_val) in enumerate(zip(snp_chrs, snp_chunk['pos'])):
        try:
            pos_int = int(pos_val)
            if chr_val in chr_to_tree and chr_to_tree[chr_val].overlaps(pos_int):
                mask[i] = True
        except (ValueError, TypeError):
            logging.debug(f"Skipping SNP with invalid position: {pos_val}")
            continue

    return mask if keep else ~mask

def filter_by_chromosomes(snp_chunk, chromosomes, keep=True):
    if not chromosomes:
        return np.ones(len(snp_chunk), dtype=bool) if keep else np.zeros(len(snp_chunk), dtype=bool)
    chr_list = [str(c).replace('chr', '') for c in chromosomes]
    chr_series = snp_chunk['chr'].astype(str).str.replace('chr', '')
    in_chr_mask = chr_series.isin(chr_list).values
    final_mask = in_chr_mask if keep else ~in_chr_mask
    return final_mask

#### transform
def random_haploidise(geno):
    """
    Randomly convert heterozygous genotypes (1) to homozygous (0 or 2).
    """
    het_mask = (geno == 1)
    if np.any(het_mask):
        random_choices = np.random.randint(0, 2, size=het_mask.sum()) * 2
        geno[het_mask] = random_choices
    return geno

def polarise(geno, snp_chunk, polarise_is_file=False, ind_data=None, config=None, sample_id=None):
    """
    Polarise SNPs to ancestral/derived orientation.
    
    Two modes:
    1. File-based: Use ancestral allele annotations from file
    2. Sample-based: Use a reference individual's alleles as ancestral
    
    Flips alleles so allele1 is always ancestral, allele2 is derived.
    Removes SNPs where ancestral state cannot be determined.
    """
    excluded_count = 0
    if polarise_is_file:
        snp_chunk = snp_chunk.copy()
        geno_filtered = geno.copy()
        ops = []
        for idx, row in snp_chunk.iterrows():
            anc = row['anc'] if 'anc' in row else np.nan
            a1 = row['allele1']
            a2 = row['allele2']
            if pd.isna(anc):
                ops.append('remove')
                continue
            comp_a1 = base_complement(a1)
            comp_a2 = base_complement(a2)
            if anc == a2:
                ops.append('flip')
            elif anc == comp_a2:
                ops.append('flip+comp')
            elif anc == comp_a1:
                ops.append('comp')
            elif anc == a1:
                ops.append('none')
            else:
                ops.append('remove')
        snp_chunk['polarise_op'] = ops

        keep_mask = snp_chunk['polarise_op'] != 'remove'
        snp_chunk = snp_chunk[keep_mask].reset_index(drop=True)
        geno_filtered = geno_filtered[keep_mask.values]
        kept_ops = snp_chunk['polarise_op'].values

        flip_indices = np.where((kept_ops == 'flip') | (kept_ops == 'flip+comp'))[0]
        comp_indices = np.where(kept_ops == 'comp')[0]
        flip_comp_indices = np.where(kept_ops == 'flip+comp')[0]
        for i, op in enumerate(kept_ops):
            if op == 'flip':
                a1, a2 = snp_chunk.at[i, 'allele1'], snp_chunk.at[i, 'allele2']
                snp_chunk.at[i, 'allele1'], snp_chunk.at[i, 'allele2'] = a2, a1
            elif op == 'flip+comp':
                a1, a2 = snp_chunk.at[i, 'allele1'], snp_chunk.at[i, 'allele2']
                snp_chunk.at[i, 'allele1'] = base_complement(a2)
                snp_chunk.at[i, 'allele2'] = base_complement(a1)
            elif op == 'comp':
                a1, a2 = snp_chunk.at[i, 'allele1'], snp_chunk.at[i, 'allele2']
                snp_chunk.at[i, 'allele1'] = base_complement(a1)
                snp_chunk.at[i, 'allele2'] = base_complement(a2)
        if len(flip_indices) > 0:
            sel_geno = geno_filtered[flip_indices]
            flipped = sel_geno.copy()
            flipped[sel_geno == 0] = 2
            flipped[sel_geno == 2] = 0
            geno_filtered[flip_indices] = flipped
        n_polarised = len(flip_indices) + len(flip_comp_indices)
        n_processed = len(kept_ops)
        excluded_count = np.sum(~keep_mask)
    else:
        sample_mask = ind_data['iid'] == sample_id
        if not sample_mask.any():
            raise ValueError(f"Sample ID '{sample_id}' not found in dataset")
        sample_idx = sample_mask.idxmax()
        if not hasattr(ind_data, '_polarisation_ref_logged'):
            logging.info(f"Using sample '{sample_id}' (index {sample_idx}) as polarisation reference")
            ind_data._polarisation_ref_logged = True
        ref_genotypes = geno[:, sample_idx]
        is_valid = (ref_genotypes != 9) & (ref_genotypes != 1)
        n_processed = np.sum(is_valid)
        current_allele1 = snp_chunk['allele1'].to_numpy()
        current_allele2 = snp_chunk['allele2'].to_numpy()
        ancs = np.where(ref_genotypes == 0, current_allele2, current_allele1)
        keep_mask = is_valid
        flip_mask = (ancs == current_allele2) & is_valid
        snp_chunk = snp_chunk[keep_mask].reset_index(drop=True)
        geno_filtered = geno[keep_mask]
        rel_flip = np.where(keep_mask)[0][flip_mask[keep_mask]]
        if rel_flip.size > 0:
            allele1_col = snp_chunk.columns.get_loc('allele1')
            allele2_col = snp_chunk.columns.get_loc('allele2')
            rel_flip = rel_flip[rel_flip < len(snp_chunk)]  # Safety check
            tmp = snp_chunk.iloc[rel_flip, allele1_col].copy()
            snp_chunk.iloc[rel_flip, allele1_col] = snp_chunk.iloc[rel_flip, allele2_col].values
            snp_chunk.iloc[rel_flip, allele2_col] = tmp.values
            sel_geno = geno_filtered[rel_flip]
            flipped = sel_geno.copy()
            flipped[sel_geno == 0] = 2
            flipped[sel_geno == 2] = 0
            geno_filtered[rel_flip] = flipped
        n_polarised = rel_flip.size
        excluded_count = np.sum(~is_valid)
    if config is not None:
        update_polarisation_counts(config, n_polarised, n_processed, excluded_count)
    return geno_filtered, snp_chunk

def base_complement(seq):
    table = str.maketrans({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'})
    return seq.translate(table) if isinstance(seq, str) and len(seq) == 1 and seq in 'ACGT' else seq

def base_complement_array(array):
    return np.array([base_complement(a) for a in array], dtype=object)

def complement_snp_alleles(snp_data, snp_id_list, logger=logging):
    snp_data = snp_data.copy()
    snp_mask = snp_data['snpID'].astype(str).isin(set(snp_id_list))
    n = snp_mask.sum()
    snp_data.loc[snp_mask, 'allele1'] = snp_data.loc[snp_mask, 'allele1'].apply(base_complement)
    snp_data.loc[snp_mask, 'allele2'] = snp_data.loc[snp_mask, 'allele2'].apply(base_complement)
    if n > 0:
        logger.info(f"{n} SNPs converted to their complement.")
    else:
        logger.info("No SNPs converted to their complement.")
    return snp_data

#### maf
def filter_by_maf(geno, snp_chunk, min_maf=None, max_maf=None, config=None, ind_data=None):
    """
    Filter SNPs based on minor allele frequency (MAF) thresholds.
    """
    if min_maf is None and max_maf is None:
        return geno, snp_chunk

    sex_chr_info = _get_sex_chr_info(snp_chunk, config)
    maf, _, _, _, _ = calculate_maf(geno, ind_data, snp_chunk, config, sex_chr_info)

    if min_maf is not None and max_maf is not None:
        mask = (maf >= min_maf) & (maf < max_maf)
    elif min_maf is not None:
        mask = maf >= min_maf
    else:
        mask = maf < max_maf

    snp_chunk_filtered = apply_mask_and_log(snp_chunk, mask, "MAF", config)
    geno_filtered = apply_mask_and_log(geno, mask, "MAF", config, snp_chunk)
    return geno_filtered, snp_chunk_filtered

def calculate_maf(geno, ind_data=None, snp_chunk=None, config=None, sex_chr_info=None):
    """
    Calculate minor allele frequency (MAF) for SNPs.
    
    Automatically detects and handles sex chromosomes.
    Returns: maf, n_missing, valid_count, missing_rate, allele1_freq.
    """
    n_ind = geno.shape[1]

    if sex_chr_info is None and snp_chunk is not None:
        sex_chr_info = _get_sex_chr_info(snp_chunk, config)

    use_autosomal_path = not (sex_chr_info and sex_chr_info['has_sex_chr'])

    if use_autosomal_path:
        return _calculate_maf_autosomal(geno, n_ind)
    else:
        return _calculate_maf_sex_chr(geno, ind_data, snp_chunk, config, n_ind, sex_chr_info)

def _calculate_maf_autosomal(geno, n_ind):
    """
    Calculate MAF for autosomal chromosomes.
    """
    is_missing = (geno == 9)
    n_missing = np.sum(is_missing, axis=1)
    valid_count = n_ind - n_missing
    missing_rate = n_missing / n_ind

    is_hom_ref = (geno == 2)
    is_het = (geno == 1)
    n_allele1 = (np.sum(is_hom_ref, axis=1) * 2) + np.sum(is_het, axis=1)
    n_chrom = valid_count * 2
    n_chrom_c = np.maximum(n_chrom, 1)
    allele1_freq = n_allele1 / n_chrom_c
    maf = np.minimum(allele1_freq, 1.0 - allele1_freq)

    return maf, n_missing, valid_count, missing_rate, allele1_freq

def _calculate_maf_sex_chr(geno, ind_data, snp_chunk, config, n_ind, sex_chr_info):
    """
    Calculate MAF accounting for sex chromosomes (X and Y).
    
    - Y chromosome: only count males (haploid)
    - X chromosome: males are haploid, females are diploid
    - Heterozygous calls on Y and male X are treated as missing
    """
    geno_processed = apply_sex_chr_missing(geno, snp_chunk, ind_data, config, sex_chr_info)

    x_indices = sex_chr_info['x_indices'] if sex_chr_info else np.array([])
    y_indices = sex_chr_info['y_indices'] if sex_chr_info else np.array([])
    non_sex_chr_indices = sex_chr_info['non_sex_chr_indices'] if sex_chr_info else np.arange(geno.shape[0])
    
    male_mask = (ind_data['sex'] == 'M').values
    female_mask = (ind_data['sex'] == 'F').values
    unknown_mask = (~male_mask) & (~female_mask)
    
    n_males = np.sum(male_mask)
    n_females = np.sum(female_mask)
    n_unknown = np.sum(unknown_mask)

    n_missing = np.zeros(geno_processed.shape[0])
    valid_count = np.zeros(geno_processed.shape[0])
    missing_rate = np.zeros(geno_processed.shape[0])
    allele1_freq = np.zeros(geno_processed.shape[0])

    if len(non_sex_chr_indices) > 0:
        non_sex_chr_geno = geno_processed[non_sex_chr_indices, :]
        is_missing_non_sex_chr = (non_sex_chr_geno == 9)
        n_missing_non_sex_chr = np.sum(is_missing_non_sex_chr, axis=1)
        n_missing[non_sex_chr_indices] = n_missing_non_sex_chr
        valid_count[non_sex_chr_indices] = n_ind - n_missing_non_sex_chr
        missing_rate[non_sex_chr_indices] = n_missing_non_sex_chr / n_ind
        
        is_hom_ref = (non_sex_chr_geno == 2)
        is_het = (non_sex_chr_geno == 1)
        n_allele1 = (np.sum(is_hom_ref, axis=1) * 2) + np.sum(is_het, axis=1)
        n_chrom = (n_ind - n_missing_non_sex_chr) * 2
        n_chrom_safe = np.maximum(n_chrom, 1)
        allele1_freq[non_sex_chr_indices] = n_allele1 / n_chrom_safe

    if len(y_indices) > 0 and n_males > 0:
        for i, snp_idx in enumerate(y_indices):
            y_geno_males = geno_processed[snp_idx, male_mask]
            y_geno_males_clean = y_geno_males.copy()
            # convert heterozygous calls to missing
            y_geno_males_clean[y_geno_males_clean == 1] = 9
            is_missing_y = (y_geno_males_clean == 9)
            n_missing_y = np.sum(is_missing_y)
            n_valid_y = n_males - n_missing_y

            n_missing[snp_idx] = n_missing_y
            valid_count[snp_idx] = n_valid_y
            missing_rate[snp_idx] = n_missing_y / n_males if n_males > 0 else 1.0
            if n_valid_y > 0:
                n_allele1_y = np.sum(y_geno_males_clean == 2)
                total_chromosomes = n_valid_y  # NOT * 2 for haploid
                allele1_freq[snp_idx] = n_allele1_y / total_chromosomes
            else:
                allele1_freq[snp_idx] = 0.0

    if len(x_indices) > 0:
        for i, snp_idx in enumerate(x_indices):
            x_geno_all = geno_processed[snp_idx, :]
            
            x_geno_males = x_geno_all[male_mask] if n_males > 0 else np.array([])
            x_geno_females = x_geno_all[female_mask] if n_females > 0 else np.array([])
            x_geno_unknown = x_geno_all[unknown_mask] if n_unknown > 0 else np.array([])
            
            x_geno_males_clean = x_geno_males.copy() if len(x_geno_males) > 0 else np.array([])
            if len(x_geno_males_clean) > 0:
                x_geno_males_clean[x_geno_males_clean == 1] = 9
            
            x_geno_females_clean = x_geno_females.copy() if len(x_geno_females) > 0 else np.array([])
            
            if config.ignore_unknown:
                x_geno_unknown_clean = np.array([])
                include_unknown = False
            else:
                x_geno_unknown_clean = x_geno_unknown.copy() if len(x_geno_unknown) > 0 else np.array([])
                include_unknown = True

            male_missing = np.sum(x_geno_males_clean == 9) if len(x_geno_males_clean) > 0 else 0
            male_valid = n_males - male_missing
            
            if male_valid > 0:
                male_allele1 = np.sum(x_geno_males_clean == 2)
                male_allele2 = np.sum(x_geno_males_clean == 0)
                male_chromosomes = male_valid
            else:
                male_allele1 = male_allele2 = male_chromosomes = 0
            
            female_missing = np.sum(x_geno_females_clean == 9) if len(x_geno_females_clean) > 0 else 0
            female_valid = n_females - female_missing
            
            if female_valid > 0:
                female_allele1 = (np.sum(x_geno_females_clean == 2) * 2) + np.sum(x_geno_females_clean == 1)
                female_allele2 = (np.sum(x_geno_females_clean == 0) * 2) + np.sum(x_geno_females_clean == 1)
                female_chromosomes = female_valid * 2  # 2 chromosomes per female
            else:
                female_allele1 = female_allele2 = female_chromosomes = 0
            
            if include_unknown and n_unknown > 0:
                unknown_missing = np.sum(x_geno_unknown_clean == 9) if len(x_geno_unknown_clean) > 0 else 0
                unknown_valid = n_unknown - unknown_missing
                
                if unknown_valid > 0:
                    unknown_allele1 = (np.sum(x_geno_unknown_clean == 2) * 2) + np.sum(x_geno_unknown_clean == 1)
                    unknown_allele2 = (np.sum(x_geno_unknown_clean == 0) * 2) + np.sum(x_geno_unknown_clean == 1)
                    unknown_chromosomes = unknown_valid * 2
                else:
                    unknown_allele1 = unknown_allele2 = unknown_chromosomes = 0
            else:
                unknown_missing = unknown_valid = 0
                unknown_allele1 = unknown_allele2 = unknown_chromosomes = 0
            
            total_valid = male_valid + female_valid + unknown_valid
            total_missing = male_missing + female_missing + unknown_missing
            total_individuals = n_males + n_females + (n_unknown if include_unknown else 0)
            
            n_missing[snp_idx] = total_missing
            valid_count[snp_idx] = total_valid
            missing_rate[snp_idx] = total_missing / total_individuals if total_individuals > 0 else 1.0
            
            total_allele1 = male_allele1 + female_allele1 + unknown_allele1
            total_chromosomes = male_chromosomes + female_chromosomes + unknown_chromosomes
            
            if total_chromosomes > 0:
                allele1_freq[snp_idx] = total_allele1 / total_chromosomes
            else:
                allele1_freq[snp_idx] = 0.0

    if sex_chr_info and sex_chr_info['has_y_chr']:
        check_and_log_female_y_genotypes(geno, snp_chunk, ind_data, config, sex_chr_info)

    maf = np.minimum(allele1_freq, 1.0 - allele1_freq)

    return maf, n_missing, valid_count, missing_rate, allele1_freq

#### geno
def filter_by_geno(geno, snp_chunk, geno_threshold, config=None, ind_data=None):
    """
    Filter SNPs based on missing genotype rate threshold.
    """
    if (config is not None and not config.ignore_sex_chr and
            ind_data is not None and snp_chunk is not None):
        sex_chr_info = _get_sex_chr_info(snp_chunk, config)
        use_autosomal_path = not (sex_chr_info and sex_chr_info['has_sex_chr'])
    else:
        use_autosomal_path = True

    if use_autosomal_path:
        return _filter_by_geno_autosomal(geno, snp_chunk, geno_threshold, config)
    else:
        return _filter_by_geno_sex_chr(geno, snp_chunk, geno_threshold, config, ind_data)

def _filter_by_geno_autosomal(geno, snp_chunk, geno_threshold, config):
    """
    Filter SNPs by missing rate for autosomal chromosomes.
    """
    n_ind = geno.shape[1]
    is_missing = (geno == 9)
    n_missing = np.sum(is_missing, axis=1)
    missing_rate = n_missing / n_ind
    mask = missing_rate <= geno_threshold
    snp_chunk_filtered = apply_mask_and_log(snp_chunk, mask, "GENO", config)
    geno_filtered = apply_mask_and_log(geno, mask, "GENO", config, snp_chunk)
    return geno_filtered, snp_chunk_filtered

def _filter_by_geno_sex_chr(geno, snp_chunk, geno_threshold, config, ind_data):
    """
    Filter SNPs by missing rate for sex chromosomes.    
    """
    n_ind = geno.shape[1]
    sex_chr_info = _get_sex_chr_info(snp_chunk, config)
    geno_processed = apply_sex_chr_missing(geno, snp_chunk, ind_data, config, sex_chr_info)

    x_indices = sex_chr_info['x_indices'] if sex_chr_info else np.array([])
    y_indices = sex_chr_info['y_indices'] if sex_chr_info else np.array([])
    non_sex_chr_indices = sex_chr_info['non_sex_chr_indices'] if sex_chr_info else np.arange(geno.shape[0])
    
    male_mask = (ind_data['sex'] == 'M').values
    female_mask = (ind_data['sex'] == 'F').values
    unknown_mask = (~male_mask) & (~female_mask)
    
    n_males = np.sum(male_mask)
    n_females = np.sum(female_mask)
    n_unknown = np.sum(unknown_mask)

    if sex_chr_info and sex_chr_info['has_y_chr']:
        check_and_log_female_y_genotypes(geno, snp_chunk, ind_data, config, sex_chr_info)

    missing_rate = np.zeros(geno_processed.shape[0])

    if len(non_sex_chr_indices) > 0:
        non_sex_chr_geno = geno_processed[non_sex_chr_indices, :]
        is_missing_non_sex_chr = (non_sex_chr_geno == 9)
        n_missing_non_sex_chr = np.sum(is_missing_non_sex_chr, axis=1)
        missing_rate[non_sex_chr_indices] = n_missing_non_sex_chr / n_ind

    if len(y_indices) > 0 and n_males > 0:
        for i, snp_idx in enumerate(y_indices):
            y_geno_males = geno_processed[snp_idx, male_mask]
            y_geno_males_clean = y_geno_males.copy()
            y_geno_males_clean[y_geno_males_clean == 1] = 9
            is_missing_y = (y_geno_males_clean == 9)
            n_missing_y = np.sum(is_missing_y)
            missing_rate[snp_idx] = n_missing_y / n_males

    if len(x_indices) > 0:
        for i, snp_idx in enumerate(x_indices):
            x_geno_all = geno_processed[snp_idx, :]
            
            x_geno_males = x_geno_all[male_mask] if n_males > 0 else np.array([])
            x_geno_females = x_geno_all[female_mask] if n_females > 0 else np.array([])
            x_geno_unknown = x_geno_all[unknown_mask] if n_unknown > 0 else np.array([])
            
            x_geno_males_clean = x_geno_males.copy() if len(x_geno_males) > 0 else np.array([])
            if len(x_geno_males_clean) > 0:
                x_geno_males_clean[x_geno_males_clean == 1] = 9
            
            x_geno_females_clean = x_geno_females.copy() if len(x_geno_females) > 0 else np.array([])
            
            if config.ignore_unknown:
                include_unknown = False
            else:
                include_unknown = True
            
            male_missing = np.sum(x_geno_males_clean == 9) if len(x_geno_males_clean) > 0 else 0
            female_missing = np.sum(x_geno_females_clean == 9) if len(x_geno_females_clean) > 0 else 0
            
            if include_unknown and n_unknown > 0:
                unknown_missing = np.sum(x_geno_unknown == 9) if len(x_geno_unknown) > 0 else 0
                total_missing = male_missing + female_missing + unknown_missing
                total_individuals = n_males + n_females + n_unknown
            else:
                total_missing = male_missing + female_missing
                total_individuals = n_males + n_females
            
            missing_rate[snp_idx] = total_missing / total_individuals if total_individuals > 0 else 1.0

    mask = missing_rate <= geno_threshold

    snp_chunk_filtered = apply_mask_and_log(snp_chunk, mask, "GENO", config)
    geno_filtered = apply_mask_and_log(geno, mask, "GENO", config, snp_chunk)
    return geno_filtered, snp_chunk_filtered

#### missing
def calculate_ind_missing(geno, ind_data=None, snp_chunk=None, config=None):
    """
    Calculate per-individual missing genotype counts.
    Excludes Y for females, counts hets as missing on Y chr.
    """
    if (config is not None and not config.ignore_sex_chr and
            ind_data is not None and snp_chunk is not None):
        sex_chr_info = _get_sex_chr_info(snp_chunk, config)
        use_autosomal_path = not (sex_chr_info and sex_chr_info['has_sex_chr'])
    else:
        use_autosomal_path = True

    if use_autosomal_path:
        return _calculate_ind_missing_autosomal(geno)
    else:
        return _calculate_ind_missing_sex_chr(geno, ind_data, snp_chunk, config)

def _calculate_ind_missing_autosomal(geno):
    """
    Calculate per-individual missing counts for autosomal chromosomes.
    """
    is_valid = (geno != 9)
    n_snps = np.sum(is_valid, axis=0)
    n_missing = geno.shape[0] - n_snps
    return n_snps, n_missing

def _calculate_ind_missing_sex_chr(geno, ind_data, snp_chunk, config):
    """
    Calculate per-individual missing counts for sex chromosomes.
    
    Excludes Y SNPs from female counts, treats hets as missing on Y chr.
    """
    original_geno = geno
    sex_chr_info = _get_sex_chr_info(snp_chunk, config)

    male_mask = (ind_data['sex'] == 'M').values
    female_mask = (ind_data['sex'] == 'F').values
    unknown_mask = (~male_mask) & (~female_mask)
    
    n_males = np.sum(male_mask)
    n_females = np.sum(female_mask)
    n_unknown = np.sum(unknown_mask)

    if sex_chr_info and sex_chr_info['has_y_chr']:
        check_and_log_female_y_genotypes(geno, snp_chunk, ind_data, config, sex_chr_info)

        x_indices = sex_chr_info['x_indices']
        y_indices = sex_chr_info['y_indices']
        non_sex_chr_indices = sex_chr_info['non_sex_chr_indices']

        n_snps = np.zeros(original_geno.shape[1], dtype=int)
        n_missing = np.zeros(original_geno.shape[1], dtype=int)

        if np.any(male_mask):
            male_geno = original_geno[:, male_mask]

            if len(y_indices) > 0:
                for y_idx in y_indices:
                    y_geno_males = male_geno[y_idx, :]
                    male_geno[y_idx, :] = np.where(y_geno_males == 1, 9, y_geno_males)

            if len(x_indices) > 0:
                for x_idx in x_indices:
                    x_geno_males = male_geno[x_idx, :]
                    male_geno[x_idx, :] = np.where(x_geno_males == 1, 9, x_geno_males)

            is_valid_males = (male_geno != 9)
            n_snps_males = np.sum(is_valid_males, axis=0)
            n_missing_males = male_geno.shape[0] - n_snps_males

            n_snps[male_mask] = n_snps_males
            n_missing[male_mask] = n_missing_males

        if np.any(female_mask):
            female_valid_indices = np.concatenate([non_sex_chr_indices, x_indices]) if len(x_indices) > 0 else non_sex_chr_indices

            if len(female_valid_indices) > 0:
                geno_females_filtered = original_geno[female_valid_indices, :]
                is_valid_females = (geno_females_filtered != 9)
                n_snps_females = np.sum(is_valid_females, axis=0)
                n_missing_females = len(female_valid_indices) - n_snps_females

                n_snps[female_mask] = n_snps_females[female_mask]
                n_missing[female_mask] = n_missing_females[female_mask]

        if np.any(unknown_mask):
            if config.ignore_unknown:
                unknown_geno = original_geno[:, unknown_mask]
                if len(y_indices) > 0:
                    for y_idx in y_indices:
                        y_geno_unknown = unknown_geno[y_idx, :]
                        unknown_geno[y_idx, :] = np.where(y_geno_unknown == 1, 9, y_geno_unknown)
                if len(x_indices) > 0:
                    for x_idx in x_indices:
                        x_geno_unknown = unknown_geno[x_idx, :]
                        unknown_geno[x_idx, :] = np.where(x_geno_unknown == 1, 9, x_geno_unknown)
                is_valid_unknown = (unknown_geno != 9)
                n_snps_unknown = np.sum(is_valid_unknown, axis=0)
                n_missing_unknown = unknown_geno.shape[0] - n_snps_unknown
                n_snps[unknown_mask] = n_snps_unknown
                n_missing[unknown_mask] = n_missing_unknown
            else:
                unknown_valid_indices = np.concatenate([non_sex_chr_indices, x_indices]) if len(x_indices) > 0 else non_sex_chr_indices
                if len(unknown_valid_indices) > 0:
                    geno_unknown_filtered = original_geno[unknown_valid_indices, :]
                    is_valid_unknown = (geno_unknown_filtered != 9)
                    n_snps_unknown = np.sum(is_valid_unknown, axis=0)
                    n_missing_unknown = len(unknown_valid_indices) - n_snps_unknown
                    n_snps[unknown_mask] = n_snps_unknown[unknown_mask]
                    n_missing[unknown_mask] = n_missing_unknown[unknown_mask]
        return n_snps, n_missing

    is_valid = (original_geno != 9)
    n_snps = np.sum(is_valid, axis=0)
    n_missing = original_geno.shape[0] - n_snps
    return n_snps, n_missing

#### freq
def calculate_snp_stats(geno, snp_chunk, ind_data=None, config=None):
    """
    Calculate SNP statistics.
    """
    if (config is not None and not config.ignore_sex_chr and
            ind_data is not None and snp_chunk is not None):
        sex_chr_info = _get_sex_chr_info(snp_chunk, config)
        use_autosomal_path = not (sex_chr_info and sex_chr_info['has_sex_chr'])
    else:
        use_autosomal_path = True

    if use_autosomal_path:
        return _calculate_snp_stats_autosomal(geno, snp_chunk)
    else:
        return _calculate_snp_stats_sex_chr(geno, snp_chunk, ind_data, config)

def _calculate_snp_stats_autosomal(geno, snp_chunk):
    maf, n_missing, valid_count, missing_rate, allele1_freq = _calculate_maf_autosomal(geno, geno.shape[1])
    allele2_freq = 1.0 - allele1_freq
    minor_allele = np.where(allele1_freq <= 0.5, snp_chunk['allele1'].values, snp_chunk['allele2'].values)
    major_allele = np.where(allele1_freq > 0.5, snp_chunk['allele1'].values, snp_chunk['allele2'].values)
    n_total_per_snp = np.full(len(snp_chunk), geno.shape[1], dtype=int)

    stats_df = pd.DataFrame({
        'chr': snp_chunk['chr'],
        'snpID': snp_chunk['snpID'],
        'n_snps': valid_count.astype(int),
        'n_missing': n_missing.astype(int),
        'n_total': n_total_per_snp,
        'missing_rate': np.round(missing_rate, 5),
        'allele1': snp_chunk['allele1'],
        'allele2': snp_chunk['allele2'],
        'allele1_freq': np.round(allele1_freq, 5),
        'allele2_freq': np.round(allele2_freq, 5),
        'major_allele': major_allele,
        'minor_allele': minor_allele,
        'maf': np.round(maf, 5)
    })
    return stats_df, maf

def _calculate_snp_stats_sex_chr(geno, snp_chunk, ind_data, config):
    sex_chr_info = _get_sex_chr_info(snp_chunk, config)
    maf, n_missing, valid_count, missing_rate, allele1_freq = calculate_maf(geno, ind_data, snp_chunk, config, sex_chr_info)
    allele2_freq = 1.0 - allele1_freq
    minor_allele = np.where(allele1_freq <= 0.5, snp_chunk['allele1'].values, snp_chunk['allele2'].values)
    major_allele = np.where(allele1_freq > 0.5, snp_chunk['allele1'].values, snp_chunk['allele2'].values)
    n_total_per_snp = np.full(len(snp_chunk), geno.shape[1], dtype=int)  # Default: all individuals

    if sex_chr_info and sex_chr_info['has_y_chr']:
        male_mask = (ind_data['sex'] == 'M').values
        female_mask = (ind_data['sex'] == 'F').values
        unknown_mask = (~male_mask) & (~female_mask)
        
        n_males = np.sum(male_mask)
        n_females = np.sum(female_mask)
        n_unknown = np.sum(unknown_mask)
        
        n_total_per_snp[sex_chr_info['y_indices']] = n_males
        
        if sex_chr_info['has_x_chr'] and config.ignore_unknown:
            n_total_per_snp[sex_chr_info['x_indices']] = n_males + n_females

    stats_df = pd.DataFrame({
        'chr': snp_chunk['chr'],
        'snpID': snp_chunk['snpID'],
        'n_snps': valid_count.astype(int),
        'n_missing': n_missing.astype(int),
        'n_total': n_total_per_snp,
        'missing_rate': np.round(missing_rate, 5),
        'allele1': snp_chunk['allele1'],
        'allele2': snp_chunk['allele2'],
        'allele1_freq': np.round(allele1_freq, 5),
        'allele2_freq': np.round(allele2_freq, 5),
        'major_allele': major_allele,
        'minor_allele': minor_allele,
        'maf': np.round(maf, 5)
    })
    return stats_df, maf

def apply_mask_and_log(df, mask, reason, config, log_df=None):
    is_array = isinstance(df, np.ndarray)
    filtered = df[mask] if is_array else df.iloc[mask].reset_index(drop=True)
    if config and hasattr(config, 'verbose') and config.verbose and hasattr(config, 'removed_snp_log') and config.removed_snp_log:
        removed_indices = np.where(~mask)[0]
        if len(removed_indices) > 0:
            source = log_df if is_array and log_df is not None else df
            try:
                with open(config.removed_snp_log, 'a') as f:
                    for idx in removed_indices:
                        if idx >= len(source): continue
                        snp = source.iloc[idx] if hasattr(source, 'iloc') else source[idx]
                        snp_id = snp['snpID'] if isinstance(snp, dict) or hasattr(snp,'__getitem__') else snp[0] if isinstance(snp, (list,tuple,np.ndarray)) else snp
                        f.write(f"{snp_id}\t{reason}\n")
            except Exception as e:
                logging.warning(f"Error writing to removed SNP log: {e}")
    return filtered

#### logging
def update_polarisation_counts(config, n_polarised, n_processed, excluded_count=0):
    if not hasattr(config, '_polarisation_counts'):
        config._polarisation_counts = {
            'total_polarised': 0,
            'total_processed': 0,
            'total_excluded': 0}
    
    config._polarisation_counts['total_polarised'] += n_polarised
    config._polarisation_counts['total_processed'] += n_processed
    config._polarisation_counts['total_excluded'] += excluded_count

def log_polarisation(config):
    if hasattr(config, '_polarisation_counts'):
        counts = config._polarisation_counts
        total_polarised = counts['total_polarised']
        total_processed = counts['total_processed']
        total_excluded = counts['total_excluded']
        
        if total_processed > 0:
            logging.info(f"Polarisation summary: {total_polarised} SNPs polarised, "
                        f"{total_processed} SNPs processed, {total_excluded} SNPs excluded due to missing/heterozygous reference")

def log_sample_sex_counts(ind_data, logger=logging):
    total_samples = len(ind_data)
    male_count = (ind_data['sex'] == 'M').sum()
    female_count = (ind_data['sex'] == 'F').sum()
    ambiguous_count = total_samples - male_count - female_count
    logger.info(f"{total_samples} samples ({male_count} males, {female_count} females, {ambiguous_count} ambiguous) loaded from .ind.")

def log_snp_counts_with_sex(snp_data, config, logger=logging):
    total_snps = len(snp_data)
    has_sex_chr = detect_sex_chromosomes(snp_data, config)
    if has_sex_chr:
        logger.info(f"{total_snps} SNPs loaded from .snp file and sex chromosomes detected in the dataset")
    else:
        logger.info(f"{total_snps} SNPs loaded from .snp file and no sex chromosomes detected in the dataset")
    return has_sex_chr

def log_missingness(miss_df, logger=logging):
    if miss_df is None or miss_df.empty or 'missing_rate' not in miss_df:
        logger.info('Missingness summary: No valid individuals for missing rate calculation')
        return
    valid_missing_rate = miss_df['missing_rate'].dropna()
    if len(valid_missing_rate) > 0:
        logger.info(f"Missingness summary:")
        logger.info(f"  Missing rate:")
        logger.info(f"    Min: {valid_missing_rate.min():.4f}")
        logger.info(f"    Max: {valid_missing_rate.max():.4f}")
        logger.info(f"    Mean: {valid_missing_rate.mean():.4f}")
    else:
        logger.info('Missingness summary: No valid individuals for missing rate calculation')

def log_chr_missingness(chr_id, missing_rate_chr, y_chr_id, logger=logging):
    valid_missing_rate_chr = missing_rate_chr[~np.isnan(missing_rate_chr)]
    if len(valid_missing_rate_chr) > 0:
        if str(chr_id) == str(y_chr_id):
            logger.info(f"Chromosome {chr_id} (Y) missingness summary (males only):")
        else:
            logger.info(f"Chromosome {chr_id} missingness summary:")
        logger.info(f"  Missing rate - Min: {valid_missing_rate_chr.min():.4f}, Max: {valid_missing_rate_chr.max():.4f}, Mean: {valid_missing_rate_chr.mean():.4f}")

def log_freq(combined_stats, logger=logging):
    if combined_stats is None or combined_stats.empty:
        logger.info("SNP statistics summary: No statistics calculated.")
        return
    logger.info(f"SNP statistics summary:")
    logger.info(f"  Total SNPs: {len(combined_stats)}")
    logger.info(f"  Missing rate:")
    logger.info(f"    Min: {combined_stats['missing_rate'].min():.4f}")
    logger.info(f"    Max: {combined_stats['missing_rate'].max():.4f}")
    logger.info(f"    Mean: {combined_stats['missing_rate'].mean():.4f}")
    logger.info(f"  MAF:")
    logger.info(f"    Min: {combined_stats['maf'].min():.4f}")
    logger.info(f"    Max: {combined_stats['maf'].max():.4f}")
    logger.info(f"    Mean: {combined_stats['maf'].mean():.4f}")
