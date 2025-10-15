import argparse
import logging
from .config import EigenConfig
from .run import EigenDataset
from .processors import (
    filter_individuals,
    filter_snps,
    random_haploidise,
    update_individual_populations,
    update_snp_ids,
    load_genetic_map,
    update_genetic_distances,
    polarise,
    complement_snp_alleles,
    apply_update)
from .io_handlers import read_ind_file, read_snp_file
import time
import os
import shutil
from .utils import get_memory_usage, log_runtime
import pandas as pd

__version__ = "1.0"

def parse_args():
    parser = argparse.ArgumentParser(description='A Python Package for EIGENSTRAT Data Processing')
    parser.add_argument('--version', action='version', version=f'pygenstrat {__version__}')
    
    # io
    parser.add_argument('--prefix', '-p',
                        help='Input EIGENSTRAT prefix for GENO/SNP/IND files')
    parser.add_argument('--out', '-o', help='Output file prefix')
    parser.add_argument('--out-type', choices=['anc', 'eig'], default='anc',
                        help='Output data type: anc (binary format) or eig (text format) (default: anc)')
    parser.add_argument('--chunk-size', type=int, default=1000,
                        help='Chunk size for processing (default: 1000)')

    # filtering/stats
    parser.add_argument('--min-maf', '-maf', type=float,
                        help='Minimum minor allele frequency')
    parser.add_argument('--max-maf', type=float,
                        help='Maximum minor allele frequency')
    parser.add_argument('--geno', type=float,
                        help='Maximum missing call rates per-variant')
    parser.add_argument('--missing', action='store_true',
                        help='Calculate missingness per-sample')
    parser.add_argument('--missing-by-chr', action='store_true',
                        help='Calculate missingness per-sample for each chromosome separately')
    parser.add_argument('--freq', action='store_true', 
                        help='Calculate allele frequencies and missingness per-variant')

    # update and/or transform
    parser.add_argument('--update-ind', type=str,
                        help='File containing individual and new population names for updating population names')
    parser.add_argument('--update-snp', type=str,
                        help='File containing old SNP IDs and new SNP IDs for updating SNP identifiers')
    parser.add_argument('--genetic-distance', type=str,
                        help='Path to genetic distance map file containing data for all chromosomes, or use "zero" to set all genetic distances to 0')
    parser.add_argument('--map-unit', choices=['cM', 'M'], default='cM',
                        help='Output unit for genetic distances: centiMorgans (cM) or Morgans (M) (default: cM, requires --genetic-distance if set)')
    parser.add_argument('--flip-strand', type=str,
                        help='File containing SNP IDs (one per line) whose alleles should be complemented (e.g., A<->T, C<->G) in the SNP file.')
    parser.add_argument('--random-haploidise', action='store_true',
                        help='Convert heterozygous genotypes to homozygous at random')
    parser.add_argument('--polarise', type=str,
                        help='Polarisation reference: either a two-column file mapping snpID to ancestral_allele, or a sample ID from the dataset to use as ancestral reference (makes allele1 ancestral and allele2 derived).')

    # sex-chromosome specific options
    parser.add_argument('--sex-chr', type=str, nargs='+', default=['23', '24'],
                        help='Sex chromosome identifiers (default: 23 24). Can also use X Y or custom values')
    parser.add_argument('--ignore-sex-chr', action='store_true',
                        help='Ignore sex chromosome information for faster processing when sex chromosomes are not included')
    parser.add_argument('--sex-chr-missing', action='store_true',
                        help='Set female genotypes to missing for Y chromosome SNPs and male heterozygous genotypes to missing for X chromosome SNPs')
    parser.add_argument('--ignore-unknown', action='store_true',
                        help='Ignore individuals with unknown sex in sex chromosome calculations (MAF, missing rate, etc.)')

    # filter indv/pop
    ind_filter_group = parser.add_mutually_exclusive_group()
    ind_filter_group.add_argument('--keep-indv', type=str,
                                  help='File containing individual IDs to keep')
    ind_filter_group.add_argument('--remove-indv', type=str,
                                  help='File containing individual IDs to remove')
    ind_filter_group.add_argument('--keep-pop', type=str,
                                  help='File containing population names to keep')
    ind_filter_group.add_argument('--remove-pop', type=str,
                                  help='File containing population names to remove')

    # filter snp/region
    snp_filter_group = parser.add_mutually_exclusive_group()
    snp_filter_group.add_argument('--keep-snps', type=str,
                                  help='File containing SNP IDs to keep')
    snp_filter_group.add_argument('--remove-snps', type=str,
                                  help='File containing SNP IDs to remove')
    snp_filter_group.add_argument('--keep-region', type=str,
                                  help='File containing genomic regions to keep (chr start_pos end_pos)')
    snp_filter_group.add_argument('--remove-region', type=str,
                                  help='File containing genomic regions to remove (chr start_pos end_pos)')

    # filter chr
    chr_filter_group = parser.add_mutually_exclusive_group()
    chr_filter_group.add_argument('--keep-chr', type=str,
                                help='File containing chromosome IDs to keep')
    chr_filter_group.add_argument('--remove-chr', type=str,
                                help='File containing chromosome IDs to remove')

    # verbose
    parser.add_argument('--verbose', action='store_true',
                        help='Output detailed logs about removed SNPs and individuals')

    args = parser.parse_args()
    return args

def main():

    print("pygenstrat running")

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s')
    
    start_time = time.time()
    start_memory = get_memory_usage()
    
    try:
        args = parse_args()
        args_dict = vars(args)
        
        has_updates = args.update_ind or args.update_snp or args.genetic_distance or args.flip_strand
        has_stats = args.missing or args.missing_by_chr or args.freq
        has_filters = any([
            args.min_maf is not None,
            args.max_maf is not None,
            args.geno is not None,
            args.keep_indv is not None,
            args.remove_indv is not None,
            args.keep_pop is not None,
            args.remove_pop is not None,
            args.keep_snps is not None,
            args.remove_snps is not None,
            args.keep_region is not None,
            args.remove_region is not None,
            args.keep_chr is not None,
            args.remove_chr is not None,
            args.random_haploidise,
            args.polarise,
            args.sex_chr_missing])
        is_simple_conversion = args.out_type is not None and not has_updates and not has_stats and not has_filters
        should_write_output = has_filters or is_simple_conversion or (has_updates and args.out is not None)
        
        if is_simple_conversion and args.out is None:
            args.out = 'pygenstrat'
            args_dict['out'] = 'pygenstrat'
        
        config = EigenConfig(**args_dict)

        if not should_write_output:
            config.out_geno = None
            config.out_snp = None
            config.out_ind = None

        interpolators = None
        if args.genetic_distance and args.genetic_distance != "zero":
            logging.info(f"Loading genetic distance map from file: {args.genetic_distance}")
            _, interpolators = load_genetic_map(args.genetic_distance)
            config.interpolators = interpolators
        
        if has_updates and not args.out and not has_stats and not has_filters:
            logging.info("Running in update-only (in-place) mode: applying update(s) only.")
            ind_data = read_ind_file(config.ind_path)
            snp_data = read_snp_file(config.snp_path)
            _, _, changes_made = apply_update(
                ind_data, snp_data, config, in_place=True)
            if changes_made:
                log_runtime(start_time, start_memory)
                return 0
            else:
                logging.warning("No updates were applied.")
                log_runtime(start_time, start_memory)
                return 1
        
        if not has_updates and not has_stats and not has_filters and not is_simple_conversion:
            logging.error("No operations specified. Please specify filtering options, update options, or stats options.")
            log_runtime(start_time, start_memory)
            return 1
        
        try:
            if not all(os.path.exists(f) for f in [config.geno_path, config.snp_path, config.ind_path]):
                missing = [f for f in [config.geno_path, config.snp_path, config.ind_path] if not os.path.exists(f)]
                raise FileNotFoundError(f"Required files not found: {', '.join(missing)}")
            
            if has_updates and args.out:
                logging.info("Updates will be applied to output files")
                
                if args.update_ind:
                    config.update_ind = args.update_ind
                if args.update_snp:
                    config.update_snp = args.update_snp
                if args.genetic_distance:
                    config.genetic_distance = args.genetic_distance
                if args.flip_strand:
                    config.flip_strand = args.flip_strand

            EigenDataset(config).process_chunks()
            
        except FileNotFoundError as e:
            logging.error(f"File error: {e}")
            log_runtime(start_time, start_memory)
            return 2
        except MemoryError:
            logging.error("Not enough memory to complete operation. Try reducing chunk size.")
            log_runtime(start_time, start_memory)
            return 3
        except Exception as e:
            logging.error(f"Processing error: {e}")
            log_runtime(start_time, start_memory)
            return 4

        log_runtime(start_time, start_memory)
        
        return 0

    except Exception as e:
        logging.error(f"Error in main: {e}")
        log_runtime(start_time, start_memory)
        raise

if __name__ == "__main__":
    main()