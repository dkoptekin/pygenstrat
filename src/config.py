from dataclasses import dataclass
import os
import logging

@dataclass
class EigenConfig:
    # inputs
    prefix: str = None
    out: str = None
    out_type: str = "anc"
    chunk_size: int = 1000
    min_maf: float = None
    max_maf: float = None
    geno: float = None
    missing: bool = False
    missing_by_chr: bool = False
    verbose: bool = False
    freq: bool = False
    keep_indv: str = None
    remove_indv: str = None
    keep_pop: str = None
    remove_pop: str = None
    keep_snps: str = None
    remove_snps: str = None
    keep_region: str = None
    remove_region: str = None
    keep_chr: str = None
    remove_chr: str = None
    update_ind: str = None
    update_snp: str = None
    genetic_distance: str = None
    map_unit: str = "cM"
    random_haploidise: bool = False
    polarise: str = None
    sex_chr: list = None
    ignore_sex_chr: bool = False
    sex_chr_missing: bool = False
    ignore_unknown: bool = False
    flip_strand: str = None

    # outputs
    out_geno: str = None
    out_snp: str = None
    out_ind: str = None
    out_missing: str = None
    out_missing_by_chr: str = None
    removed_snp_log: str = None
    removed_ind_log: str = None
    out_freq: str = None

    def __post_init__(self):
        if not self.prefix:
            raise ValueError("Input prefix must be provided (--prefix)")

        if self.sex_chr is None:
            self.sex_chr = ['23', '24']

        rename_sex_chr = []
        for chr_id in self.sex_chr:
            chr_str = str(chr_id).upper().replace('CHR', '')
            if chr_str == 'X':
                rename_sex_chr.append('23')
            elif chr_str == 'Y':
                rename_sex_chr.append('24')
            else:
                rename_sex_chr.append(chr_str)
        self.sex_chr = rename_sex_chr

        self.input_type = 'eigenstrat'
        self.geno_path = f"{self.prefix}.geno"
        self.snp_path = f"{self.prefix}.snp"
        self.ind_path = f"{self.prefix}.ind"
        
        logging.info(f"Input EIGENSTRAT files:")
        logging.info(f"  GENO: {self.geno_path}")
        logging.info(f"   SNP: {self.snp_path}")
        logging.info(f"   IND: {self.ind_path}")
        
        if self.out:
            self.out_base = self.out
        else:
            self.out_base = 'pygenstrat'
            
        self.out_geno = f"{self.out_base}.geno"
        self.out_snp = f"{self.out_base}.snp"
        self.out_ind = f"{self.out_base}.ind"
        self.out_missing = f"{self.out_base}.missing" if self.missing else None
        self.out_missing_by_chr = f"{self.out_base}" if self.missing_by_chr else None
        self.removed_snp_log = f"{self.out_base}.removed.snp" if self.verbose else None
        self.removed_ind_log = f"{self.out_base}.removed.ind" if self.verbose else None
        self.out_freq = f"{self.out_base}.frq" if self.freq else None

        if self.verbose:
            logging.info(f"Verbose logging enabled: removed SNPs will be written to {self.removed_snp_log}")
            logging.warning(
                "NOTE: SNPs are filtered in sequence (MAF → GENO → SNP_ID → REGION → CHR). "
                "Each SNP is only logged for the first filter that removes it.")

    def validate_files(self):
        _check_file_exists(self.geno_path, "GENO file")
        _check_file_exists(self.snp_path, "SNP file")
        _check_file_exists(self.ind_path, "IND file")
        _check_file_exists(self.keep_indv, "Keep individuals file")
        _check_file_exists(self.remove_indv, "Remove individuals file")
        _check_file_exists(self.keep_pop, "Keep populations file")
        _check_file_exists(self.remove_pop, "Remove populations file")
        _check_file_exists(self.keep_snps, "Keep SNPs file")
        _check_file_exists(self.remove_snps, "Remove SNPs file")
        _check_file_exists(self.update_snp, "SNP update file")
        _check_file_exists(self.update_ind, "IND update file")
        _check_file_exists(self.genetic_distance, "Genetic distance file", allow_zero=True)
        _check_file_exists(self.flip_strand, "Flip strand file")

    def validate_parameters(self):  
        if self.out_type and self.out_type not in ["anc", "eig"]:
            raise ValueError(f"out_type must be either 'anc' (binary) or 'eig' (text)")
            
        if self.min_maf is not None:
            if self.min_maf < 0 or self.min_maf > 0.5:
                raise ValueError(f"min-maf must be between 0 and 0.5")
        
        if self.max_maf is not None:
            if self.max_maf < 0 or self.max_maf > 0.5:
                raise ValueError(f"max-maf must be between 0 and 0.5")
        
        if self.min_maf is not None and self.max_maf is not None:
            if self.min_maf > self.max_maf:
                raise ValueError(f"min-maf ({self.min_maf}) cannot be greater than max-maf")

        if self.geno is not None:
            if self.geno < 0 or self.geno > 1:
                raise ValueError(f"geno must be between 0 and 1")
        
        if self.map_unit and self.map_unit not in ["cM", "M"]:
            raise ValueError(f"map_unit must be either 'cM' or 'M'")
        
        if self.map_unit and not self.genetic_distance:
            raise ValueError("--map-unit can only be used together with --genetic-distance (must specify a genetic map file or 'zero')")
        
        if self.genetic_distance and self.map_unit == "cM":
            logging.info("Genetic distance output unit: centiMorgans (cM)")
        elif self.genetic_distance and self.map_unit == "M":
            logging.info("Genetic distance output unit: Morgans (M)")

        if self.keep_region and self.remove_region:
            raise ValueError("Cannot specify both --keep-region and --remove-region")
        
        if self.keep_region and self.keep_snps:
            raise ValueError("Cannot specify both --keep-region and --keep-snps")
        
        if self.remove_region and self.remove_snps:
            raise ValueError("Cannot specify both --remove-region and --remove-snps")
        
        if self.keep_chr and self.remove_chr:
            raise ValueError("Cannot specify both --keep-chr and --remove-chr")
        
        if self.keep_indv and self.remove_indv:
            raise ValueError("Cannot specify both --keep-indv and --remove-indv")
        
        if self.keep_pop and self.remove_pop:
            raise ValueError("Cannot specify both --keep-pop and --remove-pop")

        # Sex chromosome options info/warnings
        if self.sex_chr_missing and self.ignore_sex_chr:
            logging.warning(
                "Both --sex-chr-missing and --ignore-sex-chr specified. Sex chromosome filtering will be ignored due to --ignore-sex-chr.")

        if self.sex_chr_missing:
            logging.info(
                "Sex chromosome filtering enabled: female genotypes will be set to missing for Y chromosome SNPs and male heterozygous genotypes will be set to missing for X chromosome SNPs")


def _check_file_exists(file_path, description, allow_zero=False):
    if file_path and file_path != "zero" and not os.path.exists(file_path):
        raise ValueError(f"{description} not found: {file_path}")