# pygenstrat

A Python package for handling and processing EIGENSTRAT datasets with advanced performance features and configurable options.

## Overview

`pygenstrat` allows users to efficiently manage large EIGENSTRAT datasets, providing tools for format conversion,
filtering, and statistical analysis. With optimised I/O operations, it handles large genomic datasets with high
performance.

## EIGENSTRAT File Formats

`pygenstrat` supports both **EIGEN format** (text-based) and **ANC format** (binary)

### File Structure

**Genotype file (.geno)**

```bash
000212110991...  <- SNP 1 (one character per individual)
120210099110...  <- SNP 2
...
...
```

Each line represents one SNP, with values:
- 0: Homozygous for alternate allele (zero copies of reference allele)
- 1: Heterozygous (one copy of reference allele)
- 2: Homozygous for reference allele (two copies of reference allele)
- 9: Missing data

**Individual file (.ind)**

```bash
Sample1 M Population1
Sample2 F Population1
Sample3 M Population2
...
...
```

Each line contains sample ID, sex (M/F), and population name.

**SNP file (.snp)**

```bash
rs123 1 0.000 1000000 A G
rs456 1 0.010 1050000 C T
...
```

Each line contains SNP ID, chromosome, genetic distance, physical position, and alleles.

For more detail about the Eigenstrat format, you can check [here](https://reich.hms.harvard.edu/software/InputFileFormats).

## Features

- **Format Conversion**: Convert between ANC (binary) and EIGEN (text) formats
- **Filtering Options**:
    - Individual/population-based filtering
    - SNP filtering
    - Minor allele frequency thresholds
    - Missing data rate thresholds for SNPs
    - Region-based filtering
    - Chromosome-based filtering
- **Statistical Analysis**:
    - Allele frequency calculation
    - Missing rates for SNPs
    - Per-individual missing rates
- **Updating Files**:
    - Update population information for individuals
    - Update SNP IDs
    - Update genetic distances using genetic distance maps (cM or M) or set to zero
    - Flip strand
- **Genotype Management**:
    - Random haploidisation of heterozygous genotypes
    - Polarisation using file-based or sample-based ancestral reference
- **Sex Chromosome Handling**:
    - User-configurable sex chromosome identifiers (default '23' for X, '24' for Y; accepts 'X', 'Y', or custom values)
    - Sex-aware filtering and missing-data handling:
        - Exclude females from Y chromosome SNPs
        - Mask heterozygous male genotypes on the X chromosome (set to missing)
        - Option to ignore all sex-specific logic
      - **Handling of individuals with unknown sex:** unknown-sex samples are always treated like females: 
        - they are included as diploid for X chromosome analysis and statistics, and not included in Y chromosome
          statistics 
        - If `--ignore-unknown` is used they are excluded from all sex chromosome-specific computations.
    - Comprehensive support for human/non-human datasets with custom sex chromosome schemes
    - Enables robust X and Y chromosome-specific quality control, summary statistics, and population analysis

## Installation

### Option 1: Direct Installation

1. Install the package:

```bash
git clone git@github.com:dkoptekin/pygenstrat.git

cd pygenstrat

pip install --user .
```

### Option 2: Using Virtual Environment

1. Create and activate a virtual environment:

```bash
python -m venv pygenstrat_env

source pygenstrat_env/bin/activate
```

2. Install the package:

```bash
git clone git@github.com:dkoptekin/pygenstrat.git

cd pygenstrat

pip install .
```

### Path Configuration

If you see the warning that `pygenstrat` is installed in a directory not on PATH:

```bash
# For bash users (Linux/Mac)
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# For zsh users (Linux/Mac)
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

## Command Line Options

| Parameter             | Description                                                                                                                                                                                     | Type                     | Default      |
|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------|--------------|
| `--prefix`            | Input file prefix for .geno/.snp/.ind files                                                                                                                                                     | /path/to/input/basename  | Required     |
| `--out`               | Output file prefix                                                                                                                                                                              | /path/to/output/basename | 'pygenstrat' |
| `--out-type`          | Output format type ('anc' or 'eigen')                                                                                                                                                           | choice                   | 'anc'        |
| `--chunk-size`        | Number of SNPs to process at once                                                                                                                                                               | int                      | 1000         |
| `--keep-indv`         | File containing individual IDs to keep                                                                                                                                                          | filename                 | None         |
| `--remove-indv`       | File containing individual IDs to remove                                                                                                                                                        | filename                 | None         |
| `--keep-pop`          | File containing population IDs to keep                                                                                                                                                          | filename                 | None         |
| `--remove-pop`        | File containing population IDs to remove                                                                                                                                                        | filename                 | None         |
| `--keep-snps`         | File containing SNP IDs to keep                                                                                                                                                                 | filename                 | None         |
| `--remove-snps`       | File containing SNP IDs to remove                                                                                                                                                               | filename                 | None         |
| `--keep-region`       | Regions to keep. Can be a file containing regions (chr start end) or direct specification like '1:1000-2000'                                                                                    | filename or str          | None         |
| `--remove-region`     | Regions to remove. Can be a file containing regions (chr start end) or direct specification like '1:1000-2000'                                                                                  | filename or str          | None         |
| `--keep-chr`          | Chromosomes to keep. Can be a file containing chromosome IDs or direct specification like '1,2,23' or range '1-22'                                                                              | filename or str          | None         |
| `--remove-chr`        | Chromosomes to remove. Can be a file containing chromosome IDs or direct specification like '23,24' or range '1-22'                                                                             | filename or str          | None         |
| `--min-maf`           | Minimum minor allele frequency threshold                                                                                                                                                        | float                    | None         |
| `--max-maf`           | Maximum minor allele frequency threshold                                                                                                                                                        | float                    | None         |
| `--geno`              | Maximum SNP missing rate                                                                                                                                                                        | float                    | None         |
| `--sex-chr`           | Sex chromosome identifiers (e.g., '23 24' or 'X Y')                                                                                                                                             | list                     | ['23', '24'] |
| `--ignore-sex-chr`    | Ignore sex chromosome information, for faster processing when the dataset only includes autosomal SNPs                                                                                          | flag                     | False        |
| `--sex-chr-missing`   | Set female genotypes to missing for Y chromosome SNPs and male heterozygous genotypes to missing for X chromosome SNPs                                                                          | flag                     | False        |
| `--missing`           | Calculate per-individual missing rate                                                                                                                                                           | flag                     | False        |
| `--missing-by-chr`    | Calculate per-individual missing rate for each chromosome separately (writes one file per chromosome)                                                                                           | flag                     | False        |
| `--freq`              | Calculate allele frequencies and missing rates                                                                                                                                                  | flag                     | False        |
| `--verbose`           | Enable verbose output                                                                                                                                                                           | flag                     | False        |
| `--update-ind`        | File containing individual and new population names                                                                                                                                             | filename                 | None         |
| `--update-snp`        | File containing old SNP IDs and new SNP IDs                                                                                                                                                     | filename                 | None         |
| `--genetic-distance`  | Path to genetic distance map file with chr/position/cM data OR 'zero' to set all cM to 0                                                                                                        | filename or 'zero'       | None         |
| `--map-unit`          | Unit for genetic distances ('cM' or 'M')                                                                                                                                                        | choice                   | 'cM'         |
| `--random-haploidise` | Convert heterozygous genotypes to homozygous at random                                                                                                                                          | flag                     | False        |
| `--polarise`          | Polarisation reference: either a two-column file mapping snpID to ancestral_allele, or a sample ID from the dataset to use as ancestral reference (makes allele1 ancestral and allele2 derived) | filename or str          | None         |
| `--flip-strand`       | File containing SNPs to convert to their complementary bases (swapped to the opposite DNA strand)                                                                                               | filename                 | None         |

## Processing Order

When multiple options are used, `pygenstrat` applies them in this order:

1. **Update operations**
    - Update individual populations (--update-ind)
    - Update SNP IDs (--update-snp)
    - Update genetic distances (--genetic-distance)
    - Flip-strand (--flip-strand)

2. **Sample-level filters**
    - Individual filtering (--keep-indv/--remove-indv)
    - Population filtering (--keep-pop/--remove-pop)

3. **Chunk-level processing** (applied per chunk during data processing):
    - Polarisation (--polarise)
    - Random haploidisation (--random-haploidise)
    - Sex chromosome missing filtering (--sex-chr-missing)
    - Sample filtering (--keep-indv/--remove-indv)
    - SNP ID filtering (--keep-snps/--remove-snps)
    - Region filtering (--keep-region/--remove-region)
    - Chromosome filtering (--keep-chr/--remove-chr)
    - Missing genotype rate filtering (--geno)
    - Minor allele frequency filtering (--min-maf/--max-maf)

> NOTE: 
> 
> When using **only update options** (`--update-ind`, `--update-snp`, `--genetic-distance`), the original files will be
updated in place and backups will be created with `.backup` extension.
> 
>When combining update options with filtering options and specifying an output prefix (--out), `pygenstrat` first applies the updates to the data in memory, then applies filters, and finally writes the output files. The original input files remain unchanged.

## Input File Formats

### Required Input Files

`pygenstrat` requires three input files with the same prefix (specified with `--prefix`):

- `[prefix].geno`: Genotype data file
- `[prefix].snp`: SNP information file
- `[prefix].ind`: Individual/sample information file

### Update Input Files

When using update options, `pygenstrat` expects specific file formats:

**Individual Population Update File** (`--update-ind`): Each line contains the individual ID and the new population name, separated by whitespace.

```
sample1 NewPopulationLabel
sample2 NewPopulationLabel
...
```

**SNP ID Update File** (`--update-snp`): Each line contains the old SNP ID and the new SNP ID, separated by whitespace.

```
old_snp_id1 new_snp_id1
old_snp_id2 new_snp_id2
...
```

**Genetic Map File** (`--genetic-map`): The file should contain chromosome IDs, physical positions, and genetic map
positions. The programme uses the Chromosome, Position(bp), and Map(cM) columns.

```
Chromosome Positionbp RatecM/Mb MapcM
chr1 55550 2.981822 0.000000
chr1 82571 2.082414 0.080572
...
```

**Flip-Strand File** (`--flip-strand`): Each line contains the SNP ID.

```
rs123 
rs456 
...
```

**Polarisation File** (`--polarise` with a file): Each line contains the SNP ID and the ancestral allele, separated by
whitespace or tab.

```
rs123 A
rs456 C
rs789 G
...
```

### Filter Input Files

For filtering options, `pygenstrat` accepts simple list files:

**Individual Filter Files** (`--keep-indv`, `--remove-indv`): Each line contains one individual ID.

```
sample1
sample2
...
```

**Population Filter Files** (`--keep-pop`, `--remove-pop`): Each line contains one population name.

```
Population1
Population2
...
```

**SNP Filter Files** (`--keep-snp`, `--remove-snp`): Each line contains one SNP ID.

```
rs123
rs456
...
```

**Region Filter Files** (`--keep-region`, `--remove-region`): Each line contains chromosome, start position, and end position, separated by whitespace.

```
1 1000000 2000000
2 5000000 6000000
...
```

You can also specify regions directly:

- As `chrID:start-end` (example: `--keep-region 1:1000000-2000000`).
- For multiple regions: `--keep-region "1:1000-2000,2:5000-6000"`

**Chromosome Filter Files** (`--keep-chr`, `--remove-chr`): Each line contains one chromosome ID.

```
1
2
23
...
```

You can also specify chromosomes directly:

- Use commas to separate: `--keep-chr 1,2,3`
- Use ranges: `--keep-chr 1-22`
- Combine: `--keep-chr 1,2,5-7,22`

## Usage Examples

### Format Conversion

```bash
#Convert ANC to EIG format
pygenstrat --prefix /path/to/input/basename --out-type eig --out /path/to/output/basename

#Convert EIG to ANC format
pygenstrat --prefix /path/to/input/basename --out-type anc --out /path/to/output/basename
```

### Filtering

```bash
#Filter by  MAF
pygenstrat --prefix /path/to/input/basename --min-maf 0.01 --out /path/to/output/basename

#Filter missing rate
pygenstrat --prefix /path/to/input/basename --keep-indv individuals.txt --geno 0.05 --out /path/to/output/basename

#Complex filtering with output format specification
pygenstrat --prefix /path/to/input/basename \
        --keep-pop europe.txt \
        --remove-snps exclude_snps.txt \
        --min-maf 0.05 \
        --geno 0.90 \
        --verbose \
        --out-type eig \
        --out /path/to/output/basename
```

### Statistical Analysis

```bash
#Calculate allele frequencies
pygenstrat --prefix /path/to/input/basename --freq --out /path/to/output/basename

#Calculate missing data rates per individual
pygenstrat --prefix /path/to/input/basename --missing --out /path/to/output/basename

#Calculate per-individual missing rate per chromosome (writes one file per chromosome)
pygenstrat --prefix /path/to/input/basename --missing-by-chr --out /path/to/output/basename

#Both statistics with filtering 
pygenstrat --prefix /path/to/input/basename --freq --missing --min-maf 0.01 --out /path/to/output/basename
```

### Updating Files

```bash
#Update population information for individuals
pygenstrat --prefix /path/to/input/basename --update-ind new_populations.txt

#Update SNP IDs
pygenstrat --prefix /path/to/input/basename --update-snp new_snp_ids.txt

#Update genetic distances using a genetic distance map 
pygenstrat --prefix /path/to/input/basename --genetic-distance genetic_map_GRCh37.txt

#Update genetic distances and output in Morgans M instead of centiMorgans cM
pygenstrat --prefix /path/to/input/basename --genetic-distance genetic_map_GRCh37.txt --map-unit M

#Set all genetic distances in output to zero
pygenstrat --prefix /path/to/input/basename --genetic-distance zero 

#Flip strand 
pygenstrat --prefix /path/to/input/basename --flip-strand flipstrands.txt 
```

### Genotype Management / Transformation

```bash
# Sample-based polarisation: Use a sample ID from your dataset as the polarisation reference
pygenstrat --prefix /path/to/input/basename --polarise SampleID --out /path/to/output_polarised

# File-based polarisation: Use a two-column file mapping SNP IDs to ancestral alleles
pygenstrat --prefix /path/to/input/basename --polarise ancestral_alleles.txt --out /path/to/output_polarised

# Randomly convert heterozygous genotypes to homozygous
pygenstrat --prefix /path/to/input/basename --random-haploidise --out /path/to/output/basename
```

## Output Files

| File Extension  | Description                                                                                              |
|-----------------|----------------------------------------------------------------------------------------------------------|
| `.geno`         | Processed genotype data                                                                                  |
| `.snp`          | Processed SNP information                                                                                |
| `.ind`          | Processed individual information                                                                         |
| `.missing`      | Per-individual missing data rates (if --missing used)                                                    |
| `.chr*.missing` | Per-individual missing data rates **per chromosome** (if --missing-by-chr used, one file per chromosome) |
| `.frq`          | Allele frequency statistics (if --freq used)                                                             |
| `.removed.snp`  | List of filtered SNPs with reasons (if --verbose used)                                                   |
| `.backup`       | Backup of original file (when using update options)                                                      |

