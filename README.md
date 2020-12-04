# NemaScan

GWA Mapping and Simulation with _C. elegans, C. tropicalis, and C. briggsae_

## Required Software Packages (In User's PATH)

1. [R-v3.4.1](https://www.r-project.org/)
1. [nextflow-v19.07.0](https://www.nextflow.io/docs/latest/getstarted.html)
1. [BCFtools-v1.9](https://samtools.github.io/bcftools/bcftools.html)
1. [plink-v1.9](https://www.cog-genomics.org/plink2)
1. [R-cegwas2](https://github.com/AndersenLab/cegwas2)
1. [R-tidyverse-v1.2.1](https://www.tidyverse.org/)
1. [R-correlateR](https://github.com/AEBilgrau/correlateR)
1. [R-rrBLUP-v4.6](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf)
1. [R-RSpectra-v0.13-1](https://github.com/yixuan/RSpectra)
1. [R-ggbeeswarm-v0.6](https://github.com/eclarke/ggbeeswarm)
1. [R-ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)

## Downloading NemaScan
```
git clone https://github.com/AndersenLab/NemaScan.git
```

### Mappings Profile
```
nextflow main.nf -profile mappings --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz --traitfile input_data/elegans/phenotypes/PC1.tsv --sthresh [BF] [EIGEN] [#]
```
#### Required Mapping Parameters

* `--vcf` - a VCF file with variant data. There should also be a tabix-generated index file (.tbi) in the same folder as the specified VCF file that has the same name as the VCF except for the addition of the `.tbi` extension. (generated using `tabix -p vcf vcfname.vcf.gz`). If this flag is not used, a VCF for the _C. elegans_ species will be downloaded from [CeNDR](https://elegansvariation.org/data/release/latest).

* `--traitfile` - is a tab-delimited formatted (.tsv) file that contains trait information.  Each phenotype file should be in the following format (replace trait_name with the phenotype of interest):

| strain | trait_name_1 | trait_name_2 |
| --- | --- | --- |
| JU258 | 32.73 | 19.34 |
| ECA640 | 34.065378 | 12.32 |
| ... | ... | ... | 124.33 |
| ECA250 | 34.096 | 23.1 |

#### Optional Mapping Parameters

* `--sthresh` - This determines the signficance threshold required for performing post-mapping analysis of a QTL. `BF` corresponds to Bonferroni correction, `EIGEN` corresponds to correcting for the number of independent markers in your data set, and `user-specified` corresponds to a user-defined threshold, where you replace user-specified with a number. For example `--sthresh=4` will set the threshold to a `-log10(p)` value of 4. We recommend using the strict `BF` correction as a first pass to see what the resulting data looks like.

* `--out` - A user-specified output directory name.

* `--group_qtl` - QTL within this distance of each other (bp) will be grouped as a single QTL by `Find_GCTA_Intervals_*.R`.

* `--ci_size` - The number of markers for which the detection interval will be extended past the last significant marker in the interval.


### Simulations Profile

```
nextflow main.nf -profile simulations --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz --simulate_nqtl --simulate_reps 2 --simulate_h2 input_data/all_species/simulate_h2.csv --simulate_eff input_data/all_species/simulate_effect_sizes.csv --simulate_strains input_data/all_species/simulate_strains.tsv --out example_simulation_output
module load R/3.6.3
Rscript bin/Assess_Simulated_Mappings.R example_simulation_output
```

#### Required Simulation Parameters

* `--vcf` - a VCF file with variant data. There should also be a tabix-generated index file (.tbi) in the same folder as the specified VCF file that has the same name as the VCF except for the addition of the `.tbi` extension. (generated using `tabix -p vcf vcfname.vcf.gz`). If this flag is not used, a VCF for the _C. elegans_ species will be downloaded from [CeNDR](https://elegansvariation.org/data/release/latest).

* `--simulate_maf` -  a single column CSV file that defines the minor allele frequency threshold used to filter the VCF prior to simulations (Default: 0.05).

* `--simulate_nqtl` - a single column CSV file that defines the number of QTL to simulate (format: one number per line, no column header)   (Default is provided: `input_data/all_species/simulate_nqtl.csv`).

* `--simulate_reps` - The number of replicates to simulate per number of QTL and heritability (Default: 2).

* `--simulate_h2` - A CSV file with phenotype heritability. (format: one value per line, no column header) (Default is located: input_data/all_species/simulate_h2.csv).

* `--simulate_eff` - A CSV file specifying a range of causal QTL effects. QTL effects will be drawn from a uniform distribution bound by these two values. If the user wants to specify _Gamma_ distributed effects, the value in this file can be simply specified as "gamma". (format: one value per line, no column header) (Default is located: input_data/all_species/simulate_effect_sizes.csv).

* `--simulate_strains` - A TSV file specifying the population in which to simulate GWA mappings. Multiple populations can be simulated at once, but causal QTL will be drawn independently for each population as a result of minor allele frequency and LD pruning prior to mapping. (format: one line per population; supplied population name and a comma-separated list of each strain in the population) (Default is located: input_data/all_species/simulate_strains.tsv).

#### Optional Simulation Parameters

* `--simulate_qtlloc` - A .bed file specifying genomic regions from which causal QTL are to be drawn after MAF filtering and LD pruning. (format: CHROM START END for each genomic region, with no header. NOTE: CHROM is specified as NUMERIC, not roman numerals as is convention in _C. elegans_)(Default is located: input_data/all_species/simulate_locations.bed).

* `--group_qtl` - QTL within this distance of each other (bp) will be grouped as a single QTL by `Find_GCTA_Intervals_*.R`.

* `--ci_size` - The number of markers for which the detection interval will be extended past the last significant marker in the interval.

### Annotations Profile (in development)

`nextflow main.nf --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz -profile annotations --species briggsae --wb_build WS270`

* `--species` - specifies what species information to download from WormBase (options: elegans, briggsae, tropicalis).

* `--wb_build` - specifies what WormBase build to download annotation information from (format: WSXXX, where XXX is a number greater than 270 and less than 277).

### Parameters

* `nextflow main.nf --help` - will display the help message


### Essential R Scripts

* `Get_GenoMatrix_Eigen.R` - Takes a genotype matrix and chromosome name as input and identifies the number significant eigenvalues.
* `Fix_Isotype_names_bulk.R` - Take sample names present in phenotype data and changes them to isotype names found on [CeNDR](elegansvariation.org) when the `--traitfile` flag is used.
* `create_causal_QTLs.R` - (_Simulations Profile Only_) Simulates a number of QTL effects specified by `simulate_nqtl.csv` drawn from either 1) a user-specified range or 2) a _Gamma_ distribution with shape = 0.4 and scale = 1.66. Effects are randomly assigned positive or negative directional effects. This script will also source `simulate_locations.bed` if provided. Otherwise, all variants present after MAF filtering and LD pruning are eligible to be selected as causal.
* `Find_GCTA_Intervals_Maps*.R` - (_Mapping Profile Only_) Assigns QTL "detection" intervals using the `--group_qtl` and `--ci_size` parameters if the number of significant markers does not exceed 15% of whole marker set (as this is a strong indication of high phenotypic noise due to genomic structure or low sampling). *NOTE* These regions should _not_ be treated statistical confidence intervals.
* `Find_GCTA_Intervals_*.R` - (_Simulation Profile Only_) Assigns QTL "detection" intervals for each simulated mapping using the `--group_qtl` and `--ci_size` parameters if the number of significant markers does not exceed 15% of whole marker set (as this is a strong indication of high phenotypic noise due to genomic structure or low sampling). *NOTE* These regions should _not_ be treated statistical confidence intervals.
* `Assess_Simulated_Mappings.R` - (_independent from nextflow pipeline_) Analyzes simulated mappings and gathers results and metadata for 1) all detected QTL and 2) undetected simulated QTL. Creates an .RData file within the Simulations directory that can be analyzed locally to determine performance for parameters of interest to the user.
* `pipeline.plotting.R` - (_Mapping Profile Only_) Generates 1) manhattan plots for each trait in the traitfile, 2) LD heatmaps for 

### Input Data Folder Structure

```
all_species
  ├── rename_chromosomes
  ├── simulate_effect_sizes.csv
  ├── simulate_h2.csv
  ├── simulate_maf.csv
  ├── simulate_nqtl.csv
  ├── simulate_strains.tsv
  ├── simulate_locations.bed
briggsae
  ├── annotations
      ├── GTF file
      ├── refFlat file
elegans
  ├── genotypes  
      ├── vcf
      ├── vcf_index
  ├── phenotypes
      ├── PC1.tsv
  ├── annotations
      ├── GTF file
      ├── refFlat file
tropicalis
  ├── annotations
      ├── GTF file
      ├── refFlat file
```

### Mapping Output Folder Structure

```
Genotype_Matrix
  ├── Genotype_Matrix.tsv
  ├── total_independent_tests.txt
 Mapping
  ├── Raw
      ├── traitname_lmm-exact_inbred.fastGWA
      ├── traitname_lmm-exact.loco.mlma
  ├── Processed
      ├── traitname_LMM_EXACT_INBRED_mapping.tsv
      ├── traitname_LMM_EXACT_LOCO_mapping.tsv
  ├── QTL_Regions
      ├── traitname_LMM_EXACT_INBRED_qtl_region.tsv
      ├── traitname_LMM_EXACT_LOCO_qtl_region.tsv
Plots
  ├── ManhattanPlots
      ├── traitname_manhattan.plot.png
  ├── LDPlots
      ├── traitname_LD.plot.png (if > 1 QTL detected)
  ├── EffectPlots
      ├── traitname_[QTL.INFO]_LOCO_effect.plot.png (if detected)
      ├── traitname_[QTL.INFO]_INBRED_effect.plot.png (if detected)
```

#### Genotype_Matrix
* `Genotype_Matrix.tsv` - pruned LD-pruned genotype matrix used for GWAS and construction of kinship matrix
* `total_independent_tests.txt` - number of independent tests determined through spectral decomposition of the genotype matrix

#### Mapping

##### Raw
* `traitname_lmm-exact_inbred.fastGWA` - Raw mapping results from GCTA's fastGWA program using an inbred kinship matrix
* `traitname_lmm-exact.loco.mlma` - Raw mapping results from GCTA's mlma program using a kinship matrix constructed from all chromosomes except for the chromosome containing each tested variant.
##### Processed
* `traitname_LMM_EXACT_INBRED_mapping.tsv` - Processed mapping results from lmm-exact_inbred raw mappings. Contains additional information nested such as 1) rough intervals (see parameters for calculation) and estimates of the variance explained by the detected QTL 2) phenotype information and genotype status for each strain at the detected QTL.
* `traitname_lmm-exact.loco.mlma` - Processed mapping results from lmm-exact.loco.mlma raw mappings. Contains same additional information as above.
##### QTL_Regions
* `traitname_*_qtl_region.tsv` - Contains only QTL information for each mapping. If no QTL are detected, an empty data frame is written.

#### Plots
* `traitname_manhattan.plot.png` - Standard output for GWA; association of marker differences with phenotypic variation in the population.
* `traitname_LD.plot.png` - If more than 1 QTL are detected for a trait, a plot showing the linkage disequilibrium between each QTL is generated.
* `traitname_[QTL.INFO]_INBRED_effect.plot.png` - Phenotypes for each strain are plotted against their marker genotype at the peak marker for each QTL detected for a trait. The dot representing each strain is shaded according to the percentage of the chromosome containing the QTL that is characterized as a selective sweep region.

### Simulation Output Folder Structure

```
Genotype_Matrix
  ├── [strain_set]_[MAF]_Genotype_Matrix.tsv
  ├── [strain_set]_[MAF]_total_independent_tests.txt
Simulations
  ├── NemaScan_Performance.example_simulation_output.RData
  ├── [specified effect range (simulate_effect_sizes.csv)]
      ├── [specified number of simulated QTL (simulate_nqtl.csv)]
          ├── Mappings
              ├── [nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_processed_LMM_EXACT_INBRED_mapping.tsv
              ├── [nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_processed_LMM_EXACT_LOCO_mapping.tsv
              ├── [nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_lmm-exact_inbred.fastGWA
              ├── [nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_lmm-exact.loco.mlma
          ├── Phenotypes
              ├── [nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_sims.phen
              ├── [nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_sims.par
  ├── (if applicable) [NEXT specified effect range]
      ├── ...
  ├── (if applicable) [NEXT specified effect range]
      ├── ...
```

#### Genotype_Matrix folder
* `*Genotype_Matrix.tsv` - pruned LD-pruned genotype matrix used for GWAS and construction of kinship matrix. This will be appended with the chosen minor allele frequency cutoff and strain set, as they are generated separately for each strain set.
* `*total_independent_tests.txt` - number of independent tests determined through spectral decomposition of the genotype matrix. This will be also be appended with the chosen minor allele frequency cutoff and strain set, as they are generated separately for each strain set.

#### Simulations
* `NemaScan_Performance.*.RData` - RData file containing all simulated and detected QTL from each successful simulated mapping. Contains:
1. Simulated and Detected status for each QTL.
2. Minor allele frequency and simulated or estimated effect for each QTL.
3. Detection interval according to specified grouping size and CI extension.
4. Estimated variance explained for each detected QTL.
5. Simulation parameters and the algorithm used for that particular regime.

##### Mappings
* As with the mapping profile, raw and processed mappings for each simulation regime are nested within folders corresponding each specified effect range and number of simulated QTL. QTL region files are not provided in the simulation profile; this information along with other information related to mapping performance are iteratively gathered in the generation of the performance .RData file.

##### Phenotypes
* `[nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_sims.phen` - Simulated strain phenotypes for each simulation regime.
* `[nQTL]_[rep]_[h2]_[MAF]_[effect range]_[strain_set]_sims.par` - Simulated QTL effects for each simulation regime. NOTE: Simulation regimes with identical numbers of simulated QTL, replicate indices, and simulated heritabilities should have _identical_ simulated QTL and effects.
