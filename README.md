# Migratory divide

Code related to García-Berro et al., 2025 "A north-south hemispheric migratory divide in butterflies". Deposited for reference and to allow results to be reproduced independently. Brief description of the pipelines of the study:

- **Figure1_heatmap_observations.R**: Make distribution plots of latitudinally sliced observational data from public repositories
- **Figure2_Population_genetics_with_adegenet.R**: Visualising quality of RAD loci assembled using ipyrad, filtering steps and population statistic analysis (PCA, allele contributions) 
- **Figure3_heterozygosity.sh**: Calculates per-individual genomic heterozygosity from a multi-sample VCF file (RAD-seq loci data obtained with ipyrad). For each individual, it computes the proportion of heterozygous genotypes among all non-missing genotypes
- **Figure3_heteroygosity_tests_and_plots.R**: Testing heterozygosity differences between populations and plotting against PC1 and plot differentiation genome scans obtained with Pixy 
- **mapping.sh**: WGS from reads to bam file: mapping, mark duplicates and indel realignment (for phylogeny and breakpoint id) 
- **goterms_topGO.R**: Enrichment tests of GO terms in candidate genes
- **Breakpoint_identification_Breakdancer.sh**: Run Breakdancer to identify breakpoints based on a read-paired read mapping information 
