# Migratory divide

Code related to García-Berro et al., 2025 "A north-south hemispheric migratory divide in butterflies". Deposited for reference and to allow results to be reproduced independently. Brief description of the pipelines of the study:

- Distribution plots of observational data from public repositories (Figure1_heatmap_observations.R)
- Visualising quality of RAD loci assembled using ipyrad, filtering and population statistic analysis (PCA, allele contributions) (Figure2_Population_genetics_with_adegenet.R)
- Testing heterozygosity differences according to populations and ploting against PC1 and plot differentiation genome scans obtained with Pixy (Figure3_heteroygosity_tests_and_plots.R)
- WGS from reads to bam file: mapping, mark duplicates and indel realignment (for phylogeny and breakpoint id) (mapping.sh)
- Enrichment tests of GO terms in candidate genes (goterms_topGO.R)
- Run Breakdancer to identify breakpoints based on a read-paired read mapping information (Breakpoint_identification_Breakdancer.sh)
