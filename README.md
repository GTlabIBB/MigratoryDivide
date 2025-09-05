# Migratory divide

Code related to García-Berro et al., 2025 "A north-south hemispheric migratory divide in butterflies". Deposited for reference and to allow results to be reproduced independently. Brief description of the pipelines of the study:

- **Figure1_heatmap_observations.R**: Makes distribution plots of latitudinally sliced observational data from public repositories
- **Figure2_Population_genetics_with_adegenet.R**: Visualizes quality of RAD loci assembled using ipyrad, performs filtering steps and population statistic analysis (PCA, allele contributions)
- **Figure3_Scans_visualization_and_tests.R**: calculates and visualizes genome-wide patterns of genetic diversity (π), genetic differentiation (FST), and divergence (DXY) from VCF-derived summary statistics obtained with Pixy. Includes simple outlier detection and statistical tests to explore population differences.
- **Figure3_heterozygosity.sh**: Calculates per-individual genomic heterozygosity from a multi-sample VCF file (RAD-seq loci data obtained with ipyrad). For each individual, it computes the proportion of heterozygous genotypes among all non-missing genotypes
- **Figure3_heteroygosity_tests_and_plots.R**: Tests and visualization of heterozygosity differences between populations
- **mapping.sh**: From WGS reads to bam file: performs mapping, mark duplicates and indel realignment (for phylogenetic inference and breakpoint identification) 
- **goterms_topGO.R**: Tests for enrichment of GO terms among candidate genes (e.g. FST outliers)
- **Breakpoint_identification_Breakdancer.sh**: Run Breakdancer to identify breakpoints based on a read-paired mapping information 
