
# This loop iterates over the columns of the vcf file (=individuals) 

for i in {10..300} # interval up to the total amount of individuals in vcf
do
all=$(awk -F "\t" "{print $1}" data2_md80_shortnames_autosomal_293_chr8.vcf | wc -l) # num of genotyped positions for the whole vcf file (a constant number for all indvs)
col=$(awk -F "\t" "{print \$$i}" data2_md80_shortnames_autosomal_293_chr8.vcf) # genotypes of sample $i
missing=$(echo $col | grep -o ":0,0,0,0" | wc -l) # missing genotypes of sample $i
genotyped=$(($all - $missing - 1)) # the "1" is the sample name in variable $all, discount
heterozygots1=$(echo $col | grep -o "1/0" | wc -l) # num of heterozygous genotypes of sample $i
heterozygots2=$(echo $col | grep -o "0/1" | wc -l) # num of heterozygous genotypes of sample $i
heterototal=$(($heterozygots1 + $heterozygots2)) # sum the combinations ancestral/derived and derived/ancestral
python -c "print($heterototal / $genotyped)" # proportion of heterozygote positions for sample $i
done