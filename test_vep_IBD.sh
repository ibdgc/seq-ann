#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate seq-ann
cd anno/vep_data/input

###### convert into vcf format ######
cut -f 1 test_LauraMetaCDsummaryStats5e7 > VarIDs.txt
sed -i '1d' VarIDs.txt
awk -F':' '{print $1"\t"$2"\t"$0"\t"$3"\t"$4"\t""30""\t""PASS""\t""AC=1"}' VarIDs.txt > IBDvcf_without_header.txt
echo -e "##fileformat=VCFv4.1\n##INFO=<ID=AC,Number=.,Type=String,Description="ac">\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=242193529>\n##contig=<ID=chr3,length=198295559>\n##contig=<ID=chr4,length=190214555>\n##contig=<ID=chr5,length=181538259>\n##contig=<ID=chr6,length=170805979>\n##contig=<ID=chr7,length=159345973>\n##contig=<ID=chr8,length=145138636>\n##contig=<ID=chr9,length=138394717>\n##contig=<ID=chr10,length=133797422>\n##contig=<ID=chr11,length=135086622>\n##contig=<ID=chr12,length=133275309>\n##contig=<ID=chr13,length=114364328>\n##contig=<ID=chr14,length=107043718>\n##contig=<ID=chr15,length=101991189>\n##contig=<ID=chr16,length=90338345>\n##contig=<ID=chr17,length=83257441>\n##contig=<ID=chr18,length=80373285>\n##contig=<ID=chr19,length=58617616>\n##contig=<ID=chr20,length=64444167>\n##contig=<ID=chr21,length=46709983>\n##contig=<ID=chr22,length=50818468>\n##contig=<ID=chrX,length=156040895>\n##contig=<ID=chrY,length=57227415>\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > vcfheader.txt
cat vcfheader.txt IBDvcf_without_header.txt > IBD.vcf

rm VarIDs.txt IBDvcf_without_header.txt vcfheader.txt

###### sort vcf ######
picard SortVcf I=IBD.vcf O=sortedIBD.vcf

###### VEP annotation ######
ml singularity/3.6.4
singularity exec -B ./vep_data:/opt/vep/.vep \
  vep_109.sif vep \
  --dir /opt/vep/.vep \
  --cache \
  --offline \
  --format vcf \
  --vcf \
  --force_overwrite \
  --input_file /opt/vep/.vep/input/sortedIBD.vcf \
  --output_file /opt/vep/.vep/output/annotated_IBD.vcf \
  --check_existing \
  --ccds \
  --hgvs \
  --hgvsg \
  --symbol \
  --numbers \
  --regulatory \
  --canonical \
  --protein \
  --biotype \
  --nearest symbol \
  --gene_phenotype \
  --pubmed \
  --uniprot \
  --tsl \
  --appris \
  --mane \
  --mane_select \
  --af \
  --af_1kg \
  --af_gnomade \
  --af_gnomadg \
  --max_af \
  --variant_class \
  --mirna \
  --pick_allele \
  --plugin TSSDistance \
  --plugin NearestExonJB \
  --plugin Conservation,/sc/arion/projects/Itan_lab/vep_data/Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw

###### Extract annotations ######
cd vep_data/output
conda activate seq-ann
INFILE=annotated_IBD.vcf
OUTFILE=vep_IBD.tsv
OUTHEADER=header_vep_IBD.tsv

sed -i "/^#/s/[()]//g" $INFILE
sed -i "/^#/s/-/_/g" $INFILE
sed -i "/^#/s/++/plus_plus/g" $INFILE

bgzip -@ 8 -c $INFILE > ${INFILE}.gz
tabix -p vcf ${INFILE}.gz

bcftools +split-vep ${INFILE}.gz \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\n' \
  -d \
  -A tab \
  -o $OUTFILE
bcftools +split-vep ${INFILE}.gz -l > $OUTHEADER

cut -f 2 $OUTHEADER > tmp1
tr '\n' '\t' < tmp1 > tmp2
sed 's/\t$//' tmp2 | awk '{print "CHROM""\t""POS""\t""ID""\t""REF""\t""ALT""\t"$0}' > header.txt
cat header.txt $OUTFILE > vep_annotated_IBD.txt
rm tmp1 tmp2 header.txt

conda deactivate
