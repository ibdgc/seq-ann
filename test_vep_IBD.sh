#!/bin/bash

###### convert into vcf format ######
cut -f 1 test_LauraMetaCDsummaryStats5e7 > VarIDs.txt
sed -i '1d' VarIDs.txt
awk -F':' '{print $1"\t"$2"\t"$0"\t"$3"\t"$4"\t""30""\t""PASS""\t""AC=1"}' VarIDs.txt > IBDvcf_without_header.txt
echo -e "##fileformat=VCFv4.1\n##INFO=<ID=AC,Number=.,Type=String,Description="ac">\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=242193529>\n##contig=<ID=chr3,length=198295559>\n##contig=<ID=chr4,length=190214555>\n##contig=<ID=chr5,length=181538259>\n##contig=<ID=chr6,length=170805979>\n##contig=<ID=chr7,length=159345973>\n##contig=<ID=chr8,length=145138636>\n##contig=<ID=chr9,length=138394717>\n##contig=<ID=chr10,length=133797422>\n##contig=<ID=chr11,length=135086622>\n##contig=<ID=chr12,length=133275309>\n##contig=<ID=chr13,length=114364328>\n##contig=<ID=chr14,length=107043718>\n##contig=<ID=chr15,length=101991189>\n##contig=<ID=chr16,length=90338345>\n##contig=<ID=chr17,length=83257441>\n##contig=<ID=chr18,length=80373285>\n##contig=<ID=chr19,length=58617616>\n##contig=<ID=chr20,length=64444167>\n##contig=<ID=chr21,length=46709983>\n##contig=<ID=chr22,length=50818468>\n##contig=<ID=chrX,length=156040895>\n##contig=<ID=chrY,length=57227415>\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > vcfheader.txt
cat vcfheader.txt IBDvcf_without_header.txt > IBD.vcf

rm VarIDs.txt IBDvcf_without_header.txt

###### sort vcf ######
picard SortVcf I=IBD.vcf O=sortedIBD.vcf

###### VEP annotation ######
vep \
  -i sortedIBD.vcf \
  -o annotated_IBD.vcf \
  --verbose \
  --force_overwrite \
  --format vcf \
  --vcf \
  --assembly GRCh38 \
  --cache \
  --offline \
  --sift b \
  --polyphen b \
  --ccds \
  --hgvs \
  --symbol \
  --numbers \
  --domains \
  --regulatory \
  --canonical \
  --protein \
  --biotype \
  --af \
  --af_1kg \
  --af_gnomade \
  --af_gnomadg \
  --max_af \
  --variant_class \
  --mirna \
  --pick_allele
  --plugin --plugin CADD,/path/to/plugin/data/CADD/GRCh38/whole_genome_SNVs.tsv.gz,/path/to/plugin/data/CADD/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz
  #CADD plugin is required for the MSC annotation

###### Extract annotations ######
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

###### MSC and GDI annotations ######
python MSC_GDI.py \
  vep_IBD.tsv \
  header_vep_IBD.tsv \
  msc_gdi_vep_IBD
