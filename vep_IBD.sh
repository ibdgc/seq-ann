#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate seq-ann
cd anno/vep_data/input

##### commands for parallel processing of chromosome-based files #####
N=23
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( $@; )
    printf '%.3d' $? >&3
    )&
}

###### convert summary stats into vcf format ######
cut -f 1 input_file > VarIDs.txt
echo -e "##fileformat=VCFv4.1\n##INFO=<ID=AC,Number=.,Type=String,Description="ac">\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=242193529>\n##contig=<ID=chr3,length=198295559>\n##contig=<ID=chr4,length=190214555>\n##contig=<ID=chr5,length=181538259>\n##contig=<ID=chr6,length=170805979>\n##contig=<ID=chr7,length=159345973>\n##contig=<ID=chr8,length=145138636>\n##contig=<ID=chr9,length=138394717>\n##contig=<ID=chr10,length=133797422>\n##contig=<ID=chr11,length=135086622>\n##contig=<ID=chr12,length=133275309>\n##contig=<ID=chr13,length=114364328>\n##contig=<ID=chr14,length=107043718>\n##contig=<ID=chr15,length=101991189>\n##contig=<ID=chr16,length=90338345>\n##contig=<ID=chr17,length=83257441>\n##contig=<ID=chr18,length=80373285>\n##contig=<ID=chr19,length=58617616>\n##contig=<ID=chr20,length=64444167>\n##contig=<ID=chr21,length=46709983>\n##contig=<ID=chr22,length=50818468>\n##contig=<ID=chrX,length=156040895>\n##contig=<ID=chrY,length=57227415>\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > vcfheader.txt
all_chroms=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
for i in ${all_chroms[@]}
do
(
awk -F':' -v a=$i '$1==a' VarIDs.txt > $i\_VarIDs.txt
awk -F':' '{print $1"\t"$2"\t"$0"\t"$3"\t"$4"\t""30""\t""PASS""\t""AC=1"}' $i\_VarIDs.txt > $i\_IBDvcf_without_header.txt
cat vcfheader.txt $i\_IBDvcf_without_header.txt > $i\_IBD.vcf
rm $i\_VarIDs.txt $i\_IBDvcf_without_header.txt
)
done
rm vcfheader.txt

###### sort vcf ######
sort_vcfs(){
        picard SortVcf I=$1\_IBD.vcf O=sorted$1\_IBD.vcf
}

open_sem $N
for i in ${all_chroms[@]}
do
        run_with_lock sort_vcfs $i      
done
wait

rm chr*
conda deactivate
cd ../../

###### VEP annotation ######
plugdir=/path/to/plugin/data

vep_anno(){
  singularity exec -B ./vep_data:/opt/vep/.vep \
    vep_109.sif vep \
    --dir /opt/vep/.vep \
    --cache \
    --offline \
    --format vcf \
    --vcf \
    --force_overwrite \
    --input_file /opt/vep/.vep/input/sorted$1\_IBD.vcf \
    --output_file /opt/vep/.vep/output/annotated$1\_IBD.vcf \
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
    --fasta $plugdir/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    --plugin AncestralAllele,$plugdir/AncestralAllele/homo_sapiens_ancestor_GRCh38.fa.gz \
    --plugin Conservation,$plugdir/Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
    --plugin DisGeNET,file=$plugdir/DisGeNET/all_variant_disease_pmid_associations_final.tsv.gz \
    --plugin NearestExonJB \
    --plugin TSSDistance \
    --plugin UTRAnnotator,file=$plugdir/UTRAnnotator/uORF_5UTR_GRCh38_PUBLIC.txt
}

open_sem $N
for i in ${all_chroms[@]}
do
        run_with_lock vep_anno $i
done
wait

###### Extract annotations ######
cd vep_data/output
conda activate seq-ann
extract_anno(){
  INFILE=annotated$1\_IBD.vcf
  OUTFILE=$1\_vep_IBD.tsv
  OUTHEADER=$1\_header_vep_IBD.tsv

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
  cat header.txt $OUTFILE > vep_annotated$1\_IBD.txt
  rm tmp1 tmp2 header.txt
}

open_sem $N
for i in ${all_chroms[@]}
do
        run_with_lock extract_anno  $i
done

##### extensive annotation of coding regions #####

### select variants from coding region using IMPACT ###

for i in ${all_chroms[@]}
do
awk '!/MODIFIER/' vep_annotated$i\_IBD.tsv | cut -f 3 > $i\_coding_variant_IDs.txt
bcftools view -i ID=@$i\_coding_variant_IDs.txt ../input/sorted$i\_IBD.vcf > ../input/coding_sorted$i\_IBD.vcf
done

rm chr*

conda deactivate
cd ../../

##### vep annotation #####
extensive_vep_anno(){
        singularity exec -B ./vep_data:/opt/vep/.vep \
        vep_109.sif vep \
        --dir /opt/vep/.vep \
        --cache \
        --offline \
        --format vcf \
        --vcf \
        --force_overwrite \
        --fasta $plugdir/fasta/Homo_sapiens
.GRCh38.dna.toplevel.fa.gz \
        --input_file /opt/vep/.vep/input/coding_sorted$1\_IBD.vcf \
        --output_file /opt/vep/.vep/output/coding_annotated$1\_IBD.vcf \
        --assembly GRCh38 \
        --everything \
        --check_existing \
        --hgvsg \
        --nearest symbol \
        --mane_select \
        --pick_allele \
        --plugin AncestralAllele,$plugdir/AncestralAllele/homo_sapiens_ancestor_
GRCh38.fa.gz \
        --plugin Blosum62 \
        --plugin CADD,$plugdir/CADD/whole_genome_SNVs.tsv.gz,$plugdir/CADD/gnomad.genomes.r3.0.indel.tsv.gz \
        --plugin Carol \
        --plugin Condel,$plugdir/Condel/config,b \
        --plugin Conservation,$plugdir/Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
        --plugin dbNSFP,$plugdir/dbNSFP/dbNSFP4.3a_grch38.gz,ALL \
        --plugin dbscSNV,$plugdir/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
        --plugin DisGeNET,file=$plugdir/DisGeNET/all_variant_disease_pmid_associations_final.tsv.gz \
        --plugin Downstream \
        --plugin EVE,file=$plugdir/EVE/eve_merged.vcf.gz \
        --plugin IntAct,mutation_file=$plugdir/IntAct/mutations.tsv,mapping_file=$plugdir/IntAct/mutation_gc_map.txt.gz \
        --plugin LOEUF,file=$plugdir/LOEUF/loeuf_dataset_grch38.tsv.gz,match_by=gene \
        --plugin LoFtool,$plugdir/LoFtool/LoFtool_scores.txt \
        --plugin MaxEntScan,$plugdir/MaxEntScan/ \
        --plugin mutfunc,db=$plugdir/mutfunc/mutfunc_data.db \
        --plugin neXtProt \
        --plugin NearestExonJB \
        --plugin NMD \
        --plugin pLI,$plugdir/pLI/pLI_values.txt \
        --plugin PrimateAI,$plugdir/PrimateAI/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
        --plugin ReferenceQuality,$plugdir/ReferenceQuality/sorted_GRCh38_quality_mergedfile.gff3.gz \
        --plugin REVEL,$plugdir/REVEL/new_tabbed_revel_grch38.tsv.gz \
        --plugin satMutMPRA,file=$plugdir/satMutMPRA/satMutMPRA_GRCh38_ALL.gz \
        --plugin SpliceAI,snv=$plugdir/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$plugdir/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin SpliceRegion \
        --plugin TSSDistance \
        --plugin UTRAnnotator,file=$plugdir/UTRAnnotator/uORF_5UTR_GRCh38_PUBLIC.txt
}

open_sem $N
for i in ${all_chroms[@]}
do
        run_with_lock extensive_vep_anno $i
done
wait

##### extract annotations #####
cd vep_data/output
conda activate seq-ann

extract_coding_anno(){

        INFILE=coding_annotated$1\_IBD.vcf
        OUTFILE=coding_$1\_vep_IBD.tsv
        OUTHEADER=coding_$1\_header_vep_IBD.tsv

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
}

open_sem $N
for i in ${all_chroms[@]}
do
        run_with_lock extract_coding_anno  $i
done

###### MSC, GDI and GOF-LOF annotations ######
msc_gdi_goflof_anno(){
python ../../../MSC_GDI_GoFLoF.py \
  coding_$1\_vep_IBD.tsv \
  coding_$1\_header_vep_IBD.tsv \
  msc_gdi_goflof_coding_$1\_vep_IBD
}

open_sem $N
for i in ${all_chroms[@]}
do
        run_with_lock msc_gdi_goflof_anno $i
done
wait

##### prepare final files #####
### file for all variants ####
head -1  vep_annotatedchr1_IBD.tsv > header.tsv
sed -i '1d' vep*
cat header.tsv vep* > final_annotations_IBD.tsv
bgzip final_annotations_IBD.tsv

### file for coding variants ###
head -1  msc_gdi_goflof_coding_chr1_vep_IBD.tsv > header_coding.tsv
sed -i '1d' msc_gdi_goflof*
cat header_coding.tsv msc_gdi_goflof* > final_annotations_coding_IBD.tsv
bgzip final_annotations_coding_IBD.tsv

conda deactivate
