import pandas as pd
import sys
import os

datafile = sys.argv[1]
headerfile = sys.argv[2]
outfile = sys.argv[3]

gdi_df = pd.read_csv("./GDI.csv")
msc_df = pd.read_csv("./MSC_v1.6_95.txt", usecols=['Gene', 'MSC'], sep="\t")
goflof_df = pd.read_csv("/sc/arion/scratch/karsm02/vep_plugin_resource/GOF_LOF/GoFLoF.tsv", sep="\t")

# read all column names of the vcf output
csqheader = pd.read_csv(headerfile, header=None, names=[
    "i", "name"], sep="\t")["name"].to_list()
csqheaderdedup = []
for i in range(len(csqheader)):
    if csqheader[i] in ['APPRIS', 'TSL']:
        csqheaderdedup.append(csqheader[i]+str(i))
    else:
        csqheaderdedup.append(csqheader[i])
myheader = ["CHROM", "POS", "ID", "REF", "ALT"]
myheader += csqheaderdedup  # add the columns of the input file

# compares CADD scores to MSC thresholds


def compare(x, y):
    if x != "." and not pd.isnull(y):
        if float(x) > float(y):
            return "High"
        else:
            return "Low"


chunk = 10**6
for chunk in pd.read_csv(
    datafile,
    sep="\t",
    chunksize=chunk,
    header=None,
    names=myheader
):
    chunk_gdi = pd.merge(
        chunk,
        gdi_df,
        how='left',
        left_on='SYMBOL',
        right_on='Gene',
        suffixes=('', '_gdi')
    )
    chunk_gdi.drop(
        'Gene_gdi',
        axis=1,
        inplace=True
    )

    chunk_gdi_msc = pd.merge(
        chunk_gdi,
        msc_df,
        how='left',
        left_on='SYMBOL',
        right_on='Gene',
        suffixes=('', '_msc')
    )
    chunk_gdi_msc.drop(
        'Gene_msc',
        axis=1,
        inplace=True
    )
    chunk_gdi_msc.rename(
        columns={'MSC': 'MSC_95CI'},
        inplace=True
    )

    chunk_gdi_msc['MSC_CADD_95CI'] = chunk_gdi_msc.apply(
        lambda x: compare(x.CADD_PHRED, x.MSC_95CI), axis=1)
    
    chunk_gdi_msc_goflof = pd.merge(
        chunk_gdi_msc,
        goflof_df,
        how='left',
        left_on='ID',
        right_on='ID_goflof'
    )
    
    chunk_gdi_msc_goflof.drop(
        'ID_goflof',
        axis=1,
        inplace=True
    )

    if not os.path.isfile("%s.tsv" % outfile):
        chunk_gdi_msc_goflof.to_csv(
            "%s.tsv" % outfile,
            sep="\t",
            index=False,
            na_rep='.'
        )
    else:
        chunk_gdi_msc_goflof.to_csv(
            "%s.tsv" % outfile,
            sep="\t",
            index=False,
            mode='a',
            header=False,
            na_rep='.'
        )
