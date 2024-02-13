# seq-ann
A vep annotation pipeline

## Dependencies
- Package management is handled by [conda](https://docs.conda.io/en/latest/miniconda.html).
- singularity (for VEP Docker image)
- Tools installed during first run
    - `bcftools`
    - `picard`
    - `vep`

## Setup

### First Time
```sh
# clone the repo
https://github.com/ibdgc/protocol-sharing.git
cd protocol-sharing/wgs-wes/mssm/seq-ann/setup/
./preamble
```

## Usage
- Execution of seq-ann is performed with a single shell command.

```sh
cd seq-ann
./vep_IBD.sh input.txt /path/to/plugin/dir
```
- input.txt contains one column for the variant identifier (based on GRCh38 positions) without a header
- Files required to use the VEP plugins should be downloaded separately according to the instructions at https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html. Each file should then be added to the directory with the plugin's corresponding name (e.g., /path/to/plugin/dir/dbNSFP/dbNSFP4.3a_grch38.gz).

```sh
# Example input.txt file
# This file should be placed in seq-ann/anno/vep_data/input before execution

chr12:106992760:C:T
chr12:106992822:T:TA
chr1:55039791:G:C
chr17:7661934:T:C
.
.
.

```

- There is no need to sort the variants. Chromosome-based files (chr1-22 and chrX) will be generated and processed in a parallel mode.

## Output
- The main routine will automatically generate an output directory (/seq-ann/anno/vep_data/output/).
- Two main output files will be generated: 

final_annotations_IBD.tsv
final_annotations_coding_IBD.tsv

final_annotations_IBD.tsv will contain standard VEP annotations along with annotations from `AncestralAllele`, `Conservation`, `DisGeNET`, `NearestExonJB`, `TSSDistance` and `UTRAnnotator` plugins.

final_annotations_coding_IBD.tsv will include an extensive set of annotations from multiple VEP plugins:

`Blosum62`
`CADD`
`Carol`
`Condel`
`Conservation`
`dbNSFP`
`dbscSNV`
`DisGeNET`
`Downstream`
`EVE`
`IntAct`
`LOEUF`
`LoFtool`
`MaxEntScan`
`mutfunc`
`neXtProt`
`NearestExonJB`
`NMD`
`pLI`
`PrimateAI`
`ReferenceQuality`
`REVEL`
`satMutMPRA`
`SpliceAI`
`SpliceRegion`
`TSSDistance`
`UTRAnnotator`

and the following annotations:

`Mutation Significance Cut-off` (PMID: 26820543)
`Gene Damage Index` (PMID: 26483451)
`GOF-LOF predictions` (https://doi.org/10.1101/2022.06.08.495288)
