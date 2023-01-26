# Population-scale SNP and indel calling pipeline from short read data

## Pipeline overview

1. `01_call_SNPs_indels.sh` : Call both SNPs and indels in all samples simultaneously, but in one chromosome at the time.
2. `02_concat.sh` : Concatenate VCFs for all chromosomes together, in a single large VCF file.
3. `03_filter.sh` : Filter raw SNPs and indels calls according to predefined criteria, and split SNPs and indels into separate output VCFs.

## Prerequisites

### Files
* A reference genome (`.fasta`) and its index (`.fai`) in `03_genome`

* Bam files for all samples and their indexes. These can be soft-linked/hard_linked in `04_bam` folder for easier handling : if `$BAM_PATH` is the remote path to bam files, use `for file in $(ls -1 $BAM_PATH/*); do ln -s $file ./04_bam; done`. These bam files should be named as `SAMPLEID.bam`, and can be produced using [wgs_sample_preparation pipeline](https://github.com/enormandeau/wgs_sample_preparation)

* A bam list in `02_infos`. This list can be generated with the following command, where `$BAM_DIR` is the path of the directory where bam files are located : `ls -1 $BAM_DIR/*.bam > 02_infos/bam_list.txt`

* A chromosomes list (or contigs, or sites) in `02_infos`. This list is used for parallelizing the SV calling step. It can be produced from the indexed genome file (`"$GENOME".fai`) : `less "$GENOME".fai | cut -f1 > 02_infos/chr.txt`. Warning : only add chromosomes or contigs for which variants are to be called (e.g. not unplaced contigs) 

### Software
* `bcftools` 1.12+ (done with 1.16)