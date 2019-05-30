# DeepShape
# For the estimation of Isoform-Level Ribosome Abundance and Distribution Accurately which needs Ribo-seq data only.

# Preprocess of reference

python filter_gencode_transcript.py gencode.v26.annotation.gtf gencode.v26.pc_transcripts.fa gencode.v26.pc_translations.fa
python select_rrna_from_ncrna.py Homo_sapiens.GRCh38.ncrna.fa GRCh38.rrna.fa
python build_contaminant.py GRCh38.rrna.fa hg38-tRNAs.fa Homo_sapiens human_contaminant.fa

python find_transcript_synonyms.py gencode.v26.pc_transcripts_filter.fa synonym.txt

mkdir StarIndex
mkdir StarIndex/contaminant
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir StarIndex/contaminant --genomeFastaFiles human_contaminant.fa --genomeSAindexNbases 8 --genomeChrBinNbits 11
mv Log.out StarIndex/contaminant
mkdir StarIndex/transcript
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir StarIndex/transcript/ --genomeFastaFiles gencode.v26.pc_transcripts_filter.fa --genomeSAindexNbases 11 --genomeChrBinNbits 12
mv Log.out StarIndex/transcript

# Mapping Ribo-seq data to references

gunzip -c Ref_RiboSeq_GSM546920.txt.gz > Ref_RiboSeq_GSM546920.fastq
srun -c3 cutadapt -m 5 -a TCGTATGCCGTCTTCTGCTTG --trim-n -o Ref_RiboSeq_GSM546920.adapterfree.fq Ref_RiboSeq_GSM546920.fastq
srun -c2 fastq_quality_filter -q 20 -p 50 -i Ref_RiboSeq_GSM546920.adapterfree.fq -o Ref_RiboSeq_GSM546920.HQ.fq

mkdir align
mkdir align/contaminant/
filter_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical"
STAR --runThreadN 15 --genomeDir StarIndex/contaminant/ --readFilesIn Ref_RiboSeq_GSM546920.HQ.fq --outFileNamePrefix align/contaminant/Ref_RiboSeq_GSM546920.HQ_ --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${filter_params}

mkdir align/transcript
filter_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical"
SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM"
STAR --runThreadN 15 --genomeDir StarIndex/transcript/ --readFilesIn align/contaminant/Ref_RiboSeq_GSM546920.HQ_Unmapped.out.mate1 --outFileNamePrefix align/transcript/Ref_RiboSeq_GSM546920.HQ_transcript_ ${SAM_params} ${filter_params}

samtools view align/transcript/Ref_RiboSeq_GSM546920.HQ_transcript_Aligned.out.bam | awk '{if(length($10)>=25 && length($10)<=36){print $1"\t"$3"\t"$4"\t"length($10);}}' | awk -F"|" '{printf $1"\t"$2; for(i=8;i<=NF-1;i++){if($i~/CDS:/) {gsub("CDS:","",$i); gsub("-","\t",$i); printf "\t"$i;}} printf $NF"\n"}' > Ref_RiboSeq_GSM546920.bam.reformat
python ParseBam2RiboPos.py gencode.v26.annotation.gtf gencode.v26.pc_transcripts_filter.fa Ref_RiboSeq_GSM546920.bam.reformat Ref_RiboSeq_GSM546920.bam.reformat.absolute

# DeepShape-prime

mkdir DeepShapePrimeOutputs
python DeepShape-prime.py Ref_RiboSeq_GSM546920.bam.reformat.absolute gencode.v26.pc_transcripts_filter.fa ./DeepShapePrimeOutputs 200

# DeepShape

mkdir DeepShapeOutputs_Loop199
awk '{print $1"\t"$2;}' DeepShapePrimeOutputs/runlog_199.txt | sed '$d' > TPM_DeepShapePrime_Loop199.txt
python TrueRibo_RecursiveShape_customloops.py \
 -r gencode.v26.pc_transcripts_filter.fa \
 -a WT_RPF_CT00_rep2.bam.reformat.absolute \
 -s TPM_DeepShapePrime_Loop199.txt \
 -o ./DeepShapeOutputs_Loop199 \
 -l 10
