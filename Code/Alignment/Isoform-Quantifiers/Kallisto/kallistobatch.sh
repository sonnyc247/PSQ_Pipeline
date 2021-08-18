#!/bin/bash
for fn in /external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_raw/*_R1.fastq.gz;
do
rawdir="/external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_raw"
f=`basename ${fn}`
samp="${f%_*}"
echo "Processing sample ${samp} using ${rawdir}/${samp}_R1.fastq.gz and ${rawdir}/${samp}_R2.fastq.gz"
kallisto quant -i "/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/Transcriptome_refs/kallisto_mouse_index/kallisto_mouse_transcripts.idx" \
         -o "/external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_raw/kallisto_batch_results/${samp}_quant" \
         -b 100 ${rawdir}/${samp}_R1.fastq.gz ${rawdir}/${samp}_R2.fastq.gz 
done