#!/bin/bash
for fn in /external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_raw/*_R1.fastq.gz;
do
rawdir="/external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_raw"
f=`basename ${fn}`
samp="${f%_*}"
echo "Processing sample ${samp} using ${rawdir}/${samp}_R1.fastq.gz and ${rawdir}/${samp}_R2.fastq.gz"
salmon quant -i "/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/Transcriptome_refs/salmontool_mouse_index/" -l A \
         -1 ${rawdir}/${samp}_R1.fastq.gz \
         -2 ${rawdir}/${samp}_R2.fastq.gz \
         -p 8 --validateMappings -o /external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_raw/salmon_batch_results/${samp}_quant
done