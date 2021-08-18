#!/bin/bash

##this code is adapted from https://github.com/derekhoward/he_seq, thanks to Derek Howard!
##this code ASSUMES PAIRED END READS
##this code executes RSEM quantification from reads previously aligned by STAR (STAR_pipe_for_RSEM.sh)
##this version of the script is setup to run with the *Aligned.toTranscriptome.out.bam files of interest being in a single directory
##use it as ./RSEM_pipe_post_STAR.sh /output_path/ /coord_bams_path/ /USE_THIS_RSEMDir/
##example of single-run command: rsem-calculate-expression --bam /nethome/kcni/dhoward/test1/Aligned.toTranscriptome.out.bam /nethome/kcni/dhoward/test1/processed_rsem/ /nethome/kcni/dhoward/test1/rsem_output/

output_path=$1  #desired output directory
bam_dir=$2 #directory which contains outputs from STAR (with each sample having a *Aligned.toTranscriptome.out.bam from STAR)
rsem_ref_path=$3 #path to preprocessed RSEM reference index


#initial setup/prep
cd $output_path
mkdir RSEM_scripts # will place each individual script for each sample in here
mkdir RSEM_results # will put each directory of RSEM outputs in here
mkdir RSEM_results/Gene_counts # will collect gene count tables here
mkdir RSEM_results/Isoform_counts # will collect isoform count tables here

scripts_dir="${output_path}RSEM_scripts/"
output_dir="${output_path}RSEM_results/" # path for output of RSEM processing

#generate processing script for each bam file in directory

echo "Processing files from $bam_dir"

cd $bam_dir

for bamfile in `find $bam_dir -maxdepth 1 -type f | sort`; do

    bam_fn_stem=$(basename "$bamfile" Aligned.toTranscriptome.out.bam) # extracts stem part of filename (removing extensions)
    echo "Processing $bam_fn_stem"
    echo module load RSEM >> "${scripts_dir}RSEMParamScript_${bam_fn_stem}.sh"
    echo mkdir ${output_dir}${bam_fn_stem} >> "${scripts_dir}RSEMParamScript_${bam_fn_stem}.sh"
    echo rsem-calculate-expression --paired-end --no-bam-output --bam $bamfile $rsem_ref_path "${output_dir}${bam_fn_stem}/${bam_fn_stem}" >> "${scripts_dir}RSEMParamScript_${bam_fn_stem}.sh"
    echo mv "${output_dir}${bam_fn_stem}/${bam_fn_stem}".genes.results "${output_dir}Gene_counts/" >> "${scripts_dir}RSEMParamScript_${bam_fn_stem}.sh"
    echo mv "${output_dir}${bam_fn_stem}/${bam_fn_stem}".isoforms.results "${output_dir}Isoform_counts/" >> "${scripts_dir}RSEMParamScript_${bam_fn_stem}.sh" 
    chmod +x ${scripts_dir}RSEMParamScript_${bam_fn_stem}.sh
    
done

#generate list of all parallel commands to be run (i.e: parallel running of all scripts previously generated); DO NOT NEED "parallel" at beginning if you're not doing multiple scripts per node

cd $scripts_dir

for file in $(ls *.sh); do
  echo "$scripts_dir$file" >> RSEMParamCom.txt
done

mv RSEMParamCom.txt $output_path
