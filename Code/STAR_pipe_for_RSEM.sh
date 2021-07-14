#!/bin/bash

##this code assumes data is given as STAR-acceptable files and are compressed... future to do: add argument/option for specifying data format and compression status
##it would be good if files have short names with unique identifiers, to improve prefix funcionality

datadir=$1 #directory of where all the data/fastqfiles are stored (absolute path)
outputdir=$2 #directory of where you want the output to be stored
pairend=$3 #are paired reads being given? please enter "true"/"false"; also this code assumes paired fastqfiles are arranged/ordered together (so files 1 and 2 go together, 3 and 4 go together, etc.)
genomePath=$4 # absolute path to processed genome index to be used for STAR alignment
annotationPath=$5 # abs path to genome annotation file (typically .gtf file)

echo "Processing files from $datadir ."

#initial setup/prep
cd $outputdir
mkdir STAR_scripts
mkdir STAR_results
res_fold="${outputdir}STAR_results/"
mkdir $res_fold/coord_bams
#mkdir $res_fold/pcrless_bams
mkdir $res_fold/star_logs
#mkdir $res_fold/star_logs/PCRLess
mkdir $res_fold/star_logs/Init_QCMetrics
mkdir $res_fold/star_logs/Init_Runmsgs


cd $datadir

if [ "$pairend" = false ] ; then
  echo 'Specification shows that reads are unpaired.'

  tempi=1

  # generating STAR scripts that will be run in parallel
  for file in `find $datadir -maxdepth 1 -type f | sort`; do
    # getting extension and extensionless file name
    filename=$(basename -- "$file") # gets filename withtout full path
    extension=".${filename##*.}" # extracts just extension
    filename="${filename%.*}" # extracts filename without extension
    while [[ "$filename" == *"."* ]]; do
      extension=".${filename##*.}$extension"
      filename="${filename%.*}"
    done
    echo "Script is being generated for $filename, which has the extension $extension"

    prefix=$filename

    #Writing the script to process each sample
    echo '#!/bin/bash' >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo module load STAR >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mkdir $res_fold$prefix >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo STAR --genomeDir $genomePath --runThreadN 6 --sjdbGTFfile $annotationPath --readFilesIn $file --outFileNamePrefix "$res_fold$prefix/$prefix" --quantMode TranscriptomeSAM --outSAMtype None --readFilesCommand zcat >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$prefix/$prefix"Aligned.toTranscriptome.out.bam "${res_fold}coord_bams/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$prefix/$prefix"Log.final.out "${res_fold}star_logs/Init_QCMetrics/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$prefix/$prefix"Log.out "${res_fold}star_logs/Init_Runmsgs/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"

    chmod +x $outputdir/STAR_scripts/STARParaScript$tempi.sh
    echo "Script has been generated with the name STARParaScript$tempi.sh, script output results will start with $prefix."
    ((tempi=tempi+1))
  done

elif [ "$pairend" = true ] ; then
  echo 'Specification shows that reads are paired.'

  tempi=0

  #initialize an array of odd-numbered files (i.e: the first, third, fifth, etc. file in the directory)
  forarr=()
  for file in `find $datadir -maxdepth 1 -type f | sort | awk 'NR % 2 == 1'`; do
    forarr+=("$file")
  done
  echo "The odd-numbered files are:"
  printf '%s\n' "${forarr[@]}"

  #initialize an array of even-numbered files (i.e: the second, fourth, sixth, etc. file in the directory)
  revarr=()
  for file in `find $datadir -maxdepth 1 -type f | sort | awk 'NR % 2 == 0'`; do
    revarr+=("$file")
  done
  echo "The even-numbered files are:"
  printf '%s\n' "${revarr[@]}"

  #generating STAR scripts that will be run in parallel
  for file in "${forarr[@]}"; do
    #getting extension and extensionless file name for 1st file
    filename=$(basename -- "$file")
    extension=".${filename##*.}"
    filename="${filename%.*}"
    while [[ "$filename" == *"."* ]]; do
      extension=".${filename##*.}$extension"
      filename="${filename%.*}"
    done

    file2=${revarr[$tempi]}
    #getting extension and extensionless file name for 2nd file
    filename2=$(basename -- "$file2")
    extension2=".${filename2##*.}"
    filename2="${filename2%.*}"
    while [[ "$filename2" == *"."* ]]; do
      extension2=".${filename2##*.}$extension2"
      filename2="${filename2%.*}"
    done

    prefix=$filename

    #Writing the script to process each sample
    echo '#!/bin/bash' >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo module load STAR >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mkdir $res_fold$prefix >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo STAR --genomeDir $genomePath --runThreadN 6 --sjdbGTFfile $annotationPath --readFilesIn $file $file2 --outFileNamePrefix "$res_fold$prefix/$prefix" --quantMode TranscriptomeSAM --outSAMtype None --readFilesCommand zcat >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$prefix/$prefix"Aligned.toTranscriptome.out.bam "${res_fold}coord_bams/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$prefix/$prefix"Log.final.out "${res_fold}star_logs/Init_QCMetrics/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$prefix/$prefix"Log.out "${res_fold}star_logs/Init_Runmsgs/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"

    chmod +x $outputdir/STAR_scripts/STARParaScript$tempi.sh
    echo "A script has been generated for $filename$extension and $filename2$extension2 . Script has the name STARParaScript$tempi.sh and its results will start with $prefix ."
    ((tempi=tempi+1))
  done


else

  echo "The second argument, 'pairend', was not entered as 'true' or 'false'."
  exit 1

fi

#generate list of all parallel commands to be run (i.e: parallel running of all scripts previously generated); DO NOT NEED "parallel" at beginning if you're not doing multiple scripts per node

cd $outputdir/STAR_scripts/

for file in $(ls *.sh); do
  echo "${outputdir}STAR_scripts/$file" >> STARParaCom.txt
done

mv STARParaCom.txt $outputdir

cd $outputdir

#qbatch -w 04:30:00 -b slurm -v -i -c 2 -j 2 --ppj 12 STARParaCom.txt
