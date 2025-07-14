#!/usr/bin/env bash
#
# GEMMA germline WGS alignment
#
# This script will first convert the BAM-files to (unaligned) fastq files (.fq.gz),
# Optionally runs simplify on them (reduced filesize) and then realigns them to hg38
# The script also adds headers required by the bam/sam file standard.
#
# Usage: ./run_alignment.sh batchdir *_1.fq.gz
#        ./run_alignment.sh
#
# Dir batchX should contain .bam, _1.fq.gz (or symbolic links) to all files to be
# processed. With separate batch directories you can separate the files to
# be processed to be run on multiple VMs in parallel.
#
# Argument file_filter_str should include all input files. If the filter
# contains ".bam", input files will be converted to fastq first (".fq.gz"). If the
# input filter variable contains ".fq.gz" the first step is to align them.
#
# Use fastq file naming convention '{SAMPLENAME}_{BARCODE}_{FLOWCELL}_{LANE}_1.fq.gz'
#
# When processing fastq files, all files are first aligned individually with
# their mate-pairs and then merged based on having identical {SAMPLENAME}
# to bam files with a ".merged" tag in their filename. This is done to preserve
# flowcell information.

# Note: Only for paired-end sequencing
# Note: Input bam files cannot have ".hg38.bam" suffix to be realigned
# Note: fastq files need to have ".fq.gz" suffix to be processsed
# Note: The script is currently suited to be run on a hpc3.56core VM,
#       with 56 cpus and >200GB or RAM. Modify parallelization params
#       to match the system used for the processing.


###################################################
# Adjustable params

# Dealignment parallellization params
export n_parallel_dealignments=2

# Alignment parallellization params
export bwa_mem_threads=12
export bwa_mem_sort_threads=4
export n_parallel_alignments=4

# Simplify FASTQ files? (reduces filesize ~2GB)
export do_simplify=0

# Header processing parallellization
export n_parallel_headermods=14

# Information added to bam headers
export sequencer_str="ILLUMINA"

# The parallel2 -L param can be used to connect multiple VMs to a single job
# Both the batch dir and linked dir need to be accessible to all VMs
# Leave this variable empty to not use job linking on shared network drives
export parallel_linked_base_dir=""

####################################################

#RED='\033[0;31m'
#GREEN='\033[0;32m'
#YELLOW='\033[1;33m'
#PURPLE='\033[1;35m'
#NC='\033[0m'

RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
YELLOW=$(tput setaf 3)
PURPLE=$(tput setaf 5)
NC=$(tput sgr0)

####################################################

# Environment checks
export parallel2_bin="/mnt/gemma/bin/scripts/parallel2_95.py"
export bwa_mem_bin="bwa mem"
export fasta_bin="/home/ubuntu/.cargo/bin/fasta"
export picard_jar="/mnt/gemma/apps/picard_2.25.4/picard.jar"
export ref_genome="/mnt/gemma/bin/resources/homo_sapiens/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export env_error=0

if [[ $(command -v $parallel2_bin) == "" ]] ; then echo "ERROR: 'parallel2' executable not found."; export env_error=1; fi
if [[ $(command -v $bwa_mem_bin) == "" ]] ; then echo "ERROR: 'bwa mem' binary executable not found."; export env_error=1; fi
if [[ $(command -v $fasta_bin) == "" ]] ; then echo "ERROR: 'fasta_bin' binary executable not found."; export env_error=1; fi
if [[ ! -e "$picard_jar" ]] ; then echo "ERROR: 'picard_jar' executable not found."; export env_error=1; fi
if [[ ! -e "$ref_genome" ]] ; then echo "ERROR: Reference genome (hg38) fasta file not found."; export env_error=1; fi
if [[ $env_error == 1 ]] ; then echo "ERROR: Required software not installed or errorneous paths given. Exitting."; exit 1; fi


# Prepare batch work
export batchdir=$1
export batch_filter=$2

# If no batch dir given, use current dir
[ -z "$batchdir" ] && export batchdir=$(pwd)
export batchdir=$(echo $batchdir | sed 's:/*$::')

echo "bd:$batchdir"

export input_filter=$batch_filter
export n_cpus=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
export ram=$(awk '/MemTotal/ {print $2}' /proc/meminfo)
#echo "RAM: $ram"
export processing_started=0
export n_ramgb=$(($ram / 1000000))
echo "INFO: CPUs: $n_cpus"
echo "INFO: RAM: ${n_ramgb}GB"

# Calculate required number of cpus for the alignment
export align_n_cpus=$((bwa_mem_threads*n_parallel_alignments))

if [[ $align_n_cpus > n_cpus ]] ; then
    echo "WARNING: Alignment parallelization requires more cpus than system has: ${align_n_cpus}/${n_cpus}."
fi

if [[ $batchdir == "" ]] ; then
  echo "ERROR: No batch dir specified. Usage: './run_align_batch.sh batch_dir file_filter_str'"
  exit 1
fi

# If no input filter given, look for
# ".ua.bam" and "_1.fq.gz" files
if [[ $input_filter == "" ]] ; then

   export filefilt="${batchdir}/*_1.fq.gz"
   export n_input_files=$(ls -l $filefilt 2> /dev/null | wc -l)

   if [[ $n_input_files == 0 ]] ; then
      echo "ERROR: No input files found in folder '$batchdir'."
      exit 1
   fi

   export input_filter=$(basename "$filefilt")

   echo "INFO: Using default input file filter: '${input_filter}'."
   echo "INFO: Found $n_input_files files (samples) to process."
fi

if [[ $input_filter == *"hg38.bam" ]] ; then
  echo "ERROR: input files cannot contain '.hg38.bam', because it is reserved for output files."
  exit 1
fi

if [[ $input_filter == *"hg38.rg.bam" ]] ; then
  echo "ERROR: input files cannot contain '.hg38.rg.bam', because it is reserved for output files."
  exit 1
fi

if [[ $input_filter != *".bam" && $input_filter != *".fq.gz" ]] ; then
   echo "ERROR: Unknown input file type: '${input_filter}'."
   exit 1
fi

# Remove trailing slash if it exists
export batchdir=$(echo $batchdir | sed 's:/*$::')

if [ ! -d "$batchdir" ]; then
    echo "ERROR: Batch directory '$1' not found."
    exit 1
fi

# Enter batch dir (or die)
cd "$batchdir" || exit 1
echo "INFO: Working in dir '$(pwd)'."

export n_fastq1_files=$(ls -l *_1.fq.gz 2> /dev/null | wc -l)

if [[ $n_fastq1_files -eq 0 ]] ; then
  echo "ERROR: No fastq files found for alignment processing for batch:'${batchdir}'"
  exit 1
fi

# Create directory for parallel2 Linking
export linked_dir_arg=""

export rg_tag=""
export sample_name=""
export sample_id=""

# Function for parsing bam read group information from fastq filename
# Expects file name format:
# {SAMPLENAME}_{BARCODE}_{FLOWCELL}_{LANE}_1.fq.gz
# Set values in global variables
format_read_group() {
  filename=$1

  # https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups

  #ID = Read group identifier
  #PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
  #SM = Sample
  #PL = Platform/technology used to produce the read
  #LB = DNA preparation library identifier

  #Example: D052-FFPE_FSFP202427062-1r_HCCGMDSXY_L3_1.fq.gz

  arr=(${filename//_/ })

  local sample=$(echo ${arr[0]})
  local sample_barcode=$(echo ${arr[1]})
  local flowcell=$(echo ${arr[2]})
  local lane=$(echo ${arr[3]})
  local mate=$(echo ${arr[4]})

  parse_errors=0

  if [[ -z "$sample" ]] ; then echo "ERROR: 'sample' not found."; parse_errors=1; fi
  if [[ -z "$sample_barcode" ]] ; then echo "ERROR: 'sample_barcode' not found."; parse_errors=1; fi
  if [[ -z "$flowcell" ]] ; then echo "ERROR: 'flowcell' not found."; parse_errors=1; fi
  if [[ -z "$lane" ]] ; then echo "ERROR: 'lane' not found."; parse_errors=1; fi
  if [[ -z "$mate" ]] ; then echo "ERROR: 'mate' not found."; parse_errors=1; fi
  if [[ "$mate" != "1.fq.gz" ]] ; then echo "ERROR: file '$filename' should have '_1.fq.gz' suffix. ($mate)"; parse_errors=1; fi

  if [[ $parse_errors -eq 1 ]] ; then echo "ERROR: input fastq files should be renamed to format: '{SAMPLENAME}_{BARCODE}_{FLOWCELL}_{LANE}_1.fq.gz'"; exit 1; fi

  local ID="${flowcell}.${lane}"
  local PU="${flowcell}.${lane}.${sample_barcode}"
  local PL="${sequencer_str}"
  local SM="${sample}"
  local LB="0"

  export sample_name="$sample" #sample group to merge together after alignment
  export sample_id="$SM.$ID"
  export rg_tag="@RG\tSM:${SM}\tID:${ID}\tLB:${LB}\tPL:${PL}\tPU:${PU}"
}


generate_linked_dir_name() {

  if [[ ! -z "$parallel_linked_base_dir" ]] ; then
    export linked_dir_arg=$(echo $parallel_linked_base_dir | sed 's:/*$::')
    export dt=`date '+%Y-%m-%d_%H.%M.%S'`
    export linked_dir_arg="${linked_dir_arg}/job_${dt}"
    # Do not create directory, parallel will assume it is a running job
    #mkdir $linked_dir_arg || exit 1
    echo -e "INFO: Creating linkable directory"
    echo -e "      Use '${GREEN}$parallel2_bin -n $n_parallel_alignments -L $linked_dir_arg${NC}' for linking."
    export linked_dir_arg="-L $linked_dir_arg"
    echo ""
  fi
}


export n_mates=0
export n_mates_to_align=0
export n_merge=0
export n_done=0

echo "" > sample_names.txt
echo "" > sample_groups.txt

# Run checks
for mate1 in $(echo $input_filter)
do

   export mate2=${mate1/_1.fq.gz/_2.fq.gz}

   if [[ ! -s "$mate2" ]] ; then
      echo "ERROR: Fastq mate 2 file not found for mate 1: '$mate1'."
      exit 1
   fi

   format_read_group $mate1

   echo $sample_name >> sample_names.txt

done

#cat sample_names.txt | uniq | tr '\n' ' ' > sample_groups.txt
cat sample_names.txt | uniq > sample_groups.txt

# Print
for mate1 in $(echo $input_filter)
do

   export n_mates=$((n_mates+1))

   format_read_group $mate1

   export align_outputfile=${sample_id}.hg38.rg.bam
   export merged_align_outputfile=${sample_name}.merged.hg38.rg.bam

   # Count occurances
   n_mergeable=$(egrep -c $sample_name sample_names.txt)

   if [[ $n_mergeable -lt 1 ]] ; then
        echo -e "ERROR: cannot find sample '$sample_name'." # should not happen
        exit 1
   fi

   process="ERROR"
   file_size=""

   if [[ -s "$merged_align_outputfile" ]] ; then
      process="${GREEN}DONE${NC}"
      n_done=$((n_done+1))
      file_size=$(du -D "$merged_align_outputfile" | echo "1000 * $(cut -f 1)" | bc | numfmt --to=si)

   elif [[ -s "$align_outputfile" ]] ; then
      file_size=$(du -D "$align_outputfile" | echo "1000 * $(cut -f 1)" | bc | numfmt --to=si)
      if [[ $n_mergeable -gt 1 ]] ; then
          process="${PURPLE}MERGE${NC}"
      else
          process="${GREEN}DONE${NC}"
      fi
      n_merge=$((n_merge+1))
   else
      if [[ $n_mergeable -gt 1 ]] ; then
          process="${YELLOW}ALIGN${NC}+${PURPLE}MERGE${NC}"
      else
          process="${YELLOW}ALIGN${NC}"
      fi
      n_mates_to_align=$((n_mates_to_align+1))
   fi

   if [ ! -z $file_size ] ; then file_size="$file_size"; fi

   printf "%2i: ${GREEN}%-50s${NC} %-12s [%-5s]%7s\n" $n_mates ${mate1/_1.fq.gz/} $sample_name $process $file_size

done



# Check number of files
echo ""
if [[ $n_done -eq $n_mates ]] ; then
    echo "All files are already aligned and merged. Nothing to do."
    exit 0
elif [[ $n_merge -eq $n_mates ]] ; then
    echo "INFO: All files are already aligned. Only merging files."
else
    echo "INFO: Running ${n_parallel_alignments} pairs in parallel with $bwa_mem_threads CPUs each."
    if [[ $(($bwa_mem_threads*12)) -gt ${n_ramgb} ]] ; then
        echo "WARNING: Each pipeline run uses >12GB of RAM, system has only ${n_ramgb}GB."
    fi
fi



# Test shared folder if given
if [[ ! -z $parallel_linked_base_dir ]] ; then
    echo -e "${GREEN}INFO${NC}: Linked directory for processing will be created in '$parallel_linked_base_dir'"
    echo -e "      See '$parallel2_bin --help' for more details."
    export parallel_linked_base_dir=$(echo ${parallel_linked_base_dir} | sed 's:/*$::')
    mkdir -p $parallel_linked_base_dir || exit 1
    export test_file="${parallel_linked_base_dir}/test_write"
    if [[ ! -e "$test_file" ]] ; then
      touch "$test_file"
    fi
    if [[ ! -w "$test_file" ]] ; then
       echo "ERROR: Write permission denied to '$parallel_linked_base_dir'"
       exit 1
    fi
fi


#Check for screen
if [[ $(echo $STY) == "" ]] ; then
  echo "ERROR: This script should be run inside screen only."
  exit 1
fi

read -p "INFO: Run alignment and merge for these fastq-pairs? [Y/N] " -n 1 -r

if [[ $REPLY =~ ^[Yy]$ ]]
then
  echo " [YES]"
else
  echo " [NO]"
  echo "Exitting."
  exit 0
fi



###########################
# Align FASTQ-mates (1&2) #
###########################
echo "INFO: Starting alignment for $n_fastq1_files mate pairs..."

export n_files_to_align=0

mkdir -p flagstat

if [[ $n_mates_to_align -gt 0 ]] ; then

    export parallel_run_file="do_run_bwa_mem.cmd"
    echo "INFO: Writing parallel2 commands to file: \"$batchdir/$parallel_run_file\"."
    echo "" > $parallel_run_file

    # Create parallel run file
    for mate1 in $(echo $input_filter)
    do
        echo "" >> $parallel_run_file
        # Update rg_tag, sample_id
        format_read_group $mate1
        mate2=${mate1/_1.fq.gz/_2.fq.gz}
        outfile=${sample_id}.hg38.rg.bam
	echo "#File $sample_id:" >> $parallel_run_file
        align_cmd="[ ! -s ${sample_id}.hg38.rg.bam ] && $bwa_mem_bin -t ${bwa_mem_threads} -R \"${rg_tag}\" $ref_genome ${mate1} ${mate2} | samtools fixmate -@ $bwa_mem_sort_threads -m - - | samtools sort -@ $bwa_mem_sort_threads - | samtools markdup -@ $bwa_mem_sort_threads - ${sample_id}.hg38.rg.bam;"
        echo $align_cmd >> $parallel_run_file
        echo "[ -s ${sample_id}.hg38.rg.bam ] && samtools depth $outfile | awk '{sum+=\$3} END { print \"$outfile =\",sum/NR,\"X\"}' > flagstat/${outfile/.bam/.X} &" >> $parallel_run_file
        echo "[ -s ${sample_id}.hg38.rg.bam ] && samtools flagstat -@ $bwa_mem_sort_threads $outfile  > flagstat/${outfile/.bam/.flagstat}" >> $parallel_run_file
        if [[ ! -s $outfile ]] ; then
          n_files_to_align=$((n_files_to_align+1))
        fi

    done

fi

# Read group examples:
# RGSM=${x/.bam/} RGLB=0 RGPL=illumina RGPU=None
# @RG     ID:1    LB:0    PL:illumina     SM:sample_name       PU:None
#export align_cmd="[ ! -s \${x/${fastq_file_filter/\*/}/.hg38.rg.bam} ] && $bwa_mem_bin -t ${bwa_mem_threads} -R \"@RG\tID:1\tLB:0\tPL:${sequencer_str}\tSM:\${x/_1.fq.gz/}\tPU:None\" /mnt/gemma/bin/resources/homo_sapiens/hg38 \${x} \${x/_1.fq.gz/_2.fq.gz} | samblaster | samtools view -uS - | samtools sort -@${bwa_mem_sort_threads} -o \${x/${fastq_file_filter/\*/}/.hg38.rg.bam} -;"

if [[ $n_files_to_align -eq 0 ]] ; then
    echo -e "INFO: All files already aligned. [${GREEN}OK${NC}]"
elif [[ $n_files_to_align -eq 1 ]] ; then
    echo "INFO: $n_files_to_align sample needs to be aligned."
else
    echo "INFO: $n_files_to_align samples need to be aligned."
fi

export logfile1="${batchdir}/align.log"
if [[ $n_files_to_align -gt 0 ]] ; then
  #########################
  # RUN ALIGNMENT

  generate_linked_dir_name # updates ${linked_dir_arg} 

  export processing_started=1
  start=`date +%s`

  #echo $fastq_file_filter | $parallel2_bin --job-name "align_$batchdir" ${linked_dir_arg} -n${n_parallel_alignments} "${align_cmd}"
  $parallel2_bin --job-name "align_$batchdir" -n $n_parallel_alignments -f $parallel_run_file > $logfile1 2>&1

  echo "INFO: Alignments complete"
  end=`date +%s`
  echo "INFO: Alignment runtime for $n_fastq1_files files (n${n_parallel_alignments} c${bwa_mem_threads}) was: $((end-start))"
fi

# TODO: add checks to see if everything was aligned correctly
#samtools flagstat output.bam

##########################
# Check output
export n_aligned_files=$(ls -l *.hg38.rg.bam 2> /dev/null | wc -l)
if [[ $n_aligned_files == 0 ]] ; then
  echo "ERROR: No aligned files found for further processing for batch:'${batchdir}'."
  exit 1
fi

if [[ $n_fastq1_files != $n_aligned_files ]] ; then
  export n_missing_aligned=$(( $n_fastq1_files - $n_aligned_files))
  echo "WARNING: Alignment of $n_missing_aligned files may be incomplete."
fi


####################
# Group and merge samples

if [[ ! -s sample_groups.txt ]] ; then
    echo "ERROR: Cannot find sample grouping information file (sample_groups.txt) for merging samples."
    echo "     : Unable to merge any samples without grouping information."
    exit 1
fi

export parallel_run_file2="do_run_merge_bams.cmd"
echo "INFO: Writing parallel2 commands to file: \"$batchdir/$parallel_run_file2\"."
echo "" > $parallel_run_file2

mkdir -p flagstat

n_groups_to_merge=0
n_files_to_merge=0
n_already_merged=0
open_pipes=()
output_files=()

while read group
do

    input_str=""
    n_files=0
    n_index=1
    group_error=0

    # Skip whitespace rows
    if [[ -z "${group// }" ]] ; then continue; fi

    output_file="${group}.merged.hg38.rg.bam"

    if [[ -s ${output_file} ]] ; then
      echo -e "INFO: Sample group '${group}' already merged to file '${output_file}' [${YELLOW}SKIP${NC}]"
      n_already_merged=$((n_already_merged+1))
      continue
    fi

    n_fastq1_files_for_group=$(ls -l ${group}*_1.fq.gz | wc -l)


    for file in $(echo ${group}*.hg38.rg.bam)
    do

        # Echo leaves asterisk in $file if no files are found
        if [[ $(echo $file | egrep -c '\*') -gt 0 ]] ; then
          continue
        fi

        # Skip already merged files
        if [[ $file == *".merged."* ]] ; then
            continue
        elif [[ ! -s $file ]] ; then
            echo "ERROR: File '$file' has zero filesize."
            group_error=1
            break
        fi

        #input_str="$input_str I=${file}" # MergeSamFiles
        input_str="$input_str ${file}" # samtools merge
        n_files=$((n_files + 1))
    done

    if [[ $group_error -eq 1 ]] ; then
      echo "ERROR: Sample group '${group}' could not be merged due to errors."
      continue
    fi

    if [[ $n_files -eq 0 ]] ; then
      echo -e "${RED}ERROR${NC}: No bam files to merge found for sample group '${group}'"
      continue
    fi

    # Do not merge incomplete data
    if [[ $n_fastq1_files_for_group -ne $n_files ]] ; then
      echo -e "${RED}ERROR${NC}: Sample group '${group}' has $n_files bam files but $n_fastq1_files_for_group first mate files."
      continue
    fi


    if [[ $n_files -eq 1 ]] ; then
      echo -e "INFO: Sample group '${GREEN}${group}${NC}' has only one bam file"
      # Create a symbolic link with ".merged" tag in the filename.
      single_file=$(echo ${group}*.hg38.rg.bam)
      echo -e "INFO: Creating symbolic link for merged file: '${single_file}' -> '${GREEN}${output_file}${NC}'"
      ln -s $single_file $output_file
      continue
    fi

    # Removing duplicates requires two iterations with MergeSamFiles
    #echo "[ ! -s ${output_file} ] && java -jar $picard_jar MergeSamFiles $input_str O=$output_file" >> $parallel_run_file

    if [[ ! -s ${output_file} ]] ; then

        echo -e "INFO: Merging ${GREEN}${n_files}${NC} bam files for sample group '${GREEN}${group}${NC}'"
        n_groups_to_merge=$((n_groups_to_merge + 1))
        n_files_to_merge=$((n_files + n_files_to_merge))

        echo "# Merge sample group ${group}" >> $parallel_run_file2

        # Create a special pipe file to avoid writing on disk
        pipe_file="pipe_${n_index}.bam"
        n_index=$((n_index+1))

        #open_pipes+=("$pipe_file")
        output_files+=("$output_file")

        # Make sure old pipe files are removed
        rm $pipe_file 2> /dev/null

        # Merge and remove duplicates from merged file
        #echo "mkfifo $pipe_file;" >> $parallel_run_file
        echo "samtools merge -@ 6 -u - $input_str | samtools markdup -@ $bwa_mem_sort_threads - ${output_file}; " >> $parallel_run_file2
        # Leave these quality checks running in the background
        # Calc stats
        echo "[ -s ${output_file} ] && samtools flagstat -@ $bwa_mem_sort_threads $output_file > flagstat/${output_file/.bam/.flagstat} &" >> $parallel_run_file2
        # Calc depth
        echo "[ -s ${output_file} ] && samtools depth $output_file | awk '{sum+=\$3} END { print \"$output_file =\",sum/NR,\"X\"}' > flagstat/${output_file/.bam/.X} &" >> $parallel_run_file2
        # Add dummy cmd to end task
        echo "[ -s ${output_file} ] && echo \"Alignmnet of $output_file is done.\"" >> $parallel_run_file2
    else
        echo -e "INFO: Sample group '${GREEN}${group}${NC}' already merged to file '$output_file'."
    fi

done < sample_groups.txt

# Remove temporary files
#for $pipe in "${open_pipes[@]}"
#do
#    rm $pipe 2> /dev/null
#done

if [[ $n_already_merged -gt 0 ]] ; then
    n_total_files=$((n_already_merged+n_groups_to_merge))
    echo "INFO: Found ${n_already_merged}/${n_total_files} merged bam files"
fi

if [[ $n_groups_to_merge -eq 0 ]] ; then
    echo "INFO: No bam files to merge found."
else

    echo -e "INFO: Merging a total of $n_files_to_merge files in $n_groups_to_merge groups..."

    generate_linked_dir_name # updates ${linked_dir_arg}

    #########################
    # RUN MERGER
    $parallel2_bin --job-name "merge_$batchdir" -n 4 -f $parallel_run_file2

    echo "INFO: Done merging $n_files_to_merge files to $n_groups_to_merge groups."
fi

# Check output
n_files_found=0
n_files_total=0
for outf in "${output_files[@]}"
do
    n_files_total=$((n_files_total+1))
    if [[ -s $outf ]] ; then
        n_files_found=$((n_files_found+1))
    fi
done

if [[ $n_files_found -eq $n_files_total ]] ; then
    echo -e "INFO: All $n_files_found/$n_files_total output files were generated [${GREEN}OK${NC}]"
elif [[ $n_files_found -gt 0 ]] ; then
    echo -e "ERROR: Only $n_files_found/$n_files_total output files found.  [${RED}FAIL${NC}]"
    exit 1
else
    echo -e "ERROR: No output files found.  [${RED}FAIL${NC}]"
    exit 1
fi

# RUN BQSR on the samples next!
echo ""
echo "INFO: Folder 'flagstat' contains statistics on merged bams."
echo "INFO: Filesnames ending in .X have average sequencing depth information."
echo "All done."
