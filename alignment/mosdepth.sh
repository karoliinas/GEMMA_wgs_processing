for bam in *.bam
do
  export id=$(echo $bam | awk -F"." '{print $1}'
  /data/apps/mosdepth -n -t 8 --fast-mode mosdepth/$id $bam
done

