
#create the symlinks folder if there is not one
work_dir=$PWD


#get the list of the raw bam files
file=${work_dir}/rawdata_symlinks/*_R2.fastq.gz

##loop each sample
for i in $file
do
  echo $i
  if [[ $i =~ ${work_dir}/rawdata_symlinks/(.*)_R2.fastq.gz ]] ;then
      num=${BASH_REMATCH[1]}
      #num is the inital id, like F24655
#      echo $num
     sample_id=$(grep $num ${work_dir}/id_mapping.txt | awk '{print $2}')
##sample id is the third colum of the file "id_mapping.txt"
     printf "$num\t$sample_id\n"
#rename the symlink
     mv $i ${work_dir}/rawdata_symlinks/"${sample_id}_R2.fastq.gz"
   else
     echo "not Found"
  fi

done
