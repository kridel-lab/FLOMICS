file="done_ids.txt";

for i in `cat ${file}`
do
  echo $i
  `grep $i bam_files_list.txt >>done_bam.txt`
done

`grep -vFf done_bam.txt bam_files_list.txt >unfinished_bam_list.txt`
