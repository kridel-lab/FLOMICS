
file=/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/bamQC_out_path/*/genome_results.txt
outfile='insertsize_qualimap_summary.txt'
rm insertsize_qualimap_summary.txt
echo 'sample_id	median_insert_size' >>$outfile
for i in $file
do
#   echo $i
   sample_id=$(echo $i | awk -F '/' '{print $9}')
   echo -n $sample_id >>$outfile
   echo -n -e "\t"  >>$outfile
   insert=`grep 'median insert size' $i |awk '{print $5}'`
   echo $insert >>$outfile
done
