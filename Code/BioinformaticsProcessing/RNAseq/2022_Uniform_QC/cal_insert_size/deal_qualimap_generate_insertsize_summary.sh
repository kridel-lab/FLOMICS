
file=/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/bamQC_out_path/*/genome_results.txt
outfile='insertsize_qualimap_summary.txt'
rm insertsize_qualimap_summary.txt
echo -n 'sample_id' >>$outfile
echo -n -e "\t"  >>$outfile
echo -n 'median_insert_size' >>$outfile
echo -n -e "\t"  >>$outfile
echo 'mean_insert_size' >>$outfile
for i in $file
do
#   echo $i
   sample_id=$(echo $i | awk -F '/' '{print $9}')
   echo -n $sample_id >>$outfile
   echo -n -e "\t"  >>$outfile
   median_insert=`grep 'median insert size' $i |awk '{print $5}'`
   echo -n $median_insert >>$outfile
   echo -n -e "\t"  >>$outfile
   mean_insert=`grep 'mean insert size' $i |awk '{print $5}'`
   echo $mean_insert >>$outfile
    
done
