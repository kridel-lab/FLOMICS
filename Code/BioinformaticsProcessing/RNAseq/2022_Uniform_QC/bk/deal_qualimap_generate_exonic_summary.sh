
file=/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/RNAseqQC_out_path/*/rnaseq_qc_results.txt
outfile='RNA_qualimap_summary.txt'
rm RNA_qualimap_summary.txt

echo 'sample_id exonic_perct' >>$outfile
for i in $file
do
#   echo $i
   sample_id=$(echo $i | awk -F '/' '{print $9}')
   echo -n $sample_id >>$outfile
   echo -n -e "\t"  >>$outfile
   exonic_prct=`grep 'exonic' $i |awk '{print $4}'|awk -F '[(|)]' '{print $2}'`
   echo $exonic_prct >>$outfile
done
