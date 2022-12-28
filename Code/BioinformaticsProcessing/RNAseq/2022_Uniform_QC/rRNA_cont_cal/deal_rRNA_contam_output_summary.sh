
#file=/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/rRNA_contam_cal/TGL_n_OICR2022_rrna_mapped/*_flagstat.tsv
file=$PWD/all_flagstat_tsv/*_flagstat.tsv

outfile='all_rRNA_contam_summary.txt'
rm all_rRNA_contam_summary.txt

echo 'sample_id rrna_contam_perct' >>$outfile

for i in $file
do
 sample_id=$(echo $i | awk -F '/' '{print $10}' | sed s/_flagstat.tsv//g)
   echo -n $sample_id >>$outfile
   echo -n -e "\t"  >>$outfile
   exonic_prct=`head -5 $i |tail -1|awk -F '[(|)]' '{print $2}'|awk -F '%' '{print $1}'`
   echo $exonic_prct >>$outfile
done

