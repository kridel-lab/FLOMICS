#v2 print both PF_ALIGNED_BASES ratio and PF_BASES ratio
#file=/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/rRNA_contam_cal/TGL_n_OICR2022_rrna_mapped/*_flagstat.tsv
file=$PWD/picard_RnaMetrics/*

outfile='all_picard_RnaMetrics_summary_v2.txt'
rm all_picard_RnaMetrics_summary_v2.txt

echo -n 'sample_id' >>$outfile
echo -n -e "\t"   >>$outfile
echo -n 'PF_BASES' >>$outfile
echo -n -e "\t"   >>$outfile
echo -n 'PF_ALIGNED_BASES' >>$outfile
echo -n -e "\t"   >>$outfile
echo -n 'CODING_BASES' >>$outfile
echo -n -e "\t"  >>$outfile
echo -n 'picard_RnaMetrics_perct' >>$outfile
echo -n -e "\t"   >>$outfile
echo 'picard_RnaMetrics_ALIGNED_perct' >>$outfile

for i in $file
do
 sample_id=$(echo $i | awk -F '/' '{print $9}')
   echo -n $sample_id >>$outfile
   echo -n -e "\t"  >>$outfile
   PF_BASES=`head -8 $i |tail -1 |awk '{print $1}'`
   echo -n $PF_BASES  >>$outfile
   echo -n -e "\t"  >>$outfile
   PF_ALIGNED_BASES=`head -8 $i |tail -1 |awk '{print $2}'`
   echo -n $PF_ALIGNED_BASES  >>$outfile
   echo -n -e "\t" >>$outfile
   CODING_BASES=`head -8 $i |tail -1 |awk '{print $4}'`
   echo -n $CODING_BASES >>$outfile
   echo -n -e "\t" >>$outfile
   exonic_prct=`echo "scale=5;$CODING_BASES/$PF_BASES*100"|bc`
   exonic_ALIGNED_prct=`echo "scale=5;$CODING_BASES/$PF_ALIGNED_BASES*100"|bc`
   #`head -5 $i |tail -1|awk -F '[(|)]' '{print $2}'|awk -F '%' '{print $1}'`
   echo -n $exonic_prct >>$outfile
   echo -n -e "\t" >>$outfile
   echo $exonic_ALIGNED_prct >>$outfile

done

#$1=PF_BASES	The total number of PF bases including non-aligned reads.
#$4=CODING_BASES	Number of bases in primary alignments that align to a non-UTR coding base for some gene, and not ribosomal sequence.
