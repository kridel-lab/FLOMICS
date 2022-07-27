rm outfile
file=/cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/RSeQC_analysis/*.rrna.txt

for i in $file
do
#  echo $i
  if [[ $i =~ /cluster/projects/kridelgroup/FLOMICS/2022_combined_RNA-Seq_analysis/RNA-seq_QC/RSeQC_analysis/(.*).rrna.txt ]] ;then
 num=${BASH_REMATCH[1]}
#  echo -e "$num\t\c" >>outfile
#echo -n "$num\t" >>outfile
   echo -n $num    >>outfile
   echo -n -e "\t" >>outfile
   total=`head -1 $i | awk '{print $3}'` 
   echo -n $total  >>outfile
   echo -n -e "\t" >>outfile
   rrna=`head -2 $i |tail -1| awk '{print $7}' | awk -F ':' '{print $2}'`
   echo -n $rrna >>outfile
   echo -n -e "\t" >>outfile
   not_m=`head -3 $i |tail -1| awk '{print $8}' | awk -F ':' '{print $2}'`
   echo -n $not_m >>outfile
   echo -n -e "\t" >>outfile
  # let rrna_perct=$rrna/100
 #  rrna_perct=`expr $rrna \/ $total |bc`
   rrna_perct=`echo "scale=5;$rrna/$total*100"|bc`
#  rrna_perct=awk '{print '$rrna'/'$total'}' 
  echo $rrna_perct >>outfile
 
  fi

done

