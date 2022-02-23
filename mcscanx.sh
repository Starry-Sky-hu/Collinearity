#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 13
#SBATCH --qos=queue1
##################################################################
# @Author: huyong
# @Created Time : Wed Dec 29 14:16:33 2021

# @File Name: test.sh
# @Description:
##################################################################

date
blastp -query query_n.pep.fa -db ../ref/Solanum_tuberosumDM.pep -out query_n_Solanum_tuberosumDM.blast  -max_target_seqs 5 -evalue 1e-4 -outfmt 6 -num_threads 13
date

awk '($3=="mRNA"){print}' query_n.Primary.gff3 | awk -F '\t|\tID=|;Parent=|;' '{print $1,$9,$4,$5}'  OFS='\t'|  awk '{$1="query_n_"$1}1'  | sed 's/ /\t/g'  > query_n.Primary_transcript.bed

cat ../ref/Solanum_tuberosumDM_Primary_transcript.bed query_n.Primary_transcript.bed > query_n_Solanum_tuberosumDM.gff

source activate R_Python
Rscript ~/software/Check_blastp.id_same.as.gff.R query_n_Solanum_tuberosumDM.blast query_n_Solanum_tuberosumDM.gff blast_gff_diff
conda deactivate

mkdir mcsanx
cd mcsanx

cp ../query_n_Solanum_tuberosumDM.gff .
cp ../query_n_Solanum_tuberosumDM.blast .
MCScanX -s 4 -a query_n_Solanum_tuberosumDM
perl ~/software/collinearity.pileup.pl query_n_Solanum_tuberosumDM.collinearity query_n_Solanum_tuberosumDM.gff query_n_Solanum_tuberosumDM_stat

sed -i 's/://g' query_n_Solanum_tuberosumDM_stat
sed -i 's/N=//g' query_n_Solanum_tuberosumDM_stat
sed -i 's/e_value=//g' query_n_Solanum_tuberosumDM_stat
sed -i 's/Setaria_italica//g' query_n_Solanum_tuberosumDM_stat

cat query_n_Solanum_tuberosumDM_stat | awk '{print $5,$10,$11,$1,$2,$3,$4}' OFS='\t' > query_n_Solanum_tuberosumDM_stat_infor

bedtools sort -i query_n_Solanum_tuberosumDM_stat_infor > query_n_Solanum_tuberosumDM_stat_infor.sort
bedtools merge -i query_n_Solanum_tuberosumDM_stat_infor.sort > query_n_Solanum_tuberosumDM_stat_infor.merge
bedtools merge -i query_n_Solanum_tuberosumDM_stat_infor.sort -d 100000 > query_n_Solanum_tuberosumDM_stat_infor.100kb.merge

cd ..

