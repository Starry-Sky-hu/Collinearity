#!/bin/bash
#SBATCH --partition=low,big
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J work.sh
#SBATCH --qos=queue1
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Wed Feb  9 22:04:53 2022

# @File Name: work.sh
# @Description:
##################################################################
source /public/agis/huangsanwen_group/lihongbo/software/miniconda3/bin/activate emboss

cp ../Anisodus_luridush/Solanum_tuberosumDM_Anisodus_luridush.block_pos .
orient=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $8}'))
ref_chrom=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $1}'))
ref_start=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $2}'))
ref_end=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $3}'))
query_chrom=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $4}'))
query_start=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $5}'))
query_end=($(cat Solanum_tuberosumDM_Anisodus_luridush.block_pos | awk '{print $6}'))

length=${#ref_chrom[@]}
let max=length-1
for i in `seq 0 ${max}`
do
	singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 \
		seqkit subseq /home/huyong/SolanaceaeGenomeAnalyze/reference/01_Ref/Solanum_tuberosumDM.fa \
		--chr ${ref_chrom[i]} -r ${ref_start[i]}:${ref_end[i]} -o ref_${i}.fa
	singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 \
		seqkit subseq /home/huyong/SolanaceaeGenomeAnalyze/reference/01_Ref/Anisodus_luridush.fa \
		--chr ${query_chrom[i]} -r ${query_start[i]}:${query_end[i]} -o query_${i}.fa
	if [[ ${orient[i]} == "positive" ]]; then
		diffseq -rformat excel -rusashow3 -wordsize 10 -globaldifferences ref_${i}.fa query_${i}.fa ref_query_${i}.excel ref_${i}.diff query_${i}.diff
	else
		diffseq -rformat excel -rusashow3 -wordsize 10 -globaldifferences ref_${i}.fa query_${i}.fa -sreverse2 ref_query_${i}.excel ref_${i}.diff query_${i}.diff
	fi
	sed -i 's/\[::r\]//g' ref_query_${i}.excel
	~/software/output_vcf_sv_complex_from_diffseq_x.py ref_query_${i}.excel ref_${i}.diff query_${i}.diff ref_${i}.fa query_${i}.fa ref_query_${i}.vcf_diffseq ref_query_INS_${i}.xls ref_query_DEL_${i}.xls ref_query_COMPLEX_${i}.xls
done

