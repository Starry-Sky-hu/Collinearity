#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 52
#SBATCH -J work.sh
#SBATCH --qos=queue1
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Thu Jan 27 17:38:10 2022

# @File Name: work.sh
# @Description:
##################################################################

ln -s ~/SolanaceaeGenomeAnalyze/collinearity/genetribe/Anisodus_luridush/Solanum_tuberosumDM_Anisodus_luridush.one2one .
ln -s ~/SolanaceaeGenomeAnalyze/collinearity/gffrename/Anisodus_luridush/Anisodus_luridush.all.rename.gff3 .
ln -s ~/SolanaceaeGenomeAnalyze/collinearity/genetribe/ref/Solanum_tuberosumDM.Primary.gff3 .
cut -f1 Solanum_tuberosumDM_Anisodus_luridush.one2one > ref_gene
cut -f2 Solanum_tuberosumDM_Anisodus_luridush.one2one > query_gene

mkdir colgene
ref=($(cat ref_gene))
query=($(cat query_gene))
num=${#ref[@]}
let num=num-1
for i in `seq 0 ${num}`
do
	mkdir colgene/${ref[i]}
	echo "#!/bin/bash" > colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "source activate collinerity" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	
	echo "cat <(grep \"${ref[i]}\" ../../Solanum_tuberosumDM.Primary.gff3 | grep \"CDS\" | sort -k4n | sed -n '1p' | awk '{print \$1,\$4}' ORS='\t' OFS='\t') <(grep \"${ref[i]}\" ../../Solanum_tuberosumDM.Primary.gff3 | grep \"CDS\" | sort -k4n | sed -n '\$p' | awk -F '\t|\tParent=' '{print \$5,\$9,\$8,\$7}' OFS='\t') > ${ref[i]}.bed" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "cat <(grep \"${query[i]}\" ../../Anisodus_luridush.all.rename.gff3 | grep \"CDS\" | sort -k4n | sed -n '1p' | awk '{print \$1,\$4}' ORS='\t' OFS='\t') <(grep \"${query[i]}\" ../../Anisodus_luridush.all.rename.gff3 | grep \"CDS\" | sort -k4n | sed -n '\$p' | awk -F '\t|\tParent=' '{print \$5,\$9,\$8,\$7}' OFS='\t') > ${query[i]}.bed" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	
	echo ""  >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "ref_orient=\$(cat ${ref[i]}.bed | awk '{print \$6}')" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "if [ \${ref_orient} == \"-\" ]; then" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "	awk '{ print \$1,\$2,\$3,\$4,\$5\"\\t+\"}' OFS='\t' ${ref[i]}.bed > ${ref[i]}.bed_ori" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "	singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 seqkit subseq /home/huyong/SolanaceaeGenomeAnalyze/reference/01_Ref/Solanum_tuberosumDM.fa --bed ${ref[i]}.bed_ori -u 5000 -d 5000 > ${ref[i]}.fa" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "	awk '{ if(\$6==\"-\"){print \$1,\$2,\$3,\$4,\$5\"\\t+\"} else{print \$1,\$2,\$3,\$4,\$5\"\\t-\"} }' OFS='\t' ${query[i]}.bed > ${query[i]}.bed_ori" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "	singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 seqkit subseq /home/huyong/SolanaceaeGenomeAnalyze/reference/01_Ref/Anisodus_luridush.fa --bed ${query[i]}.bed_ori -u 5000 -d 5000 > ${query[i]}.fa" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "else" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "	singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 seqkit subseq /home/huyong/SolanaceaeGenomeAnalyze/reference/01_Ref/Solanum_tuberosumDM.fa --bed ${ref[i]}.bed -u 5000 -d 5000 > ${ref[i]}.fa" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "	singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 seqkit subseq /home/huyong/SolanaceaeGenomeAnalyze/reference/01_Ref/Anisodus_luridush.fa --bed ${query[i]}.bed -u 5000 -d 5000 > ${query[i]}.fa" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "fi" >> colgene/${ref[i]}/${query[i]}_mafft.sh
    echo "" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "cat ${ref[i]}.fa ${query[i]}.fa > mafft.fa" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "mafft --thread 1 --auto mafft.fa > mafft.out" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "singularity exec ~/contianer/GenomeAssemblyContainer_v0.2 seqkit seq -w 0 mafft.out > mafft_w0.out" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	echo "python /home/huyong/software/MSA2vcf.py mafft_w0.out ${ref[i]}.bed Soltu.DM Anisodus_luridush ${query[i]}_snp_indel.vcf" >> colgene/${ref[i]}/${query[i]}_mafft.sh
	
	cd colgene/${ref[i]}
	sleep 1
	nohup bash ${query[i]}_mafft.sh &
	cd ../../
done

ls colgene/*/*vcf > vcf_list
python vcf_merge.py vcf_list Anisodus_luridush_merge.vcf Anisodus_luridush




