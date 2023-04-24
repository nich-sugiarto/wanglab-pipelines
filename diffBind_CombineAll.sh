#!/bin/bash

# This script creates a union peak set for diffBind output differential peaks
# Peaks identified as differential from the Ctrl at any timepoint are included in these peak sets

suffix1=Ctrl_UT #This includes KOs in merged sets
folder=$(pwd)


# Set up necessary files
mkdir -p PBS
mkdir -p log
mkdir -p diffBind/Combined

cd diffBind

# Loop over all folders 
for subfolder in *UT; do
factor=$(echo "$subfolder" | awk -F '_' '{print $1}' )
echo "factor: "$factor
mkdir -p Combined/$factor
mkdir -p Combined/All
cd $subfolder
cd DBA_DESEQ2
cat downregulated.bed >> ../../Combined/$factor/$factor'_DESEQ_downregulated.bed' #appends to the factor's combined downregulated bed file
cat downregulated.bed >> ../../Combined/$factor/$factor'_DESEQ_differential.bed' #appends to the factor's combined differential bed file
cat upregulated.bed >> ../../Combined/$factor/$factor'_DESEQ_upregulated.bed' #appends to the factor's combined upregulated bed file
cat upregulated.bed >> ../../Combined/$factor/$factor'_DESEQ_differential.bed' #appends to the factor's combined differential bed file
cat downregulated.bed >> ../../Combined/All/'Combined_DESEQ_downregulated.bed' #appends to the total combined downregulated bed file
cat downregulated.bed >> ../../Combined/All/'Combined_DESEQ_differential.bed' #appends to the total combined differential bed file
cat upregulated.bed >> ../../Combined/All/'Combined_DESEQ_upregulated.bed' #appends to the total combined upregulated bed file
cat upregulated.bed >> ../../Combined/All/'Combined_DESEQ_differential.bed' #appends to the total combined differential bed file
cd ..

cd DBA_EDGER
cat downregulated.bed >> ../../Combined/$factor/$factor'_EDGER_downregulated.bed' #appends to the factor's combined downregulated bed file
cat downregulated.bed >> ../../Combined/$factor/$factor'_EDGER_differential.bed' #appends to the factor's combined differential bed file
cat upregulated.bed >> ../../Combined/$factor/$factor'_EDGER_upregulated.bed' #appends to the factor's combined upregulated bed file
cat upregulated.bed >> ../../Combined/$factor/$factor'_EDGER_differential.bed' #appends to the factor's combined differential bed file
cat downregulated.bed >> ../../Combined/All/'Combined_EDGER_downregulated.bed' #appends to the total combined downregulated bed file
cat downregulated.bed >> ../../Combined/All/'Combined_EDGER_differential.bed' #appends to the total combined differential bed file
cat upregulated.bed >> ../../Combined/All/'Combined_EDGER_upregulated.bed' #appends to the total combined upregulated bed file
cat upregulated.bed >> ../../Combined/All/'Combined_EDGER_differential.bed' #appends to the total combined differential bed file
cd ..
cd ..

done

cat >${folder}/PBS/Diffbind_Combine.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=Diffbind_Combine
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/Diffbind_Combine_%j.txt -e ${folder}/log/Diffbind_Combine_%j.err.txt

#------- END OF HEADER -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

cd ${folder}
cd diffBind/Combined


EOF
cd ${folder}
cd diffBind/Combined

##Merge factor specific bed files
for f in *; do
cat >>${folder}/PBS/Diffbind_Combine.sbatch <<EOF
cd ${f}
EOF
cd $f
    for bed in *.bed; do
    cat >>${folder}/PBS/Diffbind_Combine.sbatch <<EOF
        awk '(NR>1) && (\$2 > 0 ) ' ${bed} > new_${bed}
        bedtools sort -i new_${bed} >sorted_${bed}
        rm ${bed}
        rm new_${bed}
        bedtools merge -i sorted_${bed} >${bed}
        rm sorted_${bed}
        
EOF
    done
    cd ..
    cat >>${folder}/PBS/Diffbind_Combine.sbatch <<EOF
        bedtools intersect -a ../../peakScores.bed -b ${f}_EDGER_differential.bed >${f}_Diff_EdgeR.bed
        bedtools intersect -a ../../peakScores.bed -b ${f}_DESEQ_differential.bed >${f}_Diff_DESeq.bed
        cd ..
EOF
done
cat >>${folder}/PBS/Diffbind_Combine.sbatch <<EOF
    cd All
    bedtools intersect -a ../../peakScores.bed -b Combined_EDGER_differential.bed >All_Diff_EdgeR.bed
    bedtools intersect -a ../../peakScores.bed -b Combined_DESEQ_differential.bed >All_Diff_DESeq.bed

EOF
sbatch ${folder}/PBS/Diffbind_Combine.sbatch