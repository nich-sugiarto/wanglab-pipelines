#!/bin/bash

# This script creates a union peak set for diffBind output differential peaks
# Peaks identified as differential from the Ctrl at any timepoint are included in these peak sets

suffix1=Ctrl_UT #This includes KOs in merged sets
folder=$(pwd)


# Set up necessary files
mkdir -p PBS
mkdir -p log
mkdir -p diffBind/Combined


cd diffBind/Combined
for f in *; do
rm $f
done
cd ..

# Loop over all folders 
for subfolder in *UT; do
cd $subfolder
cat downregulated.bed >> ../Combined/'Combined_downregulated.bed' #appends to the total combined downregulated bed file
cat downregulated.bed >> ../Combined/'Combined_differential.bed' #appends to the total combined differential bed file
cat upregulated.bed >> ../Combined/'Combined_upregulated.bed' #appends to the total combined upregulated bed file
cat upregulated.bed >> ../Combined/'Combined_differential.bed' #appends to the total combined differential bed file
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
for bed in *.bed; do
    cat >>${folder}/PBS/Diffbind_Combine.sbatch <<EOF
        #this removes any negative entries and sci notation(IDK why there were any...)
        awk '!(\$2 ~ /e/) && !(\$3 ~ /e/) && (\$2 > 0)' ${bed} > new_${bed}
        bedtools sort -i new_${bed} >sorted_${bed}
        rm ${bed}
        rm new_${bed}
        bedtools merge -i sorted_${bed} >${bed}
        rm sorted_${bed}
        #
        #
EOF
done
    cat >>${folder}/PBS/Diffbind_Combine.sbatch <<EOF
        awk '!(\$2 ~ /e/) && !(\$3 ~ /e/) && (NR>1) && (\$2 > 0)' ../peakScores.bed > ../ScoresOnly.bed
        #Next line gets the header for differential scores 
        awk '(NR==1)' ../peakScores.bed > Differential_Scores.bed
        bedtools intersect -a ../ScoresOnly.bed -b Combined_differential.bed >>Differential_Scores.bed
EOF

sbatch ${folder}/PBS/Diffbind_Combine.sbatch