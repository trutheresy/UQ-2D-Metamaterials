#PBS -N Han's script
#PBS -1 nodes=1:ppn=10
#PBS -1 mem=10G
cd $PBS_O_WORKDIR 
/usr/local/MATLAB/R2022a/bin/matlab -nodisplay -nodesktop -batch han_research_script