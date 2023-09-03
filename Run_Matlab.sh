#PBS -N Han_script
cd $PBS_O_WORKDIR || exit 1
/usr/local/MATLAB/R2022a/bin/matlab -nodisplay -nodesktop -batch han_script_geometry_material_properties