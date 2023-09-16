#PBS -N Han_script
cd $PBS_O_WORKDIR || exit 1
/shared/Apps/MATLAB/R2023a/bin/matlab -nodisplay -nodesktop -batch parallel_script_geometry_material_properties_K_G