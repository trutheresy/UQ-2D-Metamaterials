% clc
% clear all

% Get the full path to the currently executing script
scriptPath = mfilename('fullpath');

% Extract the path to the folder containing the script
[outputFolder, ~, ~] = fileparts(scriptPath);
%addpath('/Users/zhang/Documents/Duke/Research/2D-dispersion')
addpath('2D dispersion legacy')

save_outputs = true;
% custom_E = false;
% custom_Rho = false;
% custom_PR = false;
error_geometries = true;
mat_geo_coupled = true;

%%% DEFAULT PARAMETERS %%%
n_geometries = 1;
n_materials = 1;
% E_soft = 200e6*ones(1,n_materials);
% E_hard = 200e9*ones(1,n_materials);
% rho_soft = 1e3*ones(1,n_materials);
% rho_hard = 8e3*ones(1,n_materials);
% poisson_soft = 0*ones(1,n_materials);
% poisson_hard = 0.5*ones(1,n_materials);

%%% ALTERED PARAMETERS %%%
% importdata("scaled_default_geometry.mat")
% input_geometries = scaled_matrix_2;

%input_geometries = importdata('scaled_matrices.mat');
%input_data = pd_5_inputs;
%input_data_size = size(input_data);
%n_materials = input_data_size(2);
%n_geometries = length(input_geometries);
%E_soft = input_data(1,:);
%E_hard = input_data(1,:);
%rho_soft = input_data(1,:); %change back to correct rows later
%rho_hard = input_data(1,:);
%poisson_soft = input_data(1,:);
%poisson_hard = input_data(1,:);

%%% IMPORT DATA FILES %%%
% Specify the folder name and search string
folderName = 'gaussian 6+1 inputs quadrature rule 2nd geo study';
%folderName = 'gaussian 6+1 inputs mc study';
searchString = '_pd_2';
%searchString = '_mc_10000.mat';

% Get a list of all files in the specified folder
files = dir(fullfile(folderName, ['*' searchString '*']));

% Loop through each file
for i = 1:length(files)
    filePath = fullfile(folderName, files(i).name);

    % Load the contents of the MAT file
    load(filePath); % Variables will be loaded into the workspace
    
    % Display loaded variable names if needed
    disp(['Loaded variables from ' files(i).name ':']);
    
    % Optionally, you can further process the loaded variables here
end
%%% END IMPORT DATA FILES %%%
geometries = pd_2_geos; %change name each run
%geometries = mc_10000; %change name each run
if ndims(geometries) == 2
    n_geometries = 1;
    geometries_size = size(geometries)
elseif ndims(geometries) == 3
    geometries_size = size(geometries)
    n_geometries = geometries_size(1);
    geometries_size = geometries_size(2:end);
end

materials = pd_2_inputs'; %change name each run
%materials = mc_10000_inputs'; %change name each run
materials_size = size(materials)
n_materials = materials_size(1);

bg_size = zeros(n_materials, 1);
bg_top = zeros(n_materials, 1);
bg_bottom = zeros(n_materials, 1);

if mat_geo_coupled
    elapsedTime = 0;
    % Check if n_materials = n_geometries, as they should be in 1:1 pairs
    if n_geometries ~= n_materials
        error('Number of geometries is not equal to number of materials. Script terminated.');
    else
        disp([num2str(n_geometries) ' geometries coupled with materials'])
        E_soft = materials(:,1);
        E_hard = materials(:,2);
        rho_soft = materials(:,3);
        rho_hard = materials(:,4);
        pr_soft = materials(:,5);
        pr_hard = materials(:,6);

        tic;
        for i = 1:n_geometries
    
            [EIGENVALUE_DATA, WAVEVECTOR_DATA] = dispersion_data_custom_mat_geo_param( ...
                E_soft(i), ...
                E_hard(i), ...
                rho_soft(i), ...
                rho_hard(i), ...
                pr_soft(i), ...
                pr_hard(i), ...
                geometries(i,:,:));
                
            % Check for presence of bandgap between bands 3 and 4
            bg_bottom(i) = max(EIGENVALUE_DATA(:,3));
            bg_top(i) = min(EIGENVALUE_DATA(:,4));
            bg_size(i) = bg_top(i) - bg_bottom(i);
            if bg_size(i) > 0
                disp(['Bandgap present for input ' num2str(i)])
            else
                disp(['No bandgap for input ' num2str(i)])
            end
        end
        %toc;
        elapsedTime = toc;
        disp(['Compute time for ' num2str(n_geometries) ' input sets of geometry size ' num2str(geometries_size) ' is ' num2str(elapsedTime)])
    end

    %%% SAVE FILES - CHECK NAMES PRIOR TO RUN %%%
    if save_outputs
        cd(outputFolder);
        save(['bg_size_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(bg_size))], 'bg_size');
        save(['bg_bottom_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(bg_bottom))], 'bg_bottom');
        save(['bg_top_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(bg_top))], 'bg_top');
        save(['E_soft_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(E_soft))], 'E_soft');
        save(['E_hard_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(E_hard))], 'E_hard');
        save(['rho_soft_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(rho_soft))], 'rho_soft');
        save(['rho_hard_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(rho_hard))], 'rho_hard');
        save(['pr_soft_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(pr_soft))], 'pr_soft');
        save(['pr_hard_gaussian_' num2str(materials_size(2)) 'd_fp_5%_n' num2str(length(pr_hard))], 'pr_hard');
        save(['elapsed_time_' num2str(materials_size(2)) 'd_fp_n' num2str(n_materials)], 'elapsedTime')
    end
    %%% END SAVE FILES %%%
end

if ~mat_geo_coupled
    disp([num2str(n_geometries) ' geometries crossed with materials'])
    E_soft = materials(:,1);
    E_hard = materials(:,2);
    rho_soft = materials(:,3);
    rho_hard = materials(:,4);
    pr_soft = materials(:,5);
    pr_hard = materials(:,6);

    elapsedTime = 0;
    tic;
    for g = 1:n_geometries
        for i = 1:n_materials
            [EIGENVALUE_DATA, WAVEVECTOR_DATA] = dispersion_data_custom_mat_geo_param( ...
                materials(i,1), ...
                materials(i,2), ...
                materials(i,3), ...
                materials(i,4), ...
                materials(i,5), ...
                materials(i,6), ...
                geometries(g,:,:));
                
            % Check for presence of bandgap between bands 3 and 4
            bg_bottom(i) = max(EIGENVALUE_DATA(:,3));
            bg_top(i) = min(EIGENVALUE_DATA(:,4));
            bg_size(i) = bg_top(i) - bg_bottom(i);
            if bg_size(i) > 0
                disp(['Bandgap present for input ' num2str(i)])
            else
                disp(['No bandgap for input ' num2str(i)])
            end
        end
    end
    toc;
    elapsedTime = toc;
    disp(['Compute time for ' num2str(n_geometries) ' input sets of geometry size ' num2str(geometries_size) ' is ' num2str(elapsedTime)])

    %%% SAVE FILES - CHECK NAMES PRIOR TO RUN %%%
    if save_outputs
        cd(outputFolder);
        save(['bg_size_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(n_geometries*bg_size))], 'bg_size');
        save(['bg_bottom_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(n_geometries*bg_bottom))], 'bg_bottom');
        save(['bg_top_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(bg_top))], 'bg_top');
        save(['E_soft_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(E_soft))], 'E_soft');
        save(['E_hard_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(E_hard))], 'E_hard');
        save(['rho_soft_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(rho_soft))], 'rho_soft');
        save(['rho_hard_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(rho_hard))], 'rho_hard');
        save(['pr_soft_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(pr_soft))], 'pr_soft');
        save(['pr_hard_gaussian_' num2str(materials_size(2)) '_fp_5%_n' num2str(length(pr_hard))], 'pr_hard');
        save(['elapsed_time_' num2str(materials_size(2)) '_fp_n' num2str(n_materials)], 'elapsedTime')
    end
    %%% END SAVE FILES %%%
end


%%% GENERATE TEST FIGURE %%%
    figure
    hold on
    title(sprintf('First Four Dispersion Curves, Frequency (Hz) \n as a Function of Wavevector k (m^{-1}) Magnitude'))
    scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,1),"DisplayName",'DC1')
    scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,2),"DisplayName",'DC2')
    scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,3),"DisplayName",'DC3')
    scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,4),"DisplayName",'DC4')
    xlabel('k magnitude (m^{-1})')
    ylabel('Frequency (Hz)')
    legend('Location', 'southeast');
    hold off
%%% END GENERATE TEST FIGURE %%%
