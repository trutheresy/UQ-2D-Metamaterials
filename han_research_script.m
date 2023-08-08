%clc
%clear all
addpath('/Users/zhang/Documents/Duke/Research/2D-dispersion')
save_outputs = true;
% custom_E = false;
% custom_Rho = false;
% custom_PR = false;
% error_geometries = true;

%%% DEFAULT PARAMETERS %%%
n_materials = 1;
E_soft = 200e6*ones(1,n_materials);
E_hard = 200e9*ones(1,n_materials);
rho_soft = 1e3*ones(1,n_materials);
rho_hard = 8e3*ones(1,n_materials);
poisson_soft = 0*ones(1,n_materials);
poisson_hard = 0.5*ones(1,n_materials);

%%% ALTERED PARAMETERS %%%
% input_data = mc_10000_inputs;
% input_data_size = size(input_data);
% n_materials = input_data_size(2)
% E_soft = input_data(1,:);
% E_hard = input_data(2,:);
% rho_soft = input_data(3,:); %change back to correct rows later
% rho_hard = input_data(4,:);
% poisson_soft = input_data(5,:);
% poisson_hard = input_data(6,:);
geometries = error_geometries;
geometries_size = size(geometries);
n_geometries = geometries_size(1)

bg_size = zeros(n_materials, 1);
bg_top = zeros(n_materials, 1);
bg_bottom = zeros(n_materials, 1);

tic; % Start timer

for g = 1:n_geometries
    geometry = geometries(g,:,:);
    %disp(['geometry:' num2str(size(geometry))])
    for i = 1:n_materials
        disp(['g:' num2str(g) ', i:' num2str(i)])
        %[EIGENVALUE_DATA, WAVEVECTOR_DATA] = dispersion_data_custom_all_param(E_soft(i), E_hard(i), rho_soft(i), rho_hard(i), poisson_soft(i), poisson_hard(i));
        [EIGENVALUE_DATA, WAVEVECTOR_DATA] = dispersion_data_custom_mat_geo_param(E_soft(i), E_hard(i), rho_soft(i), rho_hard(i), poisson_soft(i), poisson_hard(i), geometry);
    
        bg_bottom(g*i) = max(EIGENVALUE_DATA(:,3));
        bg_top(g*i) = min(EIGENVALUE_DATA(:,4));
        bg_size(g*i) = bg_top(g*i) - bg_bottom(g*i);
        if bg_size(g*i) > 0
            disp(['Bandgap present for input ' num2str(i)])
        else
            disp(['No bandgap for input ' num2str(i)])
        end
    end
end
toc; % End timer and display elapsed time
elapsedTime = toc;

% Change these names before each run!
if save_outputs
    outputFolder = 'E:/Research/Projects/UQ 2D Metamaterials';
    %outputFolder = 'C:\Users\zhang\Documents\Duke\Research\UQ-2D-Metamaterials';
    cd(outputFolder);
    save(['bg_size_uniform_80p_5%_n' num2str(length(bg_size))], 'bg_size');
    save(['bg_bottom_uniform_80p_5%_n' num2str(length(bg_bottom))], 'bg_bottom');
    save(['bg_top_uniform_80p_5%_n' num2str(length(bg_top))], 'bg_top');
    save(['E_soft_uniform_80p_5%_n' num2str(length(E_soft))], 'E_soft');
    save(['E_hard_uniform_80p_5%_n' num2str(length(E_hard))], 'E_hard');
    save(['rho_soft_uniform_80p_5%_n' num2str(length(rho_soft))], 'rho_soft');
    save(['rho_hard_uniform_80p_5%_n' num2str(length(rho_hard))], 'rho_hard');
    save(['pr_soft_uniform_80p_5%_n' num2str(length(poisson_soft))], 'poisson_soft');
    save(['pr_hard_uniform_80p_5%_n' num2str(length(poisson_hard))], 'poisson_hard');
    save(['elapsed_time_80p_n' num2str(n_materials*n_geometries)], 'elapsedTime')
end

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
