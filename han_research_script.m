%clc
%clear all

save_outputs = true;
% custom_E = false;
% custom_Rho = false;
% custom_PR = false;

%%% ALTERED PARAMETERS %%%
input_data = mc_10000_inputs;
input_data_size = size(input_data);
n_samples = input_data_size(2);
%E_soft = input_data(1,:);
%E_hard = input_data(1,:);
%rho_soft = input_data(1,:); %change back to correct rows later
%rho_hard = input_data(1,:);
poisson_soft = input_data(1,:);
%poisson_hard = input_data(1,:);

%%% DEFAULT PARAMETERS %%%
% n_samples = 3;
E_soft = 200e6*ones(1,n_samples);
E_hard = 200e9*ones(1,n_samples);
rho_soft = 1e3*ones(1,n_samples);
rho_hard = 8e3*ones(1,n_samples);
%poisson_soft = 0*ones(1,n_samples);
poisson_hard = 0.5*ones(1,n_samples);

bg_size = zeros(n_samples, 1);
bg_top = zeros(n_samples, 1);
bg_bottom = zeros(n_samples, 1);

for i = 1:n_samples
    [EIGENVALUE_DATA, WAVEVECTOR_DATA] = dispersion_data_custom_all_param(E_soft(i), E_hard(i), rho_soft(i), rho_hard(i), poisson_soft(i), poisson_hard(i));

    bg_bottom(i) = max(EIGENVALUE_DATA(:,3));
    bg_top(i) = min(EIGENVALUE_DATA(:,4));
    bg_size(i) = bg_top(i) - bg_bottom(i);
    if bg_size(i) > 0
        disp(['Bandgap present for input ' num2str(i)])
    else
        disp(['No bandgap for input ' num2str(i)])
    end
end

% Change these names before each run!
if save_outputs
    outputFolder = 'E:\Research\Projects\UQ 2D Metamaterials';
    cd(outputFolder);
    save(['bg_size_uniform_mc' num2str(length(bg_size))], 'bg_size');
    save(['bg_bottom_uniform_mc' num2str(length(bg_bottom))], 'bg_bottom');
    save(['bg_top_uniform_mc' num2str(length(bg_top))], 'bg_top');
    % save(['E_soft_uniform_mc' num2str(length(E_soft))], 'E_soft');
    % save(['E_hard_uniform_mc' num2str(length(E_hard))], 'E_hard');
    % save(['rho_soft_uniform_mc_' num2str(length(rho_soft))], 'rho_soft');
    % save(['rho_hard_uniform_mc_' num2str(length(rho_hard))], 'rho_hard');
    save(['pr_soft_uniform_mc_' num2str(length(poisson_soft))], 'poisson_soft');
    %save(['pr_hard_uniform_mc_' num2str(length(poisson_hard))], 'poisson_hard');    
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
