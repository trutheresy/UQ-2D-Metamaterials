clc
clear all

save_outputs = true;
% custom_E = false;
% custom_Rho = false;
% custom_PR = false;
% error_geometries = true;

%%% DEFAULT PARAMETERS %%%
n_geometries = 1;
n_materials = 1;
E_soft = 200e6*ones(1,n_materials);
E_hard = 200e9*ones(1,n_materials);
rho_soft = 1e3*ones(1,n_materials);
rho_hard = 8e3*ones(1,n_materials);
poisson_soft = 0*ones(1,n_materials);
poisson_hard = 0.5*ones(1,n_materials);

%%% ALTERED PARAMETERS %%%
input_geometries = importdata('error_geometries.mat');
%input_data = pd_5_inputs;
%input_data_size = size(input_data);
%n_materials = input_data_size(2);
n_geometries = length(input_geometries)
%E_soft = input_data(1,:);
%E_hard = input_data(1,:);
%rho_soft = input_data(1,:); %change back to correct rows later
%rho_hard = input_data(1,:);
%poisson_soft = input_data(1,:);
%poisson_hard = input_data(1,:);

bg_size = zeros(n_materials, 1);
bg_top = zeros(n_materials, 1);
bg_bottom = zeros(n_materials, 1);

for g = 1:n_geometries
    for i = 1:n_materials
        [EIGENVALUE_DATA, WAVEVECTOR_DATA] = dispersion_data_custom_mat_geo_param(E_soft(i), E_hard(i), rho_soft(i), rho_hard(i), poisson_soft(i), poisson_hard(i), input_geometries(g,:,:));
        
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

    % Change these names before each run!
    if save_outputs
        rootFolder = 'E:\Research\Projects\UQ 2D Metamaterials';
        cd(rootFolder);
        dt = datetime('now'); % get current date and time
        dt.Format = 'yyyy-MM-dd'; % set desired format
        dt_str = char(dt); % convert to string
        mkdir(dt_str)
        cd(dt_str)
        outputFolder = [num2str(g) '_of_' num2str(n_geometries)];
        save(['bg_size_uniform_q_pd' num2str(length(bg_size)-1)], 'bg_size');
        save(['bg_bottom_uniform_q_pd' num2str(length(bg_bottom)-1)], 'bg_bottom');
        save(['bg_top_uniform_q_pd' num2str(length(bg_top)-1)], 'bg_top');
        %save(['E_soft_uniform_q_pd' num2str(length(E_soft))], 'E_soft');
        %save(['E_hard_uniform_q_pd' num2str(length(E_hard)-1)], 'E_hard');
        %save(['rho_soft_uniform_q_pd' num2str(length(rho_soft)-1)], 'rho_soft');
        %save(['rho_hard_uniform_q_pd' num2str(length(rho_hard)-1)], 'rho_hard');
        %save(['pr_soft_uniform_q_pd' num2str(length(poisson_soft)-1)], 'poisson_soft');
        %save(['pr_hard_uniform_q_pd' num2str(length(poisson_hard)-1)], 'poisson_hard');    
    end
    
    % figure
    % hold on
    % title(sprintf('First Four Dispersion Curves, Frequency (Hz) \n as a Function of Wavevector k (m^{-1}) Magnitude'))
    % scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,1),"DisplayName",'DC1')
    % scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,2),"DisplayName",'DC2')
    % scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,3),"DisplayName",'DC3')
    % scatter(sqrt(sum([WAVEVECTOR_DATA(:,1).^2, WAVEVECTOR_DATA(:,2).^2],2)), EIGENVALUE_DATA(:,4),"DisplayName",'DC4')
    % xlabel('k magnitude (m^{-1})')
    % ylabel('Frequency (Hz)')
    % legend('Location', 'southeast');
    % hold off

end
