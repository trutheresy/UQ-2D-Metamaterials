clc
clear all
close all


%% Fig 1a
% Define the matrix
matrix = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 1, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 1, 1, 1, 1, 1, 1, 1, 1, 0;
    0, 1, 0, 1, 0, 0, 1, 0, 1, 0;
    0, 1, 0, 1, 0, 0, 1, 0, 1, 0;
    0, 1, 1, 1, 1, 1, 1, 1, 1, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 1, 1, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
];

% Display the matrix as an image with a grayscale colormap
imshow(1 - matrix, 'Colormap', gray, 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(matrix);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);

fig = gcf; % Get the current figure handle
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; % [left, bottom, width, height] in inches
saveas(fig, '1st_geo.fig');

% Fig 1b
% Define the matrix
matrix = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 1, 1, 1, 1, 0, 1, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 1, 1, 1, 0, 0, 1, 1, 1, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 1, 1, 1, 0, 0, 1, 1, 1, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 1, 0, 1, 1, 1, 1, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
];
figure
% Display the matrix as an image with a grayscale colormap
imshow(1 - matrix, 'Colormap', gray, 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(matrix);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);

fig = gcf; % Get the current figure handle
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; % [left, bottom, width, height] in inches
saveas(fig, '2nd_geo.fig');

close all
%% Fig2

% Define the matrix

matrix = load('DATASETS/gamma beta 6+1 inputs quadrature rule sparse study/fp_matrices_pd_1_geos_sparse.mat');

slice = zeros (40);

for i = 1:9
    slice (:,:) = matrix.pd_1_geos(i,:,:);

    figure

    % Display the matrix as an image with a grayscale colormap
    imshow(1 - slice, 'Colormap', gray, 'InitialMagnification', 'fit');

    % Determine the size of the matrix
    [m, n] = size(slice);

    % Add a black bounding box with aligned edges
    rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);

    fig = gcf; % Get the current figure handle
    fig.Units = 'inches';
    fig.Position = [0, 0, 3, 3]; % [left, bottom, width, height] in inches
        % Convert index to alphabet
    alphabet_index = char('a' + i - 1);
    
    % Save figure with alphabet naming convention
    saveas(fig, ['defect_palette_40_', alphabet_index, '.fig']);
    
end

close all

%% Fig 3

% Define the matrix sizes and input data
matrix_sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
bg_bottoms = [1170.64993236, 1155.89183287, 1151.36607409, 1149.25125262, 1148.04873935, 1147.28184608, 1146.75424695, 1146.37122131, 1146.08174241, 1145.85602108];
bg_tops = [2635.41277912, 2554.99450755, 2536.34090556, 2528.69197521, 2524.65647776, 2522.22644528, 2520.62577626, 2519.5030904, 2518.67818653, 2518.04999572];
bg_sizes = [1464.76284676, 1399.10267468, 1384.97483147, 1379.44072259, 1376.60773841, 1374.9445992, 1373.87152931, 1373.13186909, 1372.59644412, 1372.19397465];


figure;
% Plot bg_sizes against matrix_sizes in the same plot with red color and smaller marker size
plot(matrix_sizes, bg_sizes, 'o-', 'Color', 'red', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Image resolution (px)', 'FontSize', 12);
ylabel('Bandgap size (Hz)', 'FontSize', 12);
fig = gcf; % Get the current figure handle
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; % [left, bottom, width, height] in inches
saveas(fig, 'P_FEM_convergence_afo_resolution_a.fig');

figure;
% Plot bg_center against matrix_sizes in the same plot with blue color and smaller marker size
plot(matrix_sizes, (bg_tops+bg_bottoms)./2, 'o-', 'Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Image resolution (px)', 'FontSize', 12);
ylabel('Bandgap center (Hz)', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_FEM_convergence_afo_resolution_b.fig');

% figure;
% % Plot bg_tops against matrix_sizes in the same plot with green color and smaller marker size
% plot(matrix_sizes, bg_tops, 'o-', 'Color', 'green', 'LineWidth', 1.5, 'MarkerSize', 4);
% xlabel('Image Resolution (px)', 'FontSize', 12);
% ylabel('Bandgap top frequency (Hz)', 'FontSize', 12);
% 
% figure;
% % Plot bg_bottoms against matrix_sizes in the same plot with blue color and smaller marker size
% plot(matrix_sizes, bg_bottoms, 'o-', 'Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 4);
% xlabel('Image Resolution (px)', 'FontSize', 12);
% ylabel('Bandgap bottom frequency (Hz)', 'FontSize', 12);

close all

%% Fig 4

% Assuming bg_sizes_n100_fp and bg_centers_n100_fp are MATLAB arrays

figure
p10 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_10p_5%_n100.mat');
p20 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_20p_5%_n100.mat');
p30 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_30p_5%_n100.mat');
p40 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_40p_5%_n100.mat');
p50 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_50p_5%_n100.mat');
p60 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_60p_5%_n100.mat');
p70 = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_size_uniform_70p_5%_n100.mat');

p10b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_10p_5%_n100.mat');
p20b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_20p_5%_n100.mat');
p30b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_30p_5%_n100.mat');
p40b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_40p_5%_n100.mat');
p50b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_50p_5%_n100.mat');
p60b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_60p_5%_n100.mat');
p70b = load('DATASETS/geometry_defect_effect_asymtote_study_fp_gaussian/bg_bottom_uniform_70p_5%_n100.mat');


histogram(p10.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'DisplayName', '10\times10 pixels');
hold on;
histogram(p20.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'DisplayName', '20\times20 pixels');
histogram(p30.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0 1 0], 'EdgeColor', 'none', 'DisplayName', '30\times30 pixels');
histogram(p40.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0.5 0 0.5], 'EdgeColor', 'none', 'DisplayName', '40\times40 pixels');
hold off;
xlabel('Bandgap size (Hz)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.1, 3.2]; 
saveas(fig, 'P_geo_defect_fp_convergence_40p-70p_a.fig');
% title('100 MC Samples'' Bandgap Sizes');

figure
histogram(p10b.bg_bottom+p10.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0 0 1], 'EdgeColor', 'none','DisplayName', '10\times10 pixels');
hold on;
histogram(p20b.bg_bottom+p20.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'DisplayName', '20\times20 pixels');
histogram(p30b.bg_bottom+p30.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0 1 0], 'EdgeColor', 'none', 'DisplayName', '30\times30 pixels');
histogram(p40b.bg_bottom+p40.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0.5 0 0.5], 'EdgeColor', 'none', 'DisplayName', '40\times40 pixels');
hold off;
xlabel('Bandgap center (Hz)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.1, 3.2];  
saveas(fig, 'P_geo_defect_fp_convergence_40p-70p_b.fig');
% title('100 MC Samples'' Bandgap Locations');

% suptitle('Histograms of 100 MC Samples'' Bandgap Sizes and Center Locations');

% Assuming bg_sizes_n100_fp and bg_centers_n100_fp are MATLAB arrays

figure
histogram(p40.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0.5 0 0.5], 'EdgeColor', 'none', 'DisplayName', '40\times40 pixels');
hold on;
histogram(p50.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none', 'DisplayName', '50\times50 pixels');
histogram(p60.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0 1 1], 'EdgeColor', 'none', 'DisplayName', '60\times60 pixels');
histogram(p70.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'DisplayName', '70\times70 pixels');
hold off;
xlabel('Bandgap size (Hz)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.1, 3.2];  
saveas(fig, 'P_geo_defect_fp_convergence_40p-70p_c.fig');
% title('100 MC Samples'' Bandgap Sizes');

figure
histogram(p40b.bg_bottom+p40.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0.5 0 0.5], 'EdgeColor', 'none', 'DisplayName', '40\times40 pixels');
hold on;
histogram(p50b.bg_bottom+p50.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none', 'DisplayName', '50\times50 pixels');
histogram(p60b.bg_bottom+p60.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [0 1 1], 'EdgeColor', 'none', 'DisplayName', '60\times60 pixels');
histogram(p70b.bg_bottom+p70.bg_size./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'DisplayName', '70\times70 pixels');
hold off;
xlabel('Bandgap center (Hz)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.1, 3.2]; 
saveas(fig, 'P_geo_defect_fp_convergence_40p-70p_d.fig');
% title('100 MC Samples'' Bandgap Locations');
% 
% suptitle('Histograms of 100 MC Samples'' (FP) Computed Bandgap Size and Center Locations for Different Resolutions (40-70 px)');

close all

%% Fig 5

% Assuming rho_soft_dist, rho_hard_dist, K_soft_dist, K_hard_dist, G_soft_dist,
% G_hard_dist, geo_fp_dist are MATLAB probability distributions

data = load('DATASETS/gamma beta 6+1 inputs mc study/joint_dist_mat_geo_mc_10000.mat');

figure
histogram(data.DATASETS/mc_10000_inputs(1,:), 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('\rho_{soft} (kg/m^3)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_a.fig');

figure
histogram(data.DATASETS/mc_10000_inputs(2,:), 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('\rho_{stiff} (kg/m^3)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_b.fig');

figure
histogram(data.DATASETS/mc_10000_inputs(3,:)./10^6, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('K_{soft} (MPa)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_c.fig');

figure
histogram(data.DATASETS/mc_10000_inputs(4,:)./10^9, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('K_{stiff} (GPa)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_d.fig');

figure
histogram(data.DATASETS/mc_10000_inputs(5,:)./10^6, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('G_{soft} (MPa)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_e.fig');

figure
histogram(data.DATASETS/mc_10000_inputs(6,:)./10^9, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('G_{stiff} (GPa)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_f.fig');

figure
histogram(data.DATASETS/mc_10000_inputs(7,:), 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('FP', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_gamma_input_mc_10000_g.fig');

close all

%% Fig 6

matrix = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 1, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 1, 1, 1, 1, 1, 1, 1, 1, 0;
    0, 1, 0, 1, 0, 0, 1, 0, 1, 0;
    0, 1, 0, 1, 0, 0, 1, 0, 1, 0;
    0, 1, 1, 1, 1, 1, 1, 1, 1, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 1, 1, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
];

% Display the matrix as an image with a grayscale colormap
imshow(1 - matrix, 'Colormap', gray, 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(matrix);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; 
saveas(fig, 'defect_trio_1st_geo_a.fig');

figure

EdgePixels = [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.];

imshow(1 - EdgePixels, 'Colormap',  [0 0 1; 1 1 1], 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(EdgePixels);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; 
saveas(fig, 'defect_trio_1st_geo_b.fig');

figure

DefectiveGeometry = [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 1. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 1. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0.
0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.];

% Display the matrix as an image with a grayscale colormap
imshow(1 - DefectiveGeometry, 'Colormap', gray, 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(DefectiveGeometry);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; 
saveas(fig, 'defect_trio_1st_geo_c.fig');

close all

%% Fig 7

size100 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_size_gamma_7d_fp_5%_n100.mat');
size1000 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_size_gamma_7d_fp_5%_n1000.mat');
size10000 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_size_gamma_7d_fp_5%_n10000.mat');

top100 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_top_gamma_7d_fp_5%_n100.mat');
top1000 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_top_gamma_7d_fp_5%_n1000.mat');
top10000 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_top_gamma_7d_fp_5%_n10000.mat');

bottom100 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_bottom_gamma_7d_fp_5%_n100.mat');
bottom1000 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_bottom_gamma_7d_fp_5%_n1000.mat');
bottom10000 = load('DATASETS/gamma beta 6+1 inputs mc study/bg_bottom_gamma_7d_fp_5%_n10000.mat');

% Plot for 100 MC Samples
figure;
hold on;
% histogram(top100.bg_top, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
% histogram(bottom100.bg_bottom, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram(size100.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram((top100.bg_top+bottom100.bg_bottom)./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
xlabel("Model output (Hz)", 'FontSize', 11);
ylabel("Probability density", 'FontSize', 11);
legend('Bandgap size', 'Bandgap center', 'Location', 'northwest', 'FontSize', 11);
% title('Histograms of 100 MC Samples');
ylim([0 12e-3])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 2.5];  
saveas(fig, 'bgtbs_trio_hist_7d_gamma_a.fig');

% Plot for 1000 MC Samples
figure;
hold on;
% histogram(top1000.bg_top, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
% histogram(bottom1000.bg_bottom, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram(size1000.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram((top1000.bg_top+bottom1000.bg_bottom)./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
xlabel("Model output (Hz)", 'FontSize', 11);
ylabel("Probability density", 'FontSize', 11);
legend('Bandgap size', 'Bandgap center', 'Location', 'northwest', 'FontSize', 11);
% title('Histograms of 1000 MC Samples');
ylim([0 12e-3])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 2.5]; 
saveas(fig, 'bgtbs_trio_hist_7d_gamma_b.fig');

% Plot for 10000 MC Samples
figure;
hold on;
% histogram(top10000.bg_top, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
% histogram(bottom10000.bg_bottom, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram(size10000.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram((top10000.bg_top+bottom10000.bg_bottom)./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
xlabel("Model output (Hz)", 'FontSize', 11);
ylabel("Probability density", 'FontSize', 11);
legend('Bandgap size', 'Bandgap center', 'Location', 'northwest', 'FontSize', 11);
% title('Histograms of 10000 MC Samples');
ylim([0 12e-3])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 2.5]; 
saveas(fig, 'bgtbs_trio_hist_7d_gamma_c.fig');

close all

%% Fig 8

surrogate_output_mc_pd1 = load('DATASETS/gamma beta 6+1 inputs mc study/surrogate_outputs_bgs_mc_100_pd_1.mat');
[mc_pd1_y,mc_pd1_x] = ksdensity(surrogate_output_mc_pd1.pd_1_outputs);
surrogate_output_mc_pd2 = load('DATASETS/gamma beta 6+1 inputs mc study/surrogate_outputs_bgs_mc_100_pd_2.mat');
[mc_pd2_y,mc_pd2_x] = ksdensity(surrogate_output_mc_pd2.pd_2_outputs);

figure;
histogram(size100.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('100 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 100 MC Samples');
ylim([0 8e-3])
xlim([900 2100])
hold on
plot(mc_pd1_x, mc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_pd2_x, mc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'best', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.6, 3.2]; 
saveas(fig, 'bgs_trio_hist_mc_fit_7d_n100_gamma_a.fig');

figure;
histogram(size1000.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('1000 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 1000 MC Samples');
ylim([0 8e-3])
xlim([900 2100])
hold on
plot(mc_pd1_x, mc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_pd2_x, mc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.6, 3.2]; 
saveas(fig, 'bgs_trio_hist_mc_fit_7d_n100_gamma_b.fig');

figure;
histogram(size10000.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('10000 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 10000 MC Samples');
ylim([0 8e-3])
xlim([900 2100])
hold on
plot(mc_pd1_x, mc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_pd2_x, mc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.6, 3.2];  
saveas(fig, 'bgs_trio_hist_mc_fit_7d_n100_gamma_c.fig');


surrogate_bgc_output_mc_pd1 = load('DATASETS/gamma beta 6+1 inputs mc study/surrogate_outputs_bgc_mc_100_pd_1.mat');
[mc_bgc_pd1_y,mc_bgc_pd1_x] = ksdensity(surrogate_bgc_output_mc_pd1.pd_1_outputs);
surrogate_bgc_output_mc_pd2 = load('DATASETS/gamma beta 6+1 inputs mc study/surrogate_outputs_bgc_mc_100_pd_2.mat');
[mc_bgc_pd2_y,mc_bgc_pd2_x] = ksdensity(surrogate_bgc_output_mc_pd2.pd_2_outputs);
 
figure;
histogram((top100.bg_top+bottom100.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('100 MC samples', 'Location', 'northeast');
ylim([0 13e-3])
xlim([1500 2300])
hold on
plot(mc_bgc_pd1_x, mc_bgc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_bgc_pd2_x, mc_bgc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'best', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.6, 3.2]; 
saveas(fig, 'bgc_trio_hist_mc_fit_7d_n100_gamma_a.fig');

figure;
histogram((top1000.bg_top+bottom1000.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('1000 MC samples', 'Location', 'northeast');
ylim([0 13e-3])
xlim([1500 2300])
hold on
plot(mc_bgc_pd1_x, mc_bgc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_bgc_pd2_x, mc_bgc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.6, 3.2];  
saveas(fig, 'bgc_trio_hist_mc_fit_7d_n100_gamma_b.fig');

figure;
histogram((top10000.bg_top+bottom10000.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('10000 MC samples', 'Location', 'northeast');
ylim([0 13e-3])
xlim([1500 2300])
hold on
plot(mc_bgc_pd1_x, mc_bgc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_bgc_pd2_x, mc_bgc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.6, 3.2]; 
saveas(fig, 'bgc_trio_hist_mc_fit_7d_n100_gamma_c.fig');

% Ensemble

% bandgap size
fig = figure;
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Plotting for N = 100
nexttile;
histogram(size100.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
ylim([0 8e-3])
xlim([900 2100])
hold on
plot(mc_pd1_x, mc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_pd2_x, mc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);

% Plotting for N = 1000
nexttile;
histogram(size1000.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylim([0 8e-3])
xlim([900 2100])
hold on
plot(mc_pd1_x, mc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_pd2_x, mc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

% Plotting for N = 10000
nexttile;
histogram(size10000.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylim([0 8e-3])
xlim([900 2100])
hold on
plot(mc_pd1_x, mc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_pd2_x, mc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

fig.Units = 'inches';
fig.Position = [0, 0, 10.4, 3]; 
saveas(fig, 'bgs_trio_hist_mc_fit_7d_n100_gamma.fig');


% bandgap center frequency
fig = figure;
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% Plotting for N = 100
nexttile;
histogram((top100.bg_top + bottom100.bg_bottom) ./ 2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
ylim([0 13e-3])
xlim([1500 2300])
hold on
plot(mc_bgc_pd1_x, mc_bgc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_bgc_pd2_x, mc_bgc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);

% Plotting for N = 1000
nexttile;
histogram((top1000.bg_top + bottom1000.bg_bottom) ./ 2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 12);
ylim([0 13e-3])
xlim([1500 2300])
hold on
plot(mc_bgc_pd1_x, mc_bgc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_bgc_pd2_x, mc_bgc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

% Plotting for N = 10000
nexttile;
histogram((top10000.bg_top + bottom10000.bg_bottom) ./ 2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 12);
ylim([0 13e-3])
xlim([1500 2300])
hold on
plot(mc_bgc_pd1_x, mc_bgc_pd1_y, 'r', 'LineWidth', 1.5)
plot(mc_bgc_pd2_x, mc_bgc_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=100', 'PD=2, N=100', 'Location', 'northeast', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

fig.Units = 'inches';
fig.Position = [0, 0, 10.4, 3]; 
saveas(fig, 'bgc_trio_hist_mc_fit_7d_n100_gamma.fig');

close all
%% Fig 9

% Define the colormap range and resolution
mapResolution = 100;
colormapRange = linspace(0, 1, mapResolution);

% Define the starting and ending colors for each transition
colors = [
    1, 1, 1;  % White
    1, 1, 0;  % Yellow
    1, 0.5, 0;  % Orange
    1, 0, 0;  % Red
    0.6350, 0.0780, 0.1840  % Brown
];

% Initialize the custom colormap
custom_colormap = zeros(mapResolution, 3);

% Create a custom colormap transitioning through the specified colors
for c = 1:size(colors, 1) - 1
    segmentSize = mapResolution / (size(colors, 1) - 1);
    startIdx = round((c - 1) * segmentSize) + 1;
    endIdx = round(c * segmentSize);
    for i = 1:3
        custom_colormap(startIdx:endIdx, i) = linspace(colors(c, i), colors(c+1, i), endIdx - startIdx + 1);
    end
end

% Create a 2D histogram for 10000 MC Samples
figure;
histogram2(size10000.bg_size, (top10000.bg_top+bottom10000.bg_bottom)./2, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Bandgap center (Hz)", 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
ylim([1400 2300]);
xlim([900 2000]);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3];
grid on
set(gca, 'Layer', 'top')
box on
saveas(fig, '2d_hist_7d_input_pd1_quad_gamma_a.fig');

% Create a 2D histogram for Quadrature Surrogate
figure;

surrogate_outputs_a_bgs = load('DATASETS/gamma beta 6+1 inputs quadrature rule study/surrogate_outputs_2d_bgs_qr_pd_1.mat');
surrogate_outputs_a_bgc = load('DATASETS/gamma beta 6+1 inputs quadrature rule study/surrogate_outputs_2d_bgc_qr_pd_1.mat');

histogram2(surrogate_outputs_a_bgs.pd_1_outputs, surrogate_outputs_a_bgc.pd_1_outputs, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Bandgap center (Hz)", 'FontSize', 12);

colormap(custom_colormap);
colorbar;
clim([0 100]);
ylim([1400 2300]);
xlim([900 2000]);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3];
grid on
set(gca, 'Layer', 'top')
box on
saveas(fig, '2d_hist_7d_input_pd1_quad_gamma_b.fig');


% Create a 2D histogram for MC Surrogate
figure;

surrogate_outputs_b_bgs = load('DATASETS/gamma beta 6+1 inputs mc study/surrogate_outputs_2d_bgs_mc_100_pd_1.mat');
surrogate_outputs_b_bgc = load('DATASETS/gamma beta 6+1 inputs mc study/surrogate_outputs_2d_bgc_mc_100_pd_1.mat');

histogram2(surrogate_outputs_b_bgs.pd_1_outputs, surrogate_outputs_b_bgc.pd_1_outputs, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Bandgap center (Hz)", 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
ylim([1400 2300]);
xlim([900 2000]);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3];
grid on
set(gca, 'Layer', 'top')
box on
saveas(fig, '2d_hist_7d_input_pd1_quad_gamma_c.fig');

% Create a 2D histogram for SG Surrogate
figure;

surrogate_outputs_c_bgs = load('DATASETS/gamma beta 6+1 inputs quadrature rule sparse study/surrogate_outputs_2d_bgs_sg_pd_1.mat');
surrogate_outputs_c_bgc = load('DATASETS/gamma beta 6+1 inputs quadrature rule sparse study/surrogate_outputs_2d_bgc_sg_pd_1.mat');
histogram2(surrogate_outputs_c_bgs.pd_1_outputs, surrogate_outputs_c_bgc.pd_1_outputs, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Bandgap center (Hz)", 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
ylim([1400 2300]);
xlim([900 2000]);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3];
grid on
set(gca, 'Layer', 'top')
box on
saveas(fig, '2d_hist_7d_input_pd1_quad_gamma_d.fig');

close all
%% Fig 10

% Assuming rho_soft_dist, rho_hard_dist, K_soft_dist, K_hard_dist, G_soft_dist,
% G_hard_dist, geo_fp_dist are MATLAB probability distributions

RhoSoft = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/rho_soft_gaussian_7d_fp_5%_n10000.mat');
RhoHard = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/rho_hard_gaussian_7d_fp_5%_n10000.mat');
ESoft = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/E_soft_gaussian_7d_fp_5%_n10000.mat');
EHard = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/E_hard_gaussian_7d_fp_5%_n10000.mat');
PrSoft = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/pr_soft_gaussian_7d_fp_5%_n10000.mat');
PrHard = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/pr_hard_gaussian_7d_fp_5%_n10000.mat');
FPGaussian = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/geo_fp_dist_trunc_mc_10000.mat');

figure
histogram(RhoSoft.rho_soft, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('\rho_{soft} (kg/m^3)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_a.fig');

figure
histogram(RhoHard.rho_hard, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('\rho_{stiff} (kg/m^3)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_b.fig');

figure
histogram(ESoft.E_soft./10^6, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('E_{soft} (MPa)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_c.fig');

figure
histogram(EHard.E_hard./10^9, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('E_{stiff} (GPa)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_d.fig');

figure
histogram(PrSoft.pr_soft, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('\nu_{soft}', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_e.fig');

figure
histogram(PrHard.pr_hard, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('\nu_{stiff}', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_f.fig');

figure
histogram(FPGaussian.DATASETS/mc_10000_geos, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel('FP', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 2.5, 2.3]; 
saveas(fig, 'P_7d_input_mc_10000_2nd_geo_g.fig');

close all

%% Fig 11

matrix = [
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;
    0., 1., 0., 1., 1., 1., 1., 0., 1., 0.;
    0., 0., 0., 1., 0., 0., 1., 0., 0., 0.;
    0., 1., 1., 1., 0., 0., 1., 1., 1., 0.;
    0., 1., 0., 0., 0., 0., 0., 0., 1., 0.;
    0., 1., 0., 0., 0., 0., 0., 0., 1., 0.;
    0., 1., 1., 1., 0., 0., 1., 1., 1., 0.;
    0., 0., 0., 1., 0., 0., 1., 0., 0., 0.;
    0., 1., 0., 1., 1., 1., 1., 0., 1., 0.;
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
];

% Display the matrix as an image with a grayscale colormap
imshow(1 - matrix, 'Colormap', gray, 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(matrix);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);

fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; 
saveas(fig, 'defect_trio_2nd_geo_a.fig');

figure

EdgePixels = [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
 0. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
 0. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 1. 1. 0. 0. 0.;
 0. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 0.;
 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.];

imshow(1 - EdgePixels, 'Colormap',  [0 0 1; 1 1 1], 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(EdgePixels);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);

fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; 
saveas(fig, 'defect_trio_2nd_geo_b.fig');

figure

DefectiveGeometry = [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 0. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 0. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0. 0. 1. 1. 1. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.;
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.];
% Display the matrix as an image with a grayscale colormap
imshow(1 - DefectiveGeometry, 'Colormap', gray, 'InitialMagnification', 'fit');

% Determine the size of the matrix
[m, n] = size(DefectiveGeometry);

% Add a black bounding box with aligned edges
rectangle('Position', [0.5, 0.5, n, m], 'EdgeColor', 'black', 'LineWidth', 2);

fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 3]; 
saveas(fig, 'defect_trio_2nd_geo_c.fig');

close all

%% Fig 12

size100_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_size_gaussian_7d_fp_5%_n100.mat');
size1000_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_size_gaussian_7d_fp_5%_n1000.mat');
size10000_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_size_gaussian_7d_fp_5%_n10000.mat');

top100_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_top_gaussian_7d_fp_5%_n100.mat');
top1000_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_top_gaussian_7d_fp_5%_n1000.mat');
top10000_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_top_gaussian_7d_fp_5%_n10000.mat');

bottom100_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_bottom_gaussian_7d_fp_5%_n100.mat');
bottom1000_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_bottom_gaussian_7d_fp_5%_n1000.mat');
bottom10000_2 = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/bg_bottom_gaussian_7d_fp_5%_n10000.mat');

% Plot for 100 MC Samples
figure;
hold on;

% histogram(top100_2.bg_top, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
% histogram(bottom100_2.bg_bottom, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram(size100_2.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram((top100_2.bg_top+bottom100_2.bg_bottom)./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
xlabel("Model output (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
legend('Bandgap size', 'Bandgap center', 'Location', 'northeast', 'FontSize', 12);
% legend('Bandgap top', 'Bandgap bottom', 'Bandgap size', 'Location', 'best', 'FontSize', 12);
% title('Histograms of 100 MC Samples');
box on
ylim([0 0.02])
xlim([700 1800])
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 2.5]; 
saveas(fig, 'bgtbs_trio_hist_7d_gaussian_2nd_geo_a.fig');

% Plot for 1000 MC Samples
figure;
hold on;
% histogram(top1000_2.bg_top, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
% histogram(bottom1000_2.bg_bottom, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram(size1000_2.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram((top1000_2.bg_top+bottom1000_2.bg_bottom)./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
xlabel("Model output (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
legend('Bandgap size', 'Bandgap center', 'Location', 'northeast', 'FontSize', 12);
% title('Histograms of 1000 MC Samples');
box on
ylim([0 0.02])
xlim([700 1800])
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 2.5]; 
saveas(fig, 'bgtbs_trio_hist_7d_gaussian_2nd_geo_b.fig');

% Plot for 10000 MC Samples
figure;
hold on;
% histogram(top10000_2.bg_top, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
% histogram(bottom10000_2.bg_bottom, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram(size10000_2.bg_size, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
histogram((top10000_2.bg_top+bottom10000_2.bg_bottom)./2, 50, 'Normalization', 'pdf', 'FaceAlpha', 1, 'EdgeAlpha', 0);
xlabel("Model output (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
legend('Bandgap size', 'Bandgap center', 'Location', 'northeast', 'FontSize', 12);
% title('Histograms of 10000 MC Samples');
box on
ylim([0 0.02])
xlim([700 1800])
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3, 2.5]; 
saveas(fig, 'bgtbs_trio_hist_7d_gaussian_2nd_geo_c.fig');

close all

%% Fig 13

surrogate_output_gaussian_mc_pd1 = load('DATASETS/gaussian 6+1 inputs quadrature rule 2nd geo study/surrogate_outputs_bgs_pd_1.mat');
[gaussian_pd1_y,gaussian_pd1_x] = ksdensity(surrogate_output_gaussian_mc_pd1.pd_1_outputs);
surrogate_output_gaussian_mc_pd2 = load('DATASETS/gaussian 6+1 inputs quadrature rule 2nd geo study/surrogate_outputs_bgs_pd_2.mat');
[gaussian_pd2_y,gaussian_pd2_x] = ksdensity(surrogate_output_gaussian_mc_pd2.pd_2_outputs);

figure;
histogram(size100_2.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% legend('100 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 100 MC Samples');
ylim([0 0.02])
xlim([675 975])
hold on
plot(gaussian_pd1_x, gaussian_pd1_y, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x, gaussian_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'best', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.2, 3];
saveas(fig, 'bgs_trio_hist_quad_fit_7d_2nd_geo_a.fig');

figure;
histogram(size1000_2.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
% legend('1000 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 1000 MC Samples');
ylim([0 0.02])
xlim([675 975])
hold on
plot(gaussian_pd1_x, gaussian_pd1_y, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x, gaussian_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northeast', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.2, 3];
saveas(fig, 'bgs_trio_hist_quad_fit_7d_2nd_geo_b.fig');

figure;
histogram(size10000_2.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
% legend('10000 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 10000 MC Samples');
ylim([0 0.02])
xlim([675 975])
hold on
plot(gaussian_pd1_x, gaussian_pd1_y, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x, gaussian_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northeast', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.2, 3];
saveas(fig, 'bgs_trio_hist_quad_fit_7d_2nd_geo_c.fig');


surrogate_output_gaussian_mc_pd1_bgc = load('DATASETS/gaussian 6+1 inputs quadrature rule 2nd geo study/surrogate_outputs_bgc_pd_1.mat');
[gaussian_pd1_y_bgc,gaussian_pd1_x_bgc] = ksdensity(surrogate_output_gaussian_mc_pd1_bgc.pd_1_outputs);
surrogate_output_gaussian_mc_pd2_bgc = load('DATASETS/gaussian 6+1 inputs quadrature rule 2nd geo study/surrogate_outputs_bgc_pd_2.mat');
[gaussian_pd2_y_bgc,gaussian_pd2_x_bgc] = ksdensity(surrogate_output_gaussian_mc_pd2_bgc.pd_2_outputs);

figure;
histogram((top100_2.bg_top+bottom100_2.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
% legend('100 MC samples', 'Location', 'northeast');
ylim([0 0.0125])
xlim([1300 1800])
hold on
plot(gaussian_pd1_x_bgc, gaussian_pd1_y_bgc, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x_bgc, gaussian_pd2_y_bgc, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'best', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.2, 3];
saveas(fig, 'bgc_trio_hist_quad_fit_7d_2nd_geo_a.fig');

figure;
histogram((top1000_2.bg_top+bottom1000_2.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
% legend('1000 MC samples', 'Location', 'northeast');
ylim([0 0.0125])
xlim([1300 1800])
hold on
plot(gaussian_pd1_x_bgc, gaussian_pd1_y_bgc, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x_bgc, gaussian_pd2_y_bgc, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northwest', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.2, 3];
saveas(fig, 'bgc_trio_hist_quad_fit_7d_2nd_geo_b.fig');

figure;
histogram((top10000_2.bg_top+bottom10000_2.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
% legend('10000 MC samples', 'Location', 'northeast');
ylim([0 0.0125])
xlim([1300 1800])
hold on
plot(gaussian_pd1_x_bgc, gaussian_pd1_y_bgc, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x_bgc, gaussian_pd2_y_bgc, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northwest', 'FontSize', 10);
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4.2, 3]; 
saveas(fig, 'bgc_trio_hist_quad_fit_7d_2nd_geo_c.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bandgap size
fig = figure;
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Plotting for N = 100
nexttile;
histogram(size100_2.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
% legend('100 MC Samples', 'Location', 'northeast');
% title('MC Regression Overlaid on 100 MC Samples');
ylim([0 0.02])
xlim([675 975])
hold on
plot(gaussian_pd1_x, gaussian_pd1_y, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x, gaussian_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'best', 'FontSize', 10);

% Plotting for N = 1000
nexttile;
histogram(size1000_2.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 10);
ylim([0 0.02])
xlim([675 975])
hold on
plot(gaussian_pd1_x, gaussian_pd1_y, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x, gaussian_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northeast', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

% Plotting for N = 10000
nexttile;
histogram(size10000_2.bg_size, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap size (Hz)", 'FontSize', 10);
ylim([0 0.02])
xlim([675 975])
hold on
plot(gaussian_pd1_x, gaussian_pd1_y, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x, gaussian_pd2_y, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northeast', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

fig.Units = 'inches';
fig.Position = [0, 0, 11.5, 3]; 
saveas(fig, 'bgs_trio_hist_quad_fit_7d_2nd_geo.fig');

%%%%%

% bandgap center frequency
fig = figure;
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% Plotting for N = 100
nexttile;
histogram((top100_2.bg_top+bottom100_2.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 10);
ylabel("Probability density", 'FontSize', 10);
ylim([0 0.0125])
xlim([1300 1800])
hold on
plot(gaussian_pd1_x_bgc, gaussian_pd1_y_bgc, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x_bgc, gaussian_pd2_y_bgc, 'k--', 'LineWidth', 1.5)
legend('MC, N=100', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'best', 'FontSize', 10);

% Plotting for N = 1000
nexttile;
histogram((top1000_2.bg_top+bottom1000_2.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 10);
ylim([0 0.0125])
xlim([1300 1800])
hold on
plot(gaussian_pd1_x_bgc, gaussian_pd1_y_bgc, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x_bgc, gaussian_pd2_y_bgc, 'k--', 'LineWidth', 1.5)
legend('MC, N=1000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northwest', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

% Plotting for N = 10000
nexttile;
histogram((top10000_2.bg_top+bottom10000_2.bg_bottom)./2, 'NumBins', 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
xlabel("Bandgap center (Hz)", 'FontSize', 10);
ylim([0 0.0125])
xlim([1300 1800])
hold on
plot(gaussian_pd1_x_bgc, gaussian_pd1_y_bgc, 'r', 'LineWidth', 1.5)
plot(gaussian_pd2_x_bgc, gaussian_pd2_y_bgc, 'k--', 'LineWidth', 1.5)
legend('MC, N=10000', 'PD=1, N=128', 'PD=2, N=2187', 'Location', 'northwest', 'FontSize', 10);
yticklabels({}); % Hide y-axis labels for this subplot

fig.Units = 'inches';
fig.Position = [0, 0, 11.5, 3]; 
saveas(fig, 'bgc_trio_hist_quad_fit_7d_2nd_geo.fig');

close all

%% Fig 14

% Create a 2D histogram for 10000 MC Samples
figure;
histogram2(size10000_2.bg_size, (top10000_2.bg_top+bottom10000_2.bg_bottom)./2, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel('Bandgap size (Hz)', 'FontSize', 12);
ylabel('Bandgap center (Hz)', 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
xlim([600 1050])
ylim([1300 1900])
box on
grid on;
set(gca, 'Layer', 'top');
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3]; 
saveas(fig, '2d_hist_7d_input_pd1_quad_2nd_geo_a.fig');

% Create a 2D histogram for Quadrature Surrogate
surrogate_outputs_a_bgs_Gaussian = load('DATASETS/gaussian 6+1 inputs quadrature rule 2nd geo study/surrogate_outputs_2d_q_bgs_1.mat');
surrogate_outputs_a_bgc_Gaussian = load('DATASETS/gaussian 6+1 inputs quadrature rule 2nd geo study/surrogate_outputs_2d_q_bgc_1.mat');

figure;
histogram2(surrogate_outputs_a_bgs_Gaussian.pd_1_outputs, surrogate_outputs_a_bgc_Gaussian.pd_1_outputs, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel('Bandgap size (Hz)', 'FontSize', 12);
ylabel('Bandgap center (Hz)', 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
xlim([600 1050]);
ylim([1300 1900]);
box on
grid on
set(gca, 'Layer', 'top')
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3]; 
saveas(fig, '2d_hist_7d_input_pd1_quad_2nd_geo_b.fig');

% Create a 2D histogram for MC Surrogate

surrogate_outputs_b_bgs_Gaussian = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/surrogate_outputs_2d_r_bgs_1.mat');
surrogate_outputs_b_bgc_Gaussian = load('DATASETS/gaussian 6+1 inputs mc 2nd geo study/surrogate_outputs_2d_r_bgc_1.mat');


figure;
histogram2(surrogate_outputs_b_bgs_Gaussian.pd_1_outputs, surrogate_outputs_b_bgc_Gaussian.pd_1_outputs, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel('Bandgap size (Hz)', 'FontSize', 12);
ylabel('Bandgap center (Hz)', 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
xlim([600 1050]);
ylim([1300 1900]);
box on
grid on
set(gca, 'Layer', 'top')
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3]; 
saveas(fig, '2d_hist_7d_input_pd1_quad_2nd_geo_c.fig');

% Create a 2D histogram for SG Surrogate

surrogate_outputs_c_bgs_Gaussian = load('DATASETS/gaussian 6+1 inputs sparse grid 2nd geo study/surrogate_outputs_2d_sg_bgs_1.mat');
surrogate_outputs_c_bgc_Gaussian = load('DATASETS/gaussian 6+1 inputs sparse grid 2nd geo study/surrogate_outputs_2d_sg_bgc_1.mat');


figure;
histogram2(surrogate_outputs_c_bgs_Gaussian.pd_1_outputs, surrogate_outputs_c_bgc_Gaussian.pd_1_outputs, [50, 50], 'FaceColor', 'flat', 'DisplayStyle', 'tile', 'EdgeColor', 'none', 'ShowEmptyBins', 'off');
xlabel('Bandgap size (Hz)', 'FontSize', 12);
ylabel('Bandgap center (Hz)', 'FontSize', 12);
colormap(custom_colormap);
colorbar;
clim([0 100]);
xlim([600 1050]);
ylim([1300 1900]);
box on
grid on
set(gca, 'Layer', 'top')
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 4, 3]; 
saveas(fig, '2d_hist_7d_input_pd1_quad_2nd_geo_d.fig');

close all
%% Fig 15

ESoft_Uni100 = load('DATASETS/mc_E_soft/E_soft_uniform_mc100.mat');
ESoft_Uni1000 = load('DATASETS/mc_E_soft/E_soft_uniform_mc1000.mat');
ESoft_Uni10000 = load('DATASETS/mc_E_soft/E_soft_uniform_mc10000.mat');

bgSize_Uni100 = load('DATASETS/mc_E_soft/bg_size_uniform_mc100.mat');
bgSize_Uni1000 = load('DATASETS/mc_E_soft/bg_size_uniform_mc1000.mat');
bgSize_Uni10000 = load('DATASETS/mc_E_soft/bg_size_uniform_mc10000.mat');

% Histograms for E_soft data, N=100
figure;
histogram(ESoft_Uni100.E_soft./10^6, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel("E_{soft} (MPa)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% title('E_soft, N=100');
ylim([0 1.4e-2])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_1d_input_E_soft_q_bgs_trio_hist_a.fig');

% Histograms for E_soft data, N=1000
figure;
histogram(ESoft_Uni1000.E_soft./10^6, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel("E_{soft} (MPa)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% title('E_soft, N=1000');
ylim([0 1.4e-2])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_1d_input_E_soft_q_bgs_trio_hist_b.fig');

% Histograms for E_soft data, N=10000
figure;
histogram(ESoft_Uni10000.E_soft./10^6, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel("E_{soft} (MPa)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% title('E_soft, N=10000');
ylim([0 1.4e-2])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_1d_input_E_soft_q_bgs_trio_hist_c.fig');

% Histograms for Bandgap Size data, N=100
figure;
histogram(bgSize_Uni100.bg_size, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel("bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% title('Bandgap Size, N=100');
ylim([0 2.8e-3])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_1d_input_E_soft_q_bgs_trio_hist_d.fig');

% Histograms for Bandgap Size data, N=1000
figure;
histogram(bgSize_Uni1000.bg_size, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel("bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% title('Bandgap Size, N=1000');
ylim([0 2.8e-3])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_1d_input_E_soft_q_bgs_trio_hist_e.fig');

% Histograms for Bandgap Size data, N=10000
figure;
histogram(bgSize_Uni10000.bg_size, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
xlabel("bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
% title('Bandgap Size, N=10000');
ylim([0 2.8e-3])
box on
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, 'P_1d_input_E_soft_q_bgs_trio_hist_f.fig');

close all

%% Fig 16

bgs_mc = load('DATASETS/1D_quad_pd_comparison/MC_10000.mat');
bgs_pd2 = load('DATASETS/1D_quad_pd_comparison/surrogate_PD_2_10000.mat');
bgs_pd3 = load('DATASETS/1D_quad_pd_comparison/surrogate_PD_3_10000.mat');
bgs_pd4 = load('DATASETS/1D_quad_pd_comparison/surrogate_PD_4_10000.mat');
bgs_pd5 = load('DATASETS/1D_quad_pd_comparison/surrogate_PD_5_10000.mat');

[bgs_mc_y,bgs_mc_x] = ksdensity(bgs_mc.MC_samples);
[bgs_pd2_y,bgs_pd2_x] = ksdensity(bgs_pd2.surrogate_samples);
[bgs_pd3_y,bgs_pd3_x] = ksdensity(bgs_pd3.surrogate_samples);
[bgs_pd4_y,bgs_pd4_x] = ksdensity(bgs_pd4.surrogate_samples);
[bgs_pd5_y,bgs_pd5_x] = ksdensity(bgs_pd5.surrogate_samples);

% Histograms for MC
figure;
histogram(bgs_mc.MC_samples, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
hold on
plot(bgs_mc_x, bgs_mc_y, 'r', 'LineWidth', 1.5)
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
box on
xlim([1000 1820])
ylim([0 0.002])
legend('MC', 'MC KDE', 'location', 'best', 'FontSize', 12)
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, '1d_input_E_soft_q_bgs_quinto_hist_a.fig');

% Histograms for surrogate pd 2
figure;
histogram(bgs_pd2.surrogate_samples, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
hold on
plot(bgs_pd2_x, bgs_pd2_y, 'r', 'LineWidth', 1.5)
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
box on
xlim([1000 1820])
ylim([0 0.002])
legend('PD=2', 'PD=2, KDE', 'location', 'best', 'FontSize', 12)
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, '1d_input_E_soft_q_bgs_quinto_hist_b.fig');

% Histograms for surrogate pd 3
figure;
histogram(bgs_pd3.surrogate_samples, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
hold on
plot(bgs_pd2_x, bgs_pd2_y, 'r', 'LineWidth', 1.5)
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
legend('PD=3', 'PD=3, KDE', 'location', 'best', 'FontSize', 12)
box on
xlim([1000 1820])
ylim([0 0.002])
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, '1d_input_E_soft_q_bgs_quinto_hist_c.fig');

% Histograms for surrogate pd 4
figure;
histogram(bgs_pd4.surrogate_samples, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
hold on
plot(bgs_pd2_x, bgs_pd2_y, 'r', 'LineWidth', 1.5)
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
legend('PD=4', 'PD=4, KDE', 'location', 'best', 'FontSize', 12)
box on
xlim([1000 1820])
ylim([0 0.002])
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, '1d_input_E_soft_q_bgs_quinto_hist_d.fig');

% Histograms for surrogate pd 5
figure;
histogram(bgs_pd5.surrogate_samples, 50, 'Normalization', 'pdf', 'EdgeColor', 'black');
hold on
plot(bgs_pd2_x, bgs_pd2_y, 'r', 'LineWidth', 1.5)
xlabel("Bandgap size (Hz)", 'FontSize', 12);
ylabel("Probability density", 'FontSize', 12);
legend('PD=5', 'PD=5, KDE', 'location', 'best', 'FontSize', 12)
box on
xlim([1000 1820])
ylim([0 0.002])
fig = gcf;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3]; 
saveas(fig, '1d_input_E_soft_q_bgs_quinto_hist_e.fig');

close all