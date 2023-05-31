%clc
%clear all

save_outputs = false;

% Check these names before each run!
%load('E_hard_uniform_10000.mat');
%load('E_soft_uniform_test_2000.mat');
E_soft = 100e5;
E_hard = 100e8*ones(1,length(E_soft));
rho_min = 1e3;
rho_max = 8e3;
bg_size = zeros(length(E_hard), 1);
bg_top = zeros(length(E_hard), 1);
bg_bottom = zeros(length(E_hard), 1);

for i = 1:length(E_soft)
    [EIGENVALUE_DATA, WAVEVECTOR_DATA] = Han_dispersion_data(E_hard(i), E_soft(i));
    bg_bottom(i) = max(EIGENVALUE_DATA(:,2));
    bg_top(i) = min(EIGENVALUE_DATA(:,3));
    bg_size(i) = bg_top(i) - bg_bottom(i);
    if bg_size(i) > 0
        disp(['Bandgap present for input ' num2str(i)])
    else
        disp(['No bandgap for input ' num2str(i)])
    end
end

% Change these names before each run!
if save_outputs
    save(['bg_size_uniform_test_' num2str(length(bg_size))], 'bg_size');
    save(['bg_bottom_uniform_test_' num2str(length(bg_bottom))], 'bg_bottom');
    save(['bg_top_uniform_test_' num2str(length(bg_top))], 'bg_top');
    save(['E_hard_uniform_test_' num2str(length(E_hard))], 'E_hard');
    save(['E_soft_uniform_test_' num2str(length(E_soft))], 'E_soft');
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
