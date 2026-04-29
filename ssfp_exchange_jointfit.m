clc; close all; clearvars;
% ssfp_exchange_jointfit - Multi-site Global Joint Fitting Engine
%
% Performs non-linear least squares optimization to extract chemical exchange 
% parameters from SSFP NMR data.
%
% USAGE:
%   Run script directly. Toggle 'num_sites' variable for 2 or 3 site models.
%
% INPUT DATA REQUIREMENTS:
%   Requires .mat files containing:
%   - 'peaks': Vector of complex or absolute peak intensities.
%   - 'TR_vals': Vector of repetition times (ms).
%
%
% SEE ALSO: chem_exchange_sim, lsqnonlin, MultiStart

% Processing options
num_sites = 2;      % Toggle between 2 and 3 site joint fitting

% Constants 
R1 = 1/10;          % Longitudinal relaxation rate (s^-1)
skip_points = 10;    % points to omit from fitting from the start
skip_end_points = 0; % points to omit from fitting from the end

% Load data
load('/home/mark/NMR Data/XmX processing/Data/Simulated data/mark_SSFP_test_data_offres_R2_less_than_10_A.mat','peaks','TR_vals','FLIP');
M_0_A = abs(peaks(1+skip_points:end - skip_end_points));
load('/home/mark/NMR Data/XmX processing/Data/Simulated data/mark_SSFP_test_data_offres_R2_less_than_10_B.mat','peaks','TR_vals','FLIP');
M_0_B = abs(peaks(1+skip_points:end - skip_end_points));

if num_sites == 3
    % Replace path below with actual Peak C data
    % load('/home/mark/NMR Data/XmX processing/XmX raw data/LF_Acetylacetone_24_11_25_DMSO_NaOH_ao/ACAC_peak_C_ao_FA_10.mat','peaks','TR_vals','FA');
    % M_0_C = abs(peaks(1+skip_points:end - skip_end_points));
    % For demonstration, allocating a dummy vector if the file doesn't exist yet:
    M_0_C = zeros(size(M_0_A)); 
    M_0_C = M_0_C(:);
else
    M_0_C = [];
end

% 2 site exchange boundaries
% Parameter bounds: [MA0, MB0, kex_AB, kex_BC, kex_AC, nuA, nuB, nuC, R2]
if num_sites == 2
    lb = [0.01, 0.01, 1, 0, 0, 0,    0,    0, 1];   
    ub = [0.99, 0.99, 5, 0, 0, 1000, 1000, 0, 10];

% 3 site exchange boundaries
% Parameter bounds: [MA0, MB0,MC0, kex_AB, kex_BC, kex_AC, nuA, nuB, nuC, R2]
elseif num_sites == 3
    lb = [0.01, 0.01, 0.1, 0.1, 0.1, 0, 0, 0, 0.1];   
    ub = [0.99, 0.99, 20, 20, 20, 1000, 1000, 1000, 10];
end

Nstart = 1000;  % number of random starting points within bounds
N_boot = 50;    % number of bootstrap runs
n_grid = 50;    % chi square map resolution

% Acquisition times
TR_vals = TR_vals(1+skip_points:end - skip_end_points);  % seconds
TR_vals = TR_vals(:); % Ensure column vector

M_0_A = M_0_A(:);
M_0_B = M_0_B(:);


% Objective function evaluating datasets simultaneously
function residuals = objective_function_global(params, Tacq, M_exp_A, M_exp_B, M_exp_C, R1, FLIP, num_sites)
    
    pop_raw = [params(1), params(2)];
    kex = [params(3), params(4), params(5)];
    nu  = [params(6), params(7), params(8)];
    R2  = params(9);
    
    % Population constraint
    if nu(3) == 0 && kex(2) == 0 && kex(3) == 0 
        MC0 = 0;
        total_active = sum(pop_raw);
        MA0 = params(1)/total_active;
        MB0 = params(2)/total_active;
    else
        MA0_raw = params(1);
        MB0_raw = params(2);
        MC0_raw = max(0, 1 - MA0_raw - MB0_raw);
        
        total = MA0_raw + MB0_raw + MC0_raw;
        MA0 = MA0_raw/total;
        MB0 = MB0_raw/total;
        MC0 = MC0_raw/total;
    end

    pop_constrained = [MA0, MB0, MC0];
    [M_A, M_B, M_C] = chem_exchange_sim(FLIP, Tacq, pop_constrained, nu, kex, R1, R2);

    % Concatenate signals dynamically based on site number
    if num_sites == 2
        M_sim_cat = [M_A(:); M_B(:)];
        M_exp_cat = [M_exp_A(:); M_exp_B(:)];
    else
        M_sim_cat = [M_A(:); M_B(:); M_C(:)];
        M_exp_cat = [M_exp_A(:); M_exp_B(:); M_exp_C(:)];
    end
    
    % Calculate single global scale factor
    scale = (M_sim_cat'*M_exp_cat)/(M_sim_cat'*M_sim_cat);
    
    % Residuals scaled collectively
    res_signal = (scale*M_sim_cat - M_exp_cat)/mean(abs(M_exp_cat));

    residuals = res_signal;
end

% lsqnonlin options 
options = optimoptions('lsqnonlin','Display','off','MaxIterations',1000,'TolFun',1e-9,'TolX',1e-9);

% Create a problem for MultiStart
x0 = lb + (ub-lb)/2;  
problem = createOptimProblem('lsqnonlin', ...
    'x0', x0, ...
    'objective', @(p)objective_function_global(p,TR_vals,M_0_A,M_0_B,M_0_C,R1,FLIP,num_sites), ...
    'lb', lb, 'ub', ub, ...
    'options', options);

npar = numel(lb);

% Generate Latin Hypercube samples
Xlhs = bsxfun(@plus,lb,lhsdesign(Nstart,npar).*(ub-lb));
startSet = CustomStartPointSet(Xlhs);
ms = MultiStart('Display','off','UseParallel',true);
rng default
[best_params,fval,exitflag,output,all_solutions] = run(ms, problem, Nstart);

% Extract parameter matrix
param_matrix = vertcat(all_solutions.X);  
fvals = [all_solutions.Fval];             

% Plot histograms
param_names = {'MA0','MB0','kAB','kBC','kAC','nuA','nuB','nuC','R2'};
figure
for i = 1:size(param_matrix,2)
    subplot(3,3,i)
    histogram(param_matrix(:,i),15,'FaceColor',[0.2 0.2 0.8],'EdgeColor','k')
    xlabel(param_names{i},'FontSize',12)
    ylabel('Count','FontSize',12)
    title(sprintf('%s (mean=%.3f)',param_names{i},mean(param_matrix(:,i))))
    grid on
end

if best_params(8) == 0 && best_params(4) == 0 && best_params(5) == 0 
    total_act = best_params(1) + best_params(2);
    MA_final = best_params(1)/total_act;
    MB_final = best_params(2)/total_act;
    MC_final = 0;
else
    MA_raw = max(0, best_params(1));
    MB_raw = max(0, best_params(2));
    MC_raw = max(0, 1 - MA_raw - MB_raw);
    total = MA_raw + MB_raw + MC_raw;
    MA_final = MA_raw/total;
    MB_final = MB_raw/total;
    MC_final = MC_raw/total;
end

% Extract optimized parameters and compute global scale
[M_A_unscaled, M_B_unscaled, M_C_unscaled] = chem_exchange_sim(FLIP, TR_vals, [MA_final,MB_final], ...
                                               [best_params(6),best_params(7),best_params(8)], ...
                                               [best_params(3),best_params(4),best_params(5)], R1, best_params(9));

if num_sites == 2
    M_sim_opt_cat = [M_A_unscaled(:); M_B_unscaled(:)];
    M_exp_cat = [M_0_A(:); M_0_B(:)];
else
    M_sim_opt_cat = [M_A_unscaled(:); M_B_unscaled(:); M_C_unscaled(:)];
    M_exp_cat = [M_0_A(:); M_0_B(:); M_0_C(:)];
end

global_scale = (M_sim_opt_cat'*M_exp_cat)/(M_sim_opt_cat'*M_sim_opt_cat);

M_A_opt = global_scale*M_A_unscaled(:);
M_B_opt = global_scale*M_B_unscaled(:);
M_C_opt = global_scale*M_C_unscaled(:);

% Global Goodness of fit tests
residuals_A = M_A_opt - M_0_A;
residuals_B = M_B_opt - M_0_B;
if num_sites == 3
    residuals_C = M_C_opt - M_0_C;
end

boot_params = zeros(N_boot, numel(best_params));

% Setup variables and normalization factor for Bootstrap
res_vec = objective_function_global(best_params,TR_vals,M_0_A,M_0_B,M_0_C,R1,FLIP,num_sites);
norm_fact = mean(abs(M_exp_cat));
num_pts = numel(TR_vals);

opt_fast = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 100, 'TolFun', 1e-8);
fprintf('Running Bootstrap... ');

% Split concatenated residuals dynamically
res_A = res_vec(1:num_pts);
res_B = res_vec(num_pts+1:2*num_pts);
if num_sites == 3
    res_C = res_vec(2*num_pts+1:3*num_pts);
else
    res_C = [];
end

parfor b = 1:N_boot
    shuffled_res_A = res_A(randi(num_pts, [num_pts, 1]));
    shuffled_res_B = res_B(randi(num_pts, [num_pts, 1]));
    
    noisy_M0_A = M_A_opt(:) + (shuffled_res_A(:)*norm_fact);
    noisy_M0_B = M_B_opt(:) + (shuffled_res_B(:)*norm_fact);
    
    if num_sites == 3
        shuffled_res_C = res_C(randi(num_pts, [num_pts, 1]));
        noisy_M0_C = M_C_opt(:) + (shuffled_res_C(:)*norm_fact);
    else
        noisy_M0_C = [];
    end
    
    p_start = best_params.*(1 + 0.02*randn(size(best_params)));
    p_start = max(min(p_start, ub), lb);
    
    [p_boot, ~, ~] = lsqnonlin(@(p) objective_function_global(p,TR_vals,noisy_M0_A,noisy_M0_B,noisy_M0_C,R1,FLIP,num_sites), ...
                            p_start, lb, ub, opt_fast);
    boot_params(b, :) = p_boot;
end

param_errors = std(boot_params);
inactive_mask = (ub - lb) < 1e-5; 
param_errors(inactive_mask) = 0;
fprintf('Done.');

% Map calculations
R2_range = linspace(lb(9), ub(9), n_grid); 
k_indices = [3,4,5]; 
k_label_names = {'k_{AB}', 'k_{BC}', 'k_{AC}'};

all_maps = cell(1,3);
K_mesh_list = cell(1,3);
R_mesh_list = cell(1,3);

for s = 1:3
    k_idx = k_indices(s);
    if ub(k_idx) == lb(k_idx); continue; end % Skip inactive kinetic bounds
    
    k_range = linspace(lb(k_idx), ub(k_idx), n_grid);
    [K_mesh, R_mesh] = meshgrid(k_range, R2_range);
    
    K_mesh_list{s} = K_mesh;
    R_mesh_list{s} = R_mesh;
    
    K_flat = K_mesh(:);
    R_flat = R_mesh(:);
    chi2_flat = zeros(size(K_flat));
    
    fprintf('Calculating Map %d/3 (%s)... ', s, k_label_names{s});
    
    parfor idx = 1:numel(K_flat)
        p_temp = best_params;        
        p_temp(k_idx) = K_flat(idx); 
        p_temp(9) = R_flat(idx);     
        res = objective_function_global(p_temp,TR_vals,M_0_A,M_0_B,M_0_C,R1,FLIP,num_sites);
        chi2_flat(idx) = sum(res.^2); 
    end
    
    all_maps{s} = log10(reshape(chi2_flat, n_grid, n_grid));
    fprintf('Done.\n');
end

% Map Plotting Logic
valid_maps = ~cellfun(@isempty, all_maps);
n_active_maps = sum(valid_maps);

if n_active_maps > 0
    % Find global limits across active maps only
    map_mins = cellfun(@(x) min(x(:)), all_maps(valid_maps));
    map_maxs = cellfun(@(x) max(x(:)), all_maps(valid_maps));
    global_min = min(map_mins);
    global_max = max(map_maxs);

    figure('Name', 'Unified Error Surface Analysis', 'Color', 'w', 'Position', [50, 200, 500*n_active_maps, 500]) 
    t = tiledlayout(1, n_active_maps, 'TileSpacing', 'Loose', 'Padding', 'Compact');

    for s = 1:3
        if ~valid_maps(s); continue; end
        nexttile
        k_idx = k_indices(s);
        
        contourf(K_mesh_list{s}, R_mesh_list{s}, all_maps{s}, 25, 'LineColor', 'none')
        hold on
        plot(best_params(k_idx), best_params(9), 'r*', 'MarkerSize', 12, 'LineWidth', 2)
        colormap(jet)
        clim([global_min, global_max])
        c = colorbar; c.Label.String = 'log_{10}(\chi^2)';
        xlabel(sprintf('%s (s^{-1})', k_label_names{s}), 'FontSize', 11)
        ylabel('R_2 (s^{-1})', 'FontSize', 11)
        title(sprintf('Error Surface: %s vs R_2', k_label_names{s}), 'FontSize', 13)
        grid on
    end
end

% Identify active frequencies and exchanges for dynamic plots
if num_sites == 2
    k_map_indices = [3, 3];
    nu_map_indices = [6, 7];
    map_labels = {'k_{AB} vs \nu_A', 'k_{AB} vs \nu_B'};
else
    k_map_indices = [3, 4, 5];
    nu_map_indices = [6, 7, 8];
    map_labels = {'k_{AB} vs \nu_A', 'k_{BC} vs \nu_B', 'k_{AC} vs \nu_C'};
end

n_maps = numel(k_map_indices);
all_maps_nu = cell(1, n_maps);
K_mesh_list_nu = cell(1, n_maps);
NU_mesh_list = cell(1, n_maps);

for s = 1:n_maps
    k_idx = k_map_indices(s);
    nu_idx = nu_map_indices(s);
    
    k_range = linspace(lb(k_idx), ub(k_idx), n_grid);
    nu_range = linspace(lb(nu_idx), ub(nu_idx), n_grid);
    [K_mesh, NU_mesh] = meshgrid(k_range, nu_range);
    
    K_mesh_list_nu{s} = K_mesh;
    NU_mesh_list{s} = NU_mesh;
    K_flat = K_mesh(:); NU_flat = NU_mesh(:);
    chi2_flat = zeros(size(K_flat));
    
    fprintf('Calculating Dynamic Map %d/%d (%s)... ', s, n_maps, map_labels{s});
    parfor idx = 1:numel(K_flat)
        p_temp = best_params;        
        p_temp(k_idx) = K_flat(idx); 
        p_temp(nu_idx) = NU_flat(idx);     
        res = objective_function_global(p_temp,TR_vals,M_0_A,M_0_B,M_0_C,R1,FLIP,num_sites);
        chi2_flat(idx) = sum(res.^2); 
    end
    all_maps_nu{s} = log10(reshape(chi2_flat, n_grid, n_grid));
    fprintf('Done.\n');
end

all_names = {'M_{A0}','M_{B0}','k_{AB}','k_{BC}','k_{AC}','\nu_A','\nu_B','\nu_C','R_2'};
figure('Name', 'Dynamic Frequency vs Exchange Analysis', 'Color', 'w', 'Position', [50, 50, 500*n_maps, 500]); 
t2 = tiledlayout(1, n_maps, 'TileSpacing', 'Loose', 'Padding', 'Compact');

for s = 1:n_maps
    nexttile
    k_idx = k_map_indices(s); nu_idx = nu_map_indices(s);
    contourf(K_mesh_list_nu{s}, NU_mesh_list{s}, all_maps_nu{s}, 25, 'LineColor', 'none')
    hold on; plot(best_params(k_idx), best_params(nu_idx), 'r*', 'MarkerSize', 12, 'LineWidth', 2)
    colormap(jet); clim([min(all_maps_nu{s}(:)), max(all_maps_nu{s}(:))]);
    c = colorbar; c.Label.String = 'log_{10}(\chi^2)';
    xlabel(sprintf('%s (s^{-1})', all_names{k_idx}), 'FontSize', 12, 'Interpreter', 'tex');
    ylabel(sprintf('%s (Hz)', all_names{nu_idx}), 'FontSize', 12, 'Interpreter', 'tex');
    title(map_labels{s}, 'FontSize', 13, 'Interpreter', 'tex'); grid on
end

% Fit Curves Plotting
title_str = { ...
    sprintf('Global Fit: M_{A0}=%.3f\\pm%.3f, M_{B0}=%.3f\\pm%.3f', MA_final, param_errors(1), MB_final, param_errors(2)), ...
    sprintf('k_{AB}=%.2f\\pm%.2f s^{-1}, R_{2}=%.3f\\pm%.3f s^{-1}', best_params(3), param_errors(3), best_params(9), param_errors(9)), ...
    sprintf('\\nu_{A}=%.1f\\pm%.1f Hz, \\nu_{B}=%.1f\\pm%.1f Hz', best_params(6), param_errors(6), best_params(7), param_errors(7))
};

norm_plot = max(M_exp_cat); 

figure('Name', 'Global SSFP Fit - A', 'Color', 'w', 'Position', [50, 100, 1000, 600])
plot(TR_vals*1e3, M_0_A/norm_plot, 'b-','LineWidth',2,'DisplayName','Exp Peak A'); hold on
plot(TR_vals*1e3, M_A_opt/norm_plot, 'r--','LineWidth',2,'DisplayName','Fit Peak A')
xlabel('TR [ms]', 'FontSize', 12); ylabel('Normalized Intensity', 'FontSize', 12); legend('Location', 'best')
title('Peak A ',title_str,'Interpreter','tex', 'FontSize', 13); set(gca,'FontSize',12)

figure('Name', 'Global SSFP Fit - B', 'Color', 'w', 'Position', [50, 100, 1000, 600])
plot(TR_vals*1e3, M_0_B/norm_plot, 'b-','LineWidth',2,'DisplayName','Exp Peak B'); hold on
plot(TR_vals*1e3, M_B_opt/norm_plot, 'r--','LineWidth',2,'DisplayName','Fit Peak B')
xlabel('TR [ms]', 'FontSize', 12); ylabel('Normalized Intensity', 'FontSize', 12); legend('Location', 'best')
title('Peak B ',title_str,'Interpreter','tex', 'FontSize', 13); set(gca,'FontSize',12)

if num_sites == 3
    figure('Name', 'Global SSFP Fit - C', 'Color', 'w', 'Position', [50, 100, 1000, 600])
    plot(TR_vals*1e3, M_0_C/norm_plot, 'b-','LineWidth',2,'DisplayName','Exp Peak C'); hold on
    plot(TR_vals*1e3, M_C_opt/norm_plot, 'r--','LineWidth',2,'DisplayName','Fit Peak C')
    xlabel('TR [ms]', 'FontSize', 12); ylabel('Normalized Intensity', 'FontSize', 12); legend('Location', 'best')
    title('Peak C ',title_str,'Interpreter','tex', 'FontSize', 13); set(gca,'FontSize',12)
end

