clc; close all; clearvars;

% Constants 
R1 = 1/11;          % Longitudinal relaxation rate (s^-1)
FLIP = 10;          % Flip angle
skip_points = 0;    % points to omit from fitting from the start
skip_end_points = 0; % points to omit from fitting from the end

% Load data
load('/home/mark/NMR Data/XmX processing/XmX raw data/LF_Acetylacetone_24_11_25_DMSO_NaOH_ao/ACAC_peak_A_ao_FA_10.mat','peaks','TR_vals');
M_0_A = abs(peaks(1+skip_points:end - skip_end_points));
load('/home/mark/NMR Data/XmX processing/XmX raw data/LF_Acetylacetone_24_11_25_DMSO_NaOH_ao/ACAC_peak_B_ao_FA_10.mat','peaks','TR_vals');
M_0_B = abs(peaks(1+skip_points:end - skip_end_points));
% load('/home/mark/NMR Data/XmX processing/Data/Simulated data/mark_SSFP_test_10_noise_data_offres_R2_eff_5.mat','FIDA_Tacq','FIDB_Tacq','Tacq1');

% Acquisition times
Tacq = TR_vals(1+skip_points:end - skip_end_points)*1e-3;  % seconds
Tacq = Tacq(:); % Ensure column vector

% Load both peak arrays
% M_0_A = abs(FIDA_Tacq(1+skip_points:end - skip_end_points));
% M_0_B = abs(FIDB_Tacq(1+skip_points:end - skip_end_points));

M_0_A = M_0_A(:);
M_0_B = M_0_B(:);

err_A = 0.1*M_0_A;
err_B = 0.1*M_0_B;

% Parameter bounds
% [MA0, MB0, kex_AB, kex_BC, kex_AC, nuA, nuB, nuC, R2]
lb = [0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0.1];   
ub = [0.99, 0.99, 20, 0, 0, 1000, 1000, 0, 10];
Nstart = 1000;  % number of random starting points within bounds
N_boot = 50; % number of bootstrap runs
n_grid = 50; % chi square map resolution

% Objective function evaluating both datasets simultaneously
function residuals = objective_function_global(params, Tacq, M_exp_A, M_exp_B, R1, FLIP)
    
    pop_raw = [params(1), params(2)];
    kex = [params(3), params(4), params(5)];
    nu  = [params(6), params(7), params(8)];
    R2  = params(9);
    
    % Population constraint
    if nu(3) == 0 && kex(2) == 0 && kex(3) == 0 
        % 2-site regime 
        MC0 = 0;
        total_active = sum(pop_raw);
        MA0 = params(1)/total_active;
        MB0 = params(2)/total_active;
    else
        % 3-site regime
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

    % Concatenate signals for global scale projection
    M_sim_cat = [M_A(:); M_B(:)];
    M_exp_cat = [M_exp_A(:); M_exp_B(:)];
    
    % Calculate single global scale factor
    scale = (M_sim_cat'*M_exp_cat)/(M_sim_cat'*M_sim_cat);
    
    % Residuals scaled collectively
    res_signal = (scale*M_sim_cat - M_exp_cat)/mean(abs(M_exp_cat));

    residuals = res_signal;
end

% lsqnonlin options 
options = optimoptions('lsqnonlin','Display','off','MaxIterations',1000,'TolFun',1e-9,'TolX',1e-9);

% Create a problem for MultiStart
x0 = lb + (ub-lb)/2;  % single initial guess in middle of bounds (MultiStart will generate others)
problem = createOptimProblem('lsqnonlin', ...
    'x0', x0, ...
    'objective', @(p)objective_function_global(p,Tacq,M_0_A,M_0_B,R1,FLIP), ...
    'lb', lb, 'ub', ub, ...
    'options', options);

npar = numel(lb);

% Generate Latin Hypercube samples in [lb,ub]
Xlhs = bsxfun(@plus,lb,lhsdesign(Nstart,npar).*(ub-lb));
startSet = CustomStartPointSet(Xlhs);
ms = MultiStart('Display','off','UseParallel',true);
rng default
[best_params,fval,exitflag,output,all_solutions] = run(ms, problem, Nstart);

% Extract parameter matrix from all solutions
param_matrix = vertcat(all_solutions.X);  % Each row = one local minimum
fvals = [all_solutions.Fval];             % Objective (resnorm) values


% Plot histograms for each parameter
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
    % 2-site regime
    total_act = best_params(1) + best_params(2);
    MA_final = best_params(1)/total_act;
    MB_final = best_params(2)/total_act;
    MC_final = 0;
else
    % 3-site regime
    MA_raw = max(0, best_params(1));
    MB_raw = max(0, best_params(2));
    MC_raw = max(0, 1 - MA_raw - MB_raw);
    total = MA_raw + MB_raw + MC_raw;
    MA_final = MA_raw/total;
    MB_final = MB_raw/total;
    MC_final = MC_raw/total;
end

% Extract optimized parameters and compute global scale
[M_A_unscaled, M_B_unscaled, ~] = chem_exchange_sim(FLIP, Tacq, [MA_final,MB_final], ...
                                               [best_params(6),best_params(7),best_params(8)], ...
                                               [best_params(3),best_params(4),best_params(5)], R1, best_params(9));

M_sim_opt_cat = [M_A_unscaled(:); M_B_unscaled(:)];
M_exp_cat = [M_0_A(:); M_0_B(:)];
global_scale = (M_sim_opt_cat'*M_exp_cat)/(M_sim_opt_cat'*M_sim_opt_cat);

M_A_opt = global_scale * M_A_unscaled(:);
M_B_opt = global_scale * M_B_unscaled(:);

% Global Goodness of fit tests
residuals_A = M_A_opt - M_0_A;
residuals_B = M_B_opt - M_0_B;
RMSE = rmse([M_A_opt; M_B_opt], M_exp_cat);

boot_params = zeros(N_boot, numel(best_params));

% Setup variables and normalization factor for Bootstrap
res_vec = objective_function_global(best_params,Tacq,M_0_A,M_0_B,R1,FLIP);
norm_fact = mean(abs(M_exp_cat));
num_pts = numel(Tacq);

% Define opt_fast for the bootstrap fits
opt_fast = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxIterations', 100, 'TolFun', 1e-8);

fprintf('Running Bootstrap... ');

% Split concatenated residuals back into A and B to maintain dataset integrity
res_A = res_vec(1:num_pts);
res_B = res_vec(num_pts+1:end);

parfor b = 1:N_boot
    % Shuffle independently within datasets
    shuffled_res_A = res_A(randi(num_pts, [num_pts, 1]));
    shuffled_res_B = res_B(randi(num_pts, [num_pts, 1]));
    
    noisy_M0_A = M_A_opt(:) + (shuffled_res_A(:)*norm_fact);
    noisy_M0_B = M_B_opt(:) + (shuffled_res_B(:)*norm_fact);
    
    p_start = best_params.*(1 + 0.02*randn(size(best_params)));
    p_start = max(min(p_start, ub), lb);
    
    [p_boot, ~, ~] = lsqnonlin(@(p) objective_function_global(p,Tacq,noisy_M0_A,noisy_M0_B,R1,FLIP), ...
                            p_start, lb, ub, opt_fast);
    boot_params(b, :) = p_boot;
end

% Final outputs
param_errors = std(boot_params);

% Clean up: set errors for inactive parameters to 0
inactive_mask = (ub - lb) < 1e-5; 
param_errors(inactive_mask) = 0;

fprintf('Done. RMSE = %.2f\n', RMSE);

R2_range = linspace(lb(9), ub(9), n_grid); 
k_indices = [3,4,5]; 
k_label_names = {'k_{AB}', 'k_{BC}', 'k_{AC}'};

all_maps = cell(1,3);
K_mesh_list = cell(1,3);
R_mesh_list = cell(1,3);

for s = 1:3
    k_idx = k_indices(s);
    k_range = linspace(lb(k_idx), ub(k_idx), n_grid);
    [K_mesh, R_mesh] = meshgrid(k_range, R2_range);
    
    K_mesh_list{s} = K_mesh;
    R_mesh_list{s} = R_mesh;
    
    % Flatten the grid for a single parfor pass
    K_flat = K_mesh(:);
    R_flat = R_mesh(:);
    chi2_flat = zeros(size(K_flat));
    
    fprintf('Calculating Map %d/3 (%s)... ', s, k_label_names{s});
    
    % parfor here runs on the entire 2500-point grid at once
    parfor idx = 1:numel(K_flat)
        p_temp = best_params;        
        p_temp(k_idx) = K_flat(idx); 
        p_temp(9) = R_flat(idx);     
        res = objective_function_global(p_temp,Tacq,M_0_A,M_0_B,R1,FLIP);
        chi2_flat(idx) = sum(res.^2); 
    end
    
    % Reshape back to the 2D grid
    all_maps{s} = log10(reshape(chi2_flat, n_grid, n_grid));
    fprintf('Done.\n');
end

% Determine global limits for colorbar
global_min = min(cellfun(@(x) min(x(:)), all_maps));
global_max = max(cellfun(@(x) max(x(:)), all_maps));

% Plotting loop 
figure('Name', 'Unified Error Surface Analysis', 'Color', 'w', 'Position', [50, 200, 1600, 500]) 
t = tiledlayout(1, 3, 'TileSpacing', 'Loose', 'Padding', 'Compact');

for s = 1:3
    nexttile
    k_idx = k_indices(s);
    
    contourf(K_mesh_list{s}, R_mesh_list{s}, all_maps{s}, 25, 'LineColor', 'none')
    hold on
    plot(best_params(k_idx), best_params(9), 'r*', 'MarkerSize', 12, 'LineWidth', 2)
    
    colormap(jet)
    clim([global_min, global_max])
    
    c = colorbar;
    c.Label.String = 'log_{10}(\chi^2)';
    
    xlabel(sprintf('%s (s^{-1})', k_label_names{s}), 'FontSize', 11)
    ylabel('R_2 (s^{-1})', 'FontSize', 11)
    title(sprintf('Error Surface: %s vs R_2', k_label_names{s}), 'FontSize', 13)
    grid on
end
 
% Identify active exchange rates and frequencies
active_k = find(ub(3:5) > lb(3:5)) + 2; % Indices 3, 4, or 5
active_nu = find(ub(6:8) > lb(6:8)) + 5; % Indices 6, 7, or 8

% Define map pairings dynamically
if numel(active_k) == 1 && active_k == 3 % 2-Site Case (AB only)
    % Plot k_AB vs nu_A and k_AB vs nu_B
    k_map_indices = [3, 3];
    nu_map_indices = [6, 7];
    map_labels = {'k_{AB} vs \nu_A', 'k_{AB} vs \nu_B'};
else % 3-Site Case 
    % Plot k_AB vs nu_A, k_BC vs nu_B, k_AC vs nu_C
    k_map_indices = [3, 4, 5];
    nu_map_indices = [6, 7, 8];
    map_labels = {'k_{AB} vs \nu_A', 'k_{BC} vs \nu_B', 'k_{AC} vs \nu_C'};
end

n_maps = numel(k_map_indices);
all_maps_nu = cell(1, n_maps);
K_mesh_list_nu = cell(1, n_maps);
NU_mesh_list = cell(1, n_maps);

% Dynamic data generation loop
for s = 1:n_maps
    k_idx = k_map_indices(s);
    nu_idx = nu_map_indices(s);
    
    k_range = linspace(lb(k_idx), ub(k_idx), n_grid);
    nu_range = linspace(lb(nu_idx), ub(nu_idx), n_grid);
    [K_mesh, NU_mesh] = meshgrid(k_range, nu_range);
    
    K_mesh_list_nu{s} = K_mesh;
    NU_mesh_list{s} = NU_mesh;
    
    K_flat = K_mesh(:);
    NU_flat = NU_mesh(:);
    chi2_flat = zeros(size(K_flat));
    
    fprintf('Calculating Dynamic Map %d/%d (%s)... ', s, n_maps, map_labels{s});
    
    parfor idx = 1:numel(K_flat)
        p_temp = best_params;        
        p_temp(k_idx) = K_flat(idx); 
        p_temp(nu_idx) = NU_flat(idx);     
        res = objective_function_global(p_temp,Tacq,M_0_A,M_0_B,R1,FLIP);
        chi2_flat(idx) = sum(res.^2); 
    end
    all_maps_nu{s} = log10(reshape(chi2_flat, n_grid, n_grid));
    fprintf('Done.\n');
end

% Determine global limits for colorbar
global_min_nu = min(cellfun(@(x) min(x(:)), all_maps_nu));
global_max_nu = max(cellfun(@(x) max(x(:)), all_maps_nu));

% Dynamic plotting loop 
all_names = cell(1, 9);
all_names{1} = 'M_{A0}';
all_names{2} = 'M_{B0}';
all_names{3} = 'k_{AB}';
all_names{4} = 'k_{BC}';
all_names{5} = 'k_{AC}';
all_names{6} = '\nu_A';
all_names{7} = '\nu_B';
all_names{8} = '\nu_C';
all_names{9} = 'R_2';

figure('Name', 'Dynamic Frequency vs Exchange Analysis', 'Color', 'w', 'Position', [50, 50, 500*n_maps, 500]); 
t2 = tiledlayout(1, n_maps, 'TileSpacing', 'Loose', 'Padding', 'Compact');

for s = 1:n_maps
    nexttile
    k_idx = k_map_indices(s);
    nu_idx = nu_map_indices(s);

    contourf(K_mesh_list_nu{s}, NU_mesh_list{s}, all_maps_nu{s}, 25, 'LineColor', 'none')
    hold on
    plot(best_params(k_idx), best_params(nu_idx), 'r*', 'MarkerSize', 12, 'LineWidth', 2)

    colormap(jet)

    % Set limits based on the specific map's data for better contrast
    clim([min(all_maps_nu{s}(:)), max(all_maps_nu{s}(:))]);

    c = colorbar;
    c.Label.String = 'log_{10}(\chi^2)';

    xlabel_str = sprintf('%s (s^{-1})', all_names{k_idx});
    ylabel_str = sprintf('%s (Hz)', all_names{nu_idx});

    xlabel(xlabel_str, 'FontSize', 12, 'Interpreter', 'tex');
    ylabel(ylabel_str, 'FontSize', 12, 'Interpreter', 'tex');

    title(map_labels{s}, 'FontSize', 13, 'Interpreter', 'tex')
    grid on
end

% Plot experimental and fitted data for both peaks
title_str = { ...
    sprintf('Global Fit: M_{A0}=%.3f\\pm%.3f, M_{B0}=%.3f\\pm%.3f', ...
            MA_final, param_errors(1), MB_final, param_errors(2)), ...
    sprintf('k_{AB}=%.2f\\pm%.2f s^{-1}, R_{2}=%.3f\\pm%.3f s^{-1}', ...
            best_params(3), param_errors(3), best_params(9), param_errors(9)) ...
    sprintf('\\nu_{A}=%.1f\\pm%.1f Hz, \\nu_{B}=%.1f\\pm%.1f Hz', ...
             best_params(6), param_errors(6), best_params(7), param_errors(7))
};

norm_plot = max(M_exp_cat); % Normalize plotting based on the highest point across both arrays

figure('Name', 'Global SSFP Fit', 'Color', 'w', 'Position', [50, 100, 1000, 600])
plot(Tacq*1e3, M_0_A/norm_plot, 'b-','LineWidth',2,'DisplayName','Exp Peak A')
hold on
plot(Tacq*1e3, M_A_opt/norm_plot, 'r--','LineWidth',2,'DisplayName','Fit Peak A')
xlabel('TR [ms]', 'FontSize', 12)
ylabel('Normalized Intensity', 'FontSize', 12);
legend('Location', 'best')
title('Peak A ',title_str,'Interpreter','tex', 'FontSize', 13)
set(gca,'FontSize',12)

figure('Name', 'Global SSFP Fit', 'Color', 'w', 'Position', [50, 100, 1000, 600])
plot(Tacq*1e3, M_0_B/norm_plot, 'b-','LineWidth',2,'DisplayName','Exp Peak B')
hold on
plot(Tacq*1e3, M_B_opt/norm_plot, 'r--','LineWidth',2,'DisplayName','Fit Peak B')
xlabel('TR [ms]', 'FontSize', 12)
ylabel('Normalized Intensity', 'FontSize', 12);
legend('Location', 'best')
title('Peak B ',title_str,'Interpreter','tex', 'FontSize', 13)
set(gca,'FontSize',12)

% Plot residuals globally
figure('Name', 'Global Residuals', 'Color', 'w', 'Position', [50, 100, 1000, 400])
errorbar(Tacq*1e3, residuals_A, err_A, 'bo','MarkerSize',5,'CapSize',4,'LineWidth',1, 'DisplayName','Res A')
hold on
errorbar(Tacq*1e3, residuals_B, err_B, 'ro','MarkerSize',5,'CapSize',4,'LineWidth',1, 'DisplayName','Res B')
yline(0,'--k','LineWidth',1, 'HandleVisibility', 'off')

xlabel('TR [ms]', 'FontSize', 12)
ylabel('Residuals', 'FontSize', 12)
legend('Location', 'best')
title(sprintf('Global RMSE = %.3e', RMSE), 'FontSize', 13)
set(gca,'FontSize',12)