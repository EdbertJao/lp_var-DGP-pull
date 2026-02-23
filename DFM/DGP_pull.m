%% DGP STORE (DFM): extraction-only
% Creates DGPSTORE files for:
%   dgp_type in {'G','MP'}
%   estimand_type in {'ObsShock','IV'}
% Keeps settings.simul.n_MC = 5000 and settings.specifications.random_n_spec = 100 (from Settings).
%
% Output per config (single .mat):
%   - X_all:     T x n_y x n_MC        (single)  simulated panel (data_sim_all.data_y)
%   - shock_all: T x n_MC              (single)  true shock series (data_sim_all.data_shock)
%   - Z_all:     T x nrho x nsig x n_MC (single) IV grid (only if estimand_type=='IV')
%   - var_select: n_spec x n_var       (double)  DGP variable indices
%   - var_name_short/long: 1 x n_y     (cell)    names aligned with X_all columns
%   - var_name_short_spec/long_spec: n_spec x n_var (cell) names per DGP
%   - irf_full:  IRF_hor x n_y         (double)  full-panel population IRF
%   - target_irf: IRF_hor x n_spec     (double)  truth IRF per DGP/spec
%   - plus minimal settings + identifiers

clc; clear; close all;

%% Paths (assumes this script is in the same folder as run_dfm.m)
addpath(genpath(fullfile('..', 'Auxiliary_Functions')));
addpath(genpath(fullfile('..', 'Estimation_Routines')));
addpath(genpath('Subroutines'))  % includes dgp_irfs_stats, ABCD_fun_DFM, etc.
addpath(genpath('Settings'));

%% Experiment knobs (kept consistent with run_dfm defaults)
spec_id    = 1;   % seed for drawing the 100 specs (DGP definitions)
lag_type   = 4;   % unused in extraction, but Settings/shared uses it
mode_type  = 1;   % robustness mode used by Settings/check_mode
estim_diagn = 0;  % leave 0 for batch run

%% Config lists (ONLY the ones you care about)
dgp_list      = {'G','MP'};
estimand_list = {'ObsShock','IV'};

%% Output folder
out_root = fullfile('Results', 'DGP_Stores');
if ~exist(out_root, 'dir'); mkdir(out_root); end

%% Main loop over 4 configurations
for id = 1:numel(dgp_list)
    dgp_type = dgp_list{id};

    for ie = 1:numel(estimand_list)
        estimand_type = estimand_list{ie};

        fprintf('\n====================================\n');
        fprintf('Building DGP store for: dgp_type=%s, estimand=%s\n', dgp_type, estimand_type);
        fprintf('====================================\n');

        %% SETTINGS (same structure as run_dfm.m)
        clear settings DF_model DFM_estimate

        run(fullfile('Settings', 'shared'));
        run(fullfile('Settings', dgp_type));
        run(fullfile('Settings', estimand_type));
        run(fullfile('Settings', 'check_mode'));

        % Sanity: enforce that we are *not* shrinking n_MC / n_spec here
        % (We rely on Settings values: n_MC=5000, random_n_spec=100)

        %% Estimate DFM from dataset
        DFM_estimate = DFM_est(DF_model.n_fac, DF_model.n_lags_fac, DF_model.n_lags_uar, ...
                               DF_model.reorder, DF_model.levels, DF_model.coint_rank);

        if estim_diagn == 1
            run_estim_diagn;
        end

        %% Store estimated DFM parameters into DF_model
        DF_model.Phi       = DFM_estimate.Phi;
        DF_model.Sigma_eta = DFM_estimate.Sigma_eta;

        DF_model.Lambda    = DFM_estimate.Lambda;
        DF_model.delta     = DFM_estimate.delta;
        DF_model.sigma_v   = DFM_estimate.sigma_v;

        DF_model.variable_name_code  = DFM_estimate.bpnamevec;
        DF_model.variable_name_short = DFM_estimate.bplabvec_short;
        DF_model.variable_name_long  = DFM_estimate.bplabvec_long;
        DF_model.trans_code          = DFM_estimate.bptcodevec;

        %% IV grid setup (only if IV)
        if strcmp(estimand_type, 'IV')
            DF_model.IV.rho_grid     = DF_model.IV.rho * settings.est.IV.IV_persistence_scale;
            DF_model.IV.sigma_v_grid = DF_model.IV.sigma_v * settings.est.IV.IV_strength_scale;
        end

        %% ABCD form
        [DF_model.n_y, DF_model.n_fac] = size(DF_model.Lambda);
        DF_model.ABCD = ABCD_fun_DFM(DF_model);

        %% Shock weights (same as run_dfm.m)
        shock_weight = zeros(DF_model.n_fac + DF_model.n_y, 1);
        if settings.est.estimate_shock_weight == 1
            shock_weight(1:DF_model.n_fac) = DF_model.ABCD.D(settings.est.shock_optimize_var_IRF, 1:DF_model.n_fac);
            shock_weight = shock_weight / sqrt(shock_weight' * shock_weight);
        else
            shock_weight = zeros(n_eps, 1);
            shock_weight(settings.est.manual_shock_pos) = 1;
        end
        settings.est.shock_weight = shock_weight;
        clear shock_weight;

        %% Draw the 100 DGP definitions ONCE
        settings.specifications = pick_var_fn(DF_model, settings, spec_id);

        %% Compute population IRFs / truth ONCE
        DF_model = dgp_irfs_stats(DF_model, settings, estimand_type);

        %% Prepare output file (one per configuration)
        fname = sprintf('DGPSTORE_DFM_%s_%s_specid%d_mode%d.mat', dgp_type, estimand_type, spec_id, mode_type);
        outfile = fullfile(out_root, fname);

        if exist(outfile, 'file')
            delete(outfile); % avoid accidental appends with old sizes
        end

        % Use matfile for incremental writes
        M = matfile(outfile, 'Writable', true);

        %% Save small metadata immediately
        M.spec_id       = spec_id;
        M.dgp_type      = dgp_type;
        M.estimand_type = estimand_type;
        M.lag_type      = lag_type;
        M.mode_type     = mode_type;

        % Store settings pieces that define the run
        M.settings_simul = settings.simul;
        M.settings_est   = rmfield(settings.est, intersect(fieldnames(settings.est), {'methods_name'})); % keep small
        M.settings_specifications = settings.specifications;

        % Store names and truth IRFs
        M.var_name_short = DF_model.variable_name_short;
        M.var_name_long  = DF_model.variable_name_long;
        M.var_name_code  = DF_model.variable_name_code;
        M.trans_code     = DF_model.trans_code;

        M.irf_full    = DF_model.irf;        % IRF_hor x n_y
        M.target_irf  = DF_model.target_irf; % IRF_hor x n_spec

        var_select = settings.specifications.var_select; % n_spec x n_var
        M.var_select = var_select;

        % Per-DGP name tables (n_spec x n_var)
        M.var_name_short_spec = DF_model.variable_name_short(var_select);
        M.var_name_long_spec  = DF_model.variable_name_long(var_select);

        %% Preallocate large arrays on disk WITHOUT allocating in RAM
        T    = settings.simul.T;
        n_y  = DF_model.n_y;
        n_MC = settings.simul.n_MC;

        % X_all: T x n_y x n_MC
        M.X_all(T, n_y, n_MC) = single(0);
        % shock_all: T x n_MC
        M.shock_all(T, n_MC) = single(0);

        % Z_all only for IV
        if settings.est.with_IV == 1
            nrho = length(DF_model.IV.rho_grid);
            nsig = length(DF_model.IV.sigma_v_grid);
            M.IV_rho_grid     = DF_model.IV.rho_grid;
            M.IV_sigma_v_grid = DF_model.IV.sigma_v_grid;

            % Z_all: T x nrho x nsig x n_MC
            M.Z_all(T, nrho, nsig, n_MC) = single(0);
        end

        %% Monte Carlo: simulate and write only
        fprintf('Simulating n_MC=%d draws...\n', n_MC);

        for i_MC = 1:n_MC
            if mod(i_MC, 100) == 0
                fprintf('  MC %d / %d\n', i_MC, n_MC);
            end

            rng(settings.simul.seed(i_MC), 'twister');

            data_sim_all = generate_data(DF_model, settings);

            % Write to disk (cast to single for size)
            M.X_all(:,:,i_MC)   = single(data_sim_all.data_y);
            M.shock_all(:,i_MC) = single(data_sim_all.data_shock);

            if settings.est.with_IV == 1
                M.Z_all(:,:,:,i_MC) = single(data_sim_all.data_z);
            end
        end

        fprintf('DONE. Wrote: %s\n', outfile);
    end
end

fprintf('\nAll configurations completed. Output folder:\n  %s\n', out_root);
