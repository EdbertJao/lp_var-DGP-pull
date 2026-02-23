function data_sim_select = select_data_fn(data_sim_all, settings, i_spec)
% Extract simulated data for one DGP (one specification).
% NEW: attach selected variable indices + names + truth IRFs.

var_select = settings.specifications.var_select;
with_IV    = settings.est.with_IV;

% IV grid selection indices
if with_IV == 1
    rho_select_grid_idx     = settings.specifications.rho_select_grid_idx;
    sigma_v_select_grid_idx = settings.specifications.sigma_v_select_grid_idx;

    this_rho_grid_idx   = rho_select_grid_idx(i_spec);
    this_sigma_grid_idx = sigma_v_select_grid_idx(i_spec);
else
    this_rho_grid_idx   = 1;
    this_sigma_grid_idx = 1;
end

% Core fields used by estimation code
var_idx = var_select(i_spec,:);
data_sim_select.data_y     = data_sim_all.data_y(:, var_idx);
data_sim_select.data_z     = data_sim_all.data_z(:, this_rho_grid_idx, this_sigma_grid_idx);
data_sim_select.data_shock = data_sim_all.data_shock;

% =========================
% NEW: attach metadata
% =========================
data_sim_select.meta = struct();
data_sim_select.meta.i_spec  = i_spec;
data_sim_select.meta.var_idx = var_idx;

% Selected variable names (if available)
if isfield(data_sim_all,'meta')
    if isfield(data_sim_all.meta,'variable_name_short')
        data_sim_select.meta.var_name_short = data_sim_all.meta.variable_name_short(var_idx);
    end
    if isfield(data_sim_all.meta,'variable_name_long')
        data_sim_select.meta.var_name_long  = data_sim_all.meta.variable_name_long(var_idx);
    end
end

% IV selection metadata
data_sim_select.meta.iv_rho_grid_idx   = this_rho_grid_idx;
data_sim_select.meta.iv_sigma_grid_idx = this_sigma_grid_idx;

if with_IV == 1 && isfield(data_sim_all,'meta') && isfield(data_sim_all.meta,'IV_rho_grid')
    data_sim_select.meta.iv_rho     = data_sim_all.meta.IV_rho_grid(this_rho_grid_idx);
    data_sim_select.meta.iv_sigma_v = data_sim_all.meta.IV_sigma_v_grid(this_sigma_grid_idx);
end

% =========================
% NEW: attach truth
% =========================
data_sim_select.truth = struct();

% “Truth IRF” used for evaluation: one column per DGP/spec
if isfield(data_sim_all,'truth') && isfield(data_sim_all.truth,'target_irf')
    data_sim_select.truth.target_irf = data_sim_all.truth.target_irf(:, i_spec); % 21 x 1
end

% Optional: full-panel truth IRFs for the selected variables (21 x n_var)
if isfield(data_sim_all,'truth') && isfield(data_sim_all.truth,'irf_full')
    data_sim_select.truth.irf_selected_vars = data_sim_all.truth.irf_full(:, var_idx); % 21 x 5
end

end
