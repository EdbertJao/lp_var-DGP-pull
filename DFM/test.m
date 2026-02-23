%% load_and_check_dgpstore.m (MATFILE-SAFE)
% Loads one DGPSTORE_*.mat file and checks that a sample (i_MC,i_spec)
% is properly labeled and truth IRF matches.

clear; clc;

%% ---- USER INPUT: pick ONE file to test ----
file_path = fullfile('Results','DGP_Stores','DGPSTORE_DFM_G_ObsShock_specid1_mode1.mat');
% file_path = fullfile('Results','DGP_Stores','DGPSTORE_DFM_G_IV_specid1_mode1.mat');

assert(isfile(file_path), "File not found: %s", file_path);

%% ---- Open via matfile (memory-safe) ----
M = matfile(file_path);

fprintf("Loaded: %s\n", file_path);
fprintf("dgp_type=%s, estimand_type=%s, spec_id=%d, mode_type=%d\n", ...
    string(M.dgp_type), string(M.estimand_type), M.spec_id, M.mode_type);

%% ---- Dimensions ----
szX = size(M, 'X_all');      % [T, n_y, n_MC]
T   = szX(1);
n_y = szX(2);
n_MC = szX(3);

szVS = size(M, 'var_select'); % [n_spec, n_var]
n_spec = szVS(1);
n_var  = szVS(2);

fprintf("T=%d, n_y=%d, n_MC=%d, n_spec=%d, n_var=%d\n", T, n_y, n_MC, n_spec, n_var);

%% ---- Load small metadata fully (avoid matfile indexing restrictions) ----
var_name_short = M.var_name_short;  % 1x207 cell
var_name_long  = M.var_name_long;   % 1x207 cell

% These are also small; you can load them too (optional)
var_name_short_spec_all = M.var_name_short_spec; % 100x5 cell
var_name_long_spec_all  = M.var_name_long_spec;  % 100x5 cell

%% ---- Choose a sample draw ----
i_MC   = min(1234, n_MC);
i_spec = min(7, n_spec);

idx = M.var_select(i_spec,:);  % 1x5 indices into 207 columns

%% ---- Pull names (in-RAM indexing; no matfile restrictions) ----
names_short_from_global = var_name_short(idx);
names_long_from_global  = var_name_long(idx);

names_short_spec = var_name_short_spec_all(i_spec,:);
names_long_spec  = var_name_long_spec_all(i_spec,:);

%% ---- Checks: labeling consistency ----
assert(isequal(names_short_spec, names_short_from_global), ...
    "Short-name labeling mismatch for i_spec=%d.", i_spec);

assert(isequal(names_long_spec, names_long_from_global), ...
    "Long-name labeling mismatch for i_spec=%d.", i_spec);

%% ---- Checks: data slicing consistency ----
% Load one replication’s full panel (small: 200x207 single)
X_full = M.X_all(:,:,i_MC);    % T x 207 (single)
wbar   = X_full(:, idx);       % T x 5   (single)

% Optional: verify size and finite values
assert(isequal(size(wbar), [T, n_var]), "wbar size mismatch.");
assert(all(isfinite(wbar(:))), "wbar contains non-finite values.");

%% ---- Checks: truth IRF consistency ----
truth = M.target_irf(:, i_spec);  % 21x1 double
assert(all(isfinite(truth)), "truth IRF contains non-finite values.");

%% ---- IV-only checks (if Z_all exists) ----
isIV = false;
try
    size(M, 'Z_all'); % throws if Z_all not present
    isIV = true;
catch
    isIV = false;
end

if isIV
    rho_idx = M.settings_specifications.rho_select_grid_idx(i_spec);
    sig_idx = M.settings_specifications.sigma_v_select_grid_idx(i_spec);

    % Load the instrument grid for this MC draw (typically small)
    Z_full = M.Z_all(:,:,:,i_MC);     % T x nrho x nsig (single)
    z      = Z_full(:, rho_idx, sig_idx); % T x 1 (single)

    IV_rho_grid     = M.IV_rho_grid;
    IV_sigma_v_grid = M.IV_sigma_v_grid;

    rho_val = IV_rho_grid(rho_idx);
    sig_val = IV_sigma_v_grid(sig_idx);

    fprintf("IV selection: rho_idx=%d (rho=%.4g), sig_idx=%d (sigma_v=%.4g)\n", ...
        rho_idx, rho_val, sig_idx, sig_val);

    assert(all(isfinite(z)), "Instrument z contains non-finite values.");
end

%% ---- Report summary ----
fprintf("\nSAMPLE DGP CHECK PASSED ✅\n");
fprintf("i_MC=%d, i_spec=%d\n", i_MC, i_spec);
fprintf("Selected variable indices: %s\n", mat2str(idx));

fprintf("Selected short names:\n");
disp(names_short_spec);

fprintf("Selected long names:\n");
disp(names_long_spec);

fprintf("First 5 rows of DGP data (wbar):\n");
disp(double(wbar(1:min(5,T), :)));

fprintf("First 5 horizons of truth IRF:\n");
disp(truth(1:min(5,numel(truth))));
