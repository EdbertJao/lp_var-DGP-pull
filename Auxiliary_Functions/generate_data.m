function data_sim = generate_data(model,settings)
% Function for generating simulated data
    % Use a general ABCD representation of the encompassing model (DFM, DSGE or others):
        % state transition:  s_t = A * s_{t-1} + B * e_t
        % measurement eq:    y_t = C * s_{t-1} + D * e_t

    % Warning: shock_weight' * e_t corresponds to true shock \epsilon_{1t} in our paper
    
    % external IV: z_t = \rho * z_{t-1} + \alpha * (shock_weight' * \epsilon_t) + \nu_t

% unpack settings

T      = settings.simul.T;
T_burn = settings.simul.T_burn;

A = model.ABCD.A;
B = model.ABCD.B;
C = model.ABCD.C;
D = model.ABCD.D;

[n_s,n_e] = size(B);

with_IV = settings.est.with_IV;

if with_IV == 1
    rho_grid      = model.IV.rho_grid;
    alpha         = model.IV.alpha;
    sigma_v_grid  = model.IV.sigma_v_grid;
else % meaningless placeholders in the case of no IV
    rho_grid = 0.1;
    alpha = 1;
    sigma_v_grid = 1;
end

shock_weight = settings.est.shock_weight;

% draw shocks

data_e = randn(T_burn+T,n_e);

% simulate states & measurement error

s = zeros(n_s,1);
data_s = NaN(T_burn + T, n_s);
for t = 1:(T_burn + T)
    s = A * s + B * data_e(t,:)';
    data_s(t,:) = s';
end

% simulate observables

data_y = data_s(T_burn:end-1,:)*C' + data_e(T_burn+1:end,:)*D';

% simulate IV
z = NaN(T+T_burn, length(rho_grid), length(sigma_v_grid)); % multiple IV persistence/strength setups
for idx = 1:length(rho_grid)
    for iidx = 1:length(sigma_v_grid)
        z(:,idx,iidx) = filter(1, [1 -rho_grid(idx)], alpha * data_e * shock_weight + sigma_v_grid(iidx) * randn(T+T_burn,1));
    end
end

data_z = z(T_burn+1:end,:,:);

% collect results and shift timing

data_sim.data_y     = data_y;
data_sim.data_shock = data_e((T_burn+1):(T_burn+T), :) * shock_weight;
data_sim.data_z     = data_z;

% =========================
% NEW: metadata + truth
% =========================
data_sim.meta = struct();

% Variable name metadata (aligned with columns of data_y)
if isfield(model,'variable_name_short'); data_sim.meta.variable_name_short = model.variable_name_short; end
if isfield(model,'variable_name_long');  data_sim.meta.variable_name_long  = model.variable_name_long;  end
if isfield(model,'variable_name_code');  data_sim.meta.variable_name_code  = model.variable_name_code;  end
if isfield(model,'trans_code');          data_sim.meta.trans_code          = model.trans_code;          end

% IV grid values (useful for later inspection)
data_sim.meta.with_IV = settings.est.with_IV;
if settings.est.with_IV == 1
    data_sim.meta.IV_rho_grid     = model.IV.rho_grid;
    data_sim.meta.IV_sigma_v_grid = model.IV.sigma_v_grid;
end

% Truth objects (small; convenient to carry along)
data_sim.truth = struct();
if isfield(model,'irf');        data_sim.truth.irf_full   = model.irf;        end   % 21 x 207
if isfield(model,'target_irf'); data_sim.truth.target_irf = model.target_irf; end   % 21 x n_spec

end