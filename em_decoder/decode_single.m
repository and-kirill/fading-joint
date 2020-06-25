% This file is part of the faing-joint distribution
% https://github.com/and-kirill/fading-joint.
% Copyright (c) 2019 Kirill Andreev k.andreev@skoltech.ru.
% 
% This program is free software: you can redistribute it and/or modify  
% it under the terms of the GNU General Public License as published by  
% the Free Software Foundation, version 3.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
 
% You should have received a copy of the GNU General Public License 
% along with this program. If not, see <http://www.gnu.org/licenses/>.

function data_out = decode_single(tx_scheme, sim_settings, sigma, rx, h)
K = tx_scheme.K_max;
var_msg = zeros(1, tx_scheme.code.n * K); % Start with zero input LLRs

% Prior distribution for fading coefficients
h_prior_re = repmat(struct('mu', 0, 'var', 100, 'c', 1), 1, K);
h_prior_im = repmat(struct('mu', 0, 'var', 100, 'c', 1), 1, K);

% Q messages (for fading coefficients)
h_Q_re = repmat(h_prior_re, K * tx_scheme.code.n, 1);
h_Q_im = repmat(h_prior_im, K * tx_scheme.code.n, 1);

% R messages (for fading coefficients)
h_R_re = h_Q_re;
h_R_im = h_Q_im;
% Channel estimate:
h_hat = repmat(...
    struct(...
    're', struct('mu', [], 'var', [], 'c', []),...
    'im', struct('mu', [], 'var', [], 'c', [])...
    ), sim_settings.decode.n_iter_out_max, K...
    );

n = tx_scheme.code.n; % The code length
k = tx_scheme.code.k; % Information word size
% Outer iterations
est_iwd = zeros(K, k);
est_cwd = zeros(K, n);
for it = 1:sim_settings.decode.n_iter_out_max
    sim_settings.gm.merge_dist = sim_settings.gm.merge_dist_min +  sim_settings.gm.merge_dist_decay ^ (it - 1);
    for user = 1:K
        % External iteration (M-step???):
        func_nodes_idx = randperm(tx_scheme.interleaver.T, sim_settings.decode.n_func_nodes);
        [h_R_re, rows_idx] = update_h_R(h_R_re, h_Q_re, real(rx), func_nodes_idx, var_msg, sigma, tx_scheme, sim_settings, user);
        [h_R_im, ~]        = update_h_R(h_R_im, h_Q_im, imag(rx), func_nodes_idx, var_msg, sigma, tx_scheme, sim_settings, user);
        
        [h_Q_re, h_hat_re] = update_h_Q(h_Q_re, h_R_re, h_prior_re(user), rows_idx, tx_scheme, sim_settings, user);
        [h_Q_im, h_hat_im] = update_h_Q(h_Q_im, h_R_im, h_prior_im(user), rows_idx, tx_scheme, sim_settings, user);
        % Full GM at every iteration is provided as channel estimate output
        h_hat(it, user).re = h_hat_re;
        h_hat(it, user).im = h_hat_im;
        
        h_samples_re = sample_h(h_Q_re, real(h), tx_scheme, sim_settings);
        h_samples_im = sample_h(h_Q_im, imag(h), tx_scheme, sim_settings);
        
        % Perform internal decoding iterations
        func_msg = decode_func_nodes_x_samples(...
            tx_scheme.func_nodes, var_msg, rx, sigma,...
            sim_settings.gm.n_samples, h_samples_re, h_samples_im...
            );
        assert(~sum(isnan(func_msg)));
        
        func_msg((user - 1) * n + 1:user * n) = ...
            func_msg((user - 1) * n + 1:user * n) - var_msg((user - 1) * n + 1:user * n);
        assert(~sum(isnan(func_msg)));
        
        % Decode LDPC code
        [~, ~, est_user_cwd, out_user_llr] = decode_soft(...
            1, tx_scheme.ldpc, func_msg((user - 1) * n + 1:user * n),...
            sim_settings.decode.n_iter_ldpc...
            );
        est_cwd(user, :) = est_user_cwd;
        est_iwd(user, :) = est_user_cwd(n - k + 1:end);
        
        var_msg((user - 1) * n + 1:user * n) = ...
            out_user_llr - func_msg((user - 1) * n + 1:user * n); % extrinsic information
        assert(~sum(isnan(var_msg)));
    end
    syndromes = mod(full(tx_scheme.code.H * est_cwd'), 2);
    % Check that no zero codeword during early stopping
    if it > sim_settings.decode.n_iter_out_min && sum(syndromes(:)) == 0 && sum(sum(est_cwd, 2) == 0) == 0
        break
    end
    
end
% Estimate the number of successfully detected different code-words
cwd_syndrome = sum(syndromes, 1);
words_decoded = sum(cwd_syndrome == 0);
K_detected = min(rank(est_cwd), words_decoded);

data_out.syndromes = syndromes;
data_out.n_iter_out = it;
data_out.K_detected = K_detected;
data_out.est_iwd = est_iwd;
data_out.est_cwd = est_cwd;
data_out.est_iwd_success = est_iwd(cwd_syndrome == 0, :);
data_out.h_hat = h_hat;
end


function [h_R, rows_idx] = update_h_R(h_R, h_Q, rx, func_nodes_idx, var_msg, sigma, tx_scheme, sim_settings, user)
rows_idx = zeros(1, length(func_nodes_idx));
for i = 1:length(func_nodes_idx)
    node = func_nodes_idx(i);
    pos_user = -1;
    temp_gm = struct('mu', 0, 'var', sigma^2, 'c', 1); % start with z
    for r = 1:tx_scheme.interleaver.rw(node)
        pos = tx_scheme.interleaver.S_pos(node, r) + 1;
        if floor((pos - 1) / tx_scheme.code.n) + 1 ~= user
            temp_gm = gm_convolution(temp_gm, make_llr_prod(h_Q(pos), -var_msg(pos)));
            temp_gm = gm_prune(temp_gm, sim_settings.gm);
        else
            pos_user = pos;
        end
    end
    rows_idx(i) = pos_user;
    h_R(pos_user) = temp_gm;
    h_R(pos_user).mu = h_R(pos_user).mu + rx(node);
    h_R(pos_user) = make_llr_prod(h_R(pos_user), var_msg(pos_user));
    h_R(pos_user) = gm_prune(h_R(pos_user), sim_settings.gm);
end
end


function [h_Q, h_hat] = update_h_Q(h_Q, h_R, h_prior, rows_idx, tx_scheme, sim_settings, user)
n = tx_scheme.code.n; % Code length
h_Q_update = h_prior;
for row = rows_idx
    assert(~sum(isnan(h_R(row).c)));
    h_Q_update = gm_product(h_Q_update, h_R(row));
    h_Q_update = gm_prune(h_Q_update, sim_settings.gm);
end

h_Q(((user - 1) * n + 1):1:(user * n)) = repmat(h_Q_update, n, 1);
h_hat = h_Q_update;
end


function h_samples = sample_h(h_Q, h, tx_scheme, sim_settings)
if sim_settings.gm.use_sampling
    % Sample from gaussian mixture
    h_samples = gm_sample(h_Q(1:tx_scheme.code.n:end), sim_settings.gm);
else
    % Use true channel estimation
    h_samples = repmat(h', 1, sim_settings.gm.n_samples);
end
assert(~sum(isnan(h_samples(:))));
end
