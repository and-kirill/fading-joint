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

clc;
clear;
%% General scenario parameters
codebook = 'same';
K = 4; % The number of users
snr_array = 0:20;
ldpc_filename = 'LDPC_k91_n364.alist';
min_frame_errors = 1000;

n_workers = 6;
parpool(n_workers);

%% Create simulation parameters
sim_settings.gm.n_samples =            20;
% Cumulative probability during GM truncation
sim_settings.gm.prune_cdf =            0.999;
% Maximum distance between two GM components at merge stage
% Minimum allowed distance for merge
sim_settings.gm.merge_dist_min =       1;
% Distance decay. Use exponential decay with the number of iterations
sim_settings.gm.merge_dist_decay =     0.5;
% Enable/disbale GM merge procedure
sim_settings.gm.use_merge =            true;
% Allow to merge gaussian mixture components with opposite signs
sim_settings.gm.merge_opposite_signs = true;
% Print merge parameters
sim_settings.gm.enable_debug =         false;
% Maximum GM components allowed in gaussian mixture
sim_settings.gm.prune_max_n =          500;
% Use sampled channel estimates from the gaussian mixture or use true channel estimates
sim_settings.gm.use_sampling =         true;
% Minimum variance allowed
sim_settings.gm.regularization =       0;
% Maximum user count than can be detected simultaneously
sim_settings.decode.K_max =            min(K, 4);
% The number of functional nodes.
sim_settings.decode.n_func_nodes =     50;
% The number of LDPC decoding iterations
sim_settings.decode.n_iter_ldpc =      50;
% The number of outer iterations
% Minimum iterations before early stopping
sim_settings.decode.n_iter_out_min =   5;
% Maximum outer iterations
sim_settings.decode.n_iter_out_max =   25;
% Maximum EM independent attempts
sim_settings.decode.n_try_max =        20;
% Fading types: 1) 'rayleigh' fading. Real and imaginary parts of h are
% sampled from gaussian distribution or 'phase' |h| = 1 with uniformly
% distributed phase)
sim_settings.fading =                  'rayleigh';
% Codebook: only same is supported
sim_settings.codebook =                codebook;
sim_settings.user_count =              K;
sim_settings.ldpc_filename =           preprocess_ldpc(sprintf('codes/%s', ldpc_filename));

%% Run decoder
fer = bpsk_noma_ldpc_sequential(sim_settings, snr_array, min_frame_errors);

filename = sprintf('data/out_parallel_ldpc_K=%d_code=%s_ldpc=%s_iter=%d.mat',...
    K, codebook,...
    ldpc_filename, sim_settings.decode.n_iter_out_max...
    );
data.snr_array =snr_array;
data.fer = fer;
save(filename, 'data');
