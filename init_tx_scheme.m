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

function tx_scheme = init_tx_scheme(sim_settings)
% K:             the number of users
% codebook:      interleaver sesttings ('same' or 'interleaved')
% ldpc_filename: LDPC code data
codebook = sim_settings.codebook;
ldpc_filename = sim_settings.ldpc_filename;
%% Initialize code parameters
init_ldpc = @(x) decode_soft(0, x);
[h, ~] = alist2sparse(ldpc_filename);
[H, G] = ldpc_h2g(h,2);
[k, n] = size(G);

tx_scheme.code.n = n;     % Code length
tx_scheme.code.k = k;     % The number of information bits
tx_scheme.code.R = k / n; % Code rate
tx_scheme.code.G = G;     % Code generator matrix
tx_scheme.code.H = H;     % Code parity check matrix


%% Initialize LDPC decoder
[ldpc, ~, ~] = init_ldpc(ldpc_filename);
tx_scheme.ldpc = ldpc;


%% Initialize interleaver
I = eye(n);
S = zeros(n, n * sim_settings.decode.K_max);

for u = 1:sim_settings.decode.K_max
    if strcmp(codebook, 'same')
        S(:, (1:n) + (u - 1) * n) = I;
    elseif strcmp(codebook, 'interleaved')
        S(:, (1:n) + (u - 1) * n) = I(:, randperm(n));
    else
        error('Unknown code-book');
    end
end

S = sparse(S);
[T, KT] = size(S);

S_pos = zeros(T, sim_settings.decode.K_max);
S_val = zeros(T, sim_settings.decode.K_max);
S_non_zeros = S~=0;

row_weights = transpose(sum(S_non_zeros, 2));

d_max = max(row_weights);

for ii = 1:T
    idx = find(S(ii, :));
    vals = S(ii, idx);
    dc = row_weights(ii);
    idx = idx - 1;
    if dc < d_max
        idx = idx - 1 * ones(1, d_max - dc);
        vals = vals - 1 * ones(1, d_max - dc);
    end
    S_pos(ii, :) = idx;
    S_val(ii, :) = vals;
end

tx_scheme.interleaver.T = T;            % The length of interleaved message
tx_scheme.interleaver.KT = KT;          % T x (the number of users)
tx_scheme.interleaver.S = S;            % Interleaver matrix of size T x KT
tx_scheme.interleaver.rw = row_weights; % Row weights
tx_scheme.interleaver.S_pos = S_pos;    % Positions of interleaver matrix
tx_scheme.interleaver.S_val = S_val;    % Values of interleaver matrix


%% Initialize modulator
M = 2; % Use BPSK only
hMod = modem.pskmod('M', M);

% Calculate average enery per symbol
alphabet = modulate(hMod, 0:M - 1);
Es = mean(conj(alphabet) .* alphabet);

tx_scheme.modulation.M = M;      % The size of constellation
tx_scheme.modulation.mod = hMod; % Modulator
tx_scheme.modulation.Es = Es;    % Average enrgy per symbol
% Number of information bits per channel use
tx_scheme.IBPCU = sim_settings.decode.K_max * tx_scheme.code.R * n * log2(M) / T;


%% Initialize functional nodes
tx_scheme.func_nodes = initialize_func_nodes(...
    T, KT, row_weights, S_pos, S_val,...
    real(modulate(hMod, 0:M - 1))...
    );
%% Specify maximum user count that can be decoded
tx_scheme.K_max = sim_settings.decode.K_max;
tx_scheme.fading = sim_settings.fading;
end
