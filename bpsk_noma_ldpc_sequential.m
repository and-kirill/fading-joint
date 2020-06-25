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

function fer = bpsk_noma_ldpc_sequential(sim_settings_, snr_array, min_frame_errors)
fer = zeros(1, length(snr_array));
parfor ii = 1:length(snr_array)
    sim_settings = sim_settings_;
    K = sim_settings.user_count;
    % Generate transmission scheme
    tx_scheme = init_tx_scheme(sim_settings);
    sim_settings.decode.n_func_nodes = min(sim_settings.decode.n_func_nodes, tx_scheme.code.n);

    snr = snr_array(ii);
    
    N0 = tx_scheme.modulation.Es * 10^(-snr / 10);
    sigma = sqrt(N0/2);
    tests = 0;
    wrong_dec = 0;
    while wrong_dec  < min_frame_errors
        tests = tests + 1;
        [iwd, ~, rx, h] = rx_signal(K, tx_scheme, sigma);
        min_sigma_value = 0.08;
        decoded = decode(K, tx_scheme, sim_settings, max(sigma, min_sigma_value), rx, h);
        wrong_rows = estimate_fer(decoded, K, iwd);

        if wrong_rows ~= 0
            wrong_dec = wrong_dec + wrong_rows;
            fer(ii) = wrong_dec / K / tests;
            fprintf('SNR: %2.2f\tfer = %f\t tests passed: %d\n', snr, wrong_dec / K / tests, tests);
        end
    end
end
end