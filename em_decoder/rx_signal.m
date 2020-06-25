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

function [iwd, rx_noiseless, rx, h] = rx_signal(K, tx_scheme, sigma)
% generate information word
iwd = randi(tx_scheme.modulation.M, K, tx_scheme.code.k) - 1;

% encode by LDPC code and modulate
cwd = ldpc_encode(iwd, tx_scheme.code.G, 2); % binary code
mod_cwd = modulate(tx_scheme.modulation.mod, cwd);

% Fading
if strcmp(tx_scheme.fading, 'rayleigh')
    h = sqrt(1 / 2) * (randn(1, K) + 1i * randn(1, K));
elseif strcmp(tx_scheme.fading, 'phase')
    phi = 2 * pi * rand(1, K);
    h = cos(phi) + 1i * sin(phi);
else
    error('Unknown fading type');
end

rx_noiseless = h * mod_cwd;

% add noise
noise_vector = sigma * (...
    randn(1, tx_scheme.interleaver.T) +...
    1i * randn(1, tx_scheme.interleaver.T)...
    );

rx = rx_noiseless + noise_vector;
end
