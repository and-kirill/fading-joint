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

function data_out = decode(K, tx_scheme, sim_settings, sigma, rx, h)
% Try multiple attempts to find the best approximation of EM decoder
decoded_words = [];
for attempt = 1:sim_settings.decode.n_try_max
    decoded_attempt = decode_single(tx_scheme, sim_settings, sigma, rx, h);
    n_words = size(decoded_attempt.est_iwd_success, 1);
    for i = 1:n_words
        iwd = decoded_attempt.est_iwd_success(i, :);
        if rank([decoded_words; iwd]) > rank(decoded_words)
            decoded_words = [decoded_words; iwd];
        end
    end
    % This check is added to speed-up the FER estimation. Thus, true user
    % count is used here. This does not affect decoder performance
    if rank(decoded_words) == K
        break
    end
end
data_out.est_iwd = decoded_words;
