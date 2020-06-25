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

function n_errors = estimate_fer(decoded, K, iwd)
% Exctract information and code words from decoded data
est_iwd = decoded.est_iwd;

K_total = size(est_iwd, 1);
pos = -1 * ones(1, K_total);
for ll = 1:K_total
    for qq = 1:K
        if size(intersect(est_iwd(ll, :), iwd(qq,:), 'rows'), 1) == 1
            pos(ll) = qq;
            break;
        end
    end
end

n_errors = length(setdiff(1:K, pos));

end
