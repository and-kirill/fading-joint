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

function result = gm_merge(gm, gm_settings)
if length(gm.mu) == 1
    result = gm;
    return;
end
result.mu = [];
result.var = [];
result.c = [];

N = length(gm.c);

while N > 0
    [~, i_max] = max(gm.c);
    all_idx = 1:N;
    distance_vector = gm_distance(gm.mu(i_max), gm.mu, gm.var);
    
    if gm_settings.merge_opposite_signs
        merge_idx = distance_vector < gm_settings.merge_dist;
    else
        merge_idx = (distance_vector < gm_settings.merge_dist) .* (sign(gm.mu(i_max) * gm.mu) ~= -1);
    end
    assert (sum(merge_idx) > 0);
    merge_list = all_idx(merge_idx == 1);
    [mu, var, c] = do_merge(struct(...
        'mu', gm.mu(merge_list),...
        'var', gm.var(merge_list),...
        'c', gm.c(merge_list)...
        ));
    result.mu = [result.mu, mu];
    result.var = [result.var, var];
    result.c = [result.c, c];
    % Clear merge list
    gm.mu(merge_list) = [];
    gm.var(merge_list) = [];
    gm.c(merge_list) = [];
    N = length(gm.c);
end
end

function [mu, var, c] = do_merge(gm)
c = sum(gm.c);
mu = sum(gm.c .* gm.mu) / c;
var = sum(gm.c .* ((gm.mu - mu).^2 + gm.var)) / c;
end


function d = gm_distance(mu1, mu2, var2)
% The first component is assumed to have higher weight
d = (mu1 - mu2).^2 ./ var2;
end
