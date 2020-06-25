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

function result = gm_prune(gm, gm_settings)
    [~, idx] = sort(gm.c, 'descend');
    N = 1;
    gm_cdf = cumsum(gm.c(idx));
    % TODO: binary search is possible here!
    while gm_cdf(N) < gm_settings.prune_cdf && N < length(idx)
        N = N + 1;
    end
    idx = idx(1:N);
    gm = struct(...
        'mu',  gm.mu(idx),...
        'var', gm.var(idx),...
        'c', gm.c(idx) / sum(gm.c(idx))...
        );
    if gm_settings.use_merge
        gm = gm_merge(gm, gm_settings);
    end
    % Print warning whether the gaussian mixture is significantly pruned
    if gm_settings.enable_debug && length(gm.mu) > gm_settings.prune_max_n
        prob_res = sum(gm.c(gm_settings.prune_max_n:end));
        if prob_res > 1e-2
            fprintf('W: pruned %2.3f\n', prob_res);
        end
    end
    [~, idx] = sort(gm.c, 'descend');
    idx = idx(1:min(length(idx), gm_settings.prune_max_n));
    result = struct(...
        'mu', gm.mu(idx),...
        'var', gm.var(idx),...
        'c', gm.c(idx)/sum(gm.c(idx))...
        );
    % Apply regularization
    result.var = max(result.var, gm_settings.regularization);
end
