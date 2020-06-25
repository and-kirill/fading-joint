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

function result = gm_sample(gm_array, prune_settings)
    result = zeros(length(gm_array), prune_settings.n_samples);

    for ii = 1:length(gm_array)
        qqq = gm_prune(gm_array(ii), prune_settings); %%% TODO
        h_samples = gm_sample_single(qqq, prune_settings.n_samples);
        result(ii, :) = h_samples;
    end
end
