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

function result = make_llr_prod(gm, llr)
    assert(~sum(isnan(gm.c)));
    % TODO: the following scheme works only for BPSK
    % Must be fixed for QPSK case
    result.mu = [-gm.mu, gm.mu];
    result.var = [gm.var, gm.var];

    p_plus_1 = safe_exp(llr / 2);
    p_minus_1 = safe_exp(-llr / 2);
    s = p_plus_1 + p_minus_1;
    assert(s > eps);
    p_plus_1 = p_plus_1 / s;
    p_minus_1 = p_minus_1 / s;

    result.c = [p_minus_1 * gm.c, p_plus_1 * gm.c];
    assert(~sum(isnan(result.c)));
end
