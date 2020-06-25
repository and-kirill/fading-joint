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

function result = gm_convolution(gm1, gm2)
    difference = 1e-2;
    result = struct('mu', [], 'var', [], 'c', []);
    for ii = 1:length(gm1.mu)
        for jj = 1:length(gm2.mu)
            c1 = gm1.c(ii);
            c2 = gm2.c(jj);
            if c1 * c2 < 1e-6
                continue
            end
            mu1 = gm1.mu(ii);
            var1 = gm1.var(ii);

            mu2 = gm2.mu(jj);
            var2 = gm2.var(jj);


            mu = mu1 + mu2;
            var = var1 + var2;
            c =  c1 * c2;

            ind = find(abs(result.mu - mu) < difference, 1, 'first');

            if (~isempty(ind)) && (abs(result.var(ind) - var) < difference)
                result.c(ind) = result.c(ind) + c;
            else
                result.mu = [result.mu mu];
                result.var = [result.var var];
                result.c = [result.c c];
            end
        end
    end
end
