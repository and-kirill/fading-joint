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

function result = gm_product(gm1, gm2)
    assert(~sum(isnan(gm1.c)));
    assert(~sum(isnan(gm2.c)));
    result = struct('mu', [], 'var', [], 'c', []);
    logc = [];
    for ii = 1:length(gm1.mu)
        for jj = 1:length(gm2.mu)
            mu1 = gm1.mu(ii);
            var1 = gm1.var(ii);
            c1 = gm1.c(ii);
            mu2 = gm2.mu(jj);
            var2 = gm2.var(jj);
            c2 = gm2.c(jj);

            var = var1*var2/(var1+var2);
            mu = var*( mu1/var1 + mu2/var2);
            logc = [logc, log(c1*c2/sqrt(2*pi*(var1+var2)))-(mu1-mu2)^2 / 2 / (var1+var2)];

            result.mu = [result.mu mu];
            result.var = [result.var var];
        end
    end
    result.c = exp(logc - max(logc));
    result.c = result.c./sum(result.c);
    assert(~sum(isnan(result.c)));
end
