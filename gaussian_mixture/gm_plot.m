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

function [] = gm_plot(gm, params)
    x = -10:0.0001:10;
    y = zeros(1, length(x));
    for i = 1:length(gm.mu)
        component_pdf = gauss_pdf(x, gm.mu(i), gm.var(i)) * gm.c(i);
        plot(x, component_pdf, params);
        hold on
        y = y + component_pdf;
    end
    plot(x, y, [params, '--']);
    hold on;
end

function y = gauss_pdf(x, mu, sigma2)
    y = 1/sqrt(2 * pi * sigma2) * exp(-(x - mu).^2 / (2 * sigma2));
end
