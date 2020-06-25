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

function h = gm_sample_single(gm, n_samples)
gauss_samples = randn(1, n_samples);

component_idx = discretesample(gm.c, n_samples);
h = gm.mu(component_idx) + gauss_samples .* gm.var(component_idx);
end

function x = discretesample(p, n)
assert(isfloat(p), 'discretesample:invalidarg', ...
    'p should be an array with floating-point value type.');

assert(isnumeric(n) && isscalar(n) && n >= 0 && n == fix(n), ...
    'discretesample:invalidarg', ...
    'n should be a nonnegative integer scalar.');

% process p if necessary

K = numel(p);
if ~isequal(size(p), [1, K])
    p = reshape(p, [1, K]);
end

% construct the bins

edges = [0, cumsum(p)];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

% draw bins

rv = rand(1, n);
c = histc(rv, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

% extract samples

xv = find(c);

if numel(xv) == n  % each value is sampled at most once
    x = xv;
else                % some values are sampled more than once
    xc = c(xv);
    d = zeros(1, n);
    dv = [xv(1), diff(xv)];
    dp = [1, 1 + cumsum(xc(1:end-1))];
    d(dp) = dv;
    x = cumsum(d);
end

% randomly permute the sample's order
x = x(randperm(n));
end
