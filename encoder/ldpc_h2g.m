% This file is part of the faing-joint distribution
% https://github.com/and-kirill/fading-joint.
% Copyright (c) 2019 Alexey Frolov al.frolov@skoltech.ru.
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

function [h, g] = ldpc_h2g(H,q)
% [h, g] = ldpc_h2g(H,q)
% converts tentative binary LDPC matrix H into a new matrix h
% (columns are permuted) and produces the generator matrix g
% H should be a sparse matrix in MATLAB format.
% q - Field base (power of 2) now only 2 4 8 1 32 64 128 and 256
%
% MEX file


%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   $Revision: 1.1 $  $Date: 1999/08/23 $ - implementation for GFq
