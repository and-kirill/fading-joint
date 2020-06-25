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

function [out] = ldpc_encode(in,G,q)
% function [out] = ldpc_encode(in,H,q)
% encodes data from "in" using G over GFq
% q = 2,4,8,16,32,64,128 or 256
% Please, make sure that the data you use is valid!
% Requires matlab communication toolbox for GFq operations.
% Basically, this function performs 'in*G' over GFq. 
% this function is slow, will write a C version some time

%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   $Revision: 1.0 $  $Date: 1999/11/23 $

    [k,n] = size(G); 
    gf_G = gf(G, log2(q));
    gf_in = gf(in, log2(q));
    gf_out = gf_in*gf_G;
    out = double(gf_out.x);
end
