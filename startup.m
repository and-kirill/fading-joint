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

cd cpp;
mex decode_func_nodes_x_samples.cpp -I../common
mex initialize_func_nodes.cpp -I../common
cd ..;
cd encoder;
mex ldpc_h2g.c
cd ..;

cd binary_soft_decoder;
mex decode_soft.cpp -DLOG_Q=1 -I../common
cd ..;

addpath('alist_parser');
addpath('encoder');
addpath('gaussian_mixture');
addpath('em_decoder');
% Soft decoder is stored in the binary form
addpath('binary_soft_decoder');
% Add CPP source directory
addpath('cpp');
