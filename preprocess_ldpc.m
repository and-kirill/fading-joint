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

function filename_out = preprocess_ldpc(filename_in)
filename_out = strcat(filename_in, '.ldpc');
if isfile(filename_out)
    fprintf('LDPC filename %s already exists\n', filename_out)
else
    [h, ~] = alist2sparse(filename_in);
    [H, ~] = ldpc_h2g(h,2);
    sparse2alist(2, H, filename_out);
end
end