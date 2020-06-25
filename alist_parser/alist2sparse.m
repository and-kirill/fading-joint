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
 
function [H, Q] = alist2sparse(file_name)
    bBinary = 0;
    
    file = fopen(file_name, 'r');
    
    % N M [Q]
    string = fgets(file);
    [temp, count] = sscanf(string, '%d');
    if (count ~= 2) && (count ~= 3)
        disp('Error: bad file format');
    end
    N = temp(1);
    M = temp(2);
    if count == 2
        disp('Binary matrix');
        bBinary = 1;
        Q = 2;
    else
        disp('Non-binary matrix');
        Q = temp(3);
    end
    
    % c_max r_max
    string = fgets(file);
    [temp, count] = sscanf(string, '%d');
    if (count ~= 2)
        disp('Error: bad file format');
    end
    c_max = temp(1);
    r_max = temp(2);
    
    
    % column_weights
    string = fgets(file);
    [column_weights, count] = sscanf(string, '%d');
    if (count ~= N)
        disp('Error: bad file format');
    end
    
    % row_weights
    string = fgets(file);
    [row_weights, count] = sscanf(string, '%d');
    if (count ~= M)
        disp('Error: bad file format');
    end
    
    H = sparse([], [], [], M, N);
    
    for i = 1:N
        string = fgets(file);
        [temp, count] = sscanf(string, '%d');
        if bBinary
            if count ~= c_max
                disp('Error: bad file format');
            end
            indexes = temp(1:column_weights(i));
            H(indexes, i) = 1;
        else
            if count ~= 2*c_max
                disp('Error: bad file format');
            end
            temp_matrix = reshape(temp, 2, c_max);
            indexes = temp_matrix(1, 1:column_weights(i));
            values = temp_matrix(2, 1:column_weights(i));
            H(indexes, i) = values;
        end
    end
    fclose(file);
end
