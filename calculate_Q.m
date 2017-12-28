%     Copyright 2017-2020 Dino Spiller (dinospiller@gmail.com)
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ struct_out ] = calculate_Q( struct )
%calculates the Q matrix for a composite layer
v21 = struct.v12 * struct.E22 / struct.E11;
E1= struct.E11;
E2= struct.E22;
v12= struct.v12;
G12= struct.G12;
Q = [ E1/(1-v12*v21)        v12*E2/(1-v12*v21)  0;...
            v21*E1/(1-v12*v21)     E2/(1-v12*v21)      0;...
            0                      0                   G12];
struct.Q=Q;
struct_out=struct;
end

