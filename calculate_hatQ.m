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

function [ out ] = calculate_hatQ( struct )
Q=struct.material.Q;
theta=struct.angle;

T = [ cos(theta)^2              sin(theta)^2            2*sin(theta)*cos(theta);
      sin(theta)^2              cos(theta)^2            -2*sin(theta)*cos(theta);
      -sin(theta)*cos(theta)    sin(theta)*cos(theta)   cos(theta)^2-sin(theta)^2];
R = [ 1 0 0;
      0 1 0;
      0 0 2];
Q_mod = Q*R*T*inv(R);
hatQ = inv(T)*Q_mod;

struct.Q_mod=Q_mod;
struct.Q_hat=hatQ;

out=struct;
end

