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

function [ output ] = calculate_first_ply_failure( struct )
% This function takes as input a laminate, then calculates the first-ply
% failure by applying the Tsai-Hill criterion. The hipotesys is to have
% only a longitudinal stress (Nx). The algorythm is explained at pag. 75 of
% the textbook.

% STEPS 1-2: determine Qhat, and A for the laminate (they were already
% calculated).
A= struct.laminate.A;

% auxiliary constants, as described at page 75 of the book
C3=(A(3,2)/A(3,3)*A(2,1)/A(2,2) - A(3,1)/A(3,3))/(1-A(3,2)/A(3,3)*A(2,3)/A(2,2));
C2=-A(2,3)/A(2,2)*C3 - A(2,1)/A(2,2);
C1=A(1,1)+A(1,2)*C2+A(1,3)*C3;

% STEP 3: apply unit-load Nx=1 and derive the epsilon_x_o
Nx=1;
epsilon_x_o = Nx/C1;

%now, for each lamina, calculate the main stresses
for k = 1:struct.num_layers    
    norm_stresses = struct.laminate.layup(k).Q_mod * [1 0 0]';
    C1k(k)=norm_stresses(1);
    C2k(k)=norm_stresses(2);
    C3k(k)=norm_stresses(3);
    
    sigma1=C1k(k)*epsilon_x_o;
    sigma2=C2k(k)*epsilon_x_o;
    tau12=C3k(k)*epsilon_x_o;
    
    % STEP 4: by applying the Tsai-Hill criterion, determine the
    % most-stressed layer
    sigma_LU=struct.laminate.layup(k).material.sigma1;
    sigma_TU=struct.laminate.layup(k).material.sigma2;
    tau_LTU=struct.laminate.layup(k).material.tau12;
    H(k)=(sigma1/sigma_LU)^2-(sigma1/sigma_LU)*(sigma2/sigma_LU)+(sigma2/sigma_TU)^2+(tau12/tau_LTU)^2;
end

% find the maximum value of H, thus the first layer that will brake
[max_val,first_lay_brake] = max(H);

% now find the deformation that will bring to failure, by reversing the
% Tsai-Hill formula, and imposing H=1
epsilon_x_o_failure = sqrt(1/((C1k(first_lay_brake)/sigma_LU)^2-...
                            (C1k(first_lay_brake)/sigma_LU)*(C2k(first_lay_brake)/sigma_LU)+...
                            (C2k(first_lay_brake)/sigma_TU)^2+...
                            (C3k(first_lay_brake)/tau_LTU)^2));
                        
% at this point it is easy to find the corresponding longitudinal stress:
Nx = C1*epsilon_x_o_failure;

struct.sigma_f=Nx;

output=struct;
end

