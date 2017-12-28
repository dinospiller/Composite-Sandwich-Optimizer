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

function [ out ] = calculate_laminate_ABD( struct )
new_layup = struct.layup;

%calculate medium-plane:
tot_thick = 0;
tot_dens = 0;

for ind = 1:size(new_layup)
    h(ind) = tot_thick + new_layup(ind).thick; % this vector holds the vertical position of the interlaminate,referred to the top
    tot_thick = tot_thick + new_layup(ind).thick; 
    tot_dens = tot_dens + new_layup(ind).material.dens*new_layup(ind).thick;     
end

struct.thick=tot_thick;
struct.dens = tot_dens/tot_thick;
laminate_center=tot_thick/2;

A=zeros(3);
B=zeros(3);
D=zeros(3);

%Calculate the matrices
for ind = 1:size(new_layup)
    a = -laminate_center+h(ind);% z(k)
    if(ind==1)
        b = -laminate_center;% z(k-1)
    else
        b = -laminate_center+h(ind-1);
    end
    
    A = A+new_layup(ind).Q_hat*(a-b);
    B = B+1/2*new_layup(ind).Q_hat*(a^2-b^2);
    D = D+1/3*new_layup(ind).Q_hat*(a^3-b^3);
    
end

% combine the A,B,D matrices
K=[A B ; B D];
K_hat=inv(K);

% Apparent properties of the laminate

%correct expressions
Exx=1/(tot_thick*K_hat(1,1)); %[MPa] longitudinal elastic modulus
Eyy=1/(tot_thick*K_hat(2,2)); %[MPa] transversal elastic modulus
Gxy=1/(tot_thick*K_hat(3,3)); %[MPa] shear elastic modulus

%expressions used in CLT App
A_star=inv(A);
% Exx=1/(tot_thick*A_star(1,1)); %[MPa] longitudinal elastic modulus
% Eyy=1/(tot_thick*A_star(2,2)); %[MPa] transversal elastic modulus
% Gxy=1/(tot_thick*A_star(3,3)); %[MPa] shear elastic modulus

vxy = A(1,2)/A(2,2);
vyx = A(1,2)/A(1,1);

struct.A=A;
struct.B=B;
struct.D=D;
struct.Exx=Exx;
struct.Eyy=Eyy;
struct.Gxy=Gxy;
struct.vxy=vxy;
struct.vyx=vyx;

out=struct;
end

