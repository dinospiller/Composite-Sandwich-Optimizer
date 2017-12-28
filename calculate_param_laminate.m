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

function [ output ] = calculate_param_laminate( struct , tot_thick )
% This function takes as input a param_laminate struct, which is a 
% laminate where the thickness is not previously-defined, but put as a
% parameter to this function. The layers are supposed of fixed-thickness,
% so the total thickness input will be distributed on the layers.

% calculate per-layer thickness (the whole laminate thickness will be
% calculated afterwards)
layer_num = size(struct.layup);
struct.num_layers=layer_num(1);
struct.layer_thick=tot_thick/struct.num_layers;

% keep the original laminate and waste all layers in the layup, but take
% layer 1 as "form" for the other layers.
curr_layer = struct.laminate.layup(1);
struct.laminate.layup=[curr_layer]; % delete the current layup putting an empty vector

% now the layup of the laminate will be rebuild layer-by-layer
for ind = 1:struct.num_layers
    curr_layer.material=struct.layup(ind).material;
    curr_layer.angle=struct.layup(ind).angle;
    curr_layer.thick=struct.layer_thick;
    
    struct.laminate.layup(ind,1) = calculate_hatQ(curr_layer);    
end

% calculate the whole laminate parameters
new_lamin = calculate_laminate_ABD(struct.laminate);

struct.laminate=new_lamin;

%now that the parametric laminate is defined, calculate the
%first-ply-failure stress
struct= calculate_first_ply_failure(struct);

output=struct;
end

