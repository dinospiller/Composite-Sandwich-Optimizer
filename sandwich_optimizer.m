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

%% optimizer for composite sandwich
% author: Dino Spiller
% date: october 2016
clear all;
close all;


%% Material description

% carbon unidirectional
% Skin material (UD std CF, ACPcomposites)
% create astructure of properties
% carbon = struct('E11',135000,...    %[MPa] tensile modulus
%                  'E22',10000,...    %[MPa]
%                  'G12',5000,...     %[MPa]
%                  'v12',0.3,...
%                  'sigma1',1500,...   %[MPa] tensile strength
%                  'sigma2',50,...   %[MPa]
%                  'tau12',70,...    %[MPa] shear strength
%                  'dens',1.6,...     %[kg/dm^3]
%                  'Q',[]);           %always put an empty vector

% % carbon fabric
carbon = struct('E11',75000,...    %[MPa] tensile modulus
                 'E22',75000,...    %[MPa]
                 'G12',5000,...     %[MPa]
                 'v12',0.1,...
                 'sigma1',600,...   %[MPa] tensile strength
                 'sigma2',600,...   %[MPa]
                 'tau12',90,...    %[MPa] shear strength
                 'dens',1.6,...     %[kg/dm^3]
                 'Q',[]);           %always put an empty vector

% generate the Q matrix for the material laminate
carbon = calculate_Q(carbon); 

% core material (vynil foam 5 lb. density)
% create astructure of properties

% high-quality PVC foam
% resin = struct('E',95,...     %[MPa] tensile modulus
%                 'v',0.7,...
%                 'G',27,...   %Gr= Er/(2*(1+vr)) shear modulus
%                 'sigma',2.5,... %[MPa] tensile strength
%                 'tau',1.15,... %[MPa] shear strength
%                 'dens',0.08);  %[kg/dm^3]

% % poliurethane foam
resin = struct('E',25.7,...     %[MPa] tensile modulus
                'v',0.3,...
                'G',3.93,...   %Gr= Er/(2*(1+vr)) shear modulus
                'sigma',0.8,... %[MPa] tensile strength
                'tau',0.441,... %[MPa] shear strength
                'dens',0.096);  %[kg/dm^3]

                
            
%% Laminate description
layer1 = struct('material',carbon,... [choose a material properties' struct]
                'angle',0*2*pi/360,... [angle in rads]
                'thick',0.1,...       %[mm] 
                'Q_mod',[],...        % Q*R*T*inv(R); always put an empty vector
                'Q_hat',[]);          %always put an empty vector
layer1 = calculate_hatQ(layer1);

layer2 = struct('material',carbon,... [choose a material properties' struct]
                'angle',45*2*pi/360,... [angle in rads]
                'thick',1,...       %[mm] 
                'Q_mod',[],...        % Q*R*T*inv(R); always put an empty vector
                'Q_hat',[]);          %always put an empty vector
layer2 = calculate_hatQ(layer2);

layer3 = struct('material',carbon,... [choose a material properties' struct]
                'angle',-45*2*pi/360,... [angle in rads]
                'thick',1,...       %[mm] 
                'Q_mod',[],...        % Q*R*T*inv(R); always put an empty vector
                'Q_hat',[]);          %always put an empty vector
layer3 = calculate_hatQ(layer3);


layer4 = struct('material',carbon,... [choose a material properties' struct]
                'angle',0*2*pi/360,... [angle in rads]
                'thick',0.1,...       %[mm] 
                'Q_mod',[],...        % Q*R*T*inv(R); always put an empty vector
                'Q_hat',[]);          %always put an empty vector
layer4 = calculate_hatQ(layer4);

layup=[layer1;layer2;layer3;layer4];
laminate =  struct ('layup',layup,...
                    'A',[],...  always put an empty vector: it will be calculated after
                    'B',[],...  always put an empty vector: it will be calculated after
                    'D',[],...  always put an empty vector: it will be calculated after
                    'thick',0,...always put 0: it will be calculated after
                    'dens',0,... always put 0: it will be calculated after
                    'Exx',0,...  always put 0: it will be calculated after
                    'Eyy',0,...  always put 0: it will be calculated after
                    'Gxy',0,...  always put 0: it will be calculated after
                    'vxy',0,...  always put 0: it will be calculated after
                    'vyx',0);  % always put 0: it will be calculated after
laminate= calculate_laminate_ABD(laminate);

%% PARAMETRIC LAMINATE
% differently from the previous, the parmetric laminate assumes a fixed
% thickness for every layer, in this way the apparent properties can be
% calculated by imposing the thickness of the whole laminate and also it
% can be derived the thickness of a single layer

param_laminate = struct ('layup',[struct('material',carbon,'angle',0*pi/180) ;...     layer1
                                  struct('material',carbon,'angle',90*pi/180) ;...    layer2
                                  struct('material',carbon,'angle',-90*pi/180) ;...   layer3
                                  struct('material',carbon,'angle',0*pi/180)] ,...    layer4
                           'laminate',laminate,... %it is needed a laminate struct
                           'num_layers',0,...  always put 0: it will be calculated after
                           'layer_thick',0,... %  always put 0: it will be calculated after
                           'sigma_f',0);%  [MPa] first-ply-failure always put 0: it will be calculated after

% -------- DEBUG --------
param_laminate=calculate_param_laminate(param_laminate,1);
%param_laminate=calculate_first_ply_failure(param_laminate);
% -------- MOVE AWAY FROM HERE --------


%% Calculate the minimum stiffness of the sandwich
b = 0.2;                    %[m] skate_width: this is the standard width of a skateboard
P = 4600;                   %[N] max_load: the maximum load calculated, due to a jump
eta_cent = 0.015;           %[m] max_def_at_center: maximum deformation admitted in the centre of the skateboard
L = 0.65;                   %[m] skateboard_length

%since the deformation in the center is = P*L^3/(48*EJ), it turns that:
EJ_min=P*(L^3)/(48*eta_cent); %[N*m^3/m]=[N*m^2]

% The EJ_min found guarantees the minimum stiffness of the panel


%% Calculation of failures and panel density

% The stress applied is that of a link over 2 vinculums.
% constants htat depends on this type of load
Bm = 1/4;
Bt = 1/2;
L = L*1000; % express the panel length in [mm]

% Face and core thickness intervals:
face_steps=350;
core_steps=50;
tf = linspace(0.1,3,face_steps)';  %[mm] face thickness: an interval between 0 and 3mm
tc = linspace(7,70,core_steps)'; %[mm] core thickness

% initialize the thicknesses
min_stiff_face_thick= ones(core_steps,1)*max(tf);
min_Pff_face_thick  = ones(core_steps,1)*max(tf);
min_Pcs_face_thick  = ones(core_steps,1)*max(tf);
min_Pwr_face_thick  = ones(core_steps,1)*max(tf);

prev_tf_step = 0; 
tic

for j=1:size(tf)
    %first of all, clculate laminate properties, for every face thickness
    curr_laminate(j) = calculate_param_laminate(param_laminate,tf(j));
end

for k=1:size(tc)
    
    for j=1:size(tf)
        % stiffness of the sandwich, as a function of core and face
        % thickness.
        % Thanks to semplification, the stiffness of the sandwich is equivalent to
        % that of the face, with respect to the sandwich's center:
        d=tf(j)+tc(k);
        D0(j,k)=((curr_laminate(j).laminate.Exx*b*1000*tf(j)*(d)^2)/2)*10^(-6); %([MPa*mm^4])e-6=([N*mm^2])e-6=[N*m^2]
        % specific weigth of the sandwich, as a function of core and face thickness
        sandw_dens(j,k)=((2*tf(j)*curr_laminate(j).laminate.dens) + resin.dens*tc(k))/(2*tf(j)+tc(k));%[kg/dm^3]

        % maximum load for FACE-FAILURE
        Pff(j,k)=1000*curr_laminate(j).sigma_f*b*tf(j)*d/(Bm*L);
        % maximum load for CORE-SHEAR
        Pcs(j,k)=1000*resin.tau * d / Bt;
        % maximum load for FACE-WRINKLING
        Pwr(j,k)=1000*tf(j)*d/(2*Bm*L)*(curr_laminate(j).laminate.Exx * resin.E * resin.G)^(1/3);

        % now define the functions that accomplishes the requests
        if((D0(j,k)>=EJ_min) && (tf(j)<min_stiff_face_thick(k)))
           min_stiff_face_thick(k)=tf(j);
        end
        if((Pff(j,k)>=P) && (tf(j)<min_Pff_face_thick(k)))
           min_Pff_face_thick(k)=tf(j);
        end
        if((Pcs(j,k)>=P) && (tf(j)<min_Pcs_face_thick(k)))
           min_Pcs_face_thick(k)=tf(j);
        end
        if((Pwr(j,k)>=P) && (tf(j)<min_Pwr_face_thick(k)))
           min_Pwr_face_thick(k)=tf(j);
        end
    end
    
    % now that the minimum face thickness is found, keep the worst
    % condition and plot the corresponding density of the sandwich
    min_face_thick(k) = max([min_stiff_face_thick(k),min_Pff_face_thick(k),min_Pcs_face_thick(k),min_Pwr_face_thick(k)]);
    
    % having this data, derive the corresponding density:
    panel_weigth(k)=(resin.dens * tc(k) + 2* min_face_thick(k) * curr_laminate(1).laminate.dens)*b*L/1000;
end
toc

% Now find the minimum weigth and its index
[min_w index_min_w]=min(panel_weigth);

fig_title=sprintf('Minimum weigth comes when core is: %f mm, and face is: %f mm thick \n Total panel weigth is: %3f kg',tc(index_min_w),min_face_thick(index_min_w),panel_weigth(index_min_w));
%htext = uicontrol('Style','text','String','Face material','Position', [h_pos 20 h_pos+50 v_pos+20]);
%plot(tc,[min_stiff_face_thick,min_Pff_face_thick,min_Pcs_face_thick,min_Pwr_face_thick,density']);
plot(tc,[min_stiff_face_thick,min_Pff_face_thick,min_Pcs_face_thick,min_Pwr_face_thick,panel_weigth']);
title(fig_title);
legend('stiffness constraint','face failure','core shear','face wrinkling','panel weigth');
xlabel('core thickness [mm]');
ylabel('face thickness [mm]');

%surf(tc,tf,D0);
%surf(tc,tf,sandw_dens);

            