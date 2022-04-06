%% Phase profile calculation for rectangular reflectarray with feed positioned at (x,y,z) = (x_offset,y_offset,R0)
close all
clc
clear

%Load the reflection phase profile data of the unitcell from HFSS
load ReflectionPhaseTE00.mat

%Set operating frequency as f
f = 12*10^9;
c = 3*10^8;
lamda =1000*(c/f);
k = 2*pi/lamda;

%Set phase center(feed distance)in z direction.
R0 = 198;

%Set x_offset
x_offset=-81;

%Set x_offset
y_offset=-81;

%Set aperture size consider F/D
D = R0/0.84;

%Set unit cell size(grid) as lamda/2
unit_cell = lamda/3;
resolution=20;
%resolution = D/ceil(unit_cell);
num_elements = resolution^2;



%Set a phase constant as a relative phase shift between the elements
%phi_con =0;


%Set desired beam direction as (theta_0,phi_0)
phi_0 = 0;
theta_0 = 0;
for i=-resolution/2:resolution/2-1
    for j=-resolution/2:resolution/2-1
        x(i+resolution/2+1) = (i)*unit_cell;
        y(j+resolution/2+1) = (j)*unit_cell;
        R(i+resolution/2+1,j+resolution/2+1) = sqrt((x(i+resolution/2+1)+x_offset)^2 + (y(j+resolution/2+1)+y_offset)^2+ R0^2);
        phase_profile(i+resolution/2+1,j+resolution/2+1) =rad2deg(k*R(i+resolution/2+1,j+resolution/2+1)) -rad2deg(k*(sin(theta_0)*(x(i+resolution/2+1)*cos(phi_0)+y(j+resolution/2+1)*sin(phi_0))));
        %phase_profile(i+resolution/2+1,j+resolution/2+1) = -rad2deg(k*R(i+resolution/2+1,j+resolution/2+1));
       % phase_profile(i+resolution/2+1,j+resolution/2+1) = -rad2deg(k*(sin(theta_0)*(x(i+resolution/2+1)*cos(phi_0)+y(j+resolution/2+1)*sin(phi_0))));
    end
end





%Divide the reflection phase from HFSS by 2 due to the way the parametric
%analysis is setup
ReflectionPhaseTE00(:,1) = ReflectionPhaseTE00(:,1)./2;


%Calculate the quantazation error(phase difference)between the required phase distribution and
%the one that can be closer achieved by the given unit cell. Using the
%reflection phase response response of the unicell vs size. Save the
%position of the size in idx that gives the min difference between the phase
%profiles
l=1;
for i=-resolution/2:resolution/2-1
    for j=-resolution/2:resolution/2-1
        [value,idx(l)] = min(abs(rad2deg(wrapTo2Pi(deg2rad(phase_profile(i+resolution/2+1,j+resolution/2+1)))) - rad2deg(wrapTo2Pi(deg2rad(ReflectionPhaseTE00(:,2))))));
        diff(i+resolution/2+1,j+resolution/2+1) =  min(abs(rad2deg(wrapTo2Pi(deg2rad(phase_profile(i+resolution/2+1,j+resolution/2+1)))) - rad2deg(wrapTo2Pi(deg2rad(ReflectionPhaseTE00(:,2))))));
        l=l+1;
         
    end
    
end

% Create the size database for each element
for i=1:length(idx)
    unitcell_database(i) = ReflectionPhaseTE00(idx(i),1);
end

%reshape the size to plot it in a hetamap distribution
sizes = reshape(unitcell_database, length(phase_profile),length(phase_profile));

% Calculate the max and min quantization error
Phase_Error_max = max(reshape(diff,1,length(diff)^2));
Phase_Error_min = min(reshape(diff,1,length(diff)^2));

disp(['The maximum phase error is: ', num2str(Phase_Error_max), ' degrees'])
disp(['The minimum phase error is: ', num2str(Phase_Error_min), ' degrees'])
disp(['The number of elements is :' num2str(num_elements)])

%Plot the output
figure
h1 = heatmap(rad2deg(wrapTo2Pi(deg2rad(phase_profile))));
colormap hot
h1.Title = 'Phase profile distribution';
caxis([0,360])
figure
h2 = heatmap(diff);
colormap hot
h2.Title = 'Phase difference profile distribution';
figure
h2 = heatmap(sizes);
colormap hot
h2.Title = 'Patch Width per element@ resonant patch length';

% Output the phase profile in a excel file
output = reshape(phase_profile',[length(phase_profile)^2,1]);
csvwrite('Reflectarray_Phase_Excitation.csv',output)