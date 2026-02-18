%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, reduce to target grid size
n_units = 10;                    % number of rock units contained in image
units   = imread('units.tiff');  % read in cross section with rock units
reduce  = 2;                     % reduce model resolution by factor
units   = imresize(units,1/reduce,'method','nearest');  % resize units map to model size

% set model dimensions
W       = 16e3;        % domain width (corresponds to width of geological cross section) [m]
[Nz,Nx] = size(units); % grid size according to units map
D       = W*Nz/Nx;     % domain depth (according to grid aspect ratio)
h       = W/Nx;        % grid spacing

% material properties for each rock unit (update based on your calibration)

matprop = [
% unit  kT    rho    cP     Qr    KD
   1	 1   1000   2000   0e-6   1e-9   % air/water
   2     1	 2500	1000   1e-6   1e-8   % Sa
   3     1	 2500	1000   1e-6   1e-8   % Si
   4     1	 2500	1000   1e-6   1e-8   % Gr
   5	 1	 2500	1000   1e-6   1e-8   % He1
   6	 1	 2500	1000   1e-6   1e-8   % Bg
   7	 1	 2500	1000   1e-6   1e-8   % He2
   8     1	 2500	1000   1e-6   1e-8   % Fz
   9     1	 2500	1000   1e-6   1e-8   % Ms
  10     1   2500	1000   1e-6   1e-8]; % Cm
air = units==1;

% get coefficient fields based on spatial distribution of rock units from image
kT0  = reshape(matprop(units,2),Nz,Nx);  % thermal conductivity [W/m/K]
rho0 = reshape(matprop(units,3),Nz,Nx);  % density [kg/m3]
cP0  = reshape(matprop(units,4),Nz,Nx);  % heat capacity [J/kg/K]
Qr0  = reshape(matprop(units,5),Nz,Nx);  % heat productivity [W/m3]
KD0  = reshape(matprop(units,6),Nz,Nx);  % segregation mobility [m2/Pas]
KD0  = imgaussfilt(KD0,1);               % apply some smoothing for numerical stability

% continue setting remaining model parameters, then call model routine
% ...

%*****  RUN MODEL
run('../src/main.m');