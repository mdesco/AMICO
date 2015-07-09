%
% Initialization for AMICO
%
% NB: DO NOT MODIFY THIS FILE!
%     Make a copy of it, adapt to your paths and rename it to "AMICO_Setup.m"
%

global AMICO_code_path AMICO_data_path CAMINO_path CONFIG
global niiSIGNAL niiMASK
global KERNELS bMATRIX

% Path definition: adapt these to your needs
% ==========================================
AMICO_code_path = '/Users/desm2239/Research/Source/AMICO';
AMICO_data_path = '/Users/desm2239/Research/Data';
CAMINO_path     = '/Users/desm2239/Research/Source/camino/bin';
NODDI_path      = '/Users/desm2239/Research/Source/NODDI_toolbox_v0.9/';
SPAMS_path      = '/Users/desm2239/Research/Source/spams-matlab/';

addpath( genpath(NODDI_path) )
addpath( fullfile(SPAMS_path,'build') )
addpath( fullfile(AMICO_code_path,'kernels') )
addpath( fullfile(AMICO_code_path,'models') )
addpath( fullfile(AMICO_code_path,'optimization') )
addpath( fullfile(AMICO_code_path,'other') )
addpath( fullfile(AMICO_code_path,'vendor','NIFTI') )
