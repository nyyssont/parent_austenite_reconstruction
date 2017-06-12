% Tuomo Nyyssönen 2017
% Parent austenite grain (PAG) reconstruction 4/3
% Script to export reconstruction data to a text file readable by Channel 5
% software.

clear all
close all

% Plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'Iron bcc (old)', 'color', 'light blue')};
  
% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

%% Load preprocessed data:

load('pag3_data.mat','ebsd_aus_orig_res')

%% Sort the exported data to columns:

phase_column = ebsd_aus_orig_res.phaseId;
phase_column(phase_column == 1) = 0;
phase_column(phase_column == 3) = 1;

x_column = ebsd_aus_orig_res.x;
y_column = ebsd_aus_orig_res.y;

% A bit tricky to sort out the Euler angles. There is no orientation info
% on unindexed pixels.

% Generate a pre-list consisting of only zeroes:
phi1_column = zeros(length(ebsd_aus_orig_res),1);
Phi_column = phi1_column;
phi2_column = phi1_column;

% This is just a long list;
[phi1,Phi,phi2] = Euler(ebsd_aus_orig_res('Iron fcc').orientations);

phi1_column(phase_column == 1) = phi1;
Phi_column(phase_column == 1) = Phi;
phi2_column(phase_column == 1) = phi2;

bc_column = ebsd_aus_orig_res.bc;

%% Export as a Channel 5 readable csv file:

all_columns = [phase_column,x_column,y_column,phi1_column,Phi_column,phi2_column,bc_column];

dlmwrite('EC003_50-EC003_50_05_4_austenite.txt',all_columns,';')