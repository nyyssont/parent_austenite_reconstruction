% Tuomo Nyyssönen 2017
% Parent austenite grain (PAG) reconstruction 1/3
% Script to import EBSD data for OR detection & parent austenite reconstruction. To be used with MTEX.

%% Import Script for EBSD Data

clear all
close all

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'Iron bcc (old)', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [3.6599 3.6599 3.6599], 'mineral', 'Iron fcc', 'color', 'light green')};

% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outofPlane');

%% Specify File Names

% pname = ['C:\Local\nyyssont\EBSD\Raex500_PO_ebsd\'];

% which files to be imported
% fname = [pname 'Rx1po-03um step 424x291um 9to10 bands #1.cpr'];
fname = ['EC003_50-EC003_50_05_4.cpr'];

%% Import the Data

% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,'interface','crc',...
  'convertEuler2SpatialReferenceFrame');

% Save the original ebsd data to another object:

%% Grain and boundary identification

% Assign the bcc and fcc crystal symmetries to an object:
cs_fcc = crystalSymmetry('m-3m', [3.6599 3.6599 3.6599], 'mineral', 'Iron fcc', 'color', 'light green');
cs_bcc = ebsd('Iron bcc').CS;

% For improved grain reconstruction, set nonindexed values to empty:
% ebsd('notIndexed')=[];
% ebsd('Iron fcc')=[];

%Identify grains with 3 degree boundary
[grains,ebsd('Iron bcc').grainId] = calcGrains(ebsd('Iron bcc'),'angle',3*degree);

%Remove small grains
ind = grains.grainSize < 4;
ebsd(grains(ind)).phase = 0;
ebsd(grains(ind)).grainId = 0;

%Reidentify grains with small grains removed:
[grains,ebsd('Iron bcc').grainId] = calcGrains(ebsd('Iron bcc'),'angle',3*degree);

%% Export data

save 'pag1_data.mat'
