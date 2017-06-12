% Tuomo Nyyssönen 2017
% Parent austenite grain (PAG) reconstruction 3/3
% Script to merge grains clustered with MCL software. To be used with MTEX.

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

load 'pag2_data.mat'

%% Determine script parameters

p = 1.6; %The exponent to be used in mcl clustering. Increasing the value will produce smaller clusters.

%% Determine the rotation to be used for identification

%Set the misorientation corresponding to "bcc to fcc" transformation
bcc2fcc = symrots_true * inv(fcc2bcc_true(1));
bcc_trans = orientation(bcc2fcc,cs_fcc,cs_bcc);

%% Dirty mcl:

z = sparse(mcl(:,1),mcl(:,2),mcl(:,3));
zz = z;
zzz = [zz;zeros((length(zz) - min(size(zz))),length(zz))];
zzz = zzz' + zzz;
zzzz = zzz + eye(length(zzz));
zzzzz = normc(zzzz);
zzzzzz = sparse(zzzz);

mcl_matrix = full(mcl_func(p,zzzzzz));

%%
[row_mcl,col_mcl] = find(mcl_matrix);

[a,b] = ismember(row_mcl,1:length(grains));
idx = col_mcl;
ib = accumarray(nonzeros(b), idx(a), [], @(x){x});

%# find empty cells
emptyCells = cellfun(@isempty,ib);
%# remove empty cells
ib(emptyCells) = [];

%%

ebsd_aus = ebsd;

for l = 1:length(ib);
    
    l
    
    line = ib{l};
    line = line(line>0);
    
    %Calculate potential parent austenite orientations for grains in the set:
    fcc_parent_pot = grains(line).meanOrientation*bcc_trans;
    
    %Create the orientation density function (ODF) for all potential parent
    %austenite orientations:
    odf = calcODF(fcc_parent_pot,'halfwidth',4*degree);
    
    %Determine the most intense orientation in the ODF (average parent austenite orientation):
    [center,value] = calcModes(odf,'resolution',5*degree);
    
    %Determine the amount of data points:
    n=length(fcc_parent_pot(:,1));
    
    %Determine the amount of orientation alternatives per data point:
    i=length(fcc_parent_pot(1,:));
    
    %Calculate angles of all potential orientations and the average parent 
    %austenite orientation with a single command, creates a column vector:
    misos_list=angle_outer(fcc_parent_pot,center) / degree;
    
    %Convert the list into a matrix with the previous size:
    misos=reshape(misos_list,n,i);
        
    %Find the minimum angular deviation and variant index of same on each
    %row:
    [M, I] = min(misos,[],2);

    %Turn I into proper linear indexing (to speed up script):
    I_ind = sub2ind(size(fcc_parent_pot),(1:length(I))',I);
    
    %Get the fcc parent orientations of the lowest angular deviation:
    fcc_parents = fcc_parent_pot(I_ind);
                
    %Correct the phase information of grains eligible for transformation:
    ebsd_aus(grains(line)).phase = 2;
    ebsd_aus(grains(line)).orientations = center;
    
    %Update devis:
    devis(line) = M;
end

%% Variant analysis

%Get all the misorientations from the ferritic phase to the austenitic
%phase:
fccs2bccs = inv(ebsd_aus(grains(devis>0)).orientations).*ebsd(grains(devis>0)).orientations;

%Define the fcc to bcc orientation relationship as by Morito et al.:
fcc2bcc_true = orientation(rotation(symrots_true*fcc2bcc_true(1)),cs_bcc,cs_fcc);

%Calculate the angles between the misorientations and the components of the
%full orientation relationship:
variant_angles = angle_outer(orientation(rotation(fccs2bccs),cs_fcc,1),orientation(rotation(fcc2bcc_true),cs_fcc,1)) / degree;

[M_variants, I_variants] = min(variant_angles,[],2);

%% Resolution restoration

%Get all the misorientations from the ferritic phase to the austenitic
%phase:
bccs2fccs = inv(ebsd(grains(devis>0)).orientations).*ebsd_aus(grains(devis>0)).orientations;

%Calculate the angles between the misorientations and the components of the
%full orientation relationship:
variant_angles = angle_outer(orientation(rotation(bccs2fccs),cs_fcc,1),orientation(rotation(bcc_trans),cs_fcc,1)) / degree;

[M_res, I_res] = min(variant_angles,[],2);

%Restore resolution for individual pixels in the map:
ebsd_aus_orig_res = ebsd_aus;

ebsd_aus_orig_res(grains(devis>0)).orientations = ebsd(grains(devis>0)).orientations.*bcc_trans(I_res);

%% Grain detection for reaustenitized orientation map:

%Identify grains with 5 degree boundary
[grains_aus,ebsd_aus.grainId,ebsd_aus.mis2mean] = calcGrains(ebsd_aus('Iron fcc'),'angle',10*degree);

%% Plot results

% Get interaustenitic boundaries:
gB_aus = grains_aus.boundary('Iron fcc','Iron fcc');

% Select CSL(3) grain boundaries
gB_aus_3 = gB_aus(angle_outer(gB_aus.misorientation,orientation(CSL(3,crystalSymmetry('m-3m')),cs_fcc,cs_fcc)) < 3*degree);

figure
plot(ebsd,ebsd.bc)
colormap gray
hold on;
plot(grains_aus('Iron fcc'),grains_aus('Iron fcc').meanOrientation,'faceAlpha',0.5)
hold on
plot(gB_aus_3,'lineColor','g','linewidth',2)

figure
plot(ebsd,ebsd.bc)
colormap gray
hold on;
plot(ebsd_aus_orig_res('Iron fcc'),ebsd_aus_orig_res('Iron fcc').orientations,'faceAlpha',0.5)
hold on
plot(gB_aus_3,'lineColor','g','linewidth',2)

%Comparison of grain boundary vs. reconstruction analysis
figure
plot(grains_aus('Iron fcc'),grains_aus('Iron fcc').meanOrientation,'faceAlpha',0.5)
hold on
plot(bound_bcc(not(ismember(bound_bcc_ids_sorted,OR_rows,'rows'))),'linecolor','k','linewidth',2)

%% Save data

save('pag3_data.mat')
