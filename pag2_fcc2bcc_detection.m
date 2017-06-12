% Tuomo Nyyssönen 2017
% Parent austenite grain (PAG) reconstruction 2/3
% Script to determine the orientation relationship between martensite and
% austenite from martensitic EBSD data. To be used with MTEX.

clear all
close all

% Plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'Iron bcc (old)', 'color', 'light blue')};
  
% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

%% Input script parameters

cutoff = 5; % maximum allowed angular deviation from the intermartensitic misorientation list

%% Load preprocessed data:

load 'pag1_data.mat'

%% Determine the poles to be used in fcc PF plots:

h_bcc = [Miller(1,0,0,ebsd('Iron bcc').CS),Miller(1,1,0,ebsd('Iron bcc').CS),Miller(1,1,1,ebsd('Iron bcc').CS)];
h_fcc = [Miller(1,0,0,cs_fcc),Miller(1,1,0,cs_fcc),Miller(1,1,1,cs_fcc)];

%% Parent austenite boundary detection

% Determine the initial orientation relationship as Kurdyumov-Sachs:
fcc2bcc_KS = (orientation('map',Miller(1,1,1,cs_fcc),Miller(0,1,1,cs_bcc),Miller(-1,0,1,cs_fcc),Miller(-1,-1,1,cs_bcc)));

% Symmetry order correction according to Morito:
cs_symmetries = cs_fcc.properGroup;

symrots_true(1,1) = cs_symmetries(1);
symrots_true(2,1) = cs_symmetries(21);
symrots_true(3,1) = cs_symmetries(11);
symrots_true(4,1) = cs_symmetries(16);
symrots_true(5,1) = cs_symmetries(24);
symrots_true(6,1) = cs_symmetries(8);
symrots_true(7,1) = cs_symmetries(18);
symrots_true(8,1) = cs_symmetries(10);
symrots_true(9,1) = cs_symmetries(20);
symrots_true(10,1) = cs_symmetries(3);
symrots_true(11,1) = cs_symmetries(7);
symrots_true(12,1) = cs_symmetries(23);
symrots_true(13,1) = cs_symmetries(19);
symrots_true(14,1) = cs_symmetries(2);
symrots_true(15,1) = cs_symmetries(9);
symrots_true(16,1) = cs_symmetries(22);
symrots_true(17,1) = cs_symmetries(17);
symrots_true(18,1) = cs_symmetries(12);
symrots_true(19,1) = cs_symmetries(5);
symrots_true(20,1) = cs_symmetries(15);
symrots_true(21,1) = cs_symmetries(4);
symrots_true(22,1) = cs_symmetries(14);
symrots_true(23,1) = cs_symmetries(6);
symrots_true(24,1) = cs_symmetries(13);

% Determine the full orientation relationship according to cubic symmetry
% as by Morito:
fcc2bcc_true = orientation(symrots_true*rotation(fcc2bcc_KS),cs_bcc,1);

for l = 1:20
    
    % Determine intermartensitic misorientations:
    bcc_trans =  inv(orientation(rotation(fcc2bcc_true),cs_bcc,1)) * (orientation(rotation(fcc2bcc_true(1)),cs_bcc,1));
    
    % Determine the interferritic grain boundaries from ebsd data:
    bound_bcc = grains.boundary('Iron bcc','Iron bcc');
    
    % Determine all neighboring grain pairs:
    bound_bcc_ids_sorted = sort(bound_bcc.grainId,2);
    bound_bcc_ids = unique(bound_bcc_ids_sorted,'rows');
        
    % Determine the misorientations between all grain pairs using mean orientations:
    misos_bcc_grains = inv(grains(bound_bcc_ids(:,1)).meanOrientation).*grains(bound_bcc_ids(:,2)).meanOrientation;
    
    % Determine the angles between the possible misorientation set and the
    % grain pair misorientations:
    angles = angle_outer(misos_bcc_grains,bcc_trans) / degree;
    
    % Determine the minimum angle and corresponding index between the
    % misorientation set and each grain pair misorientation:
    [omega, omega_I] = min(angles,[],2);
    
    % Limit the list of indexes of the misorientations to those whose minimum angle between the
    % misorientation set is below the cutoff value:
    omega_I_cutoff = omega_I(omega < cutoff);
    
    % Get the misorientations between intermartensitic grain pairs.
    misos_bcc_grains_omega_cutoff = misos_bcc_grains(omega < cutoff);
    
    % Get the grain pairs whose minimum angle between the misorientation set is
    % below the cutoff value (intramartensitic grain pairs).
    OR_rows = bound_bcc_ids(omega < cutoff,:);

    % Get the intramartensitic grain pairs whose index is below 7 (intrapacket
    % grain pairs):
    OR_rows_omega_I_cutoff = OR_rows(omega_I_cutoff < 7,:);
            
    %% Return all misorientations to "primary" misorientation:
    
    % Allocate object for misorientations;
    primary_misos = orientation(rotation(fcc2bcc_true(1)),cs_bcc,1);
    
    % For loop for correcting additional symmetry problem in
    % misorientations (each misorientation has 24 different possible
    % crystallographic orientations according to cubic symmetry):
    
    for k = 1:24
        
        prime_angle_miso = angle_outer(orientation(rotation(misos_bcc_grains_omega_cutoff(omega_I_cutoff == k)),cs_bcc,1),orientation(symrots_true*rotation(bcc_trans(k)),cs_bcc,1));
        
        if isempty(prime_angle_miso) == 0;
            [omega_prime,omega_I_prime] = min(prime_angle_miso,[],2);
            
            z = inv(symrots_true(omega_I_prime)).*orientation(rotation(misos_bcc_grains_omega_cutoff(omega_I_cutoff == k)),cs_bcc,1);
            
            primary_misos = [primary_misos;z];
        end
        
    end
    
    %% Handle returning the misorientations to primary orientation relationship:
    
    angle_getting_there = angle_outer(primary_misos,orientation(rotation(bcc_trans),cs_bcc,1));
    [omega_getting_there,omega_I_getting_there] = min(angle_getting_there,[],2);
    
    fcc2bcc_candidate = rotation(fcc2bcc_true(omega_I_getting_there)).*primary_misos;
    
    %Return the candidate rotation closest to the "original":
    
    candidate = mean(fcc2bcc_candidate);
    
    [M2, I2] = min(angle_outer(candidate*symrots_true,rotation(fcc2bcc_KS)));
    candidate = candidate*symrots_true(I2);    
    
    fcc2bcc_true = orientation(rotation(symrots_true*candidate),cs_bcc,1);
    
    fcc2bcc_basis_miso = orientation(rotation(fcc2bcc_true(1)),cs_bcc,cs_fcc)
    angle_planes(l) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,1,-1,cs_bcc)),Miller(0,1,1,cs_fcc))) / degree;
    angle_directions(l) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,0,1,cs_bcc,'uvw')),Miller(1,1,-1,cs_bcc,'uvw'))) / degree;

end

angle_planes(21) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,1,-1,cs_bcc)),Miller(0,1,1,cs_fcc))) / degree;
angle_directions(21) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,0,1,cs_bcc,'uvw')),Miller(1,1,-1,cs_bcc,'uvw'))) / degree;

%% Plot figures:

figure
plotPDF(orientation(rotation(fcc2bcc_true),cs_bcc,1),h_bcc,'marker','o','MarkerSize',4,'MarkerFaceColor','white','MarkerEdgeColor',[0 0 0],'LineWidth',2,'MarkerColor','white')

figure
plot(ebsd,ebsd.bc)
colormap gray
hold on;
plot(bound_bcc,'linewidth',2)
plot(bound_bcc(ismember(bound_bcc_ids_sorted,OR_rows,'rows')),'linecolor','r','linewidth',1,'faceAlpha',0.5)
plot(bound_bcc(ismember(bound_bcc_ids_sorted,OR_rows_omega_I_cutoff,'rows')),'linecolor','g','linewidth',1,'faceAlpha',0.5)
plot(bound_bcc(not(ismember(bound_bcc_ids_sorted,OR_rows,'rows'))),'linecolor','k','linewidth',2)

% %% Exclude data from MCL analysis:
% 
% [a,b] = ismember(ebsd.grainId,1:length(grains));
% idx = 1:numel(ebsd.grainId);
% ib = accumarray(nonzeros(b), idx(a), [], @(x){x});
% ebsd_bs = ebsd.bs;
% 
% mean_bs = cellfun(@(x)mean(ebsd_bs(x)), ib, 'UniformOutput', true);
% b_limit = find(mean_bs>130);
% 
% figure
% plot(ebsd_orig,ebsd_orig.bc)
% colormap gray
% hold on;
% plot(grains(b_limit),'faceAlpha',0.2)

%% Create data for MCL analysis:

surv_val = 1-cdf('Burr',omega,2,5,1);

surv_val(surv_val<0.001) = 0;

mcl = [bound_bcc_ids, surv_val];

% log1 = ismember(mcl,b_limit);
% 
% log2 = not(log1(:,1) | log1(:,2) | log1(:,3));
% 
% mcl = mcl(log2,:);
% 
% dlmwrite('C:\Local\nyyssont\Programs\cygwin64\home\nyyssont\mcl.txt',mcl,'delimiter',' ','newline','pc');
% 
% %% Perform mcl analysis:
% 
% system('C:\Local\nyyssont\Programs\cygwin64\bin\bash --login -c "C:/Local/nyyssont/Programs/cygwin64/home/nyyssont/pag.sh"');
% 
% data = dlmread('C:\Local\nyyssont\Programs\cygwin64\home\nyyssont\dump.data.mci.I16');

%% Save data:

save('pag2_data.mat')
