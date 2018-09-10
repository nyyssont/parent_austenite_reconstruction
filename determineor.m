function [fcc2bcc,omega,bound_bcc_ids,bound_bcc,OR_rows,OR_rows_omega_I_cutoff,bound_bcc_ids_sorted,hist_all_omega_I,symrots_true,bcc2fcc_odf] = determineor(symmetries,grains,cond)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

w = warning ('off','all');

symmetry1_string = symmetries{1}.mineral;

% Determine the initial orientation relationship as Kurdjumov-Sachs:
fcc2bcc_KS = (orientation('map',Miller(1,1,1,symmetries{2}),Miller(0,1,1,symmetries{1}),Miller(-1,0,1,symmetries{2}),Miller(-1,-1,1,symmetries{1})));

% Symmetry order correction according to Morito et al.:
cs_symmetries = symmetries{1}.properGroup;

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

symrots_true = rotation(symrots_true);

% Determine the full orientation relationship according to cubic symmetry
% as by Morito et al.:
fcc2bcc_true = orientation(symrots_true*rotation(fcc2bcc_KS),symmetries{1},1);

if cond
    figure
    subplot(2,1,1)
    h_1(1) = animatedline('Color','r');
    hold on
    h_1(2) = animatedline('Color','g');
    axis([0 20 0 5])
    z1 = legend(h_1,'Angle between {111}a and {011}m','Angle between <110>a and <111>m','Location','Northwest');
    set(z1,'interpreter','none');
    
    subplot(2,1,2)
    h_2 = bar(1:24,zeros(24,1));
    ylim([0 0.3]);
end
l = 1;

while  l<3 || (angle_planes(end)-angle_planes(end-1))>0.01 % Determine end condition as reasonably diminished angular deviation between close-packed planes
    
    % Determine intermartensitic misorientations:
    bcc_trans =  inv(orientation(rotation(fcc2bcc_true),symmetries{1},1)) * (orientation(rotation(fcc2bcc_true(1)),symmetries{1},1));
    
    % Determine the interferritic grain boundaries from ebsd data:
    bound_bcc = grains.boundary(symmetry1_string(1:end-1),symmetry1_string(1:end-1));
    
    % Determine all neighboring grain pairs:
    bound_bcc_ids_sorted = sort(bound_bcc.grainId,2);
    [bound_bcc_ids,~,bbcc_ic] = unique(bound_bcc_ids_sorted,'rows');
        
    % Determine the misorientations between all grain pairs using mean orientations:
    misos_bcc_grains = inv(grains(bound_bcc_ids(:,1)).meanOrientation).*grains(bound_bcc_ids(:,2)).meanOrientation;
    
    % Determine the angles between the possible misorientation set and the
    % grain pair misorientations:
    angles = angle_outer(misos_bcc_grains,bcc_trans) / degree;
    
    % Determine the minimum angle and corresponding index between the
    % misorientation set and each grain pair misorientation:
    [omega, omega_I] = min(angles,[],2);

    [hist_omega_I, bins_omega_I] = hist(omega_I,1:24);
        
    % Determine cutoff as value encompassing 90 % of misorientations
    [hist_omega, bins_omega] = hist(omega,100);
    hist_omega = cumsum(hist_omega)/sum(hist_omega);
    I_hist_omega = find(hist_omega>0.9,1,'first');
    cutoff = bins_omega(I_hist_omega);
    cutoff(cutoff>5) = 5;
    all_cutoffs(l) = cutoff;

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
    
    all_omega_I = omega_I(bbcc_ic);
    all_omega_I = all_omega_I(ismember(bound_bcc_ids_sorted,OR_rows,'rows'));
    [hist_all_omega_I, bins_all_omega_I] = hist(all_omega_I,1:24);
    hist_all_omega_I = hist_all_omega_I/sum(hist_all_omega_I);
if cond
    set(h_2, 'YData', hist_all_omega_I);
end
    %% Return all misorientations to "primary" misorientation:
    
    % Allocate object for misorientations;
    primary_misos = orientation(rotation(fcc2bcc_true(1)),symmetries{1},1);
    
    % For loop for correcting additional symmetry problem in
    % misorientations (each misorientation has 24 different possible
    % crystallographic orientations according to cubic symmetry):
    
    for k = 1:24
        
        prime_angle_miso = angle_outer(orientation(rotation(misos_bcc_grains_omega_cutoff(omega_I_cutoff == k)),symmetries{1},1),orientation(symrots_true*rotation(bcc_trans(k)),symmetries{1},1));
        
        if isempty(prime_angle_miso) == 0;
            [omega_prime,omega_I_prime] = min(prime_angle_miso,[],2);
            
            z = inv(symrots_true(omega_I_prime)).*orientation(rotation(misos_bcc_grains_omega_cutoff(omega_I_cutoff == k)),symmetries{1},1);
            
            primary_misos = [primary_misos;z];
        end
        
    end
    
    %% Handle returning the misorientations to primary orientation relationship:
    
    angle_getting_there = angle_outer(primary_misos,orientation(rotation(bcc_trans),symmetries{1},1));
    [omega_getting_there,omega_I_getting_there] = min(angle_getting_there,[],2);
    
    fcc2bcc_candidate = rotation(fcc2bcc_true(omega_I_getting_there)).*primary_misos;
    
    %Return the candidate rotation closest to the "original":
    
    candidate = mean(fcc2bcc_candidate);
    
    [M2, I2] = min(angle_outer(candidate*symrots_true,rotation(fcc2bcc_KS)));
    candidate = candidate*symrots_true(I2);    
    
    fcc2bcc_true = orientation(rotation(symrots_true*candidate),symmetries{1},1);
    
    fcc2bcc_basis_miso = orientation(rotation(fcc2bcc_true(1)),symmetries{1},symmetries{2});
    angle_planes(l) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,1,-1,symmetries{1})),Miller(0,1,1,symmetries{2}))) / degree;
    angle_directions(l) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,0,1,symmetries{1},'uvw')),Miller(1,1,-1,symmetries{1},'uvw'))) / degree;
    % Update the graph:
if cond
    addpoints(h_1(1),l,angle_planes(l));
    addpoints(h_1(2),l,angle_directions(l));
    drawnow
end   
    l = l+1;

    
end

bcc2fcc_candidate = orientation(rotation(inv(fcc2bcc_candidate.project2FundamentalRegion(fcc2bcc_KS))),symmetries{2},1);
bcc2fcc_odf = calcODF(orientation(rotation(bcc2fcc_candidate),symmetries{2},1),'halfwidth',1*degree);

angle_planes(21) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,1,-1,symmetries{1})),Miller(0,1,1,symmetries{2}))) / degree;
angle_directions(21) = min(angle(fcc2bcc_basis_miso * (symrots_true*Miller(1,0,1,symmetries{1},'uvw')),Miller(1,1,-1,symmetries{1},'uvw'))) / degree;

fcc2bcc = fcc2bcc_true;

end

