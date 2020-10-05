function [ebsd_aus,devis] = recaus(ebsd_aus,grains,symmetries,fcc2bcc,ib)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Symmetry order correction according to Morito et al.:
cs_symmetries = symmetries{1}.properGroup.rot;

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
% as by Morito et al.:
fcc2bcc_true = orientation(symrots_true*rotation(fcc2bcc),symmetries{1},1);

%Set the misorientation corresponding to "bcc to fcc" transformation
bcc2fcc = symrots_true * inv(fcc2bcc_true(1));
bcc_trans = orientation(bcc2fcc,symmetries{2},symmetries{1});

h = waitbar(0,'Reconstructing austenite orientations. Please wait...','Name','Austenite reconstruction in progress',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        
setappdata(h,'canceling',0)

for l = 1:length(ib);
    
    if getappdata(h,'canceling')
        break
    end
    
    waitbar(l/length(ib),h)
    
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
delete(h)
end

