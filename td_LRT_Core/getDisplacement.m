function [displacement_mag] = getDisplacement(ca,total_response,force_matrix,level)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:
%   ca, total_response, force_matrix, level
%
% output:
%   displacement_mag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%-------check if flag exists--------------------------
if ~exist('level','var')
    level = 'atom';
end
%-----------------------------------------------------
[AtomNum3,N_dir] = size(force_matrix);
N = AtomNum3/3;
displacement  =  total_response*force_matrix;
displacement_mag = zeros(N,N_dir);
sqr_displacement = displacement.^2;
for i=0:N-1
    displacement_mag(i+1,:) = sqrt(sum(sqr_displacement(i*3+1:i*3+3,:)));
end
displacement_mag = displacement_mag .* (10^10);

%------------------------------------------------------------------------

%%
%-----------calculate response in residue level-------------------------- 
if strcmp(level,'residue')
    AtomNumPerResidue = getAtomNumPerRes(ca);
    N_res = length(AtomNumPerResidue);
    res_displacement = zeros(N_res,N_dir);

    current_index = 0;
    for i=1:N_res
        current_res_length = AtomNumPerResidue(i);
        res_displacement(i,:) = sum(displacement_mag((current_index + 1):(current_index + current_res_length),:))/current_res_length;
        current_index = current_index + current_res_length;
    end
    displacement_mag = res_displacement;
end
%--------------------------------------------------------------------------