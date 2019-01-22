%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check varaible names of eigenvalues and eigenvectors matrix FIRST!!%%%%
%   add useful_func package
%   read PDB file select atoms for perturbation
%   (full path is advised)
%
%   time_sequence is used for calculating the trajectories
%   beta is the solvent friction for Langivin damping
%   lowest_mode_number is modes to reform covariance matrix
%   highest_mode_number is modes to reform covariance matrix
%   theta_interval and phi_interval are to tune the pertubating directions
%   select inpluse force manitude and duration time 
%
%
%   level can be 'atom' or 'residue'. 
%   Choose 'atom' to calculate response in atom level, 
%   or 'residue' to coarse-grained to residue level displacement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/mnt/Tsunami_HHD/Mike/Modules/bioStructureM/core');
addpath('/mnt/Tsunami_HHD/Mike/Modules/bioStructureM/accelerator')
addpath('/mnt/Tsunami_HHD/Mike/Modules/bioStructureM/atomselector');
addpath('/mnt/Tsunami_HHD/Mike/Modules/Linear_Response_Theory/td_LRT_Core');

ca = readPDB('ATG4B_close_vacuum.pdb');
ca = assignMass(ca);
perturbed_atoms = atomselect('(resid 74 142 258 260 278 280) and name CA',ca);

time_sequence = [0.2:0.2:4.0];
solvent_friction = 30;
lowest_mode_number = 1;
highest_mode_number = 3*length(ca)-6;

theta_interval = 30;
phi_interval = 45;
impulse_force_magnitude = 4;
impulse_force_duration_time = 0.03;

level = 'atom';
output_prefix = 'ATG4B_close';

load('eigval.mat');
load('eigvec.mat');
eigvalues = eigval(7:end);
eigvectors = eigvec(:,7:end);
