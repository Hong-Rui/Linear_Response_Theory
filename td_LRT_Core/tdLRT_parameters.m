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
addpath('C:\Users\mick3101\OneDrive\useful_fuc.git\core');
addpath('C:\Users\mick3101\OneDrive\useful_fuc.git\atomselector');
addpath('C:\Users\mick3101\OneDrive\IF-tdLRT_template\tdLRT_template');

ca = readPDB('2BNE_init.pdb');
perturbed_atoms = atomselect('resname U5P',ca);

time_sequence = [0.01:0.01:1.0,1.1:0.1:10,11:1:100,110:10:1000];
solvent_friction = 30;
lowest_mode_number = 1;
highest_mode_number = length(ca);

theta_interval = 30;
phi_interval = 45;
impulse_force_magnitude = 100.0;
impulse_force_duration_time = 0.03;

level = 'atom';
output_prefix = '2BNE';

load([output_prefix '_eigvalues.mat']);
eigvalues = eigvalues(7:end);
load([output_prefix '_eigvectors.mat']);
eigvectors = eigvectors(:,7:end);

mass_text_filename = [output_prefix '.mass'];
[mass] = convert_mass(mass_text_filename);

%%%ICC parameters%%%
time_step = 0.01;
time_bound = 5;





