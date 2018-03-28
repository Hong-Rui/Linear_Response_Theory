function [total_response] = getResponseFunction(ca,lowest_mode_number,highest_mode_number,solvent_friction,time_step,eigvalues,eigvectors,mass)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:
%   ca(PDB Structure)
%   eigenvalues matrix(top 6 invariant values removed)
%   eigenvectors matrix(top 6 invariant colums removed)
%   lowest_mode_number used to reform covariance matrix
%   highest_mode_number used to reform covariance matrix
%   beta(solvent friction for Langivin damping)
%   time_step(current time step) at picosecond unit
%
%   All values are changed into SI units.
%
% output:
%   total_response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%------get required information-----------------------
N = length(ca);
mass_3N_vector = reshape(repmat(mass,1,3)',3*N,1);
mass_matrix = sqrt(mass_3N_vector)*sqrt(mass_3N_vector)';
%-----------------------------------------------------
eigvalues = eigvalues(lowest_mode_number:highest_mode_number);
M = length(eigvalues);
eigvectors = eigvectors(:,lowest_mode_number:highest_mode_number);
%1/mobility, cm^-1 -->[rad/s], angular_frequency in unit of [rad/s]
solvent_friction = solvent_friction*(3*10^10)*2*pi;
angular_frequency = sqrt(eigvalues./(4.185*10^26));
time_step = time_step*10^-12;

%%
%--------calculate time-dependent response for each atom------------------
response = zeros(M,1);
for m=1:M
	beta_1m = sqrt(solvent_friction^2-4*angular_frequency(m)^2);
	if solvent_friction >= 2*angular_frequency(m)
		mode_response = 2*10^10/(6.02*10^23*10^3)/beta_1m*exp((-1)*solvent_friction*time_step/2)*sinh(beta_1m*time_step/2);
	else
		beta_1m=imag(beta_1m);
		mode_response = 2*10^10/(6.02*10^23*10^3)/beta_1m*exp((-1)*solvent_friction*time_step/2)*sin(beta_1m*time_step/2);
	end
	response(m) = mode_response;
end

total_response = zeros(3*N,M);
for i =1:M
    total_response(:,i) = eigvectors(:,i)*response(i);
end
total_response = total_response * eigvectors';
total_response = total_response./mass_matrix;
%-------------------------------------------------------------
