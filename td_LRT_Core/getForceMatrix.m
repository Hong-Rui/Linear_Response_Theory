function [force_matrix] = getForceMatrix(ca,perturbed_atoms,theta_interval,phi_interval,impulse_force_magnitude,impulse_force_duration_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: 
%    ca, perturbed_atoms, theta_interval, phi_interval, force_mag, duration_time
% ouput:
%   perturbation force matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------
N = length(ca);
N_perturbed_atoms = length(perturbed_atoms);
%Number of force directions x perturbed sites
N_dir = 25*N_perturbed_atoms;
F = zeros(3*N,N_dir);
i_dir = 1;
for i_pa=1:N_perturbed_atoms
	i_atom = perturbed_atoms(i_pa).atomno-1;

%define impulse forces matrix
	for theta = 0:theta_interval:90
        	for phi = 0:phi_interval:315
                if(theta == 0 && phi > 0)
                    continue;
                end

            i_x = 3*i_atom+1; i_y = 3*i_atom+2; i_z = 3*i_atom+3;
            F(i_x,i_dir) = impulse_force_magnitude*sind(theta)*cosd(phi);
            F(i_y,i_dir) = impulse_force_magnitude*sind(theta)*sind(phi);
            F(i_z,i_dir) = impulse_force_magnitude*cosd(theta);         
            i_dir = i_dir + 1;
		
        	end
	end
end
%----------------------------------------------------------------
%impulse force duration time [ps]
impulse_force_duration_time = impulse_force_duration_time*10^-12;
force_matrix = F*impulse_force_duration_time;