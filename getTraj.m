function [Traj] = getTraj(parameter_filename)
run(parameter_filename);

force_matrix = getForceMatrix(ca,perturbed_atoms,theta_interval,phi_interval,impulse_force_magnitude,impulse_force_duration_time);
[AtomNum3,N_dir] = size(force_matrix);
N = AtomNum3/3;
N_t = length(time_sequence);
Traj = zeros(N*N_t,N_dir);

for i = 0:N_t-1
    total_response = getResponseFunction(ca,lowest_mode_number,highest_mode_number,solvent_friction,time_sequence(i+1),eigvalues,eigvectors,mass);
    displacement_mag = getDisplacement(ca,total_response,force_matrix,level);
    Traj(i*N+1:(i+1)*N,:) = displacement_mag;
end

save([output_prefix '_Traj.mat'],'Traj','-v7.3');

