function[Tc] = getTcfromTraj(parameter_filename,Traj)
run(parameter_filename);

if ~exist('Traj','var')
    load([output_prefix '_Traj.mat']);
end

[N_traj,N_dir] = size(Traj);
N_t = length(time_sequence);
N = N_traj/N_t;
Tc = zeros(N,N_dir);
Tc(:) = 0.01;
dmag_old = Traj(1:N,:);
for i = 1:N_t-1
    dmag = Traj(i*N+1:(i+1)*N,:);
    delta_dmag = dmag - dmag_old;
    large_idx = find(delta_dmag>0);
    Tc(large_idx(:)) = time_sequence(i+1);
    dmag_old(large_idx(:)) = dmag(large_idx(:));
    clear large_idx;
end

save([output_prefix '_Tc.mat'],'Tc','-v7.3');
