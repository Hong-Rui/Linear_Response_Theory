function[AvgTraj] = getAvgTraj(parameter_filename,Traj)
run(parameter_filename);

if ~exist('Traj','var')
    load([output_prefix '_Traj.mat']);
end

AvgTraj = mean(Traj,2);
save([output_prefix '_AvgTraj.mat'],'AvgTraj','-v7.3');
