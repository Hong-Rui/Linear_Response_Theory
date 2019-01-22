function[AvgTc] = getTcfromAvgTraj(parameter_filename,AvgTraj)
run(parameter_filename);

if ~exist('AvgTraj','var');
    load([output_prefix '_AvgTraj.mat']);
end

N = length(ca);
N_t = length(time_sequence);
AvgTc = zeros(N,1);
AvgTc(:) = 0.2;
dmag_old = AvgTraj(1:N);
for i = 1:N_t-1
    dmag = AvgTraj(i*N+1:(i+1)*N);
    delta_dmag = dmag - dmag_old;
    large_idx = find(delta_dmag>0);
    AvgTc(large_idx(:)) = time_sequence(i+1);
    dmag_old(large_idx(:)) = dmag(large_idx(:));
    clear large_idx;
end

save([output_prefix '_AvgTc.mat'],'AvgTc','-v7.3');
