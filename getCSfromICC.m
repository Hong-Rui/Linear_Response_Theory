function [CS] = getCSfromICC(parameter_filename,ICCres)
run(parameter_filename);

if ~exist('ICCres','var');
    load([output_prefix '_ICCres.mat']);
end

[N_res,~] = size(ICCres);
CS = zeros(N_res,1);
for i = 1:N_res
	CS(i) = max([ICCres(i,:),ICCres(:,i)']);
end

save([output_prefix '_CS.mat'],'CS','-v7.3');
