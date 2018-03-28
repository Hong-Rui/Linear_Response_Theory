function writeTctoPDB(parameter_filename,pdb_filename)
run(parameter_filename);

load([output_prefix '_AvgTc.mat']);

if strcmp(level,'atom')
	for i = 1:length(ca)
    	ca(i).bval = AvgTc(i);
	end
	createPDB(ca,pdb_filename);
elseif strcmp(level,'residue')
	N_res = length(Tc);
	for i=1:N_res
		logical_array = ([ca.internalResno]==i);
		residue = ca(logical_array);
		for j=[residue.atomno]
			ca(j).bval = AvgTc(i);
		end
	end
	createPDB(ca,pdb_filename);
end
