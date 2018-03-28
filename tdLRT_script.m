addpath('/home/Mike/bioStructureM/core');
addpath('/home/Mike/bioStructureM/accelerator');
addpath('/home/Mike/bioStructureM/atomselector');
addpath('/home/Mike/tdLRT_template');

Traj = getTraj('tdLRT_parameters.m');
AvgTraj = getAvgTraj('tdLRT_parameters.m',Traj);
Tc = getTcfromTraj('tdLRT_parameters.m',Traj);
AvgTc = getTcfromAvgTraj('tdLRT_parameters.m',AvgTraj);
writeTctoPDB('tdLRT_parameters.m','3N45_test_AvgTc.pdb');
[~,ICCres] = getICCfromTc('tdLRT_parameters.m',Tc);
[CS] = getCSfromICC('tdLRT_parameters.m',ICCres);
