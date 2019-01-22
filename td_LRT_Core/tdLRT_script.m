addpath('/mnt/Tsunami_HHD/Mike/Modules/bioStructureM/core');
addpath('/mnt/Tsunami_HHD/Mike/Modules/bioStructureM/accelerator');
addpath('/mnt/Tsunami_HHD/Mike/Modules/bioStructureM/atomselector');
addpath('/mnt/Tsunami_HHD/Mike/Modules/Linear_Response_Theory/td_LRT_Core');

Traj = getTraj('tdLRT_parameters.m');
AvgTraj = getAvgTraj('tdLRT_parameters.m',Traj);
Tc = getTcfromTraj('tdLRT_parameters.m',Traj);
AvgTc = getTcfromAvgTraj('tdLRT_parameters.m',AvgTraj);
%writeTctoPDB('tdLRT_parameters.m','ATG4B_close_AvgTc.pdb');
[~,ICCres] = getICCfromTc('tdLRT_parameters.m',Tc);
[CS] = getCSfromICC('tdLRT_parameters.m',ICCres);
