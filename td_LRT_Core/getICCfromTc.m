function [ICC,ICCres] = getICCfromTc(parameter_filename,Tc)
run(parameter_filename);

if ~exist('Tc','var');
    load([output_prefix '_Tc.mat']);
end

[N,N_dir] = size(Tc);
directions = N_dir/length(perturbed_atoms);
AtomNumPerResidue = getAtomNumPerRes(ca);
N_res = length(AtomNumPerResidue);

ICC = zeros(N,N);
ICCres = zeros(N_res,N_res);
initial_time = time_sequence(1);

for i_dir = 1:N_dir
    donor_atomid = zeros(N,1);
    acceptor_atomid = zeros(N,1);
    N_donor = 0;

    for i=1:N
        if Tc(i,i_dir) == initial_time
            N_donor = N_donor + 1;
            donor_atomid(N_donor) = i;
        end
    end

    if N_donor == 0
        continue;
    end

    for time = time_sequence(2:end)
        N_acceptor = 0;
        acceptor_atomid = zeros(N,1);
        for j = 1:N
            if Tc(j,i_dir) == time
                N_acceptor = N_acceptor + 1;
                acceptor_atomid(N_acceptor) = j;
                current_acceptor = ca(j);
                for i_donor = 1:N_donor
                    current_donor = ca(donor_atomid(i_donor));
                    OD_vector = current_donor.coord;
                    OA_vector = current_acceptor.coord;
                    OP_vector = perturbed_atoms(ceil(i_dir/directions)).coord;
                    PD_vector = OD_vector - OP_vector;
                    PA_vector = OA_vector - OP_vector;
                    DA_vector = OA_vector - OD_vector;
                    
                    if (norm(DA_vector) > 3.8) & ((PD_vector'*PA_vector) < 0)
                        continue;
                    end
            
                    current_donor_atomid = current_donor.atomno;
                    current_acceptor_atomid = current_acceptor.atomno;
                    current_donor_resid = current_donor.internalResno;
                    current_acceptor_resid = current_acceptor.internalResno;
                    ICC(current_donor_atomid,current_acceptor_atomid) = ICC(current_donor_atomid,current_acceptor_atomid) + 1;
                    ICCres(current_donor_resid,current_acceptor_resid) = (ICCres(current_donor_resid,current_acceptor_resid) + 1)/sqrt(AtomNumPerResidue(current_donor_resid)*AtomNumPerResidue(current_acceptor_resid));
                end
            end
        end
        
        N_donor = N_acceptor;
        donor_atomid = acceptor_atomid;
        if N_donor == 0
            break;
        end
    end
end

ICC = ICC ./ N_dir;
ICCres = ICCres ./ N_dir;

save([output_prefix '_ICC.mat'],'ICC','-v7.3');
save([output_prefix '_ICCres.mat'],'ICCres','-v7.3');
