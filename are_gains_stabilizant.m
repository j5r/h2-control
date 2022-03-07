%
% By Junior R. Ribeiro, jrodrib@usp.br, 22-fev-2022
%
% Struct = are_gains_stabilizant(Struct)
%
%   This function computes the operator A1 Eq.3.12b
%   (https://link.springer.com/book/10.1007/b138575, p.34)
%   responsable for the MS-stability. This operator is appended to the
%   struct.
%

function Struct = are_gains_stabilizant(Struct)
%
fprintf('\n --> ARE GAINS STABILIZANT...\n')
validate_are_gains_stabilizant(Struct)
%
% lmis
Struct = are_lmi_gains_stabilizant(Struct);
%
% riccati
Struct = are_riccati_gains_stabilizant(Struct);
%
% doval
Struct = are_doval_gains_stabilizant(Struct);
fprintf('...DONE.\n');
end
%
%
%
%
function validate_are_gains_stabilizant(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
assert(isfield(Struct,'lmi_solution') || isfield(Struct,'riccati_solution'),...
    'The Structure does not have the field "K" (gains).');
end
%
%
%
%
function Struct = are_lmi_gains_stabilizant(Struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for lmi_solution
if isfield(Struct,'lmi_solution')
    n_states = size(Struct.valid_states, 1);
    n = size(Struct.A,1);
    %
    diagD = [];
    for ell = 1:n_states
        diagD = blkdiag( diagD, ...
            kron(Struct.lmi_solution.cloopA(:,:,ell), Struct.lmi_solution.cloopA(:,:,ell)) );
    end
    %
    operatorA = kron(Struct.augm_Prob', eye(n^2) ) * diagD;
    %
    Struct.lmi_solution.operatorA_Osvaldo = operatorA;
    Struct.lmi_solution.max_eigs_opA_Osvaldo = max(abs(eig( Struct.lmi_solution.operatorA_Osvaldo )));
end
end
%
%
%
%
function Struct = are_riccati_gains_stabilizant(Struct)
if isfield(Struct,'riccati_solution')
    n_states = Struct.N;
    n = size(Struct.A,1);
    %
    diagD = [];
    for ell = 1:n_states
        diagD = blkdiag( diagD, ...
            kron(Struct.riccati_solution.cloopA(:,:,ell), Struct.riccati_solution.cloopA(:,:,ell)) );
    end
    %
    operatorA = kron(Struct.Prob', eye(n^2) ) * diagD;
    %
    Struct.riccati_solution.operatorA_Osvaldo = operatorA;
    Struct.riccati_solution.max_eigs_opA_Osvaldo = max(abs(eig( Struct.riccati_solution.operatorA_Osvaldo )));
    %
end
end
%
%
%
%
function Struct = are_doval_gains_stabilizant(Struct)
if isfield(Struct,'doval_solution')
    n_states = Struct.N;
    n = size(Struct.A,1);
    %
    diagD = [];
    for ell = 1:n_states
        thetaHat = ell * (ell <= Struct.No) + (Struct.No+1)*(ell>Struct.No);
        diagD = blkdiag( diagD, ...
            kron(Struct.doval_solution.cloopA(:,:,ell,thetaHat), Struct.doval_solution.cloopA(:,:,ell,thetaHat)) );
    end
    %
    operatorA = kron(Struct.Prob', eye(n^2) ) * diagD;
    %
    Struct.doval_solution.operatorA_Osvaldo = operatorA;
    Struct.doval_solution.max_eigs_opA_Osvaldo = max(abs(eig( Struct.doval_solution.operatorA_Osvaldo )));
    %
end
end
