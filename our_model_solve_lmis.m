%
% By Junior R. Ribeiro, jrodrib@usp.br, 22-fev-2022
%
% Struct = our_model_solve_lmis(Struct)
%
%   This function solves a set of LMIs in order to obtain a set of
%   stabilizant gains K(:,:,thetaHat,rho,lambda).
%
%   Solve for R, W, G and F
%   [R(ell) - augm_pi(ell)*E(theta)*E(theta)', ...
%        A(theta)*G(thetaHat,rho,lambda) + B(theta)*F(thetaHat,rho,lambda);
%    G(thetaHat,rho,lambda)'*A(theta)'+F(thetaHat,rho,lambda)'*B(theta)',...
%        G(thetaHat,rho,lambda) + G(thetaHat,rho,lambda)' + OP_D(ell,R)]
%   > 0
%
%   [W(ell), ...
%        C(theta)*G(thetaHat,rho,lambda) + D(theta)*F(thetaHat,rho,lambda);
%    G(thetaHat,rho,lambda)'*C(theta)'+F(thetaHat,rho,lambda)'*D(theta)',...
%        G(thetaHat,rho,lambda) + G(thetaHat,rho,lambda)' + OP_D(ell,R)]
%   > 0
%
%   where OP_D(ell,R) = sum_{ell2} augm_Prob(ell2,ell) * R(ell).
%   K(:,:,thetaHat,rho,lambda) = F(thetaHat,rho,lambda)*inv(G(thetaHat,rho,lambda)).
%

function Struct = our_model_solve_lmis(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
fprintf('...OUR MODEL SOLVE LMIS...\n');
%
T = Struct.T;
N = Struct.N;
No = Struct.No;
n_states = size(Struct.valid_states,1);
%
constraints = [];
yalmip('clear');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Allocating LMI Vars &&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computing the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% objective function &&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% var G positive-def.
%
n_R = [size(Struct.A,1), size(Struct.A,1)];
n_W = [size(Struct.C,1), size(Struct.C,1)];
n_G = [size(Struct.A,1), size(Struct.A,1)];
n_F = [size(Struct.B,2), size(Struct.B,1)];
objective_function = 0;
for ell = 1:n_states
    LMI_R{ell} = sdpvar(n_R(1), n_R(2));
    LMI_W{ell} = sdpvar(n_W(1), n_W(2));
    objective_function = objective_function + trace(LMI_W{ell});
end
%
for thetaHat = 1:N
    for rho = 1:T
        for lambda = 1:No
            intersection_with_valids = intersect(Struct.valid_states(:,2:4), [thetaHat, rho-1, lambda], 'rows');
            if ~isempty(intersection_with_valids) 
                % allocating only if it is valid
                LMI_G{thetaHat,rho,lambda} = sdpvar(n_G(1), n_G(2));                
                LMI_F{thetaHat,rho,lambda} = sdpvar(n_F(1), n_F(2));
                % var G being positive-definite
                constraints = [constraints, LMI_G{thetaHat,rho,lambda} >= 0];
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraints of the LMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
%%%%%%%%%%%%% In order to turn off some LMIs, we compute a scalar to
%%%%%%%%%%%%% multiply them. If some state 'ell' is never visited, its
%%%%%%%%%%%%% propagate distribution is always zero, and so, we multiply
%%%%%%%%%%%%% the respective LMI by the sum of its probability through
%%%%%%%%%%%%% time.
%%%%%%%%%%%%%  
summation_P = zeros(n_states, n_states);
for ell = 0:n_states
    summation_P =  summation_P + Struct.augm_Prob^ell ;
end
 weights = (Struct.augm_pi(:)' * summation_P);
%%%%%%%%%%%%% 
%%%%%%%%%%%%% 
for ell = 1:n_states
    [theta,thetaHat,rho,lambda] = map1to4(ell, Struct);
    %
    % operator D_ell(R), summation
    D_R_ell = zeros(n_R(1),n_R(1));
    for i = 1:n_states
        D_R_ell = D_R_ell + Struct.augm_Prob(i,ell) * LMI_R{i};
    end
    %
    blk__22 = LMI_G{thetaHat,rho+1,lambda} + LMI_G{thetaHat,rho+1,lambda}' - D_R_ell;
    %
    % LMI blocks regarding R
    blk_R11 = LMI_R{ell} - Struct.augm_pi(ell) * Struct.E(:,:,theta) * Struct.E(:,:,theta)';
    blk_R12 = Struct.A(:,:,theta) * LMI_G{thetaHat,rho+1,lambda} + Struct.B(:,:,theta) *  LMI_F{thetaHat,rho+1,lambda};    
    % LMI blocks regarding W
    blk_W11 = LMI_W{ell};
    blk_W12 = Struct.C(:,:,theta) * LMI_G{thetaHat,rho+1,lambda} + Struct.D(:,:,theta) *  LMI_F{thetaHat,rho+1,lambda};
    %    
    %%%%%%%%%%%% adding constraints
    constraints = [constraints, weights(ell) * [ blk_R11,  blk_R12; ...
                                          blk_R12',  blk__22] >= 0]; 
    constraints = [constraints,  [ blk_W11,  blk_W12; ...
                                  blk_W12', blk__22] >= 0];    
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
solvesdp(constraints,objective_function);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Retrieving solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
trace_W = 0;
for ell = 1:n_states
    [~,thetaHat,rho,lambda] = map1to4(ell,Struct);
    R(:,:,ell) = value(LMI_R{ell});
    W(:,:,ell) = value(LMI_W{ell});
    trace_W = trace_W + trace(W(:,:,ell));
    G(:,:,thetaHat,rho+1,lambda) = value(LMI_G{thetaHat,rho+1,lambda});
    F(:,:,thetaHat,rho+1,lambda) = value(LMI_F{thetaHat,rho+1,lambda});
    K(:,:,thetaHat,rho+1,lambda) = F(:,:,thetaHat,rho+1,lambda)*inv(G(:,:,thetaHat,rho+1,lambda));
end
s.F = F;
s.G = G;
s.R = R;
s.W = W;
lmi_.FGRW = s;
lmi_.K = K;
lmi_.trace_W = trace_W;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computing closed loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
n = size(Struct.A,1);
r = size(Struct.C,1);
cloopA = zeros(n, n, n_states);
cloopC = zeros(r, n, n_states);
for ell = 1:n_states
    [theta,thetaHat,rho,lambda] = map1to4(ell, Struct);
    cloopA(:,:,ell) = Struct.A(:,:,theta) + Struct.B(:,:,theta) * K(:,:,thetaHat,rho+1,lambda);
    cloopC(:,:,ell) = Struct.C(:,:,theta) + Struct.D(:,:,theta) * K(:,:,thetaHat,rho+1,lambda);
end
lmi_.cloopA = cloopA;
lmi_.cloopC = cloopC;
Struct.lmi_solution = lmi_;
%
fprintf('...DONE.\n\n');
end
