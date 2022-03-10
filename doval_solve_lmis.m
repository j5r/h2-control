%
% By Junior R. Ribeiro, jrodrib@usp.br, 22-fev-2022
%
% Struct = doval_solve_lmis(Struct)
%
%   This function solves a set of LMIs in order to obtain a set of
%   stabilizant gains K(:,:,thetaHat), according DoVal, 2002.
%   https://doi.org/10.1016/S0005-1098(01)00210-2, Eqs 12-13.
%
%   Solve for R, W, G and F
%   [R(theta) - pi(theta)*E(theta)*E(theta)', ...
%        A(theta)*G(thetaHat) + B(theta)*F(thetaHat);
%    G(thetaHat)'*A(theta)'+F(thetaHat)'*B(theta)',...
%        G(thetaHat) + G(thetaHat)' + OP_D(theta,R)]
%   > 0
%
%   [W(theta), ...
%        C(theta)*G(thetaHat) + D(theta)*F(thetaHat);
%    G(thetaHat)'*C(theta)'+F(thetaHat)'*D(theta)',...
%        G(thetaHat) + G(thetaHat)' + OP_D(theta,R)]
%   > 0
%
%   where OP_D(i,R) = sum_{j} Prob(j,i) * R(j).
%   K(:,:,thetaHat) = F(thetaHat)*inv(G(thetaHat)).
%


function Struct = doval_solve_lmis(Struct)
%
fprintf('\n --> DO VAL: SOLVE LMIS...\n');
%
T = Struct.T;
N = Struct.N;
No = Struct.No;
%
constraints = [];
yalmip('clear');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Allocating LMI Vars &&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computing the objective function &&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% making variable G positive-definite.
%
n_R = [size(Struct.A,1), size(Struct.A,1)]; % dimension of the variables
n_W = [size(Struct.C,1), size(Struct.C,1)];
n_G = [size(Struct.A,1), size(Struct.A,1)];
n_F = [size(Struct.B,2), size(Struct.B,1)];
objective_function = 0;
for theta = 1:N
    LMI_R{theta} = sdpvar(n_R(1), n_R(2));
    LMI_W{theta} = sdpvar(n_W(1), n_W(2));
    objective_function = objective_function + trace(LMI_W{theta});
end
%
%
for thetaHat = 1:No+1
    LMI_G{thetaHat} = sdpvar(n_G(1), n_G(2));
    LMI_F{thetaHat} = sdpvar(n_F(1), n_F(2));
    % var G being positive-definite
    constraints = [constraints, LMI_G{thetaHat} >= 0];
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraints of the LMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%

for theta = 1:N
    % operator D_ell(R), summation
    D_R_theta = zeros(n_R(1),n_R(2));
    for j = 1:N
        D_R_theta = D_R_theta + Struct.Prob(j,theta) * LMI_R{j};
    end
    for thetaHat = 1:No+1
        blk__22 = LMI_G{thetaHat} + LMI_G{thetaHat}' - D_R_theta;
        %
        % LMI blocks regarding R
        blk_R11 = LMI_R{theta} - Struct.pi(theta) * Struct.E(:,:,theta) * Struct.E(:,:,theta)';
        blk_R12 = Struct.A(:,:,theta) * LMI_G{thetaHat} + Struct.B(:,:,theta) *  LMI_F{thetaHat};
        %
        % LMI blocks regarding W
        blk_W11 = LMI_W{theta};
        blk_W12 = Struct.C(:,:,theta) * LMI_G{thetaHat} + Struct.D(:,:,theta) *  LMI_F{thetaHat};
        
        %
        % adding constraints
        constraints = [constraints,  [ blk_R11,  blk_R12; ...
            blk_R12',  blk__22] >= 0];
        %
        constraints = [constraints,  [ blk_W11,  blk_W12; ...
            blk_W12', blk__22] >= 0];
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVING
%
solvesdp(constraints,objective_function);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RETRIEVING SOLUTION
%
trace_W = 0;
for theta = 1:N
    R(:,:,theta) = value(LMI_R{theta});
    W(:,:,theta) = value(LMI_W{theta});
    trace_W = trace_W + trace(W(:,:,theta));
end
%
for thetaHat = 1:No+1
    G(:,:,thetaHat) = value(LMI_G{thetaHat});
    F(:,:,thetaHat) = value(LMI_F{thetaHat});
    K(:,:,thetaHat) = F(:,:,thetaHat)*inv(G(:,:,thetaHat));
end
%
%
s.F = F;
s.G = G;
s.R = R;
s.W = W;
doval_struct.FGRW = s;
doval_struct.K = K;
doval_struct.trace_W = trace_W;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTING CLOSED LOOP MATRICES
%
n = size(Struct.A,1);
r = size(Struct.C,1);
cloopA = zeros(n, n, N, No+1);
cloopC = zeros(r, n, N, No+1);
for theta = 1:N
    for thetaHat = 1:No+1
        cloopA(:,:,theta,thetaHat) = Struct.A(:,:,theta) + Struct.B(:,:,theta) * K(:,:,thetaHat);
        cloopC(:,:,theta,thetaHat) = Struct.C(:,:,theta) + Struct.D(:,:,theta) * K(:,:,thetaHat);
    end
end
doval_struct.cloopA = cloopA;
doval_struct.cloopC = cloopC;
doval_struct.help = {'1# cloopA(:,:,theta,thetaHat) = A(:,:,theta) + B(:,:,theta) * K(:,:,thetaHat)';
    '2# cloopC(:,:,theta,thetaHat) = C(:,:,theta) + D(:,:,theta) * K(:,:,thetaHat)';
    '3# [theta] is in [1:N]; [thetaHat] is in [1:No+1]'};
Struct.doval_solution = doval_struct;
%
fprintf('...DONE.\n');
end
