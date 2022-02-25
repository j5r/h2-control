function Struct = solve_lmi_proposition(Struct,epsilon_Z)
if nargin==1
    epsilon_Z = 1e-5;
end
epsilon = 1e-5;

% -- Retrieving data from the Struct
A = Struct.A;
B = Struct.B;
C = Struct.C;
D = Struct.D;
E = Struct.E;
Prob = Struct.Prob;
No = Struct.No;
T = Struct.T;
N = Struct.N;

% Example pi0 = [ .3, .5, .2], with N = 3
init_distrib = Struct.init_distrib; 

% Example pihat0 = [0, .4, .6], with No = 1
init_distrib_tHat = Struct.init_distrib_tHat; 


% -- Computing the valid states
valid_states = []; % [theta, thetaHat, rho, lambda, counter, is_observable]
counter = 1;
for theta = 1:N
    for thetaHat = 1:N
        for rho = 0:T-1
            for lambda = 1:No
                indicator_obs = (theta<=No) * (thetaHat==theta) *...
                    (lambda<=No) * (rho==0);
                indicator_unobs = (theta>No) * (thetaHat>No) *...
                    (lambda<=No) * (rho>0);
                if (indicator_obs + indicator_unobs) > 0
                    valid_states = [valid_states;
                        theta, thetaHat, rho, lambda, counter, indicator_obs];
                    counter = counter + 1;
                end
            end
        end
    end
end
Struct.valid_states_size4 = valid_states;
n_valid_states = size(valid_states, 1);

% ---------------------------------------------------------------------
% -- DEFINING LMIs.
n = size(A,1);
m = size(B,2);

% -- LMIs' Variables R and W
for ell = 1:n_valid_states
    VAR_R{ell} = sdpvar(n,n);
    VAR_W{ell} = sdpvar(size(C,1),size(C,1));
    VAR_Z{ell} = sdpvar(n,n);
end

% -- LMIs' Variables G and F
for thetaHat = 1:N
    for rho = 0:T-1
        for lambda = 1:No
        VAR_F{thetaHat,rho+1,lambda} = sdpvar(m,n);
        VAR_G{thetaHat,rho+1,lambda} = sdpvar(n,n);
        end
    end
end

% -- LMIs' Variables Initial_distribution
%    (conditioned to theta and thetaHat in U)
for theta = No+1:N
    for thetaHat = No+1:N
        VAR_IN_DISTR{theta,thetaHat} = sdpvar(1,1);
    end
end


% -- LMIs' Variables ThetaHat_transition_distributions
%    (conditioned to theta and thetaHat in U)
for thetaHat = No+1:N
    for rho = 2:T % rho=1 means observed case
        for lambda = 1:No
            VAR_THAT_PROB{thetaHat,rho,lambda} = sdpvar(1,1);
        end
    end
end


% ---------------------------------------------------------------------
% -- BUILDING LMIS
for ell = 1:n_valid_states
    theta = valid_states(ell, 1);
    thetaHat = valid_states(ell, 2);
    rho = valid_states(ell, 3);
    lambda = valid_states(ell, 4);
    % index = valid_states(ell, 5);
    is_observed = valid_states(ell, 6);
    
    % -- Block 22 is used in LMI_R and LMI_W
    blk_22 = VAR_G{thetaHat,rho+1,lambda} + VAR_G{thetaHat,rho+1,lambda}' - VAR_Z{ell};
    
    % -- Block 12
    blk_12_R = A(:,:,theta) * VAR_G{thetaHat,rho+1,lambda} + B(:,:,theta) * VAR_F{thetaHat,rho+1,lambda};
    blk_12_W = C(:,:,theta) * VAR_G{thetaHat,rho+1,lambda} + D(:,:,theta) * VAR_F{thetaHat,rho+1,lambda};
    
    % -- Block 11
    if is_observed
        blk_11_R = VAR_R{ell} - init_distrib(theta)*...
            E(:,:,theta)*E(:,:,theta)' *...
            (rho==0); % validating initial distribution
    else
        blk_11_R = VAR_R{ell} - init_distrib(theta)*...
            VAR_IN_DISTR{theta,thetaHat} * E(:,:,theta)*E(:,:,theta)' *...
            (rho==1) * (lambda==1); % validating initial distribution
    end
    blk_11_W = VAR_W{ell};
    
    LMI_R{ell} = [blk_11_R, blk_12_R; blk_12_R', blk_22];
    LMI_W{ell} = [blk_11_W, blk_12_W; blk_12_W', blk_22];
end


% -- LMIs' regarding Z_ell and all the transition probabilities of thetaHat
%    ( Operator calligraph_D over R is replaced by Z_ell )
%    ( VAR_THAT_PROB{thetaHat,rho,lambda} epsilon_Z )
%    --
%    [Z_ell,              horizontal_block;        ]
%    [horizontal_block', 2*eye(2*n*n_valid_states) ]
I = eye(n,n);

for ell = 1:n_valid_states
    theta1 = valid_states(ell, 1);
    thetaHat1 = valid_states(ell, 2);
    rho1 = valid_states(ell, 3);
    lambda1 = valid_states(ell, 4);
    % index1 = valid_states(ell, 5);
    % is_observed1 = valid_states(ell, 6);
    
    % -- Building the horizontal_block
    horizontal_block = [];
    counter = 0;
    for row = 1:n_valid_states
        theta2 = valid_states(row, 1);
        thetaHat2 = valid_states(row, 2);
        rho2 = valid_states(row, 3);
        lambda2 = valid_states(row, 4);
        % index2 = valid_states(row, 5);
        % is_observed2 = valid_states(row, 6);
        
        % -- Indicators
        ind_1 = (theta2==thetaHat2)*(theta2==lambda2)*(theta2<=No)*(rho2==0);
        ind_2 = (theta2>No)*(thetaHat2>No)*(lambda1==lambda2)*(rho2==min(T-1,rho1+1));
        
        if ind_1
            horizontal_block = [horizontal_block, ...
                epsilon_Z * Prob(theta1,theta2) * I, VAR_R{row}/epsilon_Z];
            counter = counter + 1;
        end
        if ind_2
            horizontal_block = [horizontal_block, ...
                epsilon_Z * Prob(theta1,theta2) * ...
                VAR_THAT_PROB{thetaHat2,rho2+1,lambda2} * I,...
                VAR_R{row}/epsilon_Z];
            counter = counter + 1;
        end
    end
    LMI_Z{ell} = [VAR_Z{ell},               horizontal_block;
        horizontal_block', 2*eye(2*n*counter)];
end



% -- OBJECTIVE FUNCTION
COST = 0;
for ell = 1:n_valid_states
    COST = COST + trace(VAR_W{ell});
end

% -- CONSTRAINTS
constraints = [];

% ----------------------------------------
% -- Probabilities are set to be positive
for theta = No+1:N
    for thetaHat = No+1:N
        constraints = [constraints, VAR_IN_DISTR{theta,thetaHat} >= 0];
    end
end

% -- Probabilities are set to be positive
for thetaHat = No+1:N
    for rho = 2:T % rho=1 means observed case
        for lambda = 1:No
            constraints = [constraints, VAR_THAT_PROB{thetaHat,rho,lambda} >= 0];
        end
    end
end

% ----------------------------------------
% -- The variables G are set to be positive definite (because of its inverse)
for thetaHat = 1:N
    for rho = 1:T
        for lambda = 1:No
        constraints = [constraints, VAR_G{thetaHat,rho,lambda} >=0] ;
        end
    end
end

% -- The three LMIs
for ell = 1:n_valid_states
    constraints = [constraints, LMI_R{ell} >= 0];
    constraints = [constraints, LMI_W{ell} >= 0];
    constraints = [constraints, LMI_Z{ell} >= 0];
end


% -- The summation of the distributions VAR_IN_DISTR must sum up to 1
for rho = 2:T % rho=1 means observed case
    for lambda = 1:No
        summation = 0;
        for thetaHat = No+1:N
            summation = summation + VAR_THAT_PROB{thetaHat,rho,lambda};
        end
        constraints = [ constraints, epsilon + 1 - summation >=0];
        constraints = [ constraints, epsilon - 1 + summation >=0];
    end
end


% -- The summation of the distributions VAR_IN_DISTR must sum up to 1
for theta = No+1:N
    summation = 0;
    for thetaHat = No+1:N
        summation = summation + VAR_IN_DISTR{theta,thetaHat} ;
    end
    constraints = [ constraints, epsilon + 1 - summation >=0];
    constraints = [ constraints, epsilon - 1 + summation >=0];
end


% -- Solving the LMIs
solvesdp(constraints, COST);

% -- Recovering G and F, and computing the gains
for thetaHat = 1:N
    for rho = 1:T
        for lambda = 1:No
        G(:,:,thetaHat,rho,lambda) = value(VAR_G{thetaHat,rho,lambda});
        F(:,:,thetaHat,rho,lambda) = value(VAR_F{thetaHat,rho,lambda});
        % -- Computing the gains
        K(:,:,thetaHat,rho,lambda) = F(:,:,thetaHat,rho,lambda) * inv(G(:,:,thetaHat,rho,lambda));
        end
    end
end

% -- Recovering the initial distribution of
for theta = No+1:N
    for thetaHat = No+1:N
        init_distrib_tHat(theta,thetaHat) = value(VAR_IN_DISTR{theta,thetaHat});
    end
end
Struct.init_distrib_tHat = init_distrib_tHat;


% -- Recovering W and R
TRACE_W = 0;
for i = 1:n_valid_states
    W(:,:,i) = value(VAR_W{i});
    R(:,:,i) = value(VAR_R{i});
    Z(:,:,i) = value(VAR_Z{i});
    TRACE_W = TRACE_W + trace(W(:,:,i));
end

fprintf('\n --> TRACE_W from LMIs %f\n',TRACE_W);

for thetaHat = No+1:N
    for rho = 2:T % rho=1 means observed case
        for lambda = 1:No
            transition_that(thetaHat,rho,lambda) = value(VAR_THAT_PROB{thetaHat,rho,lambda});
        end
    end
end


% -- Building the transition probability matrix, based on transition_that
for row = 1:n_valid_states
    theta1 = valid_states(row, 1);
    thetaHat1 = valid_states(row, 2);
    rho1 = valid_states(row, 3);
    lambda1 = valid_states(row, 4);
    % index1 = valid_states(row, 5);
    % is_observed1 = valid_states(row, 6);
    
    for col = 1:n_valid_states
        theta2 = valid_states(col, 1);
        thetaHat2 = valid_states(col, 2);
        rho2 = valid_states(col, 3);
        lambda2 = valid_states(col, 4);
        % index2 = valid_states(col, 5);
        % is_observed2 = valid_states(col, 6);
        
        % -- Indicators
        ind_1 = (theta2==thetaHat2)*(theta2==lambda2)*(theta2<=No)*(rho2==0);
        ind_2 = (theta2>No)*(thetaHat2>No)*(lambda1==lambda2)*(rho2==min(T-1,rho1+1));
        
        Q(row,col) = Prob(theta1, theta2) * (ind_1 + transition_that(thetaHat2,rho2+1,lambda2)*ind_2);
    end
end
Struct.augm_Prob = Q;


for ell = 1:n_valid_states
    theta = valid_states(ell, 1);
    thetaHat = valid_states(ell, 2);
    rho = valid_states(ell, 3);
    lambda = valid_states(ell, 4);
    % index = valid_states(ell, 5);
    % is_observed = valid_states(ell, 6);
    
    % -- Indicators
    ind_1 = (theta==thetaHat)*(theta==lambda)*(theta<=No)*(rho==0);
    ind_2 = (theta>No)*(thetaHat>No)*(lambda==1)*(rho==1);
    Struct.augm_init_distrib(ell) = init_distrib(theta)*...
        (ind_1 + init_distrib_tHat(theta,thetaHat)*ind_2);
end


Struct.lmi_trace_W = TRACE_W;

sstruct.F = F;
sstruct.G = G;
sstruct.W = W;
sstruct.R = R;
sstruct.Z = Z;
sstruct.transition_that = transition_that;
sstruct.init_distrib_tHat = init_distrib_tHat;
Struct.FGRW = sstruct;
Struct.K = K;

z = 1;
w = 1;
r = 1;
for ell = 1:n_valid_states
    z = z * all(eig(Z(:,:,ell))>0);
    w = w * all(eig(W(:,:,ell))>0);
    r = r * all(eig(R(:,:,ell))>0);
end
g = 1;
for thetaHat = 1:N
    for rho = 1:T
        for lambda = 1:No
        
        g = g * all(eig(G(:,:,thetaHat,rho,lambda)));
        end
    end
end
fprintf('\nAre Positive-Definite? Z = %d  W = %d  R = %d G = %d\n',z,w,r,g);


validate_lmis(Struct);
end






function lmis_are_satisfied = validate_lmis(S)

F = S.FGRW.F;
G = S.FGRW.G;
W = S.FGRW.W;
R = S.FGRW.R;
Z = S.FGRW.Z;

eig_lmi1 = [];
eig_lmi2 = [];

for ell = 1:size(S.valid_states_size4, 1)
    %theta,thetahat,rho,lambda,index
    theta = S.valid_states_size4(ell,1);
    thetaHat = S.valid_states_size4(ell,2);
    rho = S.valid_states_size4(ell,3);
    lambda = S.valid_states_size4(ell,4);
    
    
    BLK22 = G(:,:,thetaHat,rho+1) + G(:,:,thetaHat,rho+1)' - Z(:,:,ell); % bloco22 das 2 LMIs
    
    % -- R
    BLK11R = R(:,:,ell) - S.augm_init_distrib(ell) * S.E(:,:,theta) * S.E(:,:,theta)';
    BLK12R = S.A(:,:,theta)*G(:,:,thetaHat,rho+1) + S.B(:,:,theta)*F(:,:,thetaHat,rho+1);    
       
    
    LMI1 = [BLK11R, BLK12R; BLK12R', BLK22]; % LMI-1
    eig_lmi1 = [eig_lmi1, eig(LMI1)];
    
    % -- W
    BLK11W = W(:,:,ell);
    BLK12W = S.C(:,:,theta) * G(:,:,thetaHat,rho+1) + S.D(:,:,theta) * F(:,:,thetaHat,rho+1);
    
    LMI2 = [BLK11W, BLK12W; BLK12W', BLK22]; % LMI-2
    eig_lmi2 = [eig_lmi2, eig(LMI2)];
       
end

lmis_1_are_satisfied = all(all(eig_lmi1 > 0));
lmis_2_are_satisfied = all(all(eig_lmi2 > 0));



lmis_are_satisfied = lmis_1_are_satisfied*lmis_2_are_satisfied;

if ~lmis_are_satisfied
    warning('JUNIOR >> LMIs are not satisfied!');
end

end



