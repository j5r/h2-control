function Struct = solve_lmi2_our_size4b(Struct,epsilon)

fprintf('\n. SOLVE_LMI2_OUR_SIZE4(B) ...');

if nargin==1
    epsilon = 1e-5;
end


% [G,K] = solve_riccati(Struct);
% for i = 1:size(G,3)
%     F(:,:,i) = K(:,:,i)*G(:,:,i);
% end

disp('SOLVING - DO - VAL')
SVal = solve_lmi_doval(Struct);
F = SVal.FGRW.F;
G = SVal.FGRW.G;
F(:,:,3) = F(:,:,2);
G(:,:,3) = G(:,:,2);
disp('SOLVING - DO - VAL - OK');

% -- Retrieving data from the Struct
A = Struct.A;
B = Struct.B;
C = Struct.C;
D = Struct.D;
E = Struct.E;
No = Struct.No;
T = Struct.T;
N = size(Struct.Prob,2);

Struct.init_distrib_tHat = Struct.Prob * 0;
for i=1:No
    Struct.init_distrib_tHat(i,i) = 1;
end
if sum(Struct.init_distrib(No+1:end)) > 0 % The sum is not zero
    distrib = Struct.init_distrib(:)';
    distrib(1:No) = 0;
    distrib = distrib/sum(distrib);
    for i=No+1:N
        Struct.init_distrib_tHat(i,:) = distrib;
    end
else
    for i=No+1:N % Arbitrating uniform distribution since sum is zero
        Struct.init_distrib_tHat(i,:) = ones(1,N)/N;
    end
end

% -- Limit distribution of the Markov chain needed to calculate the
% augmented transition probability matrix.
dist = Struct.init_distrib * Struct.Prob^1000;
dist(1:No) = 0;
dist = dist / sum(dist);
if sum(dist) == 0
    error('--> Limit distribution is a vector of zeros...');
end


valid_ix = [1
    14
    15
    17
    18
    23
    24
    26
    27];

% -- Augmented probability distribution matrix.
[augm_Prob,~,~,Struct] = new_observer_matrix_size4(Struct,dist);
Struct.augm_Prob = augm_Prob;

% -- Number of states of the original Markov chain.
orig_n_theta = numel(Struct.init_distrib);
% -- Number of states of the augmented Markov chain.
augm_n_theta = size(Struct.augm_Prob,1);



% -- DEFINING LMIs.
% We need to define all the LMIs which will be used during the
% execution of the algorithm outside loops because of the internal
% representation of them. Using the command "newlmi" inside
% loops cause a different internal representation and might lead
% to numerical issues.
% Another thing that makes difference is declaring all the LMIs in the
% same loop, again, because of their internal representation.
% Thus, we declared each one of them in different loops.
setlmis([])

% -- The LMI variable G must be positive definite. Thus, we define
% this constraint in the form of an LMI.
LMI_G = zeros(1,orig_n_theta*T*No);
for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    LMI_G(xi) = newlmi;
end

% -- LMI R
LMI_R = zeros(1,augm_n_theta);
for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    LMI_R(xi) = newlmi;
end

% -- LMI W
LMI_W = zeros(1,augm_n_theta);
for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    LMI_W(xi) = newlmi;
end

% -- The first LMI regards the stability restriction.
LMI1_LABELS = zeros(1,augm_n_theta);
for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    LMI1_LABELS(xi) = newlmi;
end

% -- Second LMI.
LMI2_LABELS = zeros(1,augm_n_theta);
for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    LMI2_LABELS(xi) = newlmi;
end


% -- CLEAR LMI
counting_variables = 0;

% --STATING VARIABLES
% Stating var R and W
W_trace_indexes = [];
VAR_LABEL_R = zeros(1, augm_n_theta);
VAR_LABEL_W = zeros(1, augm_n_theta);

for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    [theta,thetaHat,rho,lambda] = map1to4(xi, N, T);
    
    VAR_LABEL_R(xi)= lmivar(1,[size(A,1),1]);
    [VAR_LABEL_W(xi),~,iW] = lmivar(1,[size(C,1),1]);
    
    indicator1 = (theta<=No)*(theta==thetaHat)*(rho==0)*(theta==lambda);
    indicator2 = (theta>No)*(thetaHat>No)*(rho>0)*(lambda<=No);
    if indicator1 || indicator2
        diag_w = diag(iW);
        W_trace_indexes = [W_trace_indexes; diag_w(:)];
    end
end

% -- Stating variable "mu".
VAR_LABEL_MU = zeros(1,N);
for thetaHat = No+1:N
    VAR_LABEL_MU(thetaHat) = lmivar(1, [1,1]);
end

% -- LMIs regarding "mu".
LMI3_LABELS = zeros(1,N);
for thetaHat = No+1:N
    LMI3_LABELS(thetaHat) = newlmi;
end

% -- Stating vars F and G
VAR_LABEL_G = zeros(orig_n_theta, T, No);
VAR_LABEL_F = zeros(orig_n_theta, T, No);
for rho = 1:T
    for th=1:orig_n_theta
        for lambda=1:No
            VAR_LABEL_F(th,rho,lambda) = 0;%lmivar(2,[size(B,2),size(A,1)]);
            VAR_LABEL_G(th,rho,lambda) = 0;%lmivar(1,[size(A,1),1]);
        end
    end
end


% -- BUILDING LMIS
% LMIs to guarantee G positive-definite.
% cont = 1;
% for i=1:orig_n_theta
%     for r = 1:T
%         for lambda=1:No
%             lmiterm([-LMI_G(cont), 1, 1, VAR_LABEL_G(i,r,lambda)],1,1);
%             cont = cont + 1;
%         end
%     end
% end

for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    [theta,thetaHat,rho,lambda] = map1to4(xi, N, T);
    
    
    lmiterm([-LMI_R(xi), 1, 1, VAR_LABEL_R(xi)],1,1);
    lmiterm([-LMI_W(xi), 1, 1, VAR_LABEL_W(xi)],1,1);
    
    % AQUI ....
    
    % -- First inequality
    lmiterm([-LMI1_LABELS(xi), 1, 1, VAR_LABEL_R(xi)],1,1); % R_xi
    
    indicator1_mu = (theta==thetaHat)*(theta==lambda)*(theta<=No)*(rho==0);
    indicator2_mu = (theta>No)*(thetaHat>No)*(rho==1);
    
    % The term [ -mu * (E * E') ] is added whenever the augmented
    % distribution [ pi(0)(indicator1_mu + mu*indicator2_mu) ] is
    % not zero.
    if indicator1_mu
        % - pi_{theta}(0) * (E_{theta} * E_{theta}')
        lmiterm([-LMI1_LABELS(xi), 1, 1, 0],...
            -Struct.init_distrib(theta)*E(:,:,theta)*E(:,:,theta)');
    elseif indicator2_mu
        % - pi_{theta}(0) * mu_{thetaHat}(1)* (E_{theta} * E_{theta}')
        lmiterm([-LMI1_LABELS(xi), 1, 1, VAR_LABEL_MU(thetaHat)],...
            -Struct.init_distrib(theta)*E(:,:,theta)*E(:,:,theta)',1);
    end
    
    
    lmiterm([-LMI1_LABELS(xi), 1, 2, VAR_LABEL_G(thetaHat,rho+1,lambda)],...
        A(:,:,theta)*G(:,:,theta));                                   % AG
    
    lmiterm([-LMI1_LABELS(xi), 1, 2, VAR_LABEL_F(thetaHat,rho+1,lambda)],...
        B(:,:,theta)*F(:,:,theta));                                   % BF
    
    lmiterm([-LMI1_LABELS(xi), 2, 2, VAR_LABEL_G(thetaHat,rho+1,lambda)],G(:,:,theta)+G(:,:,theta)'); %He(G)
    
    for xi_ = 1:size(valid_ix,1)
        xi_opD = valid_ix(xi_);
        lmiterm([-LMI1_LABELS(xi), 2, 2, VAR_LABEL_R(xi_opD)],...
            -1,augm_Prob(xi_opD,xi)); % opD(R)
    end
end

for xi_ = 1:size(valid_ix,1)
    xi = valid_ix(xi_);
    [theta,thetaHat,rho,lambda] = map1to4(xi, N, T);
    % --Second inequality
    lmiterm([-LMI2_LABELS(xi), 1, 1, VAR_LABEL_W(xi)],1,1) % W_xi
    
    lmiterm([-LMI2_LABELS(xi), 1, 2, VAR_LABEL_G(thetaHat,rho+1,lambda)],...
        C(:,:,theta)*G(:,:,theta));                                   % CG
    
    lmiterm([-LMI2_LABELS(xi), 1, 2, VAR_LABEL_F(thetaHat,rho+1,lambda)],...
        D(:,:,theta)*F(:,:,theta));                                   % DF
    
    lmiterm([-LMI2_LABELS(xi), 2, 2, VAR_LABEL_G(thetaHat,rho+1,lambda)],G(:,:,theta)+G(:,:,theta)'); %He(G)
    
    for xi_ = 1:size(valid_ix,1)
        xi_opD = valid_ix(xi_); % operator D over R
        lmiterm([-LMI2_LABELS(xi), 2, 2, VAR_LABEL_R(xi_opD)],...
            -1,augm_Prob(xi_opD,xi)); % opD(R)
    end
end

% -- The three LMIs regarding "mu"
for thetaHat = No+1:N
    lmiterm([-LMI3_LABELS(thetaHat),1,1,VAR_LABEL_MU(thetaHat)],1,1); % mu > 0
end

LMI4_LABEL  = newlmi;
LMI5_LABEL  = newlmi;
lmiterm([-LMI4_LABEL ,1,1,0],epsilon-1); % epsilon-1
lmiterm([-LMI5_LABEL ,1,1,0],epsilon+1); % epsilon+1
for thetaHat = No+1:N
    lmiterm([-LMI4_LABEL ,1,1,VAR_LABEL_MU(thetaHat)],1,1); % sum(i) mu(i)
    lmiterm([-LMI5_LABEL ,1,1,VAR_LABEL_MU(thetaHat)],-1,1); % -sum(i) mu(i)
end


% -- ENDING THE SETUP
LMI_SYS = getlmis;
number_of_variables = decnbr(LMI_SYS);
weigths = zeros(1, number_of_variables);
for i = W_trace_indexes % penalties over the trace of W
    weigths(i) = 1;
end

% -- SOLVING
[cost, xopt] = mincx(LMI_SYS,weigths);

% -- RECOVERING G AND F, AND COMPUTING THE GAINS
G = zeros(size(A,1),size(A,2),size(B,3),T,No);
F = zeros(size(B,2),size(B,1),size(B,3),T,No);
K = F * 0;

for thetaHat = 1:orig_n_theta
    for rho = 1:T
        for lambda=1:No
            G(:,:,thetaHat,rho,lambda) = dec2mat(LMI_SYS, xopt, VAR_LABEL_G(thetaHat,rho,lambda));
            F(:,:,thetaHat,rho,lambda) = dec2mat(LMI_SYS, xopt, VAR_LABEL_F(thetaHat,rho,lambda));
            % COMPUTING THE GAINS
            K(:,:,thetaHat,rho,lambda) = F(:,:,thetaHat,rho,lambda) * inv(G(:,:,thetaHat,rho,lambda));
        end
    end
end

% -- RECOVERING W AND R
TRACE_W = 0;
W = zeros(size(C,1),size(C,1),augm_n_theta);
R = zeros(size(A,1),size(A,2),augm_n_theta);
for xi = 1:augm_n_theta
    [theta,thetaHat,rho,lambda] = map1to4(xi, orig_n_theta, T);
    W(:,:,xi) = dec2mat(LMI_SYS, xopt, VAR_LABEL_W(xi));
    R(:,:,xi) = dec2mat(LMI_SYS, xopt, VAR_LABEL_R(xi));
    indicator1 = (theta<=No)*(theta==thetaHat)*(rho==0)*(theta==lambda);
    indicator2 = (theta>No)*(thetaHat>No)*(rho>0)*(lambda<=No);
    if indicator1 || indicator2
        TRACE_W = TRACE_W + trace(W(:,:,xi));
    end
end

% Recovering mu
vector_mu = eye(N,N);
for theta = No+1:N
    for thetaHat = No+1:N
        vector_mu(theta,thetaHat) = dec2mat(LMI_SYS, xopt, VAR_LABEL_MU(thetaHat));
    end
    vector_mu(theta,:) = vector_mu(theta,:)/sum(vector_mu(theta,:));
end
Struct.init_distrib_tHat = vector_mu;


% Re-building the initial distribution
Struct.augm_init_distrib = zeros(1,augm_n_theta);

for xi = 1:augm_n_theta
    [theta,thetaHat,rho,lambda] = map1to4(xi,orig_n_theta,T);
    
    indicator1 = (theta<=No)*(thetaHat==theta)*(rho==0)*(lambda==theta);
    indicator2 = (theta>No)*(thetaHat>No)*(rho==1)*(lambda==1);
    
    Struct.augm_init_distrib(xi) = Struct.init_distrib(theta) * ...
        (indicator1 + indicator2*vector_mu(theta,thetaHat));
end


fprintf(' --> Cost from LMIs    %.4f\n',cost);
fprintf(' --> TRACE_W from LMIs %.4f\n',TRACE_W);

Struct.lmi_trace_W = TRACE_W;
Struct.lmi_cost = cost;
%K(isnan(K)) = 0;
Struct.K = K;
structt.F = F;
structt.G = G;
structt.R = R;
structt.W = W;
structt.mu = vector_mu;
Struct.FGRW = structt;
end
