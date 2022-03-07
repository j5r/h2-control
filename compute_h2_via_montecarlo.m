
function S = compute_h2_via_montecarlo(S)
if ~isfield(S,'montecarlo')
    error('Please, config_montecarlo first.');
end
%
fprintf('\n --> COMPUTE H2 VIA MONTECARLO...\n');
if isfield(S,'lmi_solution')
    fprintf('\n ****** COMPUTING H2 FOR OUR LMI SOLUTION...\n');
    disp(' ------ Simple system...');
    tic
    S = compute_h2_via_montecarlo_for_lmi_solution_simple_system(S);
    toc
    disp(' ------ Augmented system...');  
    tic
    S = compute_h2_via_montecarlo_for_lmi_solution_augm_system(S);
    toc   
    fprintf('\b\b DONE.\n')
end

if isfield(S,'riccati_solution')
    fprintf('\n ****** COMPUTING H2 FOR RICCATI SOLUTION...\n')
    tic
    S = compute_h2_via_montecarlo_for_riccati_solution(S);
    toc
    fprintf('\b\b DONE.\n')
end
%
if isfield(S,'doval_solution')
    fprintf('\n ****** COMPUTING H2 FOR DO VAL SOLUTION...\n')
    tic
    S = compute_h2_via_montecarlo_for_doval_solution(S);
    toc
    fprintf('\b\b DONE.\n')
end
fprintf(' ...DONE.\n');
end
%
%
%
%
function S = compute_h2_via_montecarlo_for_lmi_solution_augm_system(S)
montecarlo = S.montecarlo;
lmi_solution = S.lmi_solution;
%
H2 = zeros(1,montecarlo.repetitions);
z = zeros(size(S.C,1),1);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    percent = 0;
    for mc = 1:montecarlo.MC
        for xi0 = find(S.augm_pi)
            mkchain = get_markov_chain(S.augm_Prob,montecarlo.horizon,xi0);
            for e = 1:size(S.E,2)
                z = z*0;                
                x = x*0;
                for k = 1:numel(mkchain)
                    xi = mkchain(k);
                    [theta,~,~,~] = map1to4(xi,S);
                    x = lmi_solution.cloopA(:,:,xi)*x + S.E(:,e,theta)*(k==1);
                    z  = lmi_solution.cloopC(:,:,xi)*x;
                    h2 = h2 + z'*z;
                end
            end
            h2 = h2 * S.augm_pi(xi0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/7))==0
            percent = percent + ceil(100/7);
            fprintf('\n\tOur LMI: Monte Carlo simulation %d%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2);
    fprintf('\n');
end
S.lmi_solution.h2.via_montecarlo.augm_sys.h2 = H2;
S.lmi_solution.h2.via_montecarlo.augm_sys.mean = mean(H2);
end
%
%
%
%
function S = compute_h2_via_montecarlo_for_lmi_solution_simple_system(S)
montecarlo = S.montecarlo;
lmi_solution = S.lmi_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    percent = 0;
    for mc = 1:montecarlo.MC
        for theta0 = find(S.pi)
            mkchain = get_markov_chain(S.Prob,montecarlo.horizon,theta0);
            for e = 1:size(S.E,2)              
                x = x*0;
                rho = 0;
                lambda = 1;                
                %%%% FOR K == 1
                theta = mkchain(1);
                theta_is_obs = theta <= S.No;
                %
                rho = 0*(theta_is_obs) + min(S.T-1,rho+1)*(~theta_is_obs);
                lambda = theta*(theta_is_obs) + lambda*(~theta_is_obs);
                %
                if theta_is_obs
                    distrib_thetaHat = S.pihat(:,theta);
                else
                    distrib_thetaHat = S.pihat(:,S.No+1);
                end
                %                
                thetaHat = sum(cumsum(distrib_thetaHat) < rand) + 1;
                %
                index = map4to1(theta,thetaHat,rho,lambda,S);
                %
                x = lmi_solution.cloopA(:,:,index)*x + S.E(:,e,theta);
                z  = lmi_solution.cloopC(:,:,index)*x;
                h2 = h2 + z'*z;
                %
                %%%% FOR K  > 1                
                for k = 2:numel(mkchain)
                    theta = mkchain(k);
                    theta_is_obs = theta <= S.No;
                    %
                    rho = 0*(theta_is_obs) + min(S.T-1,rho+1)*(~theta_is_obs);
                    lambda = theta*(theta_is_obs) + lambda*(~theta_is_obs);
                    %
                    if theta_is_obs
                        distrib_thetaHat = [1:S.N] == theta;
                    else
                        distrib_thetaHat = S.mu(:,rho+1,lambda);
                    end
                    %                    
                    thetaHat = sum(cumsum(distrib_thetaHat) < rand) + 1;
                    %
                    index = map4to1(theta,thetaHat,rho,lambda,S);
                    %
                    x = lmi_solution.cloopA(:,:,index)*x;
                    z  = lmi_solution.cloopC(:,:,index)*x;
                    h2 = h2 + z'*z;
                end
            end
            h2 = h2 * S.pi(theta0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/8))==0
            percent = percent + ceil(100/8);
            fprintf('\n\tOur LMI: Monte Carlo simulation %d%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2);
    fprintf('\n');
end
S.lmi_solution.h2.via_montecarlo.simple_sys.h2 = H2;
S.lmi_solution.h2.via_montecarlo.simple_sys.mean = mean(H2);
end
%
%
%
%
function S = compute_h2_via_montecarlo_for_riccati_solution(S)
montecarlo = S.montecarlo;
riccati_solution = S.riccati_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    percent = 0;
    for mc = 1:montecarlo.MC
        for theta0 = find(S.pi)
            mkchain = get_markov_chain(S.Prob,montecarlo.horizon,theta0);
            for e = 1:size(S.E,2)
                x = x*0;
                %%%% FOR K == 1
                theta = mkchain(1);
                x = riccati_solution.cloopA(:,:,theta)*x + S.E(:,e,theta);
                z  = riccati_solution.cloopC(:,:,theta)*x;
                h2 = h2 + z'*z;
                %%%% FOR K  > 1
                for k = 1:numel(mkchain)
                    theta = mkchain(k);
                    x = riccati_solution.cloopA(:,:,theta)*x;
                    z  = riccati_solution.cloopC(:,:,theta)*x;
                    h2 = h2 + z'*z;
                end
            end
            h2 = h2 * S.pi(theta0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/7))==0
            percent = percent + ceil(100/7);
            fprintf('\n\tRiccati: Monte Carlo simulation %d%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2);
    fprintf('\n');
end
S.riccati_solution.h2.via_montecarlo.h2 = H2;
S.riccati_solution.h2.via_montecarlo.mean = mean(H2);
end
%
%
%
%
function S = compute_h2_via_montecarlo_for_doval_solution(S)
montecarlo = S.montecarlo;
doval_solution = S.doval_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    percent = 0;
    for mc = 1:montecarlo.MC
        for theta0 = find(S.pi)
            mkchain = get_markov_chain(S.Prob,montecarlo.horizon,theta0);
            for e = 1:size(S.E,2)
                x = x*0;
                %%%% FOR K == 1
                theta = mkchain(1);
                thetaHat = theta*(theta<=S.No) + (S.No+1)*(theta>S.No);
                x = doval_solution.cloopA(:,:,theta,thetaHat)*x + S.E(:,e,theta);
                z  = doval_solution.cloopC(:,:,theta,thetaHat)*x;
                h2 = h2 + z'*z;
                %%%% FOR K  > 1
                for k = 2:numel(mkchain)
                    theta = mkchain(k);
                    thetaHat = theta*(theta<=S.No) + (S.No+1)*(theta>S.No);
                    x = doval_solution.cloopA(:,:,theta,thetaHat)*x; 
                    z  = doval_solution.cloopC(:,:,theta,thetaHat)*x;
                    h2 = h2 + z'*z;
                end
            end
            h2 = h2 * S.pi(theta0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/7))==0
            percent = percent + ceil(100/7);
            fprintf('\n\tDo Val: Monte Carlo simulation %d%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2);
    fprintf('\n');
end
S.doval_solution.h2.via_montecarlo.h2 = H2;
S.doval_solution.h2.via_montecarlo.mean = mean(H2);
end