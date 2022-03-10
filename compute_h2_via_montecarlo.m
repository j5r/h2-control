%
% By Junior R. Ribeiro, jrodrib@usp.br, 10-mar-2022
%
% S = compute_h2_via_montecarlo(S)
%
%   This function performs Monte Carlo simulations for the lmi_solution,
%   doval_solution and riccati_solution in order to compute the H2 cost.
%

function S = compute_h2_via_montecarlo(S)
if ~isfield(S,'montecarlo')
    error('Please, config_montecarlo(...) first.');
end
%
fprintf('\n --> COMPUTE H2 VIA MONTECARLO...\n');
try
    % starting a paralell pool
    parpool;
catch
end
%
parfor item = 1:4
    if isfield(S,'lmi_solution') && item == 1
        fprintf('\n ****** COMPUTING H2 FOR OUR LMI SOLUTION: Simple system...\n');
        tic
        Struct_{item} = compute_h2_via_montecarlo_for_lmi_solution_simple_system(S);
        fprintf('\n DONE ---- Our Simple System: elapsed time %.2f s.\n',toc);
    end
    %
    if isfield(S,'lmi_solution') && item == 2
        fprintf('\n ****** COMPUTING H2 FOR OUR LMI SOLUTION: Augmented system...\n');
        tic
        Struct_{item} = compute_h2_via_montecarlo_for_lmi_solution_augm_system(S);
        fprintf('\n DONE ---- Our Augmented System: elapsed time %.2f s.\n',toc);
    end
    %
    if isfield(S,'riccati_solution') && item == 3
        fprintf('\n ****** COMPUTING H2 FOR RICCATI SOLUTION...\n')
        tic
        Struct_{item} = compute_h2_via_montecarlo_for_riccati_solution(S);
        fprintf('\n DONE ---- Riccati solution: elapsed time %.2f s.\n',toc);
    end
    %
    if isfield(S,'doval_solution') && item == 4
        fprintf('\n ****** COMPUTING H2 FOR DO VAL SOLUTION...\n')
        tic
        Struct_{item} = compute_h2_via_montecarlo_for_doval_solution(S);
        fprintf('\n DONE ---- Do Val solution: elapsed time %.2f s.\n',toc);
    end
end
%%%%% RETRIEVING DATA
%
S.lmi_solution.h2.via_montecarlo.simple_sys = Struct_{1}.via_montecarlo.simple_sys;
S.lmi_solution.h2.via_montecarlo.augm_sys = Struct_{2}.via_montecarlo.augm_sys;
S.riccati_solution.h2.via_montecarlo = Struct_{3}.via_montecarlo;
S.doval_solution.h2.via_montecarlo = Struct_{4}.via_montecarlo;
%
fprintf(' ...DONE.\n');
end
%
%
%
%
function Struct = compute_h2_via_montecarlo_for_lmi_solution_augm_system(S)
montecarlo = S.montecarlo;
lmi_solution = S.lmi_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    h2_sum = 0;
    percent = 0;
    for mc = 1:montecarlo.MC
        for xi0 = find(S.augm_pi)
            mkchain = get_markov_chain(S.augm_Prob,montecarlo.horizon,xi0);
            for e = 1:size(S.E,2)
                x = x*0;
                %%%% FOR K == 1
                xi = mkchain(1);
                [theta,~,~,~] = map1to4(xi,S);
                x = lmi_solution.cloopA(:,:,xi)*x + S.E(:,e,theta);
                z  = lmi_solution.cloopC(:,:,xi)*x;
                h2 = h2 + z'*z;
                %%%% FOR K  > 1
                for k = 2:numel(mkchain)
                    xi = mkchain(k);
                    x = lmi_solution.cloopA(:,:,xi)*x;
                    z  = lmi_solution.cloopC(:,:,xi)*x;
                    h2 = h2 + z'*z;
                end
            end
            h2_sum = h2_sum + h2 * S.augm_pi(xi0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/3))==0
            percent = percent + 100/3;
            fprintf('\n\tOur LMI: augm_sys: Monte Carlo simulation %.2f%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2_sum);
    fprintf('\n');
end
%
S.lmi_solution.h2.via_montecarlo.augm_sys.h2 = H2;
S.lmi_solution.h2.via_montecarlo.augm_sys.mean = mean(H2);
Struct = S.lmi_solution.h2;
end
%
%
%
%
function Struct = compute_h2_via_montecarlo_for_lmi_solution_simple_system(S)
montecarlo = S.montecarlo;
lmi_solution = S.lmi_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    h2_sum = 0;
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
            h2_sum = h2_sum + h2 * S.pi(theta0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/3))==0
            percent = percent + 100/3;
            fprintf('\n\tOur LMI: simple_sys: Monte Carlo simulation %.2f%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2_sum);
    fprintf('\n');
end
%
S.lmi_solution.h2.via_montecarlo.simple_sys.h2 = H2;
S.lmi_solution.h2.via_montecarlo.simple_sys.mean = mean(H2);
Struct = S.lmi_solution.h2;
end
%
%
%
%
function Struct = compute_h2_via_montecarlo_for_riccati_solution(S)
montecarlo = S.montecarlo;
riccati_solution = S.riccati_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    h2_sum = 0;
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
            h2_sum = h2_sum + h2 * S.pi(theta0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/3))==0
            percent = percent + 100/3;
            fprintf('\n\tRiccati: Monte Carlo simulation %.2f%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2_sum);
    fprintf('\n');
end
%
S.riccati_solution.h2.via_montecarlo.h2 = H2;
S.riccati_solution.h2.via_montecarlo.mean = mean(H2);
Struct = S.riccati_solution.h2;
end
%
%
%
%
function Struct = compute_h2_via_montecarlo_for_doval_solution(S)
montecarlo = S.montecarlo;
doval_solution = S.doval_solution;
%
H2 = zeros(1,montecarlo.repetitions);
x = zeros(size(S.A,1),1);
for rep = 1:montecarlo.repetitions
    h2 = 0;
    h2_sum = 0;
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
            h2_sum = h2_sum + h2 * S.pi(theta0)/montecarlo.MC;
        end
        %
        if mod(mc,ceil(montecarlo.MC/3))==0
            percent = percent + 100/3;
            fprintf('\n\tDo Val: Monte Carlo simulation %.2f%% done...',percent);
        end
    end
    H2(rep) = sqrt(h2_sum);
    fprintf('\n');
end
%
S.doval_solution.h2.via_montecarlo.h2 = H2;
S.doval_solution.h2.via_montecarlo.mean = mean(H2);
Struct = S.doval_solution.h2;
end