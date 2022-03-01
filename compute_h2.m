%
% By Junior R. Ribeiro, jrodrib@usp.br, 01-mar-2022
%
% Struct = compute_h2(Struct)
%
%   This function computes the H2 cost by controllability and observability
%   gramian, for the lmi_solution and riccati_solution.
%

function S = compute_h2(S)
fprintf('\n --> COMPUTE H2...\n')
validate_compute_h2(S);
%
S = compute_h2_for_lmi_solution(S);
%
S = compute_h2_for_ricatti_solution(S);
fprintf('...DONE.\n');
end
%
%
%
%
function validate_compute_h2(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
assert(isfield(Struct,'lmi_solution') || isfield(Struct,'riccati_solution'),...
    'The Structure does not have the field lmi_solution nor riccati_solution.');
%
if isfield(Struct,'lmi_solution')
    lmi_solution = Struct.lmi_solution;
    assert(isfield(lmi_solution,'ctrl_gramian') || isfield(lmi_solution,'obsv_gramian'), ...
        sprintf('There is neither [*.lmi_solution.obsv_gramian] nor \n[*.lmi_solution.ctrl_gramian] field in the struct.\nPlease, compute at least one of them.'));
end
%
if isfield(Struct,'riccati_solution')
    riccati_solution = Struct.riccati_solution;
    assert(isfield(riccati_solution,'ctrl_gramian') || isfield(riccati_solution,'obsv_gramian'), ...
        sprintf('There is neither [*.riccati_solution.obsv_gramian] nor \n[*.riccati_solution.ctrl_gramian] field in the struct.\nPlease, compute at least one of them.'));
end
%
end
%
%
%
%
function S = compute_h2_for_lmi_solution(S)
%
n_states = size(S.valid_states,1);
lmi_solution = S.lmi_solution;
if isfield(lmi_solution,'obsv_gramian')
    %
    h2 = 0;
    for i = 1:n_states
        [theta,~,~,~] = map1to4(i,S);
        %
        OPER_E_i = zeros(size(S.A,1), size(S.A,1));
        for j = 1:n_states
            OPER_E_i = OPER_E_i + S.augm_Prob(i,j) * S.lmi_solution.ctrl_gramian(:,:,j);
        end
        %
        h2 = h2 + S.augm_pi(i) * trace(S.E(:,:,theta)' * OPER_E_i * S.E(:,:,theta));
    end
    %
    S.lmi_solution.h2.via_obsv_gramian = h2;
end
%
if isfield(lmi_solution,'ctrl_gramian')
    h2 = 0;
    for j = 1:n_states
        [theta,thetaHat,rho,lambda] = map1to4(j,S);
        %
        factor_1 = blkdiag(  S.C(:,:,theta)' * S.C(:,:,theta), ...
            S.D(:,:,theta)' * S.D(:,:,theta)  );
        %
        gain = S.lmi_solution.K(:,:,thetaHat, rho+1, lambda);
        gramian = S.lmi_solution.ctrl_gramian(:,:,j);
        factor_2 = [gramian,            gramian * gain';...
            gain * gramian,     gain * gramian * gain'];
        %
        
        h2 = h2 + trace(factor_1 * factor_2);
    end
    S.lmi_solution.h2.via_ctrl_gramian = h2;
end
%
end
%
%
%
%
function S = compute_h2_for_ricatti_solution(S)
%
n_states = S.N;
riccati_solution = S.riccati_solution;
if isfield(riccati_solution,'obsv_gramian')
    %
    h2 = 0;
    for i = 1:n_states
        [theta,~,~,~] = map1to4(i,S);
        %
        OPER_E_i = zeros(size(S.A,1), size(S.A,1));
        for j = 1:n_states
            OPER_E_i = OPER_E_i + S.Prob(i,j) * S.riccati_solution.ctrl_gramian(:,:,j);
        end
        %
        h2 = h2 + S.augm_pi(i) * trace(S.E(:,:,theta)' * OPER_E_i * S.E(:,:,theta));
    end
    %
    S.riccati_solution.h2.via_obsv_gramian = h2;
end
%
if isfield(riccati_solution,'ctrl_gramian')
    h2 = 0;
    for j = 1:n_states
        [theta,thetaHat,rho,lambda] = map1to4(j,S);
        %
        factor_1 = blkdiag(  S.C(:,:,theta)' * S.C(:,:,theta), ...
            S.D(:,:,theta)' * S.D(:,:,theta)  );
        %
        gain = S.lmi_solution.K(:,:,thetaHat, rho+1, lambda);
        gramian = S.riccati_solution.ctrl_gramian(:,:,j);
        %
        factor_2 = [gramian,          gramian * gain' ;...
            gain * gramian ,   gain * gramian * gain'];
        %
        h2 = h2 + trace(factor_1 * factor_2);
    end
    S.riccati_solution.h2.via_ctrl_gramian = h2;
end
%
end