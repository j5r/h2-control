%
% By Junior R. Ribeiro, jrodrib@usp.br, 01-mar-2022
%
% Struct = ogura_model_create_phisical_states(Struct)
%
%   This function creates a new field on Struct, called 'valid_states'.
%   These valid_states stores the tuples (theta,thetaHat,rho) which
%   does make sense according to their definitions in 
%   (OGURA, https://doi.org/10.1016/j.automatica.2017.11.022)
%

function Struct = ogura_model_create_phisical_states(Struct)
fprintf('\n --> OGURA MODEL: CREATE PHISICAL STATES...\n');
valid_states = zeros((Struct.N)*Struct.T*Struct.No,3);
%
% validating the compatible states
row = 0;
for theta = 1:Struct.N
    for thetaHat = 1:Struct.No
        for rho = 0:Struct.T-1
            indicator_obs = (theta==thetaHat);
            indicator_unobs = (theta>Struct.No);
            if indicator_obs + indicator_unobs > 0
                row = row + 1;
                valid_states(row,:)= [theta,thetaHat,rho];
            end            
        end
    end
end
valid_states(row+1:end,:) = [];
Struct.valid_states = valid_states;
fprintf('...DONE.\n');
end