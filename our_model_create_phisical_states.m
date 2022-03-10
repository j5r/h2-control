%
% By Junior R. Ribeiro, jrodrib@usp.br, 01-mar-2022
%
% Struct = our_model_create_phisical_states(Struct)
%
%   This function creates a new field on Struct, called 'valid_states'.
%   These valid_states stores the tuples (theta,thetaHat,rho,lambda) which
%   does make sense according to their definitions.
%

function Struct = our_model_create_phisical_states(Struct)
fprintf('\n --> OUR MODEL: CREATE PHISICAL STATES...\n');
valid_states = zeros((Struct.N^2)*Struct.T*Struct.No,4);
%
% validating the compatible states
row = 0;
for theta = 1:Struct.N
    for thetaHat = 1:Struct.N
        for rho = 0:Struct.T-1
            for lambda = 1:Struct.No
                indicator_obs = (theta==thetaHat)*(theta==lambda)*...
                    (theta<=Struct.No)*(rho==0);
                indicator_unobs = (theta>Struct.No)*(thetaHat>Struct.No)*(rho>0);
                if indicator_obs + indicator_unobs > 0
                    row = row + 1;
                    valid_states(row,:)= [theta,thetaHat,rho,lambda];
                end
            end
        end
    end
end
valid_states(row+1:end,:) = [];
Struct.valid_states = valid_states;
fprintf('...DONE.\n');
end