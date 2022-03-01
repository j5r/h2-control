% test_for_distributions
clear,clc
N = randi([4,6]);
T = 50;randi([3,6]);
No = randi([2,N-2]);
fprintf('N=%d  T=%d  No=%d',N,T,No);
%%%
P = rand(N);
P = P./sum(P,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*****% pre-allocating q_old (the road so far) and q_new (the new formulation)
q_old = zeros(N,T,No);
q_new  = zeros(N,T,No);
for lambda = 1:No
    distribution = (1:N==lambda) * P; % 1 time step
    distribution(1:No) = 0;
    distribution = distribution/sum(distribution);
    for rho = 1
        % starting the first distribution 
        q_old(:,rho,lambda) = distribution;
        q_new(:,rho,lambda) = distribution;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**********% the road so far
for lambda = 1:No %...................lambda
    distribution = q_old(:,1,lambda)';
    for rho = 2:T %...................rho
        distribution = distribution * P; % updating
        distribution(1:No) = 0;
        distribution = distribution/sum(distribution);
        q_old(:,rho,lambda) = distribution;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**********% the new formulation (Eduardo's text)
for i = 1:No %........................lambda
    for k = 2:T %.....................rho
        for m = No+1:N %.................thetaHat
            %
            %
            sum_n = 0;
            for n = No+1:N
                numerator = P(n,m) * q_new(n,k-1,i);
                denominator_sum = 0;
                for r = No+1:N
                    for s = No+1:N
                        denominator_sum = denominator_sum + P(s,r) * q_new(s,k-1,i);
                    end
                end
                sum_n = sum_n +  numerator/denominator_sum;
            end
            %
            %
            q_new(m,k,i) = sum_n;
        end
    end
end
%
%
%
q_old,
q_new,
norm_of_difference = norm(q_old(:) - q_new(:))







