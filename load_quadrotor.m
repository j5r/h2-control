function S = load_quadrotor() 
load quadrotor_3modes.mat 

%%%%%%%%%%%
A = A;
B = B;
E = Bw;
C = Cz;
D = Du;
%%%%%%%%%%%%
P = [ .5 .25 .25 0; .25 .5 0 .25 ; .01 0 .99 0 ; .01 0 0 .99 ];
pi00 = [0 .5 .5 0];
A = 1*A;
A(:,:,2) = 1*A(:,:,1);
A(:,:,3) = 1*A(:,:,1);
A(:,:,4) = A(:,:,1);
%
B = B;
B(:,:,2) = .5 * B(:,:,1);
B(:,:,3) = 1 * B(:,:,1);
B(:,:,4) = .1 * B(:,:,1);
%
E(:,:,2) =   E(:,:,1);
E(:,:,3) =   E(:,:,1);
E(:,:,4) =   E(:,:,1);
%
C(:,:,2) =   C(:,:,1);
C(:,:,3) =   C(:,:,1);
C(:,:,4) =   C(:,:,1);
%
D(:,:,2) =   D(:,:,1);
D(:,:,3) =   D(:,:,1);
D(:,:,4) =   D(:,:,1);
%%%%%%%%%%%%
Prob = P;
pi = pi00;
T = 2;
No = 2;
N = numel(pi);
pihat = eye(N,No+1);
% pihat(:,No+1) = pi(:);
% pihat(1:No,No+1) = 0;
% if sum(pihat(N,No+1)) > 0
%    pihat(:,No+1) = pihat(:,No+1)/sum(pihat(:,No+1));
% else
%     pihat(No+1:N,No+1) = 1/(N-No);
% end


S = parse_data(A,B,E,C,D,Prob,pi,pihat,T,No);
end