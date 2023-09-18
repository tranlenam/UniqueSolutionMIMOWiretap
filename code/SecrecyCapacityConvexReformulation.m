function [SecrecyCap,optsig] = SecrecyCapacityConvexReformulation(nA,nE,He,Delta_bar,P0)

% This function computes the secrecy capacity for a degraded channel using
% Lemma 1 of the following paper: A. Mukherjee, B. Ottersten and L. -N. Tran,
% "On the Secrecy Capacity of MIMO Wiretap Channels: Convex Reformulation and Efficient Numerical Methods,"
% in IEEE Transactions on Communications, vol. 69, no. 10, pp. 6865-6878, Oct. 2021, 

%% YALMIP implementation
X = sdpvar(nA,nA,'hermitian','complex');
Y = sdpvar(nA,nA,'hermitian','complex');
F = [X>=0,Y>=0,real(trace(X)) <= P0];
F = [F,[eye(nA)+Delta_bar*X*Delta_bar'-Y  Delta_bar*X*He';
    He*X'*Delta_bar'   eye(nE)+He*X*He'] >= 0];
ops = sdpsettings('solver','mosek','verbose',0);
obj = logdet(Y); % this is equivalent to maximizing the determinant of (Y)
diagnotics = optimize(F,-obj,ops);
if diagnotics.problem==0
    optsig = value(X);
    SecrecyCap = real(log(det(value(Y))));
end
%% A corresponding CVX implementation is provided below
%{
cvx_begin quiet
variable X(nA,nA) complex hermitian;
variable Y(nA,nA) complex hermitian;
obj = det_rootn(Y);


maximize (obj)
subject to

real(trace(X)) <= P0;
X == hermitian_semidefinite( nA);
Y == hermitian_semidefinite( nA );
[eye(nA)+Delta_bar*X*Delta_bar'-Y  Delta_bar*X*He';
    He*X'*Delta_bar'   eye(nE)+He*X*He'] == hermitian_semidefinite( nA + nE );
cvx_end
optsig = X;
SecrecyCap = real(log(det(Y)));
%}

end