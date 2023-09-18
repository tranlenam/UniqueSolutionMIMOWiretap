function [flag] = IsFeasible(X,P)
% flag =1 if X is feasible
% flag =0 if X is infeasible
if(real(trace(X))==P && min(real(eig(X)))>=0)
    flag = 1;
else
    flag=0;

end

