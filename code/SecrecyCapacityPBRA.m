function [SecrecyCapacity,X,SecrecyRateSeq] = SecrecyCapacityPBRA(H,Hb,He,nA,nB,nE,maxIter,P0)
% This function implements Algorithm 3 in the following paper: 
% A. Mukherjee, B. Ottersten and L. -N. Tran, "On the Secrecy Capacity of MIMO Wiretap Channels: 
% Convex Reformulation and Efficient Numerical Methods," in IEEE Transactions on Communications, 
% vol. 69, no. 10, pp. 6865-6878, Oct. 2021, 
SecrecyRateSeq = zeros(maxIter,1);
K = eye(nB+nE);
for iIter = 1:maxIter
    % find X
    H1 = K^(-0.5)*H;
    Delta = H1'*H1-He'*He;
    Delta_bar = sqrtm(Delta);
    [~,X] = SecrecyCapacityConvexReformulation(nA,nE,He,Delta_bar,P0);
    
    SecrecyRateSeq(iIter) = ComputeSecrecyRateMinMaxObj(Hb,He,X,K);
    
    % Finding K using Lemma 3
    mypsi = inv(K+H*X*H');
    
    myphi12 = mypsi(1:nB,nB+1:nB+nE);
    [U,D] = eig(myphi12*myphi12');
    d = real(diag(D));
    K_bar1 = -2*U*diag(1./(1+sqrt(1+4*d)))*U'*myphi12;
    
    K = [eye(nB) K_bar1;K_bar1' eye(nE)];
   
    if (iIter>10) % let the iterative process run for at least 10 iterations then check convergence
        if(abs(SecrecyRateSeq(iIter)-SecrecyRateSeq(iIter-5))<=1e-5)
            break
        end
    end
    
end
SecrecyRateSeq(iIter+1:end)=[];
SecrecyCapacity = SecrecyRateSeq(end);
end

