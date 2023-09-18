function [S] = projectSPC(Y,P0)
% Projection onto Y onto the set {X:PSD,trace(X)<=0}
Nt = size(Y,1);
[U, Sig]=eig(Y); % compute EVD of Y
Sig_y = real(diag(Sig)); % get the diagonal elements
mu_max = max(Sig_y);
mu_min = 0;
errortol = 1e-6;
if(mu_max>0) % we carry out the bisection to find optimal \mu
    while((mu_max-mu_min)>errortol) % the bisection starts
        mu_mid = (mu_max+mu_min)/2; % a given \mu
        Sig_X=max((Sig_y-mu_mid),0); % this follows the derivation
        if sum(Sig_X)<P0
            mu_max=mu_mid;
        else
            mu_min=mu_mid;
        end
    end % finish the bisection
    % after the bisection method, we calculate the projection using
    S=U*diag(Sig_X)*U';
else % in this case S is zero matrix
    S=zeros(Nt);
end

end

