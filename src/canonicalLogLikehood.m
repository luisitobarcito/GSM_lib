function log_lik = canonicalLogLikehood(X, param)

COVcs = param.COVcs;
Bcs = chol(COVcs, 'upper');
Kcs = size(COVcs, 1);
%% Compute log-likelihood

% compute log-likelihood
Scs = Bcs;
Acs = X*pinv(Scs);
a = sqrt(sum(Acs.*Acs, 2));
p_xcs = (a.^(1 - Kcs/2)).*besselk(Kcs/2 - 1,a)/abs(det(Bcs));

log_lik = log(p_xcs); % plus a constant

end