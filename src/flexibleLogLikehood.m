function log_lik = flexibleLogLikehood(X, param)

COVcs = param.COVcs;
COVc = param.COVc;
COVs = param.COVs;
Pind = param.Pind;
indc = param.idxc;
inds = param.idxs;

Kcs = size(COVcs, 1);
Kc = size(COVc, 1);

diagCOVs = diag(COVs);
Bcs = chol(COVcs, 'upper');
Bc = chol(COVc, 'upper');

%% Compute log-likelihood

Scs = Bcs;
Acs = X*pinv(Scs);
a = sqrt(sum(Acs.*Acs, 2));
p_xcs = (a.^(1 - Kcs/2)).*besselk(Kcs/2 - 1,a)/abs(det(Bcs));

Sc = Bc;
Ac = X(:, indc)*pinv(Sc);
a = sqrt(sum(Ac.*Ac, 2));
p_xc = (a.^(1 - Kc/2)).*besselk(Kc/2 - 1,a)/abs(det(Bc));

Ss = sqrt(diagCOVs);
As = bsxfun(@rdivide, X(:, inds), Ss(:)');
a = abs(As);
p_xs = prod(sqrt(pi/2)*bsxfun(@rdivide, exp(-a), Ss(:)'), 2);

log_lik = log(p_xcs*Pind + p_xc.*p_xs*(1 - Pind)); % plus a constant

end