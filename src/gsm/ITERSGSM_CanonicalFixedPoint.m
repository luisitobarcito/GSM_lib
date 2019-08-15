function [COVcs,INVCOVcs, f_lklhd] =  ITERSGSM_CanonicalFixedPoint(Kc, Ks, data, theeps, numiter)
%ITERGSM_EM_CANONICALFIXEDPOINT - learn parameters of single GSM normalization
%This function adjusts the parameters of a Gaussian Scale Mixtures to model
%spatially uniform normalization operation.
%Learning is based on the Expectation-Maximization algorithm, Where at each
%maximization step a series of fixed point updates that increase the likelihood
%are performed.
%
% Syntax:  [out_args] = ITERSGSM_EM_CanonicalFixedPoint(in_args) 
% out_args = COVcs, INVCOVcs, lklhd
% in_args = Kc, Ks, data, theeps, numiter
% Inputs:
%    Kc - number of center components
%    Ks - number of surround components
%    data - array of size N x (Kc + Ks) containing the data for which we want to
%           maximize the likelihood
%    theeps - (optional. default 1e-10) regularization 
%    numiter - (optional. default 50) number of EM iterations 
%
% Outputs:
%    COVcs - covariance matrix of the center-surround dependent component of the model. 
%            Array of size (Kc + Ks) x (Kc + Ks).
%    INVCOVcs - inverse covariance matrix of the center-surround dependent component.
%    lklhd - likelihood values for the parameters given data at each EM iteration 
%
%
% See also: ITERGSM_EM_CANONICALFACTORIZEDFIXEDPOINT

% Author: Luis Gonzalo Sanchez Giraldo
% June 2019; Last revision: Aug-12-2019


%%% constants
if(~exist('theeps', 'var'))   % reg param for intial guess
    theeps = 1e-10; % default 1^-10
end                                            

if(~exist('numiter', 'var'))   % # EM loops
    numiter = 50; % default 50 CEM iterations
end                                            

Kcs = Kc+Ks;                                                  % # filters
N = size(data,1);                                           % # data samples
% generate random responsibilities for the component
X = data;
X_cs = X; % initial guess
COVcs = (X_cs'*X_cs/N + theeps*eye(Kcs))/2;
iter = 0;
f_lklhd = zeros(numiter, 1);
while (iter<numiter)
    iter = iter + 1;
    fprintf('Iteration %d\n', iter)
    %%%%% M-STEP center-surround
    %%% partial M-step for Scs
    % Run fixed point iterations
    S = COVcs;
    n_dim = Kcs;
    Ssqrttmp = S;
    try
        A = X*pinv(Ssqrttmp);
    catch
        keyboard()
    end
    atmp = sqrt(sum(X.*A, 2));
    r_a = (besselk(n_dim/2, atmp)./besselk(n_dim/2 - 1, atmp));
    r_a(isnan(r_a)) = 1;
    %g_a = ((besselk(n_dim/2, atmp)./besselk(n_dim/2 - 1, atmp))./atmp).*(p_cs(:)/sum(p_cs));
    g_a = (r_a./atmp)/N;
    valid_g_a = ~isnan(g_a);
    S = X(valid_g_a, :)'*bsxfun(@times, X(valid_g_a, :), g_a(valid_g_a));
    COVcs = S;
    % compute log-likelihood 
    Scs = COVcs;
    Acs = X*pinv(Scs);
    a = sqrt(sum(X.*Acs, 2));
    p_xcs = (a.^(1 - Kcs/2)).*besselk(Kcs/2 - 1,a)/sqrt(det(COVcs));
    
    f_lklhd(iter) = mean(log(p_xcs)); % plus a constant
      
end % EM iterations
%% Just for compatibility purposes
INVCOVcs = inv(COVcs);

