function [COVcs, INVCOVcs, lklhd] =  ITERSGSM_CanonicalFactorizedFixedPoint(Kc, Ks, data, theeps, numiter)
%ITERGSM_EM_CANONICALFACTORIZEDFIXEDPOINT - learn parameters of single GSM normalization
%This function adjusts the parameters of a Gaussian Scale Mixtures to model
%spatially uniform normalization operation.
%Learning is based on the Expectation-Maximization algorithm, Where at each
%maximization step a series of fixed point updates that increase the likelihood
%are performed. The fixed point updates are made on the factorized covariance
%matrix rather than the covariance matrix directly. That is, if COV is the
%covariance matrix then we can express it as the product of a matrix A and its transpose 
%                           COV = A' * A
%The fixed point updates are performed to A instead of COV. However, the function does 
%return the covariance COV.
%
% Syntax:  [out_args] = ITERSGSM_EM_CanonicalFactorizedFixedPoint(in_args) 
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
% See also: ITERGSM_EM_INDFACTORIZEDFIXEDPOINT

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
X = data;
X_cs = X; % initial guess
COVcs = (X_cs'*X_cs/N + theeps*eye(Kcs))/2;
Bcs = chol(COVcs, 'upper');
iter = 0;
lklhd = zeros(numiter, 1);
while (iter<numiter)
    iter = iter + 1;
    fprintf('Iteration %d\n', iter)
    %%%%% center-surround
    S = Bcs;
    n_dim = Kcs;
    Ssqrttmp = S;
    try
      A = X*pinv(Ssqrttmp);
    catch
      keyboard()
    end
    atmp = sqrt(sum(A.*A, 2));
    r_a = (besselk(n_dim/2, atmp)./besselk(n_dim/2 - 1, atmp));
    r_a(isnan(r_a)) = 1;
    g_a = (r_a./atmp)/N;
    valid_g_a = ~isnan(g_a);
    S = A(valid_g_a, :)'*bsxfun(@times, X(valid_g_a, :), g_a(valid_g_a));
    Bcs = S;
    COVcs = S'*S;
    % compute log-likelihood 
    Scs = Bcs;
    Acs = X*pinv(Scs);
    a = sqrt(sum(Acs.*Acs, 2));
    p_xcs = (a.^(1 - Kcs/2)).*besselk(Kcs/2 - 1,a)/abs(det(Bcs));
    lklhd(iter) = mean(log(p_xcs)); % plus a constant
      
end % EM iterations
%% Just for compatibility purposes
INVCOVcs = inv(COVcs);

