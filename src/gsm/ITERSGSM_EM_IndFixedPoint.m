function [COVcs,INVCOVcs,COVc,INVCOVc,COVs,INVCOVs,Pcs,lklhd] =  ITERSGSM_EM_IndFixedPoint(Kc, Ks, indc, inds, data, theeps, numiter, collapse, fixediter)
%ITERGSM_EM_INDFIXEDPOINT - learn parameters of flexible normalization
%This function adjust the parameters of a two-component mixture of Gaussian
%Scale Mixtures to model nonuniform normalization operation.
%Learning is based on the Expectation-Maximization algorithm, Where at each
%maximization step a series of fixed point updates that increase the likelihood
%are performed. 
%
% Syntax:  [out_args] = ITERSGSM_EM_IndFixedPoint(in_args) 
% out_args = COVcs, INVCOVcs, COVc, INVCOVc, COVs, INVCOVs, Pcs, lklhd
% in_args = Kc, Ks, indc, inds, data, theeps, maxiter, collapse
% Inputs:
%    Kc - number of center components
%    Ks - number of surround components
%    idxc - indices of the center components
%    idxs - indices of the surround components
%    data - array of size N x (Kc + Ks) containing the data for which we want to
%           maximize the likelihood
%    theeps - (optional. default 1e-10) regularization 
%    numiter - (optional. default 50) number of EM iterations 
%    collapse - (optional. default 0) threshold value for the estimated Prior Pcs to avoid 
%               numerical instabilities. Usual range of threshold is between 0 and 0.05.
%    fixediter - (optional. default 1) number of fixed point iterations per EM update
%
% Outputs:
%    COVcs - covariance matrix of the center-surround dependent component of the model. 
%            Array of size (Kc + Ks) x (Kc + Ks).
%    INVCOVcs - inverse covariance matrix of the center-surround dependent component.
%    COVc - covariance matrix of the center units in the center-surround
%           independent component. Array of size Kc x Kc
%    INVCOVc - inverse covariance of the center units in the center-surround
%              independent component
%    COVs - covariance matrix of the surround units in the center-surround
%           independent component. Array of size Kc x Kc
%    INVCOVs - Description
%    Pcs - values of the estiamted prior probability of center-surroudn
%          dependence at each EM iteration
%    lklhd - likelihood values for the parameters given data at each EM iteration 
%
%
% See also: ITERGSM_EM_INDFACTORIZEDFIXEDPOINT

% Author: Luis Gonzalo Sanchez Giraldo
% June 2019; Last revision: Aug-12-2019

if(~exist('theeps', 'var'))
	theeps=1e-10;  % to avoid Inf and NaN
elseif isempty(theeps)
    theeps = 0;
end

% collapse value either Pcs > (1 - collapse) || Pcs < collapse
% to avoid model components becoming numerically unstable
if(~exist('collapse', 'var'))
    collapse = 0.0;  
end

%%% constants
if(~exist('maxiter', 'var'))   % # EM loops
    numiter = 50; % default 50 CEM iterations
end

if(~exist('fixediter', 'var'))   % # EM loops
    fixediter = 1; % default 50 CEM iterations
end


Kcs = Kc+Ks;   % total number of units comprising center and surround
[N, n_dim] = size(data); % N is number of data samples
assert(n_dim == Kcs, 'Number of dimensions in data is inconsistent with arguments Kc and Ks');

% generate random responsibilities for the component
resp = rand(N, 1);
X = data;
COVcs = data' * bsxfun(@times, data, resp/sum(resp)); % initial guess
COVc = data(:, indc)' * bsxfun(@times, data(:, indc), (1 - resp)/(N-sum(resp))); % initial guess
COVs = data(:, inds)' * bsxfun(@times, data(:, inds), (1 - resp)/(N-sum(resp))); % initial guess
%%% Initial guess of parameters
Pold = 0.5; 
COVcs = (COVcs + theeps * eye(Kcs)) / 2;
COVc = (COVc + theeps * eye(Kc)) / 2;
COVs = (COVs + theeps * eye(Ks)) / 2;
diagCOVs = diag(COVs);

iter = 0;
lklhd = zeros(numiter, 1);
Pcs = zeros(numiter, 1);
while (iter < numiter)
    iter = iter + 1;
    fprintf('EM Iteration %d\n', iter)
    %%%%% E-STEP center-surround
    % Compute posteriors
    Acs = X / COVcs;
    a = sqrt(sum(Acs .* X, 2));
    p_xcs = (a.^(1 - Kcs/2)) .* besselk(Kcs/2 - 1, a);
    %     p_xcs = (a.^(1 - Kcs/2)).*besselk(Kcs/2 - 1,a)/sqrt(det(COVcs));
    
    Ac = X(:, indc)/COVc;
    a = sqrt(sum(Ac .* X(:, indc), 2));
    p_xc = (a.^(1 - Kc/2)) .* besselk(Kc/2 - 1, a);
%     p_xc = (a.^(1 - Kc/2)).*besselk(Kc/2 - 1,a)/sqrt(det(COVc));
    
    Ss = sqrt(diagCOVs);
    As = bsxfun(@rdivide, X(:, inds), Ss(:)');
    a = abs(As);
    p_xs = prod(sqrt(pi/2) * exp(-a), 2);
%     p_xs = prod(sqrt(pi/2)*bsxfun(@rdivide, exp(-a), Ss(:)'), 2);
    p_cs = Pold * (p_xcs ./ (Pold*p_xcs + (1 - Pold) * (p_xc .* p_xs) * sqrt(det(COVcs)/(det(COVc)*det(COVs)))));
    Pold = mean(p_cs);
 
    %%%%% M-STEP center-surround
    %%% partial M-step for Scs
    % Run fixed point iterations
    if Pold > collapse
        S = COVcs;
        n_dim = Kcs;
        for iFx = 1 : fixediter
            try
                A = X / S;
            catch
                keyboard()
            end
            atmp = sqrt(sum(A .* X, 2));
            g_a = ((besselk(n_dim / 2, atmp) ./ besselk(n_dim / 2 - 1, atmp)) ./ atmp) .* (p_cs(:) / sum(p_cs));
            valid_g_a = ~isnan(g_a);
            S = X(valid_g_a, :)' * bsxfun(@times, X(valid_g_a, :), g_a(valid_g_a));
        end
        COVcs = S;
    end
    %%%%% E-STEP center-only
    % Compute posteriors
    Acs = X / COVcs;
    a = sqrt(sum(Acs .* X, 2));
    p_xcs = (a.^(1 - Kcs/2)) .* besselk(Kcs/2 - 1,a);
    
    p_cs = Pold * (p_xcs ./ (Pold*p_xcs + (1 - Pold) * (p_xc .* p_xs) * sqrt(det(COVcs) / (det(COVc) * det(COVs)))));
    Pold = mean(p_cs);

    %%%%% M-STEP center-only
    %%% partial M-step for Sc
    % Run fixed point iterations
    if Pold < (1 - collapse)
        S = COVc;
        n_dim = Kc;
        for iFx = 1:fixediter
            try
                A = X(:, indc) / S;
            catch
                keyboard()
            end
            atmp = sqrt(sum(A .* X(:, indc), 2));
            g_a = ((besselk(n_dim / 2, atmp) ./ besselk(n_dim / 2 - 1, atmp)) ./ atmp) .* ((1 - p_cs(:)) / sum(1 - p_cs));
            valid_g_a = ~isnan(g_a);
            S = X(valid_g_a, indc)' * bsxfun(@times, X(valid_g_a, indc), g_a(valid_g_a));
        end
        COVc = S;
    end
    %%%%% E-STEP surround-only
    % Compute posteriors
    
    Ac = X(:, indc) / COVc;
    a = sqrt(sum(Ac .* X(:, indc), 2));
    p_xc = (a.^(1 - Kc/2)).*besselk(Kc/2 - 1,a);
    
    p_cs = Pold * (p_xcs ./ (Pold * p_xcs + (1 - Pold) * (p_xc .* p_xs) * sqrt(det(COVcs) / (det(COVc) * det(COVs)))));
    Pold = mean(p_cs);
    
    %%%%% M-STEP surround-only
    %%% partial M-step for Ss
    % Run fixed point iterations
    if Pold < (1 - collapse)
        S = diagCOVs;
        for iFx = 1:fixediter
            Ssqrttmp = sqrt(S);
            try
                A = bsxfun(@rdivide, X(:, inds), Ssqrttmp(:)');
            catch
                keyboard()
            end
            g_a = ((1 - p_cs(:))/sum(1-p_cs));
            valid_g_a = ~isnan(g_a);
            Stmp = sum(bsxfun(@times, abs(X(valid_g_a, inds)), g_a(valid_g_a,:)), 1);
            S = Stmp(:).*Ssqrttmp(:);
        end
        diagCOVs = S;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute log-likelihood 
    
    Ss = sqrt(diagCOVs);
    As = bsxfun(@rdivide, X(:, inds), Ss(:)');
    a = abs(As);
    p_xs = prod(sqrt(pi/2)*bsxfun(@rdivide, exp(-a), Ss(:)'), 2);
    
    lklhd(iter) = mean(log(p_xcs * Pold + (p_xc .* p_xs) * (1 - Pold))); % plus a constant
    Pcs(iter) = Pold;
    
    
    
end % EM iterations
%% Just for compatibility purposes
COVs = diag(diagCOVs);
INVCOVcs = inv(COVcs);
INVCOVc = inv(COVc);
INVCOVs = inv(COVs);

