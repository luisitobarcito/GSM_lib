function [G, Ecs_gain] = canonicalNormalization(X, model, whiten)
% CANONICALNORMALIZATION
% This function computes a weighted normalization on the data matrix X.
% The normalization model is based on a Gaussian scale mixture distribution.
% The model represents center and context dependencies throught a multiplicative
% mixer
%
% Syntax:  [G, Ecs_gain] = canonicalNormalization(X, param)
%
% Inputs:
%    X - data array of size n x u, where n is the number of sample points
%       and u is the total number of units center plus surround context 
%       considered in the normalization
%    param - struct with the follwing fields 
%       COVcs : center-surround covariance
%       idxm : index of the main unit 
%       idxc : indexes of columns from data array corresponding to center units.
%       idxs : indexes of the surround units.
%       whiten : (optional) applies ZCA whitening transformation to Gaussian 
%       variable.
%
% Outputs:
%    G - Normalized values of the main unit X(:, idxm)
%    Ecs_gain - gain for conditional expectation of G given X and cs dependent
%    model
% 
%
% See also: FLEXIBLENORMALIZATION

% Author: Luis Gonzalo Sanchez Giraldo
% Oct 2018; Last revision: Aug-12-2019

%% Parse option arguments
if ~exist('whiten', 'var')
    whiten = false;
end

%% read model parameters
COVcs = model.COVcs; % center-surround covariance for the cs dependent component
idxm = model.idxm; % indexes of the main units to return normalized responses
idxc = model.idxc; % indexes of center units in data matrix
idxs = model.idxs; % indexes of the surround units data matrix
Kc = length(idxc);
Ks = length(idxs);
Kcs = Kc + Ks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization parameter is set internally 
reg_param = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute conditional expectation of G given X and model
Scs = real(sqrtm(COVcs));
Acs = X*pinv(Scs);
acs = sqrt(sum(Acs.*Acs, 2));
% Rcs = (besselk((Kcs - 1)/2, acs)./besselk(Kcs/2 - 1, acs)); % ORIGINAL NO REGULARIZATION
Rcs = (besselk((Kcs - 1)/2, reg_param + acs)./besselk(Kcs/2 - 1, reg_param + acs)); % REGULARIZATION ADDED (WARNING)
Rcs(isnan(Rcs)) = 1;
% Ecs_gain = Rcs./sqrt(reg_param + acs); % ORIGINAL
Ecs_gain = Rcs./sqrt(reg_param + acs); % REGULARIZATION ADDED (WARNING)
if whiten == true
    % compute ZCA matrix to whiten GSM covariance
    [eigvec_cs, eigval_cs] = eig(COVcs);
    W_cs = eigvec_cs*diag(1./sqrt(diag(eigval_cs)))*eigvec_cs';
    X_cs = X*W_cs(:, idxm);
    G = bsxfun(@times, X_cs, Ecs_gain);
else
    G = bsxfun(@times, X(:, idxm), Ecs_gain(:));
end
