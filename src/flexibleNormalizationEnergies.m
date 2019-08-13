function [G, p_cs, Ecs, Ec] = flexibleNormalizationEnergies(X, model, whiten, matchvar)
%FLEXIBLENORMALIZATIONENERGIES
%This function computes flexible normalization on the data matrix X.
% The model is a mixture of two Gaussian scale mixture distributions.
% One component represents center and context dependencies throught 
% multiplicative mixer. The second component represent center and context
% independence. The model allows context variables to be independent from each
% other as well. In this case each variable has its own scaling factor.
%
% Syntax:  [G, p_cs, p_gain, Ecs_gain, Ec_gain] = flexibleNormalization(X, model)
%
% Inputs:
%    X - data array of size n x u, where n is the number of sample points
%       and u is the total number of units center plus surround context 
%       considered in the normalization
%    model - struct with the follwing fields 
%       COVcs : center-surround covariance
%       COVc : center independent covariance
%       COVs : surround independent covariance
%       Pcs : probability of center surround dependent
%       idxm : index of the main unit 
%       idxc : indexes of columns from data array corresponding to center units.
%       idxs : indexes of the surround units.
%    matchvar - (optional) reescales the center surround dependent and center
%       surround independent components to have tha same variance in the main
%       unit.
%    whiten - (optional) rotates the center surround independent components
%       to match the center surround dependent component covariance and then
%       applies the same whitening transformation.
%
% Outputs:
%    G - Normalized values of the main unit X(:, idxm)
%    p_cs - posterior probability for center surround dependent component given
%    the observed center and surround values.
%    Ecs - energy of the center-surround pool
%    Ec - energy of the center-only pool
% 
%
% See also: FLEXIBLENORMALIZATION

% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: Aug-12-2019

%% Parse option arguments
if ~exist('whiten', 'var')
    whiten = false;
end

if ~exist('matchvar', 'var')
    matchvar = false;
end

if ~isfield(model, 's_ind')
    s_ind = false;
else
    s_ind = model.s_ind;
end

if ~isfield(model, 'collapse')
    collapse = 0;
else
    collapse = model.collapse;
end



%% read model parameters 
COVcs = model.COVcs; % center-surround covariance for the cs dependent component
COVc = model.COVc; % center covariance for the cs independent component
COVs = model.COVs; % surround covariace for the cs indepdendent component
Pcs = model.Pcs; % Prior probability of cs depenedent
idxm = model.idxm; % indexes of the main units to return normalized responses
idxc = model.idxc; % indexes of center units in data matrix
idxs = model.idxs; % indexes of the surround units data matrix



% Scales the center-only covarinace to match the center-surround variance of the
% unit to be normalized
if matchvar
    COVc = COVc*COVcs(idxm, idxm)/COVc(idxm, idxm);
end

Kc = length(idxc);
Ks = length(idxs);
Kcs = Kc + Ks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization parameter is set internally 
reg_param = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We work with factorized version of the covariance (reparameterization)
Scs = real(sqrtm(COVcs));
Sc = real(sqrtm(COVc));
if s_ind
    Ss = sqrt(diag(COVs));
else
    Ss = real(sqrtm(COVs));
end

% Compute conditional likelihood
if Pcs > collapse 
    Acs = real(X / Scs);
    acs = sqrt(sum(Acs.*Acs, 2));
    p_xcs = (acs.^(1 - Kcs/2)).*besselk(Kcs/2 - 1, reg_param + acs)/sqrt(det(COVcs));
else
    p_xcs = zeros(size(X, 1),1);
end
    
if Pcs < 1 - collapse
    Ac = real(X(:, idxc) / Sc);
    ac = sqrt(sum(Ac.*Ac, 2));
    p_xc = (ac.^(1 - Kc/2)).*besselk(Kc/2 - 1, reg_param + ac)/sqrt(det(COVc));

    if s_ind
        As = bsxfun(@rdivide, X(:, idxs), Ss(:)');
        as = abs(As);
        p_xs = prod(sqrt(pi/2)*bsxfun(@rdivide, exp(-as), Ss(:)'), 2);
    else
        As = X(:, idxs) / Ss;
        as = sqrt(sum(As.*As, 2));
        p_xs = (as.^(1 - Ks/2)).*besselk(Ks/2 - 1,as)/sqrt(det(COVs));
    end
else
    p_xc = zeros(size(X, 1),1);
    p_xs = zeros(size(X, 1),1);
end

% Compute posterior probabilities
p_cs = Pcs * p_xcs ./ (Pcs*p_xcs + (1 - Pcs) * p_xc .* p_xs);
p_cs(isnan(p_cs)) = 1;

if Pcs > collapse 
    Rcs = (besselk((Kcs - 1) / 2, reg_param + acs) ./ besselk(Kcs / 2 - 1, reg_param + acs));
    Rcs(isnan(Rcs)) = 1;
%     Ecs_gain = Rcs ./ sqrt(acs); %% ORIGINAL NO REGULARIZATION
    Ecs_gain = Rcs ./ (sqrt(reg_param + acs));%% SMALL REGULARIZATION TWEAK (WARNING)
else
    Ecs_gain = zeros(size(X, 1),1);
end

if Pcs < 1 - collapse
    Rc = (besselk((Kc - 1) / 2, reg_param + ac) ./ besselk(Kc / 2 - 1, reg_param + ac));
    Rc(isnan(Rc)) = 1;
%     Ec_gain = Rc ./ (sqrt(ac)); %% ORIGINAL NO REGULARIZATION
    Ec_gain = Rc ./ (sqrt(reg_param + ac));%% SMALL REGULARIZATION TWEAK (WARNING)
else
    Ec_gain = zeros(size(X, 1),1);
end


if whiten == true
    % compute ZCA matrices for each component of the mixture of GSM
    % Here the idea is to relate the two covariances by a transformation
    % and then whiten both using the same matrix
    if Pcs > collapse 
        [eigvec_cs, eigval_cs] = eig(COVcs);
        W_cs = eigvec_cs*diag(1./sqrt(diag(eigval_cs)))*eigvec_cs';
    else
        W_cs = eye(Kcs);
    end
    X_cs = X*W_cs(:, idxm);
    Sc_s = zeros(size(COVcs));
    Sc_s(idxc, idxc) = Sc;
    Sc_s(idxs, idxs) = sqrtm(COVs);
    if Pcs < 1 - collapse &&  Pcs > collapse
        W_c_s = pinv(Sc_s)*Scs*W_cs;
    elseif Pcs < collapse
        COVc_s = Sc_s*Sc_s;
        [eigvec_c_s, eigval_c_s] = eig(COVc_s);
        W_c_s = eigvec_c_s*diag(1./sqrt(diag(eigval_c_s)))*eigvec_c_s';
    elseif Pcs >= (1 - collapse)
        W_c_s = W_cs;
    end
    X_c = X*W_c_s(:, idxm);
    G_cs = bsxfun(@times, X_cs, p_cs.*Ecs_gain);
    G_c = bsxfun(@times, X_c, (1 - p_cs).*Ec_gain);
    G = G_cs + G_c;
    p_gain = zeros(size(G));
else
    p_gain = p_cs.*Ecs_gain + (1 - p_cs).*Ec_gain;
    G = bsxfun(@times, X(:, idxm), p_gain(:));
end
Ecs = acs;
Ec = ac;
end



