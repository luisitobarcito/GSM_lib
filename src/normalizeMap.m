function [varargout] = normalizeMap(map, norm_func, params)
%NORMALIZEMAP - apply normalization model to a map blob from convolutional layer
%Apply normalization model specified in params and the normalziation function
%to the "map" array. The fucntion can output more than one argument since we may
%want to retrieve components of the normalization as well. 
%
% Syntax:  [norm_map, varargout] = normalizeMap(map, params, norm_func)
%
% Inputs:
%    map - 5-dimensional array from a convolution map the dims are 
%          [W, H, C, G, N] W stands for width, H for height, C for channels per
%          group, G for groups, and N for number of samples.
%    params - Struct with the following fields
%       norm_func - handle to normalization function
%       data_mean - cell array {C x G}
%       factor - cell array {C x G}
%       whiten - logical
%       model - cell array {C x G} each cell contains a struct with the
%       parameters of the data model
%
%
% Outputs:
%    norm_map - array. Same size as map with the normalized values 
%
% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 10-Oct-2018
 
if ~isfield(params, 'whiten')
    params.whiten = false;
end
if ~isfield(params, 'padding')
    params.padding = 'full';
end



assert(ndims(map) == 5, 'Map must be W x H x C x G x N');
% map array dimensions
% h: height
% w: width
% ch: number of channesl per group
% grp: number of groups
% n: number of samples to normalize
[w, h, ch, grp, n] = size(map);

switch(params.padding)
    case 'full'
    % unfold map to perform apply normalization model
        unfld_map = computeFullCenterSurround(map, params.c_neigh, params.s_neigh);
    case 'valid'
        [unfld_map, eff_w_h] = computeValidCenterSurround(map, params.c_neigh, params.s_neigh);
        % change w and h tot the valid sizes
        w = eff_w_h(1);
        h = eff_w_h(2);
    otherwise
        error('Only full and valid options for padding');
end
% Sometimes we would like to retrieve the components that made the inferred
% normalization (for example in flexible normalization)
% preallocate to store the maps
nout = max(1, nargout);
[G_map{1:nout}] = deal(zeros(w, h, ch, grp, n));

% Loop through all channels
% each channel has its own normalization model
for  iGrp = 1 : grp
    for iCh = 1 : ch
        X = unfld_map(:,iCh, iGrp, :, :);
        X = reshape(X, [], size(unfld_map, 5));
        factor = params.factor{iCh, iGrp};
        X_norm = bsxfun(@minus, X, params.data_mean{iCh, iGrp});
        X_norm = bsxfun(@times, X_norm, factor);   

        [G{1:nout}]= norm_func(X_norm, params.model{iCh, iGrp});
        for iLst = 1:numel(G)    
            G_map{iLst}(:, :, iCh, iGrp, :)  = reshape(G{iLst}, [w, h, 1, 1, n]);
        end
    end
end
varargout = G_map;
end
