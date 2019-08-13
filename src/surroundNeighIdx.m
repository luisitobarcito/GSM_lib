function [neigh] = surroundNeighIdx(dist, conf)
%SURROUNDNEIGHIDX - Relative indexes for surround neighborhood
%Compute relative positions for spatial neighbors in the normalization model.
%At the time the neighborhood size is fixed to 8 positions:
%
% Syntax:  [s_neigh] = surroundNeighIdx(dist, conf)
%
% Inputs:
%    dist - distance from center 
%    conf - spatial configuration of neighbors. Two configurations are available
%    'square' and 'circular'
%
% Outputs:
%    neigh - relative spatial location of surround neighbors rounded to the
%    closest integer. This is a 2d array where first dims are row positions and
%    second dimension column postions relative to the center postion.
%
% See also: centerNeighIdx

% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 04-October-2018

%% Compute all pairs of relative positions
x = -dist : dist : dist;
y = x;
neigh = [reshape(repmat(y(:), 1, length(y)), [], 1) reshape(repmat(x(:)', length(x), 1), [],1)];
neigh(sum(abs(neigh), 2) == 0, :) = []; % remove center location from list
%% Correct postions if circular neighborhood 
switch(conf)
    case 'circular'
        neigh = round(bsxfun(@rdivide, neigh, sqrt(sum(neigh.^2, 2)))*dist);
    case 'square'
        return;
    otherwise
        error('Spatial neighbourhood distribution type %s is not defined', conf);
end
end