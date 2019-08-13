function [neigh] = centerNeighIdx(n_channel, n_groups, n_neighs)
%CENTERNEIGHIDX - indexes for center neighborhood
%Compute idexes for cenetr neighbors in the normalization model.
%
% Syntax:  [s_neigh] = surroundNeighIdx(dist, conf)
%
% Inputs:
%    dist - distance from center 
%    conf - spatial configuration of neighbors. Two configurations are available
%    'square' and 'circular'
%
% Outputs:
%    neigh - indexes of center neighborhood for each channel. Note each
%    neighborhood includes the channel
%
% See also: centerNeighIdx

% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 08-November-2018

%% Compute channel neighbors (wrap around the channels in the group)
all_channels = repmat((1:n_channel)', 1, n_neighs + 1);
all_rel_neigh = -floor(n_neighs/2) : ceil(n_neighs/2);
n_channel_per_group = n_channel / n_groups;
neigh = mod(bsxfun(@plus, all_channels, all_rel_neigh) - 1, n_channel_per_group) + 1;
neigh = permute(reshape(neigh, [n_channel_per_group, n_groups, n_neighs + 1]), [1 3 2]);
end