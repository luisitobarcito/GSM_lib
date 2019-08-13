function [cs_feat, eff_neigh] = computeFullCenterSurround(map, c_neigh, s_neigh)
%COMPUTEFULLCENTERSURROUND - compute center surround from a map array
%
% Syntax:  [cs_feat, eff_neigh] = computeFullCenterSurround(map, c_neig, s_neigh)
%
% Inputs:
%    map - 5-dimensional array [width, height, n_channels, n_groups, n_samples]
%    c_neigh - 3-dimensioanl array [n_channels, n_center_neigh, n_groups]
%    s_neigh - 2-dimensional array [n_surround_neigh, 2] [vertical, horizontal] 
%              dimensions
%
% Outputs:
%    cs_feat - 5-dimensional array [width * height, n_channels, n_groups, n_samples, n_feat]
%    eff_neigh - indicator variable with the effective neighbours
%
% Author: Luis Gonzalo Sanchez Giraldo
% Oct 2017; Last revision: 09-Oct-2018

%% 
[w, h, ch, grp, n] = size(map);
ctr = size(c_neigh, 2); % number of center channels (including the one to be normalized)
srd = size(s_neigh, 1); % number of surround channels 
% calulate size of pool (number of features)
feat = ctr + srd;
cs_feat = zeros(w, h, ch, grp, n, feat);
eff_neigh = false(w, h, srd, ch, grp, n); % use logical indexing

for iGp = 1 : grp
    for iCh = 1 : ch
        for iFt = 1 : ctr
            cs_feat(:, :, iCh, iGp, :, iFt) = reshape(map(:, :, c_neigh(iCh, iFt, iGp), iGp, :), [w, h, 1, 1, n, 1]);
        end
        for iFt = (1 + ctr) : feat
            iSr = iFt - ctr;
            h_feat_idx = max(1, 1 - s_neigh(iSr, 1)) : min(h, h - s_neigh(iSr, 1));
            w_feat_idx = max(1, 1 - s_neigh(iSr, 2)) : min(w, w - s_neigh(iSr, 2));
            h_map_idx = max(1, 1 + s_neigh(iSr, 1)) : min(h, h + s_neigh(iSr, 1));
            w_map_idx = max(1, 1 + s_neigh(iSr, 2)) : min(w, w + s_neigh(iSr, 2));
            cs_feat(w_feat_idx, h_feat_idx, iCh, iGp, :, iFt) = reshape(map(w_map_idx, h_map_idx, iCh, iGp, :), [length(w_map_idx), length(h_map_idx), n]);
            eff_neigh(w_feat_idx, h_feat_idx, iSr, iCh, iGp, :) = true;
        end
    end
end
cs_feat = reshape(cs_feat, [], ch, grp, n, feat);
eff_neigh = reshape(eff_neigh, [], srd, ch, grp, n);
end
