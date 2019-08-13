function [cs_feat, eff_dim] = computeValidCenterSurround(map, c_neigh, s_neigh)
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
%
% Author: Luis Gonzalo Sanchez Giraldo
% Oct 2017; Last revision: 08-Novemeber-2018

%% 
[w, h, ch, grp, n] = size(map);
ctr = size(c_neigh, 2); % number of center channels (including the one to be normalized)
srd = size(s_neigh, 1); % number of surround channels 
% calulate size of pool (number of features)
feat = ctr + srd;
eff_w = w + min(s_neigh(:,2)) - max(s_neigh(:,2));
eff_h = h + min(s_neigh(:,1)) - max(s_neigh(:,1));
eff_dim = [eff_w, eff_h];
eff_w_str = max(1, 1 - min(s_neigh(:,2)));
eff_w_end = min(w, w - max(s_neigh(:,2)));
eff_h_str = max(1, 1 - min(s_neigh(:,1)));
eff_h_end = min(h, h - max(s_neigh(:,1)));
cs_feat = zeros(eff_w, eff_h, ch, grp, n, feat);
for iGp = 1 : grp
    for iCh = 1 : ch
        for iFt = 1 : ctr
            cs_feat(:, :, iCh, iGp, :, iFt) = reshape(map(eff_w_str:eff_w_end, eff_h_str:eff_h_end, c_neigh(iCh, iFt, iGp), iGp, :), [eff_w, eff_h, 1, 1, n, 1]);
        end
        for iFt = (1 + ctr) : feat
            iSr = iFt - ctr;
            h_map_idx = max(1, eff_h_str + s_neigh(iSr, 1)) : min(h, eff_h_end + s_neigh(iSr, 1));
            w_map_idx = max(1, eff_w_str + s_neigh(iSr, 2)) : min(w, eff_w_end + s_neigh(iSr, 2));
            cs_feat(:, :, iCh, iGp, :, iFt) = reshape(map(w_map_idx, h_map_idx, iCh, iGp, :), [length(w_map_idx), length(h_map_idx), n]);
        end
    end
end
cs_feat = reshape(cs_feat, [], ch, grp, n, feat);
end
