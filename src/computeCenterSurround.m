function [cs_feat, eff_sz, eff_idx] = computeCenterSurround(map, c_neigh, s_neigh, n_patch, patch_sz)
% COMPUTECENTERSURROUND - compute center surround from small locations in image
%
% Syntax:  [cs_feat, eff_neigh] = computeFullCenterSurround(map, c_neig, s_neigh)
%
% Inputs:
%    map : 5-dimensional array [width, height, n_channels, n_groups, n_samples] 
%    c_neigh : 3-dimensioanl array [n_channels, n_center_neigh, n_groups]
%    s_neigh : 2-dimensional array [n_surround_neigh, 2] [vertical, horizontal]
%              dimensions
%    n_patch : number of patches to sample
%    patch_sz : size of the patches to be sampled from input map [width, height]
%
% Outputs:
%    cs_feat : 5-dimensional array [n_samp, n_channels, n_groups, n]
%    eff_sz : effective spatial size of the center surround features []
%    eff_idx : location of the start of the subsample patch in map
%
% Author: Luis Gonzalo Sanchez Giraldo
% October 2016; last revision: 23-October-2018

if ~exist('patch_sz', 'var')
   patch_sz = ones(1, 2);
end
% compute borders of the effective size
min_h = min(s_neigh(:,1));
max_h = max(s_neigh(:,1));
min_w = min(s_neigh(:,2));
max_w = max(s_neigh(:,2));
cs_sz_w = max_w - min_w + 1 + (patch_sz(1) - 1);
cs_sz_h = max_h - min_h + 1 + (patch_sz(2) - 1);
% cs_center = floor(([cs_sz_r, cs_sz_c] + 1)/2);
cs_center =  1 - [min_w, min_h];
% keep in mind we use [width, height] configuration
[w, h, ch, grp, n] = size(map);

eff_sz = [w -  cs_sz_w, h - cs_sz_h] + 1;
if exist('n_patch', 'var')
    eff_idx = zeros(n_patch, 2, n);
    for iSmp = 1 : n
        rnd_idx = randperm(eff_sz(1)*eff_sz(2));
        eff_idx(:, 1, iSmp) = cs_center(1) + mod(rnd_idx(1 : n_patch) - 1, eff_sz(1))';
        eff_idx(:, 2, iSmp) = cs_center(2) + floor((rnd_idx(1 : n_patch) - 1) / eff_sz(1))';
    end
else
% if no number of patches is given, we consider the entire effective size as
% a single patch. Valid center postions where all surround neighbors are within the 
% size of the image.
    eff_idx(:, 1, :) = repmat(cs_center(1), 1, n);
    eff_idx(:, 2, :) = repmat(cs_center(2), 1, n);
    patch_sz(1) = eff_sz(1);
    patch_sz(2) = eff_sz(2);
    n_patch = 1;
end
n_feat = size(c_neigh, 2) + size(s_neigh, 1);
patch_numel = patch_sz(1) * patch_sz(2);
cs_feat = zeros(size(eff_idx, 1) * patch_numel, ch, grp, n, n_feat);
n_s_neigh = size(s_neigh, 1);
n_c_neigh = size(c_neigh, 2);

% TODO vectorize over channels and groups
for iSmp = 1:n
    for iGrp = 1:grp
        for iCh = 1:ch
            for iCngh = 1 : n_c_neigh
                for iSPtch = 1 : n_patch
                    idx = (iSPtch - 1) * patch_numel + (1 : patch_numel);
                    w_feat_idx = eff_idx(iSPtch, 1, iSmp) + (0:(patch_sz(1)-1));
                    h_feat_idx = eff_idx(iSPtch, 2, iSmp) + (0:(patch_sz(2)-1));
                    cs_feat(idx, iCh, iGrp, iSmp, iCngh) = reshape((map(w_feat_idx, h_feat_idx, c_neigh(iCh, iCngh), iGrp, iSmp)), [], 1);
                end
            end
            for iSngh = 1 : n_s_neigh
                for iSPtch = 1 : n_patch
                    w_feat_idx = eff_idx(iSPtch, 1, iSmp) + s_neigh(iSngh, 1) + (0:(patch_sz(1)-1));
                    h_feat_idx = eff_idx(iSPtch, 2, iSmp) + s_neigh(iSngh, 2) + (0:(patch_sz(2)-1));
                    idx = (iSPtch - 1) * patch_numel + (1 : patch_numel);
                    iFt = n_c_neigh + iSngh;
                    try
                        cs_feat(idx, iCh, iGrp, iSmp, iFt) = reshape((map(w_feat_idx, h_feat_idx, iCh, iGrp, iSmp)), [], 1);
                    catch
                        w_feat_idx
                        h_feat_idx
                        keyboard();
                    end
                end
            end
        end
    end
end

