function [c_neigh, s_neigh] = centerSurroundNeighIdx(p_dist, n_center_units, spatial_conf)
x = -p_dist:p_dist:p_dist;
y = x;
s_neigh = [reshape(repmat(y(:), 1, length(y)), [], 1) reshape(repmat(x(:)', length(x), 1), [],1)];
s_neigh(sum(abs(s_neigh), 2) == 0, :) = [];
c_neigh = 0 : n_center_units - 1;

switch(spatial_conf)
    case 'circular'
        s_neigh = fix(bsxfun(@rdivide, s_neigh, sqrt(sum(s_neigh.^2, 2)))*p_dist);
    case 'squared'
        return;
    otherwise
        error('Spatial neighbourhood distribution type %s is not defined', spatial_neigh_dist);
end
