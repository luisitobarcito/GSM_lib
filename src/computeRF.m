function RF = computeRF(map_pos, map_list)
padding = map_list(end).padding;
k_size = map_list(end).size;
stride = map_list(end).stride;
RF = zeros(0,2);
for i = 1 : size(map_pos,1)
    [RF_X, RF_Y] = meshgrid(1:k_size(1), 1:k_size(2));
    RF_X = RF_X + (map_pos(i, 1) - 1)*stride - padding;
    RF_Y = RF_Y + (map_pos(i, 2) - 1)*stride - padding; 
    RF_tmp = [RF_X(:), RF_Y(:)];
    RF = union(RF, RF_tmp, 'rows');
end
if length(map_list) > 1
    RF = computeRF(RF, map_list(1:end-1));
end
end
