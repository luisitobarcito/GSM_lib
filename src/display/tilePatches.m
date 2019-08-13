function [tiled_patches] = tilePatches(patches, n_cols, n_rows)

patch_shape_dims = ndims(patches) - 1;
n_patches = size(patches, patch_shape_dims + 1);
switch( patch_shape_dims)
    case 3
        patch_depth = size(patches, patch_shape_dims);
    case 2
        patch_depth = 1;
        patches = reshape(patches, [size(patches, 1), size(patches, 2), 1, n_patches]);
    otherwise 
        error('patches have too many dimensions')
end

if ~exist('n_rows', 'var')
    n_rows = ceil(n_patches/n_cols);
end
patch_h = size(patches, 1);
patch_w = size(patches, 2);
grid_spacing = ceil([patch_h, patch_w]*0.15);
tiled_h = n_rows*patch_h + (n_rows - 1)*grid_spacing(1);  
tiled_w = n_cols*patch_w + (n_cols - 1)*grid_spacing(2);

tiled_patches = 0.5*(max(patches(:))-min(patches(:)))*ones(tiled_h, tiled_w, patch_depth);
iPtch = 1;
for iRow = 1 : n_rows
    r_idx = (iRow - 1)*(patch_h + grid_spacing(1)) + 1;
    for iCol = 1 : n_cols
        if(iPtch > n_patches)
            break;
        end
        c_idx = (iCol - 1)*(patch_w + grid_spacing(2)) + 1;
        tiled_patches(r_idx + (0:patch_h-1), c_idx + (0:patch_w-1), :) = patches(:,:,:,iPtch); 
        iPtch = iPtch + 1;
    end
end
    
end 