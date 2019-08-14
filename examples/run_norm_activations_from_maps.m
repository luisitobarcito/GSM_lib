% This script takes a list of .mat files each containing 4 -dimensional array
% of size [N, H, W, C] and normalizes the maps using the model specified in
% "model_file" variable
%   -- A .txt file where each line is the path to a training data file
%   -- A "model_file" string with the name of the file contianing a trained
%      normalization model


% Luis Gonzalo Sanchez Giraldo 2019

%%
% Load data set
if ~exist('data_files', 'var')
    error('data_files  must be specified');
end

if ~exist('data_folder', 'var')
    error('data_folder  must be specified');
end

if ~exist('model_file', 'var')
    error('model_file  must be specified');
end

if ~exist('batch_size', 'var')
    batch_size = 100;
end
% Define center and surround neighborhoods
if ~exist('n_cross_neigh', 'var')
    error('n_cross_neigh  must be specified');
end

if ~exist('n_groups', 'var')
    error('n_groups  must be specified');
end

if ~exist('p_dist', 'var')
    error('p_dist  must be specified');
end

if ~exist('spatial_conf', 'var')
    spatial_conf = 'circular';
end

all_map_files = listLinesFromText(data_files, false);
datastr = load(fullfile(data_folder, all_map_files{1}), 'data');
data = datastr.data;
n_train_maps = size(data, 1);
conv_shape = size(data);
n_conv = conv_shape(end);
    

%% Create datasets with normalized v1 activations
s_neigh = surroundNeighIdx(p_dist, spatial_conf);
% define center neighbohouds based on cifar cross normalization pool topology
c_neigh = mod(bsxfun(@plus, repmat((1:n_conv)', 1, n_cross_neigh + 1), -ceil(n_cross_neigh/2):floor(n_cross_neigh/2)) - 1, n_conv/n_groups) + 1;
c_neigh = permute(reshape(c_neigh, [n_conv/n_groups, n_groups, n_cross_neigh + 1]), [1 3 2]);
n_surround = size(s_neigh, 1);
n_dim = (n_cross_neigh + 1  + n_surround);
whiten_G = true; % whiten inference

%% Compute conv1 outputs before rectification and normalization
% Loop through the dataset select only n_train_images_per_folder at random
gsm_params = load(model_file);

for iFl = 1 : length(all_map_files)
    datastr = load(fullfile(data_folder, all_map_files{iFl}));
    data = datastr.data;
    n_train_maps = size(data, 1);
    conv_shape = size(data);
    n_conv = conv_shape(end);
    crt_batch_size =  min(batch_size, conv_shape(1));
    n_batches = floor((n_train_maps -crt_batch_size) / batch_size) + 1;
    conv_norm_data = zeros(conv_shape([2,3,4,1]));
    for iBtch = 1 : n_batches
        iImg = (iBtch - 1)*crt_batch_size + 1;
        % rescale image to square image of size sz
        conv1_im = data(iImg:(iBtch*crt_batch_size), :, :, :);
        %%% extract conv1 features using AlexNet
        tmp_map = permute(conv1_im, [2, 3, 4, 1]);
        tmp_map = reshape(tmp_map, size(tmp_map, 1), size(tmp_map, 2), size(tmp_map,3)/n_groups, n_groups, []);
        % transform conv1 data into v1 map format
        % subsample the maps
        % transform conv1 data into v1 map format
        G = normalizeMap(tmp_map, @flexibleNormalization, gsm_params);
        % save the maps as flattened
        batch_idx = iImg : min(iImg + crt_batch_size - 1, conv_shape(1)); 
        conv_norm_data(:,:,:, batch_idx) = reshape(G, size(tmp_map, 1), size(tmp_map, 2), conv_shape(end), []);
        fprintf(' . ');
        if mod(iBtch, 25) == 0
            fprintf('\n');
        end
    end % train images per class
    fprintf('\n');
    fprintf('Done normalizing data \n');
    norm_data = permute(conv_norm_data, [4, 1, 2, 3]);
    norm_data = single(real(norm_data));
    norm_data_size = size(norm_data);
    disp(sprintf('Saving normalized array size %d, %d, %d, %d', norm_data_size(1),...
        norm_data_size(2), norm_data_size(3), norm_data_size(4)));
    
    [~, filename, fileext] = fileparts(all_map_files{iFl});
    save(fullfile(data_folder, strcat(filename, '_normalized', fileext)), ...
        'norm_data', '-v7.3');
end
