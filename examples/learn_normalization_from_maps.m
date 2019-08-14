% Train center-surround covariance model from a given list of conv map array
% the script requires:
%   -- A .txt file where each line is the path to a training data file.
%   -- "model_file" string with the name to be given to the trained
%      normalization model.

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

if ~exist('samples_per_map', 'var')
    samples_per_map = 10;
end

if ~exist('n_train_maps', 'var')
    n_train_maps = 50000;
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

%% preallocate training data array 

%% Create data set for training flexible normalization covariance
% Define center and surround neighborhoods

s_neigh = surroundNeighIdx(p_dist, spatial_conf);
% define center neighborhoods 
c_neigh = mod(bsxfun(@plus, repmat((1:n_conv)', 1, n_cross_neighh + 1), -ceil(n_cross_neighh/2):floor(n_cross_neighh/2)) - 1, n_conv/n_groups) + 1;
c_neigh = permute(reshape(c_neigh, [n_conv/n_groups, n_groups, n_cross_neighh + 1]), [1 3 2]);
n_surround = size(s_neigh, 1);
n_dim = (n_cross_neighh + 1  + n_surround);

%% Training data set creation
% pre-allocate array to fill with values from deep network

data_sz = samples_per_map*n_train_maps; % pre-allocate data by adding chunks of this size to memory 
train_data = zeros(data_sz, n_dim, n_conv/n_groups, n_groups);


data_ptr = 1;
file_ptr = 1;
%% Compute conv outputs before rectification and normalization
% Loop through the dataset select only n_train_images_per_folder at random
all_map_files = listLinesFromText(data_files, false);

% select n_train_per_folder images at random from each folder
rnd_idx = randperm(length(all_map_files));

while data_ptr <= data_sz
    %% read a map file gather data from it
    map_filename = fullfile(data_folder, all_map_files{rnd_idx(file_ptr)});
    data = load(map_filename);
    data_shape = size(data);
    n_conv = data_shape(end);
    n_maps = data_shape(1);
    %%% extract conv1 features using AlexNet
    tmp_map = permute(data, [2, 3, 4, 1]);
    tmp_map = reshape(tmp_map, size(tmp_map, 1), size(tmp_map, 2), size(tmp_map,3)/n_groups, n_groups, []);
    % transform conv1 data into v1 map format
%     tmp = computeCenterSurround(tmp_map, c_neigh, s_neigh, V1RF_sz, samples_per_map);
    tmp = computeCenterSurround(tmp_map, c_neigh, s_neigh, samples_per_map);
    
    % subsample the maps
    train_data(data_ptr: min(data_ptr + samples_per_map*n_maps - 1, datat_sz),:,:,:) = reshape(permute(tmp, [1, 5, 2, 3, 4]), samples_per_map*n_maps, size(tmp,2), size(tmp, 3), size(tmp, 4));
    data_ptr = data_ptr + samples_per_map*n_maps;
    fprintf(' . ');
    if mod(iBtch, 25) == 0
        fprintf('\n');
    end
    file_ptr = file_ptr + 1;
end
fprintf('\n')
%% Train the center-surround Mixture of GSMs
Kcs = size(train_data, 2);                                     % # center-surround units per channel 
Kc = size(c_neigh, 2);                    % # center units
Ks = Kcs - Kc;                                                 % # surround units
indc = 1:Kc;      % indices of center units, phase 1
inds = length(indc) + (1 : Ks);      % indices of surround units, phase 1
USE_PARFOR = true;
STORE_RESULTS = false;
%% run CEM algorithm for each channel(orientation) and group (pyramid level) 
if USE_PARFOR
    parpool(8)
end
mgsm  = cell( [n_conv/n_groups, n_groups]);
nbrsred = cell( [n_conv/n_groups, n_groups]);
factor = cell( [n_conv/n_groups, n_groups]);
data_mean = cell( [n_conv/n_groups, n_groups]);
scale_data = 1;
for iGrp = 1:n_groups
    if USE_PARFOR
        parfor iCh = 1: n_conv/n_groups
            %%% center and rescale data
            data_mean{iCh, iGrp} = mean(train_data(:, :, iCh, iGrp), 1);
            data = bsxfun(@minus, train_data(:, :, iCh, iGrp), data_mean{iCh, iGrp});
            factor{iCh, iGrp}(indc) = std(data(:, indc));
            factor{iCh, iGrp}(inds) = factor{iCh, iGrp}(indc(n_cross_neighh/2 + 1));
            factor{iCh, iGrp} = scale_data./factor{iCh, iGrp}
            data = data*diag(factor{iCh, iGrp});
            %%% create neighbours array (for the symmetry constraint)
            nbrsred{iCh, iGrp} = zeros(Kcs, 4);
            nbrsred{iCh, iGrp}(:, 1) = iGrp - 1;
            nbrsred{iCh, iGrp}(indc, 2) = c_neigh(iCh, :, iGrp)' - 1;
            nbrsred{iCh, iGrp}(inds, 2) = iCh - 1;
            nbrsred{iCh, iGrp}(inds, 3:4) = s_neigh;
            name = 'younameit';
            theeps = 1e-10;
            centcov=5;
            %%% START CEM
            [COVcs, INVCOVcs, COVc, INVCOVc, COVs, INVCOVs, passign, f] = ITERSGSM_EM_IndFactorizedFixedPoint(Kc, Ks, indc, inds, data, theeps, 100);
            %%% store MGSM for each orientation and level
            mgsm{iCh, iGrp} .COVcs = COVcs;
            mgsm{iCh, iGrp}.COVc = COVc;
            mgsm{iCh, iGrp}.COVs = COVs;
            mgsm{iCh, iGrp}.INVCOVcs = INVCOVcs;
            mgsm{iCh, iGrp}.INVCOVc = INVCOVc;
            mgsm{iCh, iGrp}.INVCOVs = INVCOVs;
            mgsm{iCh, iGrp}.Pind = passign(end);
            mgsm{iCh, iGrp}.passing = passign;
            mgsm{iCh, iGrp}.f = f;
        end
    else
        for iCh = 1: n_conv/n_groups
            fprintf('Running channel %d from group %d\n', iCh, iGrp);
            %%% center and rescale data
            data_mean{iCh, iGrp} = mean(train_data(:, :, iCh, iGrp), 1);
            data = bsxfun(@minus, train_data(:, :, iCh, iGrp), data_mean{iCh, iGrp});
            factor{iCh, iGrp}(indc) = std(data(:, indc));
            factor{iCh, iGrp}(inds) = factor{iCh, iGrp}(indc(n_cross_neighh/2 + 1));
            factor{iCh, iGrp} = scale_data./factor{iCh, iGrp}
            data = data*diag(factor{iCh, iGrp});
            %%% create neighbours array (for the symmetry constraint)
            nbrsred{iCh, iGrp} = zeros(Kcs, 4);
            nbrsred{iCh, iGrp}(:, 1) = iGrp - 1;
            nbrsred{iCh, iGrp}(indc, 2) = c_neigh(iCh, :, iGrp)' - 1;
            nbrsred{iCh, iGrp}(inds, 2) = iCh - 1;
            nbrsred{iCh, iGrp}(inds, 3:4) = s_neigh;
            name = 'younameit';
            theeps = 1e-10;
            centcov=5;
            %%% START CEM
            [COVcs, INVCOVcs, COVc, INVCOVc, COVs, INVCOVs, passign, f] = ITERSGSM_EM_IndFactorizedFixedPoint(Kc, Ks, indc, inds, data, theeps, 100);
            %%% store MGSM for each orientation and level
            mgsm{iCh, iGrp} .COVcs = COVcs;
            mgsm{iCh, iGrp}.COVc = COVc;
            mgsm{iCh, iGrp}.COVs = COVs;
            mgsm{iCh, iGrp}.INVCOVcs = INVCOVcs;
            mgsm{iCh, iGrp}.INVCOVc = INVCOVc;
            mgsm{iCh, iGrp}.INVCOVs = INVCOVs;
            mgsm{iCh, iGrp}.Pind = passign(end);
            mgsm{iCh, iGrp}.passing = passign;
            mgsm{iCh, iGrp}.f = f;
        end
    end
end
if USE_PARFOR
    delete(gcp) 
end

%% Save model 
sind = true;
if exist(model_file, 'file')
    warning(sprintf('model %s already exists results will overwrite existing model', model_file));
end
save(model_file, 'mgsm', 'Ks', 'Kc', 'c_neigh', 's_neigh', 'nbrsred', 'sind', 'indc', 'inds', 'factor', 'data_mean', 'p_dist');
