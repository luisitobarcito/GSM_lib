%% This is an example script that normalizes the tensors contained in the files as specified in
% the "data_files" string.

% a txt file where each line has the name of a file to be processed
data_files = './example_maplist.txt'
% location of the files to be processed
data_folder = './'
% name of the normalization model to either process data or for training
model_file = './alexnet_ilsvcr2012_flex_norm_conv2_d4.mat'
% size of a batch for running the normalization function
batch_size = 2
% number of center neighbors
n_cross_neigh = 4
% number of groups in the map
n_groups = 2
% distance from center of the surround units
p_dist = 4
% spatial configuration of the surround neighborhood
spatial_conf = 'circular'

run ./run_norm_activations_from_maps.m;
