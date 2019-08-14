%% This is an example script that normalizes the tensors contained in the files as specified in
% the "data_files" string.

% a txt file where each line has the name of a file to be processed
data_files = './example_maplist.txt'
% location of the files to be processed
data_folder = './'
% name of the normalization model to either process data or for training
model_file = './example_trained_normalization_model.mat'
% number of samples per image 
samples_per_map = 10
% number of train maps to be considered. "samples_per_map" are extratcted per map
n_train_maps = 50000
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

run ./learn_normalization_from_maps.m;
