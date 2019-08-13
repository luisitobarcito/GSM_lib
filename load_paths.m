% custom path according to your own installation
gsm_lib_path = fileparts(which('load_paths.m'));
src_path = fullfile(gsm_lib_path, 'src');
addpath(genpath(src_path));
