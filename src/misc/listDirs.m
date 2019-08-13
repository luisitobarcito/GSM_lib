function [dir_list] = listDirs(root_path, relative)
%LISTDIRS - Recursively lists all the leaf dirs from root_path
% List directories that do not contain any directories inside them
%
% Syntax:  [dir_list] = listDirs(root_path)
%
% Inputs:
%    root_path - starting path 
%    relative - (Optional) List dirs relative to root_path if set true (default)
%    Otherwise it will append the root_path to each element of the list
%
% Outputs:
%    dir_list - list of all directories inside root_path
%
% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 03-October-2018

if ~exist('relative', 'var')
    relative = true;
end

% recursive search is wrapped to avoid appending root_path
if relative
    current_dir = pwd();
    cd(root_path);
    dir_list = listDirsAbs();
    cd(current_dir);
else
    dir_list = listDirsAbs();
end

if isempty(dir_list)
    dir_list{1} = './';
end

end

function [dir_list] = listDirsAbs(root_path)
if exist('root_path', 'var')
    dir_struct = dir(root_path);
else
    dir_struct = dir();
    root_path = [];
end
dir_list = {};
has_dirs = false;
for iCat = 1: length(dir_struct)
    is_dot = any(strcmp(dir_struct(iCat).name, {'.', '..'}));
    if dir_struct(iCat).isdir && ~is_dot
        new_root = fullfile(root_path, dir_struct(iCat).name, filesep);
        dir_list = cat(1, dir_list, listDirsAbs(new_root));
        has_dirs = true;
    end
end
if ~has_dirs
    dir_list = root_path;
end
end
