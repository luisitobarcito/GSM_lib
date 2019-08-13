function success = checkDir(dir_names, create_dir, verbose)
%CHECKDIR - checks if dir_name is an existing directory or creates one
%
% Syntax:  [success] = checkDir(dir_name, create_dir)
%
% Inputs:
%    dir_names - string or cell array with the names of the directories to be 
%               checked
%    create_dir - (optional) bool. If true checkDir creates a directory
%               with name 'dir_name' if directory did not exist.
%
% Outputs:
%    success - bool. True if directory exist
%
% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 15-Oct-2018 

if ~exist('create_dir', 'var')
    create_dir = false;
end
if ~exist('verbose', 'var')
    verbose = false;
end

if ischar(dir_names)
    dir_names = {dir_names};
elseif ~iscell(dir_names)
    error('dir_names must be a string or a cell array of strings');
end
for iDr = 1 : numel(dir_names)
    dir_name = dir_names{iDr};
    if verbose == true
        fprintf('Checking : %s\n', dir_name);
    end
    if isdir(dir_name)
        success = true;
    elseif create_dir
        if verbose == true
            fprintf('Directory did not exist. Creating directory : %s\n', dir_name);
        end
        success = mkdir(dir_name);
    else
        if verbose == true
            fprintf('Directory did not exist\n');
        end   
        success = false;
    end
    if false == success
        break
    end
end

end