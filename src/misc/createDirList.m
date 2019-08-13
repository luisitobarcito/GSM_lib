function dir_list = createDirList(name_list, dir_pattern)
%CREATEDIRLIST - Creates a list of directories from a list of names
%For each entry in name_list a directory the matches the dir pattern is
%created. if no dir_pattern is specified the function creates a list of 
%direcotries using the only name_list   
%
% Syntax:  [dir_list] = createDirList(name_list, dir_pattern)
%
% Inputs:
%    name_list - cell array of strings. Each entry contains one name of the
%    list 
%    dir_pattern - (optional) Is a prefixed pattern to be completed with
%    the name information. the place where the name should be added is
%    denoted by the string '<name>' in the pattern.
%
% Outputs:
%    dir_list - cell array with the same size as name_list. Each entry of
%    the array is the result of replacing the occurrences of '<name>' with
%    the corresponding name in the list. 
%
% Example:
%    name_list = {'name1'; 'name2'; 'name3'};
%    dir_pattern = './somedir/<name>/someothedir';
%    dir_list = createDirList(name_list, dir_pattern);
%    disp(dir_list)
%    './somedir/name1/somedir'    
%    './somedir/name2/somedir'    
%    './somedir/name3/somedir'
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 03-October-2018
if ischar(name_list)
    name_list = {name_list};
elseif ~iscell(name_list)
    error('name_list must be a string or a cell array of strings');
end

if ~exist('dir_pattern', 'var')
    dir_pattern = strcat('.', filesep, '<name>');
end
dir_list = cell(size(name_list));
for i_name = 1 : numel(name_list)
    dir_list{i_name} = regexprep(dir_pattern, '<name>', name_list{i_name});
end
end