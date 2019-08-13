function llist = listLinesFromText(filename, fullpath)
%LISTLINESFROMTEXT - Retruns a cell array with all the lines line from txt file
% Reads a txt file lien by line and returns a one dimensional cell array, where
% each cell contains a string with the corresponding line of text from txt file 
%
% Syntax:  [llist] = listLinesFromText(filename)
%
% Inputs:
%    filename - .txt file with '\n' newline charater
%
% Outputs:
%    llist - cell array with N entries. i-th cell contains the string
%    corresponding to the i-th line from the txt file. '\n' character is not 
%    included in the string.
%
% Author: Luis Gonzalo Sanchez Giraldo
% October 2018; Last revision: 22-October-2018

if ~exist('fullpath', 'var')
    fullpath = true;
end

fid = fopen(filename);
tline = fgetl(fid);
llist = cell(0,1);
while ischar(tline)
    if fullpath
        llist{end+1,1} = tline;
    else
        [~, name, ext] = fileparts(tline);
        llist{end+1,1} = strcat(name, ext);
    end
    tline = fgetl(fid);
end
fclose(fid);
end