function [img_out] = rescaleToSquareImage(img_in, sz, rel_dim)
%RESCALETOSQUAREIMAGE
% Image is reescaled and cropped to squared Image of size 'sz'
% The reescaling is chosen accoriding to rel_dim which can be 'shortest' or
% 'largest' 
% - 'shortest': rescales with respect to the shortest dimension
%   and crops the largest dimension symmetrically from the center of the image
% - 'largest': rescales with respect to the largest dimension
%   and zero pads the shorter dimension simetrically from the center of the 
%   image
%
% Syntax:  [img_out] = rescaleToSquareImage(img_in, sz, rel_dim)
%
% Inputs:
%    img_in - input image
%    sz - desired size of the squre image
%    rel_dim - dimension to which the image is rescaled
%
% Outputs:
%    img_out - square image of size sz
%


% Author: Luis Gonzalo Sanchez Giraldo
% August 2019; Last revision: Aug-12-2018

r = size(img_in, 1);
c = size(img_in, 2);
img_out = zeros(sz, sz, size(img_in, 3));
switch(rel_dim)
    case 'largest'
        if max(r,c) == r
           img_temp = imresize(img_in, [sz, NaN]);
           lpad_sz = floor((sz - size(img_temp, 2))/2);
           img_out(:, lpad_sz + (1:size(img_temp,2)), :) = img_temp;
        else
           img_temp = imresize(img_in, [NaN, sz]);
           upad_sz = floor((sz - size(img_temp, 1))/2);
           img_out(upad_sz + (1:size(img_temp,1)), :, :) = img_temp;
        end
    case 'shortest'
        if min(r,c) == r
           img_temp = imresize(img_in, [sz, NaN]);
           lcrop_sz = floor((size(img_temp, 2) - sz)/2);
           img_out = img_temp(:, lcrop_sz + (1:sz), :);
        else
           img_temp = imresize(img_in, [NaN, sz]);
           ucrop_sz = floor((size(img_temp, 1) - sz)/2);
           img_out = img_temp(ucrop_sz + (1:sz), :, :);
        end
    otherwise
        error('Method %s not implemented', rel_dim);
end
        
