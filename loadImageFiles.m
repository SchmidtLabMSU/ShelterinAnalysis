function data = loadImageFiles()
%--------------------------------------------------------------------------
% Calling this function will pop-up window for selecting image files that
% will be read and stored to the structure variable data, where the first
% row will contain pathname, second row filename and third row the actual
% image.
% 
% 
% Created by: Tomas Janovic, 2022
%-------------------------------------------------------------------------- 

[filename, pathname] = uigetfile( ...
{  '*.jpg;*.jpeg;*.tif;*.png;*.gif;*.bmp','All Image Files';}, ...
   'Pick images', ...
   'MultiSelect', 'on');

if isequal(filename, 0) %if no file is selected
    error('No files are loaded.');
end

data = struct([]); %memory allocation

if iscell(filename) == 0 %cell check
    data{1,1} = pathname; %path name
    data{2,1} = filename; %file name
    data{3,1} = imread([pathname filename]); %load image
else
    n = length(filename);
    for i = 1:n
        data{1,i} = pathname; %path name
        data{2,i} = filename{i}; %file name
        data{3,i} = imread([pathname filename{i}]); %load image
    end
end