function data = loadTiffMovies()

[filename, pathname] = uigetfile( ...
{  '*.tif';}, ...
   'Pick TIFF movies', ...
   'MultiSelect', 'on');

if isequal(filename, 0)
    error('No files are loaded.');
end

data = struct([]); %memory allocation

if iscell(filename) == 0 %cell check
    data{1,1} = pathname; %path name
    data{2,1} = filename; %file name
    data{3,1} = tiffreadVolume([pathname filename]); %load image
else
    n = length(filename);
    for i = 1:n
        data{1,i} = pathname; %path name
        data{2,i} = filename{i}; %file name
        data{3,i} = tiffreadVolume([pathname filename{i}]); %load image
    end
end