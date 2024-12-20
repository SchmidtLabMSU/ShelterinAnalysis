function data = analyzeFoci(localizedFociPath, imagesPath, nuclearLabelsPath, sizeThreshold)
%--------------------------------------------------------------------------
% Function for foci analysis. Input is path as a string to: foci
% localizations from ComDet (FiJi plugin) as a *.csv file, then path to
% images (*.tif files) and segmented nuclei labels from cellpose (*.mat
% files). Function will go through all images and through all cells and
% count the number of foci per cell with their respective intensities. It
% will store everything in the structure data, where first row contains
% image filename and second row analyzed values (foci number and intensity
% for each foci with respective cell index).
% 
% 
% Created by: Tomas Janovic, 2024
%-------------------------------------------------------------------------- 

%Get all filenames and allocate final structure
fociDir = dir([localizedFociPath, '*.csv']);
imagesDir = dir([imagesPath, '*.tif']);
nucLabelsDir = dir([nuclearLabelsPath, '*.mat']);
data = struct([]);

for i = 1:length(imagesDir) %for each image
    nucLabels = load([nuclearLabelsPath nucLabelsDir(i).name]).labels; %load labels from cellpose segmentation
    localizations = readtable([localizedFociPath fociDir(i).name]); %load localization of foci
    img = imread([imagesPath imagesDir(i).name]); %load an image

    labelsNew = zeros(size(nucLabels)); %allocation of new matrix for labels
    
    for j = 1:max(nucLabels(:)) %go through all segmented cells
        tempLabel = nucLabels;
        tempLabel(tempLabel ~= j) = 0; %suppress all cells except the current one
        tempLabel = imclearborder(tempLabel); %suppress label if it is connected to image border
        if max(tempLabel(:)) > 0 %if the label is still present
            if nnz(tempLabel) >= sizeThreshold %if the label has more than thresholdSize pixels
                labelsNew = labelsNew + tempLabel; %add this label to the output matrix

                % Background subtraction from a current cell
                intensityLabel = img;
                intensityLabel(tempLabel ~= j) = 0;
                labelImgIntensities = intensityLabel(intensityLabel > 0);
                h = histogram(labelImgIntensities, 1024);
                [~, whichbin] = max(h.Values);
                background = h.BinEdges(whichbin);
                img(intensityLabel > 0) = img(intensityLabel > 0) - background;
            end
        end
    end
    
    %Allocation
    values = struct([]);
    fociNr = [];
    fociIntensity = [];

    uniqueIndex = unique(labelsNew);
    tempIndex = 1;
    for j = 2:length(uniqueIndex) %for each segmented cell
        label = uniqueIndex(j);
        tempLabel = labelsNew;
        tempLabel(tempLabel~=label) = 0; %suppress all cells except the current one
        fociCount = 0;
        for k = 1:height(localizations) %for each foci
            xmin = localizations{k,8};
            xmax = localizations{k,10};
            ymin = localizations{k,9};
            ymax = localizations{k,11};
            labelROI = tempLabel(ymin:ymax, xmin:xmax);
            if sum(labelROI(:)) > 0 %if the foci position is within the current cell
                fociCount = fociCount + 1; %increment foci number
                imgROI = img(ymin:ymax, xmin:xmax); %get ROI of this foci position from original image
                fociIntensity(1, tempIndex) = j-1;
                fociIntensity(2, tempIndex) = max(imgROI(:)); %store maximum intentensity
                tempIndex = tempIndex + 1;
            end
        end
        fociNr(1,j-1) = fociCount; %store foci number for the cell
    end

    %Store all processed data
    values{1,1} = fociNr;
    values{2,1} = fociIntensity;
    data{1,i} = imagesDir(i).name;
    data{2,i} = values;
end
end