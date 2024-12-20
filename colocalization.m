function data = colocalization(localizedFociPath, localizedFociPath_CH2, nuclearLabelsPath, sizeThreshold)
%--------------------------------------------------------------------------
% Function for colocalization analysis Input is path as a string to: foci
% localizations from ComDet (FiJi plugin) as a *.csv file, then another
% localizations file for second channel, segmented nuclei labels from
% cellpose (*.mat files) and size threshold for the cells to analyze.
% Analysis will count the number of co-localization per each cell. The
% co-localization is determined that two detected foci are closer than 3px.
% It will store everything in the structure data, where first row contains
% image filename and second row analyzed values (foci number and intensity
% for each foci with respective cell index).
% 
% 
% Created by: Tomas Janovic, 2024
%-------------------------------------------------------------------------- 

%Get all filenames and allocate final structure
fociDir = dir([localizedFociPath, '*.csv']);
fociDir_CH2 = dir([localizedFociPath_CH2, '*.csv']);
nucLabelsDir = dir([nuclearLabelsPath, '*.mat']);
data = struct([]);


for i = 1:length(nucLabelsDir) %for each image
    disp(['Processing image ' num2str(i) '/' num2str(length(nucLabelsDir)) '.']);
    nucLabels = load([nuclearLabelsPath nucLabelsDir(i).name]).labels; %load labels from cellpose segmentation
    localizations = readtable([localizedFociPath fociDir(i).name]); %load localization of foci
    localizations_CH2 = readtable([localizedFociPath_CH2 fociDir_CH2(i).name]); %load localization of foci for CH2

    labelsNew = zeros(size(nucLabels)); %allocation of new matrix for labels

    for j = 1:max(nucLabels(:)) %go through all segmented cells
        tempLabel = nucLabels;
        tempLabel(tempLabel ~= j) = 0; %suppress all cells except the current one
        tempLabel = imclearborder(tempLabel); %suppress label if it is connected to image border
        if max(tempLabel(:)) > 0 %if the label is still present
            if nnz(tempLabel) >= sizeThreshold %if the label has more than thresholdSize pixels
                labelsNew = labelsNew + tempLabel; %add this label to the output matrix
            end
        end
    end
    
    %Allocation
    coLocFociNr = [];

    uniqueIndex = unique(labelsNew);
    for j = 2:length(uniqueIndex) %for each segmented cell
        label = uniqueIndex(j);
        tempLabel = labelsNew;
        tempLabel(tempLabel~=label) = 0; %suppress all cells except the current one
        coLocFociCount = 0;
        for k = 1:height(localizations)
            xmin1 = localizations{k,8};
            xmax1 = localizations{k,10};
            ymin1 = localizations{k,9};
            ymax1 = localizations{k,11};
            labelROI = tempLabel(ymin1:ymax1, xmin1:xmax1);
            if sum(labelROI(:)) > 0 %if the foci position is within the current cell
                x1Center = localizations{k,3};
                y1Center = localizations{k,4};
                for l = 1:height(localizations_CH2)
                    xmin2 = localizations_CH2{l,8};
                    xmax2 = localizations_CH2{l,10};
                    ymin2 = localizations_CH2{l,9};
                    ymax2 = localizations_CH2{l,11};
                    x2Center = localizations_CH2{l,3};
                    y2Center = localizations_CH2{l,4};
                    distance = sqrt((x1Center-x2Center)^2 + (y1Center-y2Center)^2);
                    if distance <= 3 %if the foci_CH2 position is within the current foci from CH1 (distance is less than 3 px)
                        coLocFociCount = coLocFociCount + 1; %increment colocalization
                    end
                end
            end
        end
        coLocFociNr(1,j-1) = coLocFociCount; %store colocalization number for current cell
    end

    %Store all processed data
    data{1,i} = nucLabelsDir(i).name;
    data{2,i} = coLocFociNr;
end
disp('Colocalization analysis is finished.');
end