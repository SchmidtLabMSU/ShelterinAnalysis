%% ---------------------- NUCLEI SEGMENTATION -----------------------------
% Script for nuclei segmentation using cellpose AI algorithm. In order to
% use it, you need to have installed cellpose interface and download cyto2
% or nuclei training model. As an input parameter the average cell diameter
% is needed. Results (labeled masks) will be saved to a folder as a *.mat
% files.
% 
% 
% Created by: Tomas Janovic, 2024
%-------------------------------------------------------------------------- 

clear all; close all; clc;

%% ---------- USER INPUT ----------
averageCellDiameter = 125; %average diameter of nuclei in px for segmentation (125 default)
plotResults = 1; % 0 - do not plot results, 1 - plot results
saveAsNuclearMaskTiff = 0; % assign the same label for all segmented cells

%% ---------- PROCESSING ----------
data = loadImageFiles(); %load Nuclei files
n = length(data(1,:)); %number of files loaded

%Saving
pathname = [data{1,1} 'Nuclear masks/']; %pathname for saving
if ~exist(pathname,'dir') %if folder does not exist
    mkdir(pathname); %create a folder
end

disp('Loading complete.');

%Processing files
cp = cellpose(Model = "cyto2"); %initialize cellpose model cyto2 or nuclei
for i = 1:n
    img = im2double(data{3,i}); %change to double
    
    %Normalization [0-1]
    imgNorm = img-min(img(:));
    imgNorm = imgNorm./max(imgNorm(:));
    
    %Segmentation with cellpose
    labels = segmentCells2D(cp, imgNorm, ImageCellDiameter = averageCellDiameter);

    %Saving
    fileName = data{2,i};
    fileName = fileName(1:end-4);
    if saveAsNuclearMaskTiff == 1 %Saving labels as binary tiff files
        labels(labels > 0) = 1;
        imwrite(im2uint8(labels), [pathname fileName '.tiff'], 'tiff');
    else %Saving labels as .mat files
        save([pathname fileName '.mat'], "labels");
    end

    %Plot results
    if plotResults == 1
        imgOverlay = labeloverlay(img,labels);
        figure(i)
        subplot(121)
        imshow(img, []);
        hold on
        subplot(122)
        imshow(imgOverlay);
        hold off
    end
    disp(['File ' num2str(i) ' out of ' num2str(n) ' segmented and saved.']);
end

