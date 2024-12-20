%% ------------------------- FOCI ANALYSIS --------------------------------
% Script for foci analysis.
% 
% 
% Created by: Tomas Janovic, 2024
%-------------------------------------------------------------------------- 

clear all; close all; clc;

%% --------------------------- USER INPUT ---------------------------------
% Type of analysis:
% 1 - single channel foci analysis
% 2 - co-localization analysis
% 3 - single channel foci analysis with additional masks from imagesCH2 (eg. analyze only nuclear foci of GFP-positive cells)
analysisType = 1;

% Path to input files
localizedFociPath = '/Volumes/Tom SSD/Research/Imaging/2024/2024-08-14 - HeLa13-photobleaching abundance-Sphase - 2 replica/Halo-TIN2 c33/Localized foci/'; % localized foci .xlsx files
imagesPath = '/Volumes/Tom SSD/Research/Imaging/2024/2024-08-14 - HeLa13-photobleaching abundance-Sphase - 2 replica/Halo-TIN2 c33/Halo/'; % your main images where foci were localized .tif files
nuclearLabelsPath = '/Users/admin/Documents/Research/PosDoc/Imaging/Segmented files/'; % segmented cells from cellpose .mat files

localizedFociPath_CH2 = '...';
imagesCH2Path = '...'; % your CH2 images that can be used either for colocalization or as a marker .tif files

% Name of the output .xlsx file
xlsxName = 'HeLa13-Halo-TIN2-c33- 2.xlsx';

% Size threshold (number of pixels) for segmented label to call it a nuclei - filtering out micronuclei etc.
sizeThreshold = 12000;

% Intensity threshold for colocalization - foci in CH2 images will be
% called as colocalization if the foci intensity in CH2 is above threshold
CoLocThreshold = 200;

% Intesity threshold for CH2 images - only for type 3 analysis
IntCH2Threshold = 1500;

%Photobleaching step
photobleachingStep = 108.1226667;

%% ----------------------- PROCESSING  TYPE 1 -----------------------------
if analysisType == 1
    data = analyzeFoci(localizedFociPath, imagesPath, nuclearLabelsPath, sizeThreshold); % analyze foci  
    [fociNr, fociStats, fociIntensities] = fociStatistics(data); % statistics of analyzed foci

    % ---------- Save/export ----------
    writecell(fociNr, xlsxName, 'Sheet', 'Foci Number', 'Range', 'A1');
    writecell(fociStats, xlsxName, 'Sheet', 'Foci intensity stats', 'Range', 'A1');
    writecell(fociIntensities, xlsxName, 'Sheet', 'Raw foci intensities', 'Range', 'A1');

    % ---------- Visualization ----------
    fociNr = cell2mat(fociNr(2:end, 3));
    figure(1)
    boxplot(fociNr)
    ylabel('Foci number/cell');
    
    fociInt = cell2mat(fociStats(2:end, 3));
    figure(2)
    boxplot(fociInt)
    ylabel('Median foci intensity/cell');

    allIntensities = cell2mat(fociIntensities(2:end,4));
    allIntensities = allIntensities./photobleachingStep;
    allIntensitiesFiltered = allIntensities(allIntensities>=1);
    medAbundance = median(allIntensitiesFiltered);
    hCounts = histcounts(allIntensitiesFiltered ,'BinWidth', 1, 'Normalization', 'probability');
    x = 1:length(hCounts);
    
    disp(['Median protein abundance = ' num2str(medAbundance)]);

    figure(3) % plotting the histogram with a median
    bar(x, hCounts, 1)
    hold on
    plot([medAbundance medAbundance], [0 0.028], 'LineWidth', 4, 'Color', 'black');
    hold off
    xlim([0 140]);
    ylim([0 0.03]);
    ylabel('Probability (%)');
    xlabel('Protein abundance');
    set(gca, 'FontSize', 17, 'FontWeight', 'bold', 'YMinorTick', 'off', 'XMinorTick', 'off', 'box', 'on', LineWidth = 1.5);
    set(gcf, 'position', [800, 800, 640, 300]);
    
%% ----------------------- PROCESSING  TYPE 2 -----------------------------
elseif analysisType == 2
    data = colocalization(localizedFociPath, localizedFociPath_CH2, nuclearLabelsPath, sizeThreshold);
    coLocalizationNr = {'Filename', 'Cell Nr.', 'Number of colocalizations'};

    for i = 1:length(data(1,:))
        filename = data{1,i};
        values = data{2,i}';

        for j = 1:length(values)
            coLocalizationNr{end+1,1} = filename;
            coLocalizationNr{end,2} = j;
            coLocalizationNr{end,3} = values(j);
        end
    end

    % ---------- Save/export ----------
    writecell(coLocalizationNr, xlsxName, 'Sheet', 'Colocalization', 'Range', 'A1');

    % ---------- Visualization ----------
    colocalizations = cell2mat(coLocalizationNr(2:end, 3));
    figure(1)
    boxplot(colocalizations)
    ylabel('Colocalizations/cell');

%% ----------------------- PROCESSING  TYPE 3 -----------------------------
elseif analysisType == 3
    data = analyzeFociForCH2Positive(localizedFociPath, imagesPath, nuclearLabelsPath, imagesCH2Path, sizeThreshold, IntCH2Threshold); %analyze foci
    [fociNr, fociStats, fociIntensities] = fociStatistics(data); % statistics of analyzed foci

    % ---------- Save/export ----------
    writecell(fociNr, xlsxName, 'Sheet', 'Foci Number', 'Range', 'A1');
    writecell(fociStats, xlsxName, 'Sheet', 'Foci intensity stats', 'Range', 'A1');
    writecell(fociIntensities, xlsxName, 'Sheet', 'Raw foci intensities', 'Range', 'A1');

    % ---------- Visualization ----------
    fociNr = cell2mat(fociNr(2:end, 3));
    figure(1)
    boxplot(fociNr)
    ylabel('Foci number/cell');
    
    fociInt = cell2mat(fociStats(2:end, 3));
    figure(2)
    boxplot(fociInt)
    ylabel('Median foci intensity/cell');
end