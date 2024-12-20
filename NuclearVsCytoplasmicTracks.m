%% ----------------- NUCLEAR VS CYTOPLASMIC TRACKS ------------------------
% Script that will segment single-particle trajectories (trackPar) as
% either nuclear or cytoplasmic tracks based on provided masks (*.tiff
% files)
% 
% 
% Created by: Tomas Janovic, 2024
%--------------------------------------------------------------------------

clear all, close all, clc;

%% Input data
%path to the folders where all tracks (.mat) and masks (.tiff) are stored
input_pathname_tracks = 'G:\Research\Imaging\2024\2024-12-02 - HeLa-Residence time - 3 replica\500 ms\Halo-TRF2 - RAP1-KO c32\Tracking\Raw tracks 1500 ms\';
input_pathname_masks = 'G:\Research\Imaging\2024\2024-12-02 - HeLa-Residence time - 3 replica\500 ms\Halo-TRF2 - RAP1-KO c32\Nuclear masks\';
pixel_size = 0.1533; %effective pixel size BT - 0.1083
plot_tracks = 0; %1-plot results or 0-skip plotting the results (much faster!)

%% Output folder
%specify the path and name of the folder where you want to save your data
output_pathname = 'G:\Research\Imaging\2024\2024-12-02 - HeLa-Residence time - 3 replica\500 ms\Halo-TRF2 - RAP1-KO c32\Tracking\Nuclear tracks 1500 ms\';
if ~exist(output_pathname,'dir') %if folder does not exist
    mkdir(output_pathname); %create a folder
end

%% ----- Processing ----- %%
tracking_data = dir([input_pathname_tracks,'*.mat']);
masks = dir([input_pathname_masks,'*.tiff']);
%check that number of tracking files match the number of masks
if length(tracking_data) ~= length(masks)
    error("Number of tracking files does not match the number of masks.");
end

for i = 1:length(tracking_data) %for every file
    disp(['Processing ', num2str(i) '/' num2str(length(tracking_data)) ' file.']);
    mask = im2double(imread([input_pathname_masks masks(i).name])); %load mask and convert it to double
    tracks = load([input_pathname_tracks tracking_data(i).name]);  %load respective tracks
    tracks_nuclear = struct();
    tracks_cytoplasmic = struct();
    nuc_counter = 0;
    cyt_counter = 0;
    for j = 1:length(tracks.trackedPar)
        localizations = round(tracks.trackedPar(j).xy ./ pixel_size); %convert x,y localizations to pixel positions
        if sum(sum(mask(localizations(:,2), localizations(:,1)))) >= 1 %if any of the localization is within the nuclear mask, label this track as nuclear
            nuc_counter = nuc_counter + 1;
            tracks_nuclear(nuc_counter).xy = tracks.trackedPar(j).xy;
            tracks_nuclear(nuc_counter).Frame = tracks.trackedPar(j).Frame;
            tracks_nuclear(nuc_counter).TimeStamp = tracks.trackedPar(j).TimeStamp;
        else
            cyt_counter = cyt_counter + 1;
            tracks_cytoplasmic(cyt_counter).xy = tracks.trackedPar(j).xy; %otherwise, label it as cytoplasmic
            tracks_cytoplasmic(cyt_counter).Frame = tracks.trackedPar(j).Frame;
            tracks_cytoplasmic(cyt_counter).TimeStamp = tracks.trackedPar(j).TimeStamp;
        end
    end
    settings = tracks.settings; %imaging settings
    trackedPar = tracks_nuclear;
    save([output_pathname [tracking_data(i).name(1:end-4) '_nuclear.mat']], 'trackedPar', "settings"); %save nuclear tracks
    trackedPar = tracks_cytoplasmic;
    save([output_pathname [tracking_data(i).name(1:end-4) '_cytoplasmic.mat']], 'trackedPar', "settings"); %save cytoplasmic tracks

%% ----- Plot ----- %%
    if plot_tracks == 1
        figure(i)
        subplot(121)
        if length(tracks_nuclear) > 1 %if there is at least one track
            hold on
            for j = 1:length(tracks_nuclear)
                x = tracks_nuclear(j).xy(:,1);
                y = tracks_nuclear(j).xy(:,2);
                plot(x, y, 'LineWidth', 1.5);
            end
            axis([0 settings.Width*pixel_size 0 settings.Height*pixel_size]);
            title('Nuclear tracks');
            hold off
        end
        subplot(122)
        if length(tracks_cytoplasmic) > 1 %if there is at least one track
            hold on
            for j = 1:length(tracks_cytoplasmic)
                x = tracks_cytoplasmic(j).xy(:,1);
                y = tracks_cytoplasmic(j).xy(:,2);
                plot(x, y, 'LineWidth', 1.5);
            end
            axis([0 settings.Width*pixel_size 0 settings.Height*pixel_size]);
            title('Cytoplasmic tracks');
            hold off
        end
    end
end

disp('Track sorting completed.')