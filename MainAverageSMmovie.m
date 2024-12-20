%% ------------------------ MOVIE AVERAGING -------------------------------
% Script that will average user-defined consecutive frames of a selected
% movies (*.tiff).
% 
% 
% Created by: Tomas Janovic, 2024
%--------------------------------------------------------------------------

close all; clear all; clc;

%% Input
framesToAverage = 3; %how many consecutive frames to average.

%% Loading files
data = loadTiffMovies(); %load single-molecule tiff movies

%% Saving folder
pathname = [data{1,1} 'Processed movies/']; %pathname
if ~exist(pathname,'dir') %if folder does not exist
    mkdir(pathname); %create a folder
end

%% Processing
disp('Loading finished.')

n = length(data(1,:)); %number of files loaded
for i = 1:n
    movie = data{3,i};
    averagedMovie = []; % allocation of averaged movie
    z = 1;
    for j = 1:framesToAverage:length(movie(1,1,:))-framesToAverage+1
        tempImg = movie(:,:,j:j+framesToAverage-1);
        tempImg = sum(tempImg,3)./framesToAverage;
        averagedMovie(:,:,z) = tempImg;
        z = z + 1;
    end
    %Saving movie
    options.overwrite = true;
    saveastiff(uint16(averagedMovie), [pathname data{2, i}], options);
end