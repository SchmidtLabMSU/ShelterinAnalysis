%% ----------------------- PHOTOBLEACHIN STEP -----------------------------
% Script for extracting photobleaching step size from photobleaching
% traces. Input is excel sheet where each column represents photobleaching
% trace and second excel sheet that represents background signal for
% respective movie.
% 
% 
% Created by: Tomas Janovic, 2024
%--------------------------------------------------------------------------

clear all; close all; clc;

%% --------------- Loading files & Input ---------------

trace_id = 1; % which trace (column) from .xlsx file you want to analyze
background_id = 1; % which trace (column) from .xlsx file is the background trace

exposure_time = 0.04; % exposure time in seconds

% Loading a .xlsx file
[filename, pathname] = uigetfile( ...
{  '*.xlsx';}, ...
   'Pick a photobleaching Halo traces .xlsx file', ...
   'MultiSelect', 'off');

if isequal(filename, 0) %if no file is selected
    error('No files are loaded.');
end

data = readtable([pathname filename]); %load photobleaching Halo-TRF1 traces

% Loading a .xlsx file
[filename, pathname] = uigetfile( ...
{  '*.xlsx';}, ...
   'Pick a photobleaching background traces .xlsx file', ...
   'MultiSelect', 'off');

if isequal(filename, 0) %if no file is selected
    error('No files are loaded.');
end

background = readtable([pathname filename]); %load photobleaching background traces

%% --------------- Processing ---------------
background_bleach  = table2array(background(:,background_id)); % average background photobleaching
data = table2array(data);
data = data(:,trace_id);

% Data median filtration
data = medfilt1(data, 3);
background_bleach = medfilt1(background_bleach, 3);
data = data - background_bleach;

% % extend background
% for i = 1:50
%     data(end:end+50,1) = data(end-50:end,1);
%     i = i+1;
% end

n_orig = length(data(:,1));
data_trim = data(data <= 5000); % trimming data
n_new = length(data_trim(:,1));

nBins = round(500 * (max(data_trim)/5000));
h = histogram(data_trim, nBins); % calculate histogram with n bins

counts = h.Values; % get values from histogram
edges = h.BinEdges; % get intensity values for histogram counts

% create variable x for intensity values of histogram counts
x = [];
for i = 1:length(edges)-1
    x = [x (edges(i+1)+edges(i))/2];
end

[maxCounts, maxCountsPosition] = max(counts); % get position of maximum counts

mygauss = @(x,xdata) x(1)*exp(-((xdata-x(2)).^2/(2*x(3).^2))); % Gauss function x(1) = amplitude (a), x(2) = mean (u), x(3) = standard deviation (sigma)

% initial guess of gauss function
amplitude = maxCounts;
meanValue = x(maxCountsPosition);
sigma = 20;
gaussCurve = mygauss([amplitude meanValue sigma], x);

% Optimization for gauss function fitting to histogram
sig = Inf; % initial criterial value
s = 0.5:0.5:25; % vector of possible sigma values for optimization
u = meanValue-25:0.5:meanValue+25; %u = meanValue-50:0.5:meanValue+50;
a = min(counts)+1:1:amplitude;
for i = 1:length(s)
    for j = 1:length(u)
        for k = 1:length(a)
            tempFun = mygauss([a(k) u(j) s(i)], x); % temporary gauss function
            distValue = sum((counts - tempFun).^2) / length(x); % what is the criterial distance
            if distValue < sig % if it is lower
               sig = distValue; % overwrite it
               opt_s1 = s(i); % and store this sigma value
               opt_u1 = u(j);
               opt_a1 = a(k);
            end
        end
    end
end
curveFit = mygauss([opt_a1 opt_u1 opt_s1], x);


temp_x = x;
temp_x((temp_x < (meanValue + 4*opt_s1)) | (temp_x > (meanValue + 200))) = 0;
temp_x(temp_x > 0) = 1;
temp_counts = counts.*temp_x;
[maxCountsPeak2, maxCountsPositionPeak2] = max(temp_counts); % get position of maximum counts for a second peak

% initial guess of gauss function
amplitude2 = maxCountsPeak2;
meanValue2 = x(maxCountsPositionPeak2(1));
sigma2 = 20;
gaussCurve2 = mygauss([amplitude2 meanValue2 sigma2], x);

% Optimization for gauss function fitting to histogram
sig = Inf; % initial criterial value
s = 0.5:0.5:70; % vector of possible sigma values for optimization
u = meanValue2-25:0.5:meanValue2+25; %u = meanValue2-150:0.5:meanValue2+100;
a = min(temp_counts)+1:1:amplitude2;
for i = 1:length(s)
    for j = 1:length(u)
        for k = 1:length(a)
            tempFun = mygauss([a(k) u(j) s(i)], x); % temporary gauss function
            distValue = sum((temp_counts - tempFun).^2) / length(x); % what is the criterial distance 
            if distValue < sig % if it is lower
               sig = distValue; % overwrite it
               opt_s2 = s(i); % and store this sigma value
               opt_u2 = u(j); % and mean value
               opt_a2 = a(k); % and amplitude value
            end
        end
    end
end
curveFit2 = mygauss([opt_a2 opt_u2 opt_s2], x);

singleBleachingStep = opt_u2 - opt_u1; % single photobleaching step
disp(['Single photobleaching step: ' num2str(singleBleachingStep)]);


%% --------------- Plotting ---------------
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');

time_axis = linspace((n_orig - n_new) * exposure_time, n_orig * exposure_time, n_new); % x - time axis
x_resampled = resample(x, 2, 1);
curveFit2 = resample(curveFit2, 2, 1);

figure(1) % plotting the phobleaching curve
plot(time_axis, data_trim, 'LineWidth', 1.5)
xlim([((n_orig - n_new) * exposure_time) ((n_orig * exposure_time) + 1)]);
ylabel('Fluorescence intensity (A.U.)');
xlabel('Time (s)');
% xlim([0 200]);
% ylim([0 2500]);
set(gca, 'FontSize', 17, 'FontWeight', 'bold', 'YMinorTick', 'off', 'XMinorTick', 'off', 'box', 'on', LineWidth = 1.5);
% handaxes2 = axes('position', [0.52 0.52 0.3 0.3]);
% plot(time_axis, data_trim, 'LineWidth', 1.5)
% xlim([80 160]);
% ylim([50 550]);
% set(gca, 'FontSize', 17, 'FontWeight', 'bold', 'YMinorTick', 'off', 'XMinorTick', 'off', 'box', 'on', LineWidth = 1.5);

figure(2) % plotting the histogram
bar(x, counts, 0.8)
hold on
plot(x, curveFit, 'LineWidth', 3.5)
hold on
plot(x_resampled, curveFit2, 'LineWidth', 3.5)
hold off
xlim([0 500]);
%ylim([0 750]);
ylabel('Counts');
xlabel('Fluorescence intensity (A.U.)');
set(gca, 'FontSize', 17, 'FontWeight', 'bold', 'YMinorTick', 'off', 'XMinorTick', 'off', 'box', 'on', LineWidth = 1.5);