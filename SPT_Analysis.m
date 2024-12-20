% SPT_Analysis.m Version: 1.0 Copyright (C) Jens Schmidt
clear all; close all; clc;
%%%% Overview %%%%
% This is the master single particle tracking analysis script of the
% Schmidt lab at Michigan State University it has various feature for the
% analysis of single particle tracking data generated using the MTT
% algorithm implemented in the SPT code published with the Spot On paper by
% Anders Hansen. Single particle tracking data in the trackedPar structure
% variable are necessary for this code to operate correctly. There are a
% number of different analysis modalities that can be implemented which
% will be described in detail below
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% DEFINE INPUT AND OUTPUT PATHS %%%%%%%%%%%%%%%%%%%%%%%%
% specify input path with .tracked files:
input_path=('/Users/admin/Documents/MATLAB/Raw tracks/');
output_path=('/Users/admin/Documents/MATLAB/Raw tracks/');
%%%%% make output folder if it does not already exist
if exist(output_path,'dir')
    % OK, output_path exists and is a directory (== 7).
    disp('The given output folder exists. MATLAB workspaces will be save to:');
    disp(output_path);
else
    mkdir(output_path);
    disp('The given output folder did not exist, but was just created. MATLAB workspaces will be save to:');
    disp(output_path);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Specify cores for parallel processing %%%%%
% workers = 10;
% parpool(workers)

%%%%% Specify Experiment properties %%%%%
pixelsize = 0.1533 ;%pixel size in µm
timedelay = 1.5 ;%time per imaging frame
minlength = 3 ;%minimal track length to analyze

%%%%% Specify Analyses to run SET to 1 to run %%%%%
plot_tracks = 0;%generate a plot of all trajectories, with track inforamtion saved in UserData
plot_tracksCC = 0 ;%generate a plot of all trajectories color coded by Diffusion Coefficient, with track inforamtion saved in UserData
Coloc1 = 0; % Will filter tracks based on proximity to a marker file order has to match track file order
Coloc2 = 0; % Will filter tracks based on proximity to two marker file order has to match track file order
ColocDUAL = 0; %If set to one will do frame by frame colocalizations require complete tracking of marker
PlotColoc = 0; %If set to 1 it will plot all tracks L1 tracks and NOT L1 tracks
MSD = 0; % If set to one will calculate D for each track and generate
DVS = 0; % Will generate a plot of the distance to the nearest marker and the step size of the trajectory
DHIST = 0; %Will generate a histogram of Ds
RESI = 1; %Will generate a survival plot of the track length for all tracks and separated by location
INTENSITY = 0; %Will generate a histogram of localization intensities for all particles longer than minlength
MEANSTEP = 0; %Will calculate the mean steps size of all tracks in file
LOCcount = 0; %Will calculate the number of localizations in each trackedPar file

%%%%% IF RUNNING COLOCALIZATION SPECIFY MARKER TRACK FILE LOCATION AND NAME %%%%%
coloc1_path=('E:\Research\Imaging\2024\2024-04-17 - HeLa-Halo-shelterin-SMimaging-sparse-dense\Halo-POT1 c16\Telomere trajectories\');
coloc2_path=('...');
ColocDIST1 = 0.48;%Distance in µm below which a track is considered colocalized with a marker 1
ColocDIST2 = 0.48;%Distance in µm below which a track is considered colocalized with a marker 2

%%%%% Choose DataSet %%%%%
input=input_path(1:end-1); % This code adds all the file names of the tracks in the input folder to a structure to cycle through
mat_files=dir([input_path,'*.mat']);
data_struct = struct([]);
data_struct(1).path = input_path;
data_struct(1).Include = [];
for iter = 1:length(mat_files)
    data_struct(1).workspaces{iter} = mat_files(iter).name;
    data_struct(1).Include = [data_struct(1).Include iter];
end

%%%%% Load Marker Track Files %%%%%
if Coloc1 == 1 % This code adds all the file names of the tracks of Marker 1 in the Coloc1 folder to a structure to cycle through
    input=coloc1_path(1:end-1);
    mat_files=dir([coloc1_path,'*.mat']);
    data_struct(2).path = coloc1_path;
    data_struct(2).Include = [];
    for iter = 1:length(mat_files)
        data_struct(2).workspaces{iter} = mat_files(iter).name;
        data_struct(2).Include = [data_struct(2).Include iter];
    end
end

if Coloc2 == 1 % This code adds all the file names of the tracks of Marker 1 and Marker 2 in the Coloc1 folder to a structure to cycle through
    input=coloc1_path(1:end-1);
    mat_files=dir([coloc1_path,'*.mat']);
    data_struct(2).path = coloc1_path;
    data_struct(2).Include = [];
    for iter = 1:length(mat_files)
        data_struct(2).workspaces{iter} = mat_files(iter).name;
        data_struct(2).Include = [data_struct(2).Include iter];
    end
    input=coloc2_path(1:end-1);
    mat_files=dir([coloc2_path,'*.mat']);
    data_struct(3).path = coloc2_path;
    data_struct(3).Include = [];
    for iter = 1:length(mat_files)
        data_struct(3).workspaces{iter} = mat_files(iter).name;
        data_struct(3).Include = [data_struct(3).Include iter];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Initiate Track Structure %%%%%
%This structure will contain all trackedPar track structures for each input
%file this is done for the tracks and marker files
TRACK_struct = struct([]);
TRACKL1_struct = struct([]);
TRACKL2_struct = struct([]);
for i = 1: length(data_struct(1).workspaces)
    file = data_struct(1).workspaces{i};
    filepath = ([input_path,file]);
    open(filepath);
    data1 = ans.trackedPar;
    name = data_struct(1).workspaces{i};
    TRACK_struct(i).trackedPar = data1;
    TRACK_struct(i).name = name;
end

if Coloc1 == 1
for i = 1: length(data_struct(2).workspaces)
    fileLoc1 = data_struct(2).workspaces{i};
    filepathLoc1 = ([coloc1_path,fileLoc1]);
    open(filepathLoc1);
    dataLoc1 = ans.trackedPar;
    name = data_struct(2).workspaces{i};
    TRACKL1_struct(i).trackedPar = dataLoc1;
    TRACKL1_struct(i).name = name;
end
end

if Coloc2 == 1
for i = 1: length(data_struct(3).workspaces)
    fileLoc2 = data_struct(3).workspaces{i};
    filepathLoc2 = ([coloc2_path,fileLoc2]);
    open(filepathLoc2);
    dataLoc2 = ans.trackedPar;
    name = data_struct(3).workspaces{i};
    TRACKL2_struct(i).trackedPar = dataLoc2;
    TRACKL2_struct(i).name = name;
end
end
%%%%% END Initiate Track Structure %%%%%

%%%%%%%%%%% MODULE CALCULATE MSD %%%%%%%%%%%
if MSD == 1
    disp('CALCULATING MSD D');
    SPACE_UNITS = 'um';
    TIME_UNITS = 's';
    for i = 1: length(TRACK_struct)
        %%% LOAD TRACKS %%%
        All = {};
        MSDidx=[];
        tracks = TRACK_struct(i).trackedPar;
        for k = 1: length(tracks)
            if length(tracks(k).xy)>= minlength
                temp = zeros(length(tracks(k).Frame),3);
                temp(:,1)= tracks(k).TimeStamp;
                temp(:,2)= tracks(k).xy(:,1);
                temp(:,3)= tracks(k).xy(:,2);
                All = [All, temp];
                MSDidx=[MSDidx,k];
            end
        end
        ma_all = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
        ma_all = ma_all.addAll(All);
        ma_all.labelPlotTracks;
        ma_all = ma_all.computeMSD;
        ma_all = ma_all.fitMSD(0.8);
        for j = 1:length(MSDidx)
            idx = MSDidx(j);
            TRACK_struct(i).trackedPar(idx).D = ma_all.lfit.a(j)/4;
        end
    end
end
%%%%%%%%%%% END MODULE CALCULATE MSD %%%%%%%%%%%

%%%%%%%%%%% MODULE Marker Colocalization %%%%%%%%%%%
% This code runs through each localization of a track and measures the
% distance to marker localizations and adds two new elements to the
% structure the coordinates of the closest marker position (data.neigb) for each frame
% of the track and the [MinimalDistance Stepsize] (data.distVstep). If
% colocDUAL is being used the marker is a track file and each localization
% is compared to the marker locations in the correscponding frame of the
% movie. The data will be saved as a new variable with the same structure
% as TrackedPar with the additional elemtents described above.
if Coloc1 == 1 & ColocDUAL == 0
    disp('RUNNING COLOCALIZE WITH SINGLE MARKER COORDINATES');
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        loc1 = TRACKL1_struct(i).trackedPar;
        TRACK_struct(i).trackedParL1 = distandstepSINGLE(minlength,tracks,loc1);
    end
end

if Coloc1 == 1 & ColocDUAL == 1
    disp('RUNNING COLOCALIZE WITH SINGLE MARKER TRACKS');
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        loc1 = TRACKL1_struct(i).trackedPar;
        TRACK_struct(i).trackedParL1 = distandstepDUAL(minlength,tracks,loc1);
    end
end

%Sorting based on the distance from the marker
if Coloc1 == 1
for i = 1: length(TRACK_struct)
    counterL1 = 1;
    counterNL1 = 1;
    trackedParL1 = struct;
    trackedParNL1 = struct;
    name = TRACK_struct(i).name;
    for j = 1:length(TRACK_struct(i).trackedParL1)
        len = length(TRACK_struct(i).trackedParL1(j).xy);
        if isempty(TRACK_struct(i).trackedParL1(j).distVstep) == 1
            trackedParNL1(counterNL1).xy = TRACK_struct(i).trackedParL1(j).xy;
            trackedParNL1(counterNL1).Frame = TRACK_struct(i).trackedParL1(j).Frame;
            trackedParNL1(counterNL1).TimeStamp = TRACK_struct(i).trackedParL1(j).TimeStamp;
            %trackedParNL1(counterNL1).Intensity = TRACK_struct(i).trackedParL1(j).Intensity;
            trackedParNL1(counterNL1).distVstep = TRACK_struct(i).trackedParL1(j).distVstep;
            trackedParNL1(counterNL1).neighb = TRACK_struct(i).trackedParL1(j).neighb;
            counterNL1 = counterNL1 + 1;
        elseif isempty(TRACK_struct(i).trackedParL1(j).distVstep) ~= 1 & len > minlength &  min(TRACK_struct(i).trackedParL1(j).distVstep(:,1)) < ColocDIST1
            trackedParL1(counterL1).xy = TRACK_struct(i).trackedParL1(j).xy;
            trackedParL1(counterL1).Frame = TRACK_struct(i).trackedParL1(j).Frame;
            trackedParL1(counterL1).TimeStamp = TRACK_struct(i).trackedParL1(j).TimeStamp;
            %trackedParL1(counterL1).Intensity = TRACK_struct(i).trackedParL1(j).Intensity;
            trackedParL1(counterL1).distVstep = TRACK_struct(i).trackedParL1(j).distVstep;
            trackedParL1(counterL1).neighb = TRACK_struct(i).trackedParL1(j).neighb;
            counterL1 = counterL1 + 1;
        elseif isempty(TRACK_struct(i).trackedParL1(j).distVstep) ~= 1 & len > minlength &  min(TRACK_struct(i).trackedParL1(j).distVstep(:,1)) > ColocDIST1
            trackedParNL1(counterNL1).xy = TRACK_struct(i).trackedParL1(j).xy;
            trackedParNL1(counterNL1).Frame = TRACK_struct(i).trackedParL1(j).Frame;
            trackedParNL1(counterNL1).TimeStamp = TRACK_struct(i).trackedParL1(j).TimeStamp;
            %trackedParNL1(counterNL1).Intensity = TRACK_struct(i).trackedParL1(j).Intensity;
            trackedParNL1(counterNL1).distVstep = TRACK_struct(i).trackedParL1(j).distVstep;
            trackedParNL1(counterNL1).neighb = TRACK_struct(i).trackedParL1(j).neighb;
            counterNL1 = counterNL1 + 1;
        end
    end
        trackedPar = trackedParL1;
        
        save([output_path, name(1:end-4), '_L1.mat'], 'trackedPar');
        trackedPar = trackedParNL1;
        name = TRACKL1_struct(i).name;
        save([output_path, name(1:end-4), '_NL1.mat'], 'trackedPar');
    if PlotColoc == 1
        disp('RUNNING PLOT COLOC TRACKS');
        tracks = TRACK_struct(i).trackedPar;
        name = TRACK_struct(i).name;
        h=figure('Name',name);
        subplot(2,2,1) 
        hold on
        if isfield(tracks,'xy')
            tracktotal = length(tracks);
            for j = 1:tracktotal %cycles through all particle tracks
    
            temp = size(tracks(j).xy); %obtain number of tracks for a given particle
            temp=temp(1,1);
            if temp(1)>= minlength
                plot(tracks(j).xy(:,1),tracks(j).xy(:,2),'-x');
                A=gca;
                B=A.Children(1);
                B.UserData.xy=tracks(j).xy;
                B.UserData.timestamp=tracks(j).TimeStamp;
                B.UserData.Frame=tracks(j).Frame;
            end
            end
        end
        a = gca;
        xaxis = a.XLim;
        yaxis = a.YLim;
        hold off
        subplot(2,2,2) % plot L1
        a = gca;
        a.XLim = xaxis;
        a.YLim = yaxis;
        hold on
        tracks = trackedParL1;
        if isfield(tracks,'xy')
            tracktotal = length(tracks);
            for j = 1:tracktotal %cycles through all particle tracks
    
            temp = size(tracks(j).xy); %obtain number of tracks for a given particle
            temp=temp(1,1);
            if temp(1)>= minlength
                plot(tracks(j).xy(:,1),tracks(j).xy(:,2),'-x');
                A=gca;
                B=A.Children(1);
                B.UserData.xy=tracks(j).xy;
                B.UserData.timestamp=tracks(j).TimeStamp;
                B.UserData.Frame=tracks(j).Frame;
            end
     
            end
        end
        hold off
        subplot(2,2,3) % plot L1
        a = gca;
        a.XLim = xaxis;
        a.YLim = yaxis;
        hold on
        tracks = trackedParNL1;
        if isfield(tracks,'xy')
            tracktotal = length(tracks);
            for j = 1:tracktotal %cycles through all particle tracks
    
            temp = size(tracks(j).xy); %obtain number of tracks for a given particle
            temp=temp(1,1);
            if temp(1)>= minlength
                plot(tracks(j).xy(:,1),tracks(j).xy(:,2),'-x');
                A=gca;
                B=A.Children(1);
                B.UserData.xy=tracks(j).xy;
                B.UserData.timestamp=tracks(j).TimeStamp;
                B.UserData.Frame=tracks(j).Frame;
            end
     
            end
        end
        hold off
        subplot(2,2,4) % plot L1
        a = gca;
        a.XLim = xaxis;
        a.YLim = yaxis;
        hold on
        tracks = TRACKL1_struct(i).trackedPar;
        if isfield(tracks,'xy')
            tracktotal = length(tracks);
            for j = 1:tracktotal %cycles through all particle tracks
    
            temp = size(tracks(j).xy); %obtain number of tracks for a given particle
            temp=temp(1,1);
            if temp(1)>= minlength
                plot(tracks(j).xy(:,1),tracks(j).xy(:,2),'-x');
                A=gca;
                B=A.Children(1);
                B.UserData.xy=tracks(j).xy;
                B.UserData.timestamp=tracks(j).TimeStamp;
                B.UserData.Frame=tracks(j).Frame;
            end
     
            end
        end
        hold off
    end
end
end

if Coloc2 == 1 & ColocDUAL == 0
    disp('RUNNING COLOCALIZE WITH A TWO MARKERS WITH A SINGLE COORDINATES');
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        loc1 = TRACKL1_struct(i).trackedPar;
        loc2 = TRACKL2_struct(i).trackedPar;
        TRACK_struct(i).trackedParL2 = distandstepSINGLE(minlength,tracks,loc1,loc2);
    end
end

if Coloc2 == 1 & ColocDUAL == 1
    disp('RUNNING COLOCALIZE WITH A TWO MARKERS WITH A SINGLE COORDINATES');
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        loc1 = TRACKL1_struct(i).trackedPar;
        loc2 = TRACKL2_struct(i).trackedPar;
        TRACK_struct(i).trackedParL2 = distandstepDUAL(minlength,tracks,loc1,loc2);
    end
end



%%%%%%%%%%%  END MODULE Marker Colocalization %%%%%%%%%%%

%%%%%%%%%%% MODULE PLOT TRACKS %%%%%%%%%%%
%Plots trajectories and saves a figure for each trackedPar file. This can
%be used to check the quality of the trajectories generated
if plot_tracks == 1
    disp('RUNNING PLOT TRACKS');
    for i = 1: length(TRACK_struct)
    tracks = TRACK_struct(i).trackedPar;
    name = TRACK_struct(i).name;
    if isfield(tracks,'xy')
    h=figure('Name',name(1:end-4));
    tracktotal = length(tracks)
    hold on
    
    for j = 1:tracktotal %cycles through all particle tracks
    
        temp = size(tracks(j).xy); %obtain number of tracks for a given particle
        temp=temp(1,1);
        if temp(1)>= minlength
            plot(tracks(j).xy(:,1),tracks(j).xy(:,2),'-x');
            A=gca;
            B=A.Children(1);
            B.UserData.xy=tracks(j).xy;
            B.UserData.timestamp=tracks(j).TimeStamp;
            B.UserData.Frame=tracks(j).Frame;
        end
     
    end
   
    A=gca;
    A.UserData.secperframe = timedelay;
    hold off
    savefig(h,[output_path, name(1:end-4),'.fig']);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% MODULE PLOT TRACKS CC %%%%%%%%%%%
% Plots tracks color coded by D to determine of there are any spatial
% differences in the diffusion coefficient
if plot_tracksCC == 1
    disp('RUNNING PLOT TRACKS CC');
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        name = TRACK_struct(i).name;
        h=figure('Name',name(1:end-4));
        tracktotal = length(tracks);
        hold on
  
    for j = 1:tracktotal %cycles through all particle tracks
        temp = size(tracks(j).xy); %obtain number of tracks for a given particle
        temp=temp(1,1);
        if temp(1)>= minlength
               if tracks(j).D < 0.01
               C = [0 0 0.5000];
               elseif tracks(j).D >= 0.01 & tracks(j).D < 0.05
               C = [0 0.1667 1];
               elseif tracks(j).D >= 0.05 & tracks(j).D < 0.5
               C = [0 0.833 1];
               elseif tracks(j).D >= 0.5 & tracks(j).D < 1
               C = [0.5 1 0.5];
               elseif tracks(j).D >= 1 & tracks(j).D < 2
               C = [1 0.833 0];
               elseif tracks(j).D >= 2 & tracks(j).D < 3
               C = [1 0.1667 0];
               elseif  tracks(j).D >= 3
               C = [0.5 0 0];
               end
            plot(tracks(j).xy(:,1),tracks(j).xy(:,2),'-x','Color',C);
            A=gca;
            B=A.Children(1);
            B.UserData.xy=tracks(j).xy;
            B.UserData.timestamp=tracks(j).TimeStamp;
            B.UserData.Frame=tracks(j).Frame;
        end
    end
    A=gca;
    A.UserData.secperframe = timedelay;
    hold off
    %savefig(h,[output_path, name(1:end-4)]);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% MODULE PLOT DVS %%%%%%%%%%%
% Generates a distance from marker vs step size plot. Binding events should
% have a small stepsize and short distance from the respective marker
if DVS == 1
    disp('RUNNING PLOT DVS');
    All = [];
    for i = 1: length(TRACK_struct)
        data = TRACK_struct(i).trackedPar;
        name = TRACK_struct(i).name;
        if isfield(data,'distVstep')
            for j = 1:tracktotal %cycles through all particle tracks
                All = [All ; data(j).distVstep];
            end
        end
    end
    Allunique = unique(All,'rows','stable');
    XData = All(:,1);
    YData = All(:,2);
    YData=YData+0.0001;
    var1=log10(YData);
    var2=log10(XData);
    h = figure ();
    hold on
    scatplot(var1,var2,'circles');
    xlim([-2.5 0])
    ylim([-1.5 1.0])
    yline(log10(10^0.5),'--','LineWidth',1.0);
    yline(log10(1),'--','LineWidth',1.0);
    yline(log10(10^(-0.5)),'--','LineWidth',1.0);
    yline(log10(0.1),'--','LineWidth',1.0);
    xline(log10(10^0.5),'--','LineWidth',1.0);
    xline(log10(1),'--','LineWidth',1.0);
    xline(log10(10^(-0.5)),'--','LineWidth',1.0);
    xline(log10(0.1),'--','LineWidth',1.0);
    xline(log10(10^(-1.5)),'--','LineWidth',1.0);
    xline(log10(0.01),'--','LineWidth',1.0);
    xlabel('Step-size (log10, \mum)', 'FontSize', 14,'FontWeight','bold')
    ylabel('Distance from Telomere (log10, \mum)', 'FontSize', 14,'FontWeight','bold')
    set(gca, 'FontSize', 12, 'YMinorTick','on','XMinorTick','on')
    hold off
    savefig(h,[output_path, name(1:end-4), '_dvsLOC1.fig'])   
end
%%%%%%%%%%% END MODULE PLOT DVS %%%%%%%%%%%

%%%%%%%%%%% MODULE DHIST %%%%%%%%%%%
%Generates a histogram of the diffusion coefficients for all tracks and
%tracks co-localized with  markers
if DHIST == 1
    disp('RUNNING PLOT DHIST');
    D_all = [];
    D_loc1 = [];
    D_loc2 = [];
    D_NOloc = [];
    name = TRACK_struct(1).name;
    if Coloc1 == 0 & Coloc2 == 0
    for i = 1: length(TRACK_struct)
        for j = 1:length(TRACK_struct(i).trackedPar)
            if TRACK_struct(i).trackedPar(j).D > 0
                D_all = [D_all; TRACK_struct(i).trackedPar(j).D];
            end
        end
    end
    end
    if Coloc1 == 1
        for i = 1: length(TRACK_struct)
            for j = 1:length(TRACK_struct(i).trackedPar)
                if TRACK_struct(i).trackedPar(j).D > 0
                    D_all = [D_all; TRACK_struct(i).trackedPar(j).D];
                    dists = TRACK_struct(i).trackedParL1(j).distVstep(:,1);
                    if min(dists) < ColocDIST1
                        D_loc1 = [D_loc1, TRACK_struct(i).trackedPar(j).D];
                    else
                        D_NOloc = [D_NOloc, TRACK_struct(i).trackedPar(j).D];
                    end
                end
            end
        end
    end
    if Coloc2 == 1
        for i = 1: length(TRACK_struct)
            for j = 1:length(TRACK_struct(i).trackedPar)
                if TRACK_struct(i).trackedPar(j).D > 0
                    D_all = [D_all; TRACK_struct(i).trackedPar(j).D];
                    dists = TRACK_struct(i).trackedParL1(j).distVstep(:,1);
                    dists2 = TRACK_struct(i).trackedParL2(j).distVstep2(:,1);
                    if min(dists2) < ColocDIST2
                        D_loc2 = [D_loc2, TRACK_struct(i).trackedPar(j).D];
                    elseif min(dists) < ColocDIST1
                        D_loc1 = [D_loc1, TRACK_struct(i).trackedPar(j).D];
                    else
                        D_NOloc = [D_NOloc, TRACK_struct(i).trackedPar(j).D];
                    end
                end
            end
        end
    end
    f = figure('Name',name);
    % hold on
    % subplot(3,1,1)
    h = histogram(log10(D_all), 'BinWidth',0.1);
    h.Normalization = 'probability';
    % subplot(3,1,2)
    % histogram(log10(D_loc1),'BinWidth',0.1)
    % subplot(3,1,3)
    % histogram(log10(D_NOloc),'BinWidth',0.1)
    % hold off
    savefig(f,[output_path, name(1:end-4),'_DHIST.fig']);
    if Coloc2 == 1
        f = figure('Name',name);
        hold on
        subplot(4,1,1)
        histogram(log10(D_all),'BinWidth',0.2)
        subplot(4,1,2)
        histogram(log10(D_loc1),'BinWidth',0.2)
        subplot(4,1,3)
        histogram(log10(D_loc2),'BinWidth',0.2)
        subplot(4,1,4)
        histogram(log10(D_NOloc),'BinWidth',0.2)
        hold off
        savefig(f,[output_path, name(1:end-4),'_DHIST2.fig']);
    end
end
%%%%%%%%%%% MODULE DHIST %%%%%%%%%%%

%%%%%%%%%%% RESIDENCE TIME MODULE %%%%%%%%%%%%%
%Plots a cumulative density function of track lengths. If there is no
%co-localization with markers it will simply plot the length of all
%trajectories longer than the minimum track length defined in the settings
%above. If co-localization is used it will plot the CDF for tracks
%co-localized with each marker and the CDF for tracks that do not
%co-localize with a marker. This analysis can be used with tracks generated
%for static particles by constraining the D during tracking to limit
%tracking to static particles.
if RESI == 1
    disp('RUNNING RESIDENCE TIME SURVIVAL PLOT');
    ALLresi = [];
    LOC1resi = [];
    LOC2resi = [];
    NOLOCresi = [];
    name = TRACK_struct(1).name;
    if Coloc1 == 0 & Coloc2 == 0
        for i = 1: length(TRACK_struct)
            for j = 1:length(TRACK_struct(i).trackedPar) %cycles through all particle tracks
                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength
                    resi = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    ALLresi = [ALLresi; resi];
                end
            end
        end
        u=figure('Name',[name,'_CDF']);
        cdfGraph = cdfplot(ALLresi);
        h = findobj(gca,'Type','line');
        x = get(h,'Xdata');
        y = get(h,'Ydata');
        yy = 1-y;
        r = figure(2);
        plot(x,yy, 'LineWidth', 2);
        set(gca, 'YScale', 'log');
        xlim([0 200]);
        ylim([0.01 1]);
        xlabel('Time (s)');
        ylabel('1-CDF (%)');
        title('Survival function');
        savefig(r,[output_path, name(1:end-4), ' - survival.fig']);
    end
    if Coloc1 == 1
        for i = 1: length(TRACK_struct)
            for j = 1:length(TRACK_struct(i).trackedPar) %cycles through all particle tracks

                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength
                    dists = TRACK_struct(i).trackedParL1(j).distVstep(:,1);
                end

                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength
                    resi = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    ALLresi = [ALLresi; resi];
                end
                %check if any coordinate of the trajectory is below the
                %coloralziation cutoff for marker 1 add residence time to
                %colocalized structure if not add it to not colocalized
                %structure
                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength & min(dists) <= ColocDIST1;
                    residence = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    LOC1resi = [LOC1resi residence];
                else 
                    resi2 = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    NOLOCresi = [NOLOCresi resi2];
                end
            end
        end
        u=figure('Name',[name,'_CDF']);
        hold on
        cdfplot(LOC1resi);
        cdfplot(NOLOCresi);
        hold off
        savefig(u,[output_path, name(1:end-4), 'res.fig']);
    end
    if Coloc1 == 2
        for i = 1: length(TRACK_struct)
            for j = 1:length(TRACK_struct(i).trackedPar) %cycles through all particle tracks

                if length(TRACK_struct(i).trackedPar(j).xy)> minlength
                    dists = TRACK_struct(i).trackedParL2(j).distVstep2(:,1);
                end

                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength
                    resi = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    ALLresi = [ALLresi; resi];
                end
                %check if any coordinate of the trajectory is below the
                %coloralziation cutoff for marker 1 add residence time to
                %colocalized structure
                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength & min(dists) <= ColocDIST1;
                    residence = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    LOC1resi = [LOC1resi residence];
                end
                %check if any coordinate of the trajectory is below the
                %coloralziation cutoff for marker 2 add residence time to
                %colocalized structure for marker 2
                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength & min(dists) <= ColocDIST2;
                    residence2 = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    LOC2resi = [LOC1resi residence2];
                end
                %Check if coordinates are greater than the cutoff away
                %from both markers add no NOLOC structure
                if length(TRACK_struct(i).trackedPar(j).xy)>= minlength && min(dists) > ColocDIST1 & min(dists) > ColocDIST2;
                    resi2 = TRACK_struct(i).trackedPar(j).TimeStamp(length(TRACK_struct(i).trackedPar(j).TimeStamp))-TRACK_struct(i).trackedPar(j).TimeStamp(1);
                    NOLOCresi = [NOLOCresi resi2];
                end
            end
        end
        u=figure('Name',[name,'_CDF'])
        hold on
        cdfplot(LOC1resi);
        cdfplot(LOC2resi);
        cdfplot(NOLOCresi);
        hold off
        savefig(u,[output_path, name(1:end-4), 'res.fig']);
   end
end

%%%%%%%%%%% END RESIDENCE TIME MODULE %%%%%%%%%%%%%

%%%%%%%%%%% MODULE INTENSITY %%%%%%%%%%%
% % Generates Intensity histrogram for all localizations
if INTENSITY == 1
    disp('RUNNING INTENSITY');
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        name = TRACK_struct(i).name;
        tracktotal = length(tracks);
        intensities = [];
        for i = 1:tracktotal
            if length(tracks(i).Intensity) >= minlength
                temp = tracks(i).Intensity;
                intensities = [intensities;temp];
            end
        end
    end
    histogram(intensities)
    savefig([output_path, name(1:end-4), 'INT.fig']);
end
%%%%%%%%%% END INTENSITY MODUL  %%%%%%%%%
        
%%%%%%%%%%% MODULE MEAN STEP SIZE %%%%%%%%%%%
% % Generates a scatter plot of mean step sizes for each track

if MEANSTEP == 1
    disp('RUNNING MEAN STEP SIZE');
    means = [];
    name = TRACK_struct(1).name;
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        name = TRACK_struct(i).name;
        tracktotal = length(tracks);
        stepsizes = [];
        for j = 1:tracktotal
            cords = TRACK_struct(i).trackedPar(j).xy;
            if length(cords) > minlength
                for k = 1:length(cords(:,1))-1; %Loop over all localizations for a track
                    XY = [cords(k,1) cords(k,2)];
                    XY2 = [cords(k+1,1) cords(k+1,2)];
                    step = pdist2(XY, XY2);
                    stepsizes = [stepsizes step];
                end
            end
        end
        meanstep = mean(stepsizes)
        means = [means meanstep];
    end
    u=figure('Name',[name,'_MeanSteps'])
    hold on
    scatter (1,means)
    hold off
    savefig(u,[output_path, name(1:end-4), 'MeanSteps.fig']);
end

%%%%%%%%%% END MEAN STEP SIZE  %%%%%%%%%

%%%%%%%%%%% LOCcount STEP SIZE %%%%%%%%%%%
% % Generates a structure with the file name and the total number of
% localizations in each file

if LOCcount == 1
    disp('RUNNING Localization Count');
    loccountOUT = struct();
    for i = 1: length(TRACK_struct)
        tracks = TRACK_struct(i).trackedPar;
        name = TRACK_struct(i).name;
        tracktotal = length(tracks);
        loccount = 0;
        for j = 1:tracktotal;
            loccount = loccount + length(tracks(j).Frame);
        end
        loccountOUT(i).Name = name;
        loccountOUT(i).locs = loccount;
    end
    save([output_path, name(1:end-4), '_LOCCOUNT.mat'], 'loccountOUT');
end

%%%%%%%%%% END LOCcount STEP SIZE  %%%%%%%%%




    












