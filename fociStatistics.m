function [fociNr, fociStats, fociIntensities] = fociStatistics(data)
%--------------------------------------------------------------------------
% Function that will get the statistics from analyzed foci that are
% ready to save.
% 
% 
% Created by: Tomas Janovic, 2024
%--------------------------------------------------------------------------

fociNr = {'Filename', 'Cell Nr.', 'Number of foci'};
fociStats = {'Filename', 'Cell Nr.', 'Median foci intensity'};
fociIntensities = {'Filename', 'CellNr.', 'Foci Nr.', 'Foci intensity'};

for i = 1:length(data)
    filename = data{1,i};
    values = data{2,i};
    fociCount = values{1,1}';
    intensities = values{2,1}';

    for j = 1:length(fociCount)
        fociNr{end+1,1} = filename;
        fociNr{end,2} = j;
        fociNr{end,3} = fociCount(j);
    end

    tempArray = [];
    for j = 1:length(fociCount)
        tempArray = [tempArray ones.*(1:fociCount(j))];
    end

    for j = 1:length(intensities)
        fociIntensities{end+1, 1} = filename;
        fociIntensities{end,2} = intensities(j,1);
        fociIntensities{end,3} = tempArray(j);
        fociIntensities{end,4} = intensities(j,2);
    end
end


tempInt = cell2mat(fociIntensities(2, 4));
for i = 3:length(fociIntensities)
    if isequal(fociIntensities{i,2}, fociIntensities{i-1,2}) && (fociIntensities{i,3} > fociIntensities{i-1,3})
        tempInt = [tempInt cell2mat(fociIntensities(i,4))];
    else
        fociStats{end+1,1} = fociIntensities{i-1,1};
        fociStats{end,2} = fociIntensities{i-1,2};
        fociStats{end,3} = median(tempInt);
        tempInt = [];
    end
end
fociStats{end+1,1} = fociIntensities{i-1,1};
fociStats{end,2} = fociIntensities{i-1,2};
fociStats{end,3} = median(tempInt);
end