%% average_bursts
% Trying to iterate through a directory and get a moving window averaging
% bursts
% INCOMPLETE/WIP
% George Lu, July 2022

% First let's just try averaging three bursts
% Change this directory to the SD card (has all the dated subfolders inside)
myFolder = '/Users/georgelu/Documents/ApRES burn-in tests/';

% Specify attenuation setting: 1 or 2
atten_setting = 1;

% Specify number of bursts you want to average
to_avg = 3;

% Get list of subfolders with the automated tests
subfolderID = fullfile(myFolder,'DIR*');
subfolders = dir(subfolderID);
% initialize counter of bursts acccessed
% initialize dummy vdat to store sum
vdat_temp.count = 0;
vdat_array
% Iterate through each subfolder
for i=1:length(subfolders)
    subfolder = subfolders(i).name;
    filePattern = fullfile(strcat(myFolder,subfolder),'*.DAT');
    % Get .DAT files in each folder
    fileList = dir(filePattern);
    % Iterate through all the files in each subfolder
    for j=1:length(fileList)
        filename = fileList(j).name;
        
        % iterate through bursts
        for k = 1:100 % arbitrary max of 100
            vdat = Field_load(strcat(myFolder,subfolder,'/',filename),k);
            if vdat.Code == -4 % burst not found in file
                break;
            elseif vdat.Code == -5
                break;
            else
                vdats = Field_burst_split_by_att(vdat);
                vdat_temp.count = vdat_temp.count + 1;
                vdat = Field_burst_mean(vdats(atten_setting));
                if vdat_temp.count == 3;
                end

            end
            
        end

    end
end