%% Check dates and/or times of measurements
%
% Run this script on the SD card folder to identify potential data gaps
% George Lu, July 2022

%% User parameters
% Change this directory to the SD card (has all the dated subfolders inside)
myFolder = '/Users/georgelu/Downloads/S30_201808/';

% Determine if you want to check every burst (option 0 - every 15 minutes) 
% or date of every file, which takes multiple bursts (option 1 - approximately 1 day). 
% The first takes longer whereas the second will only help identify
% much larger gaps
resolution = 0; % 0=fine, 1 = coarse

%% Rest of script
% Check if the folder exists
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

% Get list of subfolders with the automated tests
subfolderID = fullfile(myFolder,'DIR*');
subfolders = dir(subfolderID);

% Initialize array of dates and times
dates = [];
times = [];
tic
% Iterate through each subfolder
for i=1:length(subfolders)
    subfolder = subfolders(i).name;
    filePattern = fullfile(strcat(myFolder,subfolder),'*.DAT');
    % Get .DAT files in each folder
    fileList = dir(filePattern);
    % Iterate through all the files in each subfolder
    for j=1:length(fileList)
        filename = fileList(j).name;
        disp(strcat('Opening file: ',filename))
        if resolution==1 % coarse resolution - only plot dates
            % Clean up names so we get datetimes
            date_string = erase(filename,'DATA');
            date_string = erase(date_string,'.DAT');
            datetime_val = datetime(date_string,'InputFormat','yyyy-MM-dd-HHmm');
            dates = vertcat(dates,datetime_val);
        elseif resolution == 0 % fine resolution - plot dates vs time
            filecontents = fileread(strcat(myFolder,subfolder,'/',filename));
            % Use regexp to find the datetimes in a file
            expression = 'Time stamp=....-..-.. ..:..:..';
            datetimes = regexp(filecontents,expression,"match");
            for k=1:length(datetimes)
                date_string = erase(datetimes(k),'Time stamp=');
                datetime_val = datetime(date_string,'InputFormat','yyyy-MM-dd HH:mm:ss');
                time = timeofday(datetime_val);
                date = datetime(datetime_val,'Format','dd-MMM-yyyy');
                times = vertcat(times,time);
                dates = vertcat(dates,date);
            end
        end
    end
end
toc

%% Plotting
if resolution == 1
    plot(dates,0,'o');
    xlabel("Date");title("Date of each .DAT file created");set(gca,'YTick',[]);
elseif resolution == 0
    plot(dates,times,'o');
    xlabel("Date");ylabel("Time of Day");title("Time and date of each burst");
end