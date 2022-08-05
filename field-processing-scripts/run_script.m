function run_script(site_name,mode)
%% Main script to run all the checkscripts
% Make parameters optional
if ~exist('site_name','var')
     % parameter does not exist, so default it to something
      site_name = 'unspecified-site';
end
if ~exist('mode','var')
     % parameter does not exist, so default it to something
      mode = 0; % Default is unattended mode
elseif mode == 1
    subfolder = 'Survey*';
end
logfile_name = strcat('log-',site_name,'-',datestr(now,'mm-dd-yy-HH-MM'),'.txt');
diary(logfile_name)

%% close figures
close all

% add helper scripts to the path
addpath('nicholls_utils')

% Select SD card path
SD_card_directory = uigetdir(); % Can manually type the path as well
%SD_card_directory = 'D:\';      % this is the path to the SD card, assuming no other hard drive has been plugged in previously.
myFolder = strcat(SD_card_directory,'/');

disp(['Loading data from ',  SD_card_directory]);  % so that it goes in the log file

%% Data gap check
resolution = 1; % Default 1 - coarse resolution. 0 is for fine resolution
expected_gap_hours = 4;
percent_error = 0.05;
%if mode == 0
    [gap_starts,gap_durations] = check_dates_vs_time(myFolder, resolution, expected_gap_hours, percent_error);
%else
%    [gap_starts,gap_durations] = check_dates_vs_time(myFolder, resolution, expected_gap_hours, percent_error, subfolder);
%end
%% Clipping and attenuation check
% Set a maximum number of plots to avoid too many popping up
max_plots = 10;

% Set how many files and/or bursts to skip over each iteration (for efficiency)
% 1 would mean no skipping, 2 would mean every other one
file_spacing = 2;
burst_spacing = 2; % Recommend 20 for in-field 

% Choose whether to detect clipping (1) or too much attenuation (0)
clipping = 1;

% Pick a max amplitude deemed to be too attenuated (clipping is at
% amplitude of 1.25)
amplitude = 0.25;

% clipping
if mode == 0
    [pct_clipped,flagged_clipped] = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude);
else
    [pct_clipped,flagged_clipped] = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude,subfolder);
end
% Switch to attenuation
clipping = 0;
if mode == 0
    [pct_attenuated, flagged_attenuated] = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude);
else
    [pct_attenuated, flagged_attenuated] = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude,subfolder);
end
disp(strcat('Percentage of checked files flagged for clipping: ',int2str(pct_clipped),'%'));
disp(strcat('Percentage of checked files flagged for overattenuation: ',int2str(pct_attenuated),'%'));

%% Vertical velocity checker 
if mode == 0
    [c1,c2] = check_vertical_velocity(myFolder);
else
    [c1,c2] = check_vertical_velocity(myFolder, subfolder);
end

%% Generate summary report
disp([newline,newline,'%%%%% SUMMARY REPORT %%%%%',newline]);
disp('Potential Data Gaps (Duration, start time): ');
if isempty(gap_starts) %== 0     I commented out the ==0 because I think it makes it the wrong way round, i.e. if isempty is not true that means its not empty, which is not what we want. 
    data_message = 'No data gaps identified.';
    disp(data_message);
else
    for i=1:length(gap_starts)
        disp([' ', num2str(gap_durations(i)), ' hours:   ',datestr(gap_starts(i))])
    end
    data_message = ['See ', logfile_name, ' for list of ', int2str(length(gap_starts)), ' flagged data gaps.'];  
end

disp([newline,'Files flagged for clipping'])
if size(flagged_clipped,1) == 0
    disp('No files flagged for clipping');
    clipping_message = 'No clipping detected.';
else
    for i = 1:size(flagged_clipped,1)
        disp([' ', flagged_clipped(i,:)])
    end
    clipping_message = ['See ', logfile_name, ' for list of ', int2str(size(flagged_clipped,1)), ' flagged clipped files (', int2str(pct_clipped),'%).'];
end


disp([newline,'Files flagged for overattenuation'])
if size(flagged_attenuated,1) == 0
    disp('No files flagged for overattenuation');
    attenuation_message = 'No overattenuation detected.';
else
    for i = 1:size(flagged_attenuated,1)
        disp([' ', flagged_attenuated(i,:)])
    end
    attenuation_message = ['See ', logfile_name, ' for list of ', int2str(size(flagged_attenuated,1)), ' flagged overattenuated files: (', int2str(pct_attenuated),'%).'];
 end

disp(newline)
disp(strcat('Percentage of checked files flagged for clipping: ',int2str(pct_clipped),'%'));
disp(strcat('Percentage of checked files flagged for overattenuation: ',int2str(pct_attenuated),'%'));

disp([newline,'Strain rates:'])
vv_msg_1 = 'Estimated strain rates from the first and last files are:';
disp(vv_msg_1);
vv_msg_2 = [num2str(c1(1),'%.3f'),' /day and ',num2str(c2(1),'%.3f'), ' /day, respectively'];
disp(vv_msg_2)

diary off;

%% Display final message
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msg = msgbox(['\fontsize{20}', data_message,newline, clipping_message,newline, attenuation_message, newline, vv_msg_1, newline, vv_msg_2],"Results", CreateStruct);