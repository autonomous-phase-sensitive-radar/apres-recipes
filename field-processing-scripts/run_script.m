%% Main script to run all the checkscripts
diary check_log_file.txt
% Select SD card path
SD_card_directory = uigetdir(); % Can manually type the path as well
myFolder = strcat(SD_card_directory,'/');
%% Data gap check
resolution = 1; % Default 1 - coarse resolution. 0 is for fine resolution
check_dates_vs_time(myFolder, resolution);

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
pct_clipped = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude);

% Switch to attenuation
clipping = 0;
pct_attenuated = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude);

disp(strcat('Percentage of checked bursts clipped:',int2str(pct_clipped),'%'));
disp(strcat('Percentage of checked bursts with high attenuation:',int2str(pct_attenuated),'%'));

%% Strain rate checker 
check_strain_rate(myFolder);
diary off;

