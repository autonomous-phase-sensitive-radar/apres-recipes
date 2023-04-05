function [pct_flagged, flagged_files] = check_clipping_attenuation(myFolder,max_plots,file_spacing, burst_spacing, clipping ,amplitude,folder)
%% Check clipping/attenuation
% This script, which adapts Elizabeth Case's field processing scripts,
% flags data that potentially contains clipping. If one burst within a data
% file is clipped, the entire file is flagged, and the script moves on to
% examining the next file. This is for the sake of avoiding too many plots
% popping up. 
%
% By George Lu, July 2022

%% User parameters
% Change this directory to the SD card (has all the dated subfolders inside)
%myFolder = '/Users/georgelu/Downloads/S30_201808/';

% Set a maximum number of plots to avoid too many popping up
%max_plots = 20; 

% Set how many files and/or bursts to skip over each iteration (for efficiency)
% 1 would mean no skipping, 2 would mean every other one
%file_spacing = 2;
%burst_spacing = 2;

% Choose whether to detect clipping (1) or too much attenuation (0)
%clipping = 0; 

% Pick a max amplitude deemed to be too attenuated (clipping is at
% amplitude of 1.25)
%amplitude = 0.25; 
resolution = 0; 


% Additional settings that aren't as important
maxchirps = 100;
depthset = 1500;
pad = 2;
win = @blackman;
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
if ~exist('folder','var')
    % parameter does not exist, so default it to something
    subfolderID = fullfile(myFolder,'DIR*');
    extension = '*.DAT';
else
    % Get list of subfolders with the automated tests
    subfolderID = fullfile(myFolder,folder);
    extension = '*.dat';
end
subfolders = dir(subfolderID);
num_plots = 0;
total_flagged_files = 0;
total_files_checked = 0; 
flagged_files = [];
for i=1:length(subfolders)
    subfolder = subfolders(i).name;
    filePattern = fullfile(strcat(myFolder,subfolder),extension);
    % Get .DAT files in each folder
    fileList = dir(filePattern);
    % Iterate through all the files in each subfolder
    for j=1:file_spacing:length(fileList)
        filename = fileList(j).name;
        total_files_checked = total_files_checked + 1;
        getBurst = 1;
        % Coarse resolution scenario
        % Assumes 2 attenuator settings
        if resolution == 1 && getBurst && num_plots<max_plots
            for AttNum=1:2
                vdat = mean_burst_file(strcat(myFolder,subfolder,'/',filename),AttNum);
                if vdat.Code == -4 % burst not found in file
           
                    disp([num2str(BurstNo-1) ' burst(s) found in file'])
                    %return
                    getBurst = 0;
            
                elseif vdat.Code == -5
                    
                    disp(['No chirp starts found in file ' filename ' burst ' int2str(thisburst) ': - Corrupted data?'])
                    getBurst = 0;
                else
                    if any(round(vdat.vif,2)==0) || any(round(vdat.vif,2)==2.5) && clipping ==1
                        drawPlot = 1;
                        disp(['!!!!!!! The signal might be clipped in ' filename ' avg' int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ' '!!!!!!!!']); 
                        disp("Moving to look at next file");
                        getBurst = 0;

                    elseif max(vdat.vif)<1.25+amplitude && min(vdat.vif)>1.25-amplitude && clipping == 0
                        drawPlot = 1;
                        disp(['!!!!!!! The signal might be too attenuated in ' filename ' avg' int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ' '!!!!!!!!']); 
                        disp("Moving to look at next file");
                        getBurst = 0;
                    end
                    if drawPlot==1
                        num_plots = num_plots + 1;
                        total_flagged_files = total_flagged_files + 1;
                        if num_plots == max_plots
                            disp("Max number of plots reached. Stopping...")
                        end
                        % Figure 1a & b: % plot voltage vs. time series, histogram of the avg. chirp,
                        % and note if there is any cutoff
                        
                        [tax,hax,aax,pax] = open_plot(vdat);
                        
                        axes(tax), hold on
                        ht = plot(vdat.t,vdat.vif); % signal
        
                        axes(hax), hold on
                        hh = histogram(vdat.vif,'Orientation','horizontal');
                        plot(linspace(1,max(hh.Values())+100,size(vdat.t,2)),repmat(0,size(vdat.t)),'r','HandleVisibility','off'); 
                        plot(linspace(1,max(hh.Values())+100,size(vdat.t,2)),repmat(2.5,size(vdat.t)),'r','HandleVisibility','off');
                        xlim([0 max(hh.Values())+100])
        
                    % Figure 1c & d: plot amplitude/phase profile
                    
                        % phase process data
                        [rc,~,~,su] = fmcw_range(vdat,pad,depthset,win);
        
                        axes(pax), hold on
                        plot(rc,angle(su));
                        xlim([0 depthset])
        
                        axes(aax), hold on
                        plot(rc,20*log10(abs(su)));
                        xlim([0 depthset])   
                    end
                end

            end
        end
        % else fine resolution - loop through bursts in file
        burstlist = 1:100; 
        chirplist = 1:maxchirps;
        BurstNo = 0;

        while BurstNo<length(burstlist) && getBurst && num_plots<max_plots && resolution == 0

            %sets thisburst to current burst
            BurstNo = BurstNo + burst_spacing;
            thisburst = burstlist(BurstNo);
            disp(strcat('Opening file: ',filename,"-Burst",int2str(thisburst)));

            vdat = Field_load(strcat(myFolder,subfolder,'/',filename),thisburst); % load data
            
            if vdat.Code == -4 % burst not found in file
                
                disp([num2str(BurstNo-burst_spacing) ' burst(s) opened in file'])
                %return
                getBurst = 0;
            
            elseif vdat.Code == -5
                
                disp(['No chirp starts found in file ' filename ' burst ' int2str(thisburst) ': - Corrupted data?'])
                getBurst = 0;
            
            else %data is good
            
             % Split burst into various attenuator settings
            
                vdats = Field_burst_split_by_att(vdat);
                for AttSetNum = 1:length(vdats) % attenuator setting number
                
                    % Generate labels for each chirp to plot
                    shotNamePrefix = [strrep(filename,'_','-') ' b:' int2str(thisburst) ' c:'];
                    
                    %average chirps from burst
                    vdat = Field_burst_mean(vdats(AttSetNum));
                    drawPlot = 0; % Default dont draw plot
                    if any(round(vdat.vif,2)==0) || any(round(vdat.vif,2)==2.5) && clipping ==1
                        vdat.chirpname = [shotNamePrefix ' avg' int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB '];
                        drawPlot = 1;
                        disp(['!!!!!!! The signal might be clipped in ' vdat.chirpname '!!!!!!!!']); 
                        disp("Moving to look at next file");
                        getBurst = 0;
    
                    elseif max(vdat.vif)<1.25+amplitude && min(vdat.vif)>1.25-amplitude && clipping == 0
                        vdat.chirpname = [shotNamePrefix ' avg' int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB '];
                        drawPlot = 1;
                        disp(['!!!!!!! The signal might be too attenuated in ' vdat.chirpname '!!!!!!!!']); 
                        disp("Moving to look at next file");
                        getBurst = 0;
                        
                    end
                    if drawPlot==1
                        total_flagged_files = total_flagged_files + 1;
                        num_plots = num_plots + 1;
                        flagged_files = vertcat(flagged_files,vdat.chirpname);
                        if num_plots == max_plots
                            disp("Max number of plots reached. Stopping...")
                        end
                        % Figure 1a & b: % plot voltage vs. time series, histogram of the avg. chirp,
                        % and note if there is any cutoff
                        
                        [tax,hax,aax] = open_plot(vdat);
                        
                        axes(tax), hold on
                        ht = plot(vdat.t,vdat.vif); % signal
        
                        axes(hax), hold on
                        hh = histogram(vdat.vif,'Orientation','horizontal');
                        plot(linspace(1,max(hh.Values())+100,size(vdat.t,2)),repmat(0,size(vdat.t)),'r','HandleVisibility','off'); 
                        plot(linspace(1,max(hh.Values())+100,size(vdat.t,2)),repmat(2.5,size(vdat.t)),'r','HandleVisibility','off');
                        xlim([0 max(hh.Values())+100])
        
                    % Figure 1c & d: plot amplitude/phase profile
                    
                        % phase process data
                        [rc,~,~,su] = fmcw_range(vdat,pad,depthset,win);

                        axes(aax), hold on
                        plot(rc,20*log10(abs(su)));
                        xlim([0 depthset])        
                    end    
                   
                end
                
            end
            
        end
    end
end

pct_flagged = 100*total_flagged_files/total_files_checked; 


%% 
function [tax,hax,aax] = open_plot(vdat)
    figure('Position',[0.1557    0.1530    1.0413    0.4686]*1e3);
    t=tiledlayout(3,4);

    tax = nexttile(1,[2,2]);
    set(tax,'tag','tax');
    title('Voltage')
    hold on
    box on
    xlabel('Time (s)')
    ylabel('Voltage')
    ylim([-0.25 2.75])
    plot(vdat.t,repmat(0,size(vdat.t)),'r','HandleVisibility','off'); % ADC saturation level
    plot(vdat.t,repmat(2.5,size(vdat.t)),'r','HandleVisibility','off'); % ADC saturation level
    
    hax = nexttile(3,[2,2]);
    set(hax,'tag','hax');
    title('Histogram of voltage')
    xlabel('Count')
    ylabel('Voltage')
    hold on
    box on
    ylim([-0.25 2.75])
    

    % Amp subplot
    aax = nexttile(9,[1,4]);
    set(aax,'tag','aax')
    title('Amplitude (dB)')
    box on
    xlabel('Range (m)');
    ylabel('amplitude (dB Vrms)')

    title(t,vdat.chirpname);
    drawnow
