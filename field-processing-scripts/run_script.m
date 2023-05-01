function run_script(site_name)
if nargin <1
    site_name = 'unknown_site';
end
close all;
addpath 'C:\Users\apres\Desktop\ApRES\scripts\processing\field-processing-scripts'
addpath 'C:\Users\apres\Desktop\ApRES\scripts\processing\field-processing-scripts\nicholls_utils'
logfile_name = strcat('log-',site_name,'-',datestr(now,'mm-dd-yy-HH-MM'),'.txt');
diary(logfile_name);
% Select SD card path
%SD_card_directory = uigetdir(); % Can manually type the path as well below
SD_card_directory = 'D:\';      % this is the path to the SD card, assuming no other hard drive has been plugged in previously.
disp(['Site ',site_name,' , data in: ', SD_card_directory])
disp(' ')
dat_info = dir(strcat(SD_card_directory,'/**/*.DAT')); 
figure(1);
plot(datetime([dat_info.datenum],'ConvertFrom','datenum'),[dat_info.bytes]./1e6,'.');
xlabel('date of measurement');
ylabel('file size (Mbytes)');
total_data_str = ['Total Data Size: ',num2str(sum([dat_info.bytes])/1e9,3),' Gbytes'];
title(total_data_str);
disp(['Total .DAT File Count: ',num2str(size(dat_info,1))])
disp(total_data_str);
% Identify Anomalous Filesizes
outlier_indices = find([dat_info.bytes] >305e6 | [dat_info.bytes] < 295e6);
disp(' ')
disp('Irregular file sizes on:');
for i=outlier_indices
    disp(dat_info(i).date)
end
disp(' ')
disp('Quick plots from across season')
ok_indices = setdiff(1:length(dat_info),outlier_indices);

%% Plotting bursts and histograms and profiles from first, middle, last day
to_plot_indices = [ok_indices(1), ok_indices(ceil(end/2)), ok_indices(end)];
depthset = 1500;
pad = 2;
win = @blackman;
for i=to_plot_indices
    filename = strcat(dat_info(i).folder,'/',dat_info(i).name);
    disp(['Plotting file:' dat_info(i).name])
    vdat = Field_load(filename,1); 
    if vdat.Code == -4 % burst not found in file
            
            disp('Burst not found in file');

    elseif vdat.Code == -5
            
            disp(['No chirp starts found in file ' dat_info(i).name]);
        
    else %data is good
        
         % Split burst into various attenuator settings
        
         vdats = Field_burst_split_by_att(vdat);

        %% Plot

        vdat = Field_burst_mean(vdats(1));
        chirpname = dat_info(i).name;% int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB '];
         [tax,hax,aax] = open_plot(vdat,chirpname);
            
            for j = 1:length(vdats)
                vdat = Field_burst_mean(vdats(j));
                axes(tax), hold on;
                ht = plot(vdat.t,vdat.vif,'DisplayName',[int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ']); % signal
                if any(round(vdat.vif,2)==0) || any(round(vdat.vif,2)==2.5)
                    disp(['!!!!!!! The signal might be clipped in ' chirpname '!!!!!!!!']); 
                end

                axes(hax), hold on;
                hh = histogram(vdat.vif,'Orientation','horizontal','DisplayName',[int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ']);

                % phase process data
                [rc,~,~,su] = fmcw_range(vdat,pad,depthset,win);
                axes(aax), hold on
                plot(rc,20*log10(abs(su)),'DisplayName',[int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ']);
                
            end
            axes(tax);
            legend;
            axes(hax);
            legend;
            plot(linspace(1,max(hh.Values())+100,size(vdat.t,2)),repmat(0,size(vdat.t)),'r','HandleVisibility','off'); 
            plot(linspace(1,max(hh.Values())+100,size(vdat.t,2)),repmat(2.5,size(vdat.t)),'r','HandleVisibility','off');
            xlim([0 max(hh.Values())+100])
            axes(aax);
            xlim([0 depthset]);
            legend;
            
           
            
        
    end
end
end

function [tax,hax,aax] = open_plot(vdat,chirpname)
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
    title(t,chirpname);

end
