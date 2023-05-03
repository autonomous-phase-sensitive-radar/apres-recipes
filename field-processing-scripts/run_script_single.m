function run_script_single(numplots)
    if nargin==0
        numplots = 1;
    end
    close all;
    [file,path] = uigetfile('*.DAT');
    filename = fullfile(path, file);
    addpath 'C:\Users\apres\Desktop\ApRES\scripts\processing\field-processing-scripts'
    addpath 'C:\Users\apres\Desktop\ApRES\scripts\processing\field-processing-scripts\nicholls_utils'
    logfile_name = strcat('log-',file,'-',datestr(now,'mm-dd-yy-HH-MM'),'.txt');
    %diary(logfile_name);
    filecontents = fileread(filename);
    % Use regexp to find the datetimes in a file
    expression = 'Time stamp=....-..-.. ..:..:..';
    datetimes = regexp(filecontents,expression,"match");
    max_bursts = length(datetimes);
    disp(['Total bursts: ' num2str(max_bursts)])
    if numplots == 1
        file_spaced = [1];
    else
        file_spaced = floor(linspace(1,max_bursts,numplots));
    end
    depthset = 1500;
    pad = 2;
    win = @blackman;
    for i=file_spaced
        disp(['Plotting file: ' char(file) ', ' char(datetimes(i))])
        vdat = Field_load(filename,i); 
        if vdat.Code == -4 % burst not found in file
                
                disp('Burst not found in file');
    
        elseif vdat.Code == -5
                
                disp(['No chirp starts found in file ' file]);
            
        else %data is good
            
             % Split burst into various attenuator settings
            
             vdats = Field_burst_split_by_att(vdat);
    
            %% Plot
    
            vdat = Field_burst_mean(vdats(1));
            chirpname = [char(file) ', ' char(datetimes(i))];% int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB '];
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
