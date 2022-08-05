function [c1,c2] = check_vertical_velocity(myFolder,folder)
%% check_vertical_velocity
% Plots vertical velocities between two bursts
%
% George Lu, July 2022

global cfg

% Detect first and last file
% Main folder
%myFolder = '/Users/georgelu/Downloads/S30_201808/';

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

% Get list of subfolders with the unattended tests
if ~exist('folder','var')
    % parameter does not exist, so default it to something
    subfolderID = fullfile(myFolder,'DIR*');
else
    subfolderID = fullfile(myFolder,folder);
end
subfolders = dir(subfolderID);
% First file
first_folder = subfolders(1).name;
filePattern = fullfile(strcat(myFolder,first_folder),'*.DAT');
fileList = dir(filePattern);
first_file = fileList(1).name;
filename1 = strcat(myFolder,first_folder,'/',first_file); % can also modify to specific file (include full path)

% Last file
last_folder = subfolders(end).name;
filePattern = fullfile(strcat(myFolder,last_folder),'*.DAT');
fileList = dir(filePattern);
last_file = fileList(end).name;
filename2 = strcat(myFolder,last_folder,'/',last_file); % can also modify to specific file (include full path)

% The commands below do the actual calculation of vertical velocity. 
% You can mimic this format to do custom calculations
% First file - between first 2 bursts
[range,dh,dhe,dt,c1] = fmcw_melt(filename1,filename1,1,2,0);

% Last file - 2 consecutive bursts
[range,dh,dhe,dt,c2] = fmcw_melt(filename2,filename2,1,2,0);

% Long comparison between first and last file
%[range,dh,dhe,dt] = fmcw_melt(filename1,filename2,1,1,0);
%[range,dh,dhe,dt] = fmcw_melt(filename1,filename2,2,2,0);

%% Functions
function [range,dh,dhe,dt,c] = fmcw_melt(filename1,filename2,b1,b2,avg)
%% Load processing settings
global cfg
cfg = fmcw_process_config_vsr; % Load default configuration
% if nargin >2
%     % overwrite defaults
%     fieldnames = fields(cfgi);
%     fieldlist = [];
%     for ii = 1:length(fieldnames);
%         thisfield = fieldnames{ii};
%         cfg = setfield(cfg,thisfield,getfield(cfgi,thisfield));
%         fieldlist = [fieldlist ', ' thisfield];
%     end
%     fieldlist = fieldlist(3:end);
%     cfg.notes = [cfg.notes '. User modified fields: ' fieldlist '.'];
% end

%% Define Constants
daysPerYear = 365.25;

%% Select and load data
if nargin == 0
    if cfg.useTestFiles
        % test data
        filename1 = cfg.filename1;
        filename2 = cfg.filename2;
    else
        [file,path] = uigetfile({'*.dat;*.DAT;*.000;*.mat','Radar files: .dat, .DAT, .000 and .mat'},'Choose first (or both) radar files to calculate strain rate','multiselect','on');
        if isa(file,'double') % no files chosen
            return
        end
        if ischar(file) % then we've only selected one file as second is in a different dir. repeat
            filename1 = [path file];
            [file,path] = uigetfile({'*.dat;*.DAT;*.000;*.mat','Radar files: .dat, .DAT, .000 and .mat'},'Choose second radar files to calculate strain rate');
            filename2 = [path file];
        else
            filename1 = [path file{1}];
            filename2 = [path file{2}];
        end
    end
end
[path1,name1,ext1] = fileparts(filename1);
[path2,name2,ext2] = fileparts(filename2);

% Load
if avg == 0
    f = Field_load(filename1,b1); 
    g = Field_load(filename2,b2);
elseif avg == 1
    f = mean_burst_file(filename1,1);
    g = mean_burst_file(filename2,1);
end
%% Print Output
Disp(' ')
Disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
Disp(['Radar processor ' mfilename])
Disp(' ')
Disp('Using files:')
Disp('Filename                                  Date/Time')
Disp([name1 '   ' datestr(f.TimeStamp)])
Disp([name2 '   ' datestr(g.TimeStamp)])
dt = g.TimeStamp-f.TimeStamp; % time difference days
Disp(['Time interval: ' num2str(dt) ' days'])
Disp(' ')
%Disp('Processing settings:')
%Disp(cfg)
%Disp(' ')
% Display comand line to repeat this test
% Disp('To repeat this processing paste the following command from the clipboard')
% cmd = ['fmcw_melt(''' name1 ext1 ''',''' name2 ext2 ''')'];
% clipboard('copy',cmd)
% Disp(cmd)


%% Clean shots
% Chirps
if cfg.doManualChirpSelect
    Disp('Using chirp subset defined in config')
    f = fmcw_burst_subset(f,cfg.fchirplist);
    g = fmcw_burst_subset(g,cfg.gchirplist);
end

% Frequency range
if cfg.fRange(1) > 2e8 || cfg.fRange(2) < 4e8
    % Cull f.specCor range
    Disp(['Cropping frequency range to: ' num2str(cfg.fRange(1)/1e6) '-' num2str(cfg.fRange(2)/1e6) ' MHz'])
    Disp(' ')
    f = fmcw_cull_freq(f,cfg.fRange);
    g = fmcw_cull_freq(g,cfg.fRange);
end

% Keep only one attenuator setting (otherwise phase standard deviation and
% error estimate is wrong. (as attenuators have different phase delays)
attSetList = unique(f.chirpAtt,'stable');
if numel(attSetList)>1
    disp(['Warning: multiple attenuator settings found file ' name1 '. Keeping first set only.'])
    fs = fmcw_burst_split_by_att(f);
    f = fs(1); % keeing results from first attenuator setting only.
end
attSetList = unique(g.chirpAtt,'stable');
if numel(attSetList)>1
    disp(['Warning: multiple attenuator settings found file ' name2 '. Keeping first set only.'])
    gs = fmcw_burst_split_by_att(g);
    g = gs(1);
end

if cfg.doClean
    % Cull gross contaminated chirps
    [f,nBad] = fmcw_cull_bad(f,cfg.noisePowerLimit);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' name1])
    [g,nBad] = fmcw_cull_bad(g,cfg.noisePowerLimit);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' name2])
    
    % Cull noisey remaining
    Disp(['Culling noisest ' int2str(cfg.nNoisest) ' shots'])
    f = fmcw_cull_noisey(f,cfg.nNoisest);
    g = fmcw_cull_noisey(g,cfg.nNoisest);
end

%% Phase sensitive processing
% Average all shots in burst and phase process this average shot
fm = fmcw_burst_mean(f);
%[f.rangeCoarse,f.rangeFine,f.specCor,f.specRaw] = fmcw_range(fm,cfg.p,cfg.maxRange,cfg.winFun); % note - should this be changed to a weighted mean for cases where the attenuation changes within the burst
[f.rangeCoarse,f.rangeFine,f.specCor,f.specRaw] = fmcw_range(fm,cfg.p,cfg.maxRange,cfg.winFun); % note - should this be changed to a weighted mean for cases where the attenuation changes within the burst
gm = fmcw_burst_mean(g);
%[g.rangeCoarse,g.rangeFine,g.specCor,g.specRaw] = fmcw_range(gm,cfg.p,cfg.maxRange,cfg.winFun);
[g.rangeCoarse,g.rangeFine,g.specCor,g.specRaw] = fmcw_range(gm,cfg.p,cfg.maxRange,cfg.winFun);
if f.rangeCoarse~=g.rangeCoarse
    error('ranges are different')
end
dr = mean(diff(f.rangeCoarse));

%% ALIGN BULK: Bulk Align upper internals (to keep the search window smaller for the fine scale correlation)
if cfg.doBulkAllignment
    % Allign internal layers to account for changes in cable length surface
    % accumulation and firn compaction using xcor. Note this only offsets to
    % the closest integer range bin.
    Disp(['Co-registering profiles using amplitude cross-correlation'])
    Disp(['> depth range: ' mat2str(cfg.bulkAlignRange)])
    Disp(['> max offset : ' mat2str(cfg.maxOffsetM)])
    fi = find((f.rangeCoarse>=min(cfg.bulkAlignRange) & f.rangeCoarse<max(cfg.bulkAlignRange))); % depth bins to use (f)
    maxlag = ceil(cfg.maxOffsetM/dr); % max bin lags
    [~,AB.ampCor,~,AB.lags] = fmcw_xcorr(f.specCor,g.specCor,fi,maxlag);
    [AB.maxCor,ii] = max(AB.ampCor); % get index (mci) of best amplitude correlation
    AB.n = AB.lags(ii); % n is the number of steps a2 should be shifed right to match a1
    Disp(['correlation =  ' num2str(AB.maxCor) ])
    if AB.maxCor < 0.8
        Disp(['Warning: poor correlation in bulk internal allignment - check files'])
    end
    % Apply the offset to shot 2 to make this match shot 1
    if AB.n==0
        Disp('Internals match - no offset required')
        Disp(' ')
    else
        Disp(['Shifting profile 2, ' int2str(AB.n) ' steps left to align internals. (' num2str(AB.n*dr) 'm)'])
        Disp(' ')
        g.specRawUnshifted = g.specCor; % keep a copy of g.specCor for plotting
        g.specCorUnshifted = g.specCor; % keep a copy of g.specCor for plotting
        g.specRaw = circshift(g.specRaw,[0 -AB.n]); % lagg offset
        g.specCor = circshift(g.specCor,[0 -AB.n]); % lagg offset
        g.rangeFine = circshift(g.rangeFine,[0 -AB.n]); % lagg offset
        
        if cfg.doPlotAlignBulk || cfg.doPlotAll
            % plot before and after offset
            figure
            plot(f.rangeCoarse,dB(abs(f.specCor)),'r');
            hold on
            plot(g.rangeCoarse,dB(abs(g.specCorUnshifted)),'b');
            plot(g.rangeCoarse,dB(abs(g.specCor)),'c');
            legend('shot1','shot2','shot2 shifted')
            title(['Profile bulk co-registration, depth range: ' mat2str(cfg.bulkAlignRange) ' m'])
            set(gca,'xlim',cfg.bulkAlignRange)
%            set(gcf,'pos',[232 554 560 420])
        end
    end
else
    AB = nan;
end

%% Find bed and select max depth to match
switch cfg.maxDepthMethod
    case 'auto'
        % Use auto detected bed
        Disp(['Searching for bed using method: ' cfg.bedMethod])
        f.bn = fmcw_findbed(f.rangeCoarse,abs(f.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
        f.bedDepth = f.rangeCoarse(f.bn) + f.rangeFine(f.bn);
        g.bn = fmcw_findbed(g.rangeCoarse,abs(g.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
        g.bedDepth = g.rangeCoarse(g.bn) + g.rangeFine(g.bn);
        fit.maxDepth = min([f.bedDepth g.bedDepth]) - cfg.bedBuffer;
        Disp(['Shot 1: bed found at ' num2str(f.bedDepth) ' m'])
        Disp(['Shot 2: bed found at ' num2str(g.bedDepth) ' m'])
        %bed.dh = g.bedDepth - f.bedDepth;
        %Disp(['Range difference: ' num2str(bed.dh) ' m'])
        
    case 'config'
        % Use maxDepth from config
        fit.maxDepth = cfg.maxDepthConfig;
        Disp(['Using user defined max fit depth ' num2str(cfg.maxDepthConfig)])
        
    case 'manual'
        % Manually define max depth (graphical)

        % Plot Amplitudes vs range
        figure
        s1ahl = plot(f.rangeCoarse,20*log10(abs(f.specCor)),'r');
        hold on
        s2ahl = plot(g.rangeCoarse,20*log10(abs(g.specCor)),'b');
        ylabel('Vrms (dB)')
        xlabel('range (m)')
        legend([s1ahl(1) s2ahl(1)],{name1,name2},'Location','SouthWest','interpreter','none')
        
        % Graphically select range
        [fit.maxDepth,~] = ginput(1);
        y = get(gca,'ylim');
        plot([fit.maxDepth fit.maxDepth],y,'k');
        set(gca,'xlim',[0 1.5*fit.maxDepth]);
        h = fillXRange([cfg.firnDepth fit.maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
        pause(1)
        close
        Disp(['Using user defined max fit depth ' num2str(fit.maxDepth)])
end

%% ALIGN COARSE: Co-register profile segments (amplitude xcorr)
% Cross correlate internals gi segments to get vertical shift as a function of range
% xcor in big chunks to get minimum chance of wrapping errors.
AC.maxOffset = cfg.maxStrain*(fit.maxDepth-cfg.minDepth); %
AC.stepSizeM = 5; % cfg.coarseChunkWidth/2; %cfg.coarseChunkWidth/2;
binStart = [cfg.minDepth:AC.stepSizeM:fit.maxDepth-cfg.coarseChunkWidth]; % measure offset over a wider range to plot
[AC.range,AC.dh] = deal(zeros(size(binStart)));
for ii = 1:numel(binStart);
    depthRange = [binStart(ii) binStart(ii)+cfg.coarseChunkWidth];
    fi = find((f.rangeCoarse>=min(depthRange) & f.rangeCoarse<max(depthRange))); % depth bins to use (f)
    maxlag = ceil(AC.maxOffset/dr); % max bin lags
    [AC.RANGEIND(ii,:),AC.AMPCOR(ii,:),~,AC.LAGS(ii,:)] = fmcw_xcorr(f.specCor,g.specCor,fi,maxlag);
    AC.RANGE(ii,:) = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,AC.RANGEIND(ii,:));
    [~,mci] = max(AC.AMPCOR(ii,:));
    AC.range(ii) = AC.RANGE(ii,mci); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
    AC.lags(ii) = AC.LAGS(ii,mci);
    AC.ampCor(ii) = AC.AMPCOR(ii,mci);
    AC.dh(ii) = dr*AC.lags(ii); % Range offset (m) between segments
    
    % Quality checks on best correlation
    % Check whether correlation is limited by maxlag
    if mci == 1 || mci == size(AC.LAGS,2) 
        AC.ampCor(ii) = 0;
    end
    % Check prominence of peak (how much better than the next match)
    [cpk,~] = findpeaks(AC.AMPCOR(ii,:),'sortstr','descend','npeaks',2);
    if isempty(cpk) % no peaks!
        AC.ampCorProm(ii) = 0;
    elseif numel(cpk)==1
        AC.ampCorProm(ii) = 1; % this is the only maximum
    else
        AC.ampCorProm(ii) = cpk(1) - cpk(2); % Absolute prominence
    end
end
AC.isGood = AC.ampCor>=cfg.minAmpCor & AC.ampCorProm>=cfg.minAmpCorProm;

% Now fit a polnomial through the lags
AC.P = polyfit(AC.range(AC.isGood),AC.lags(AC.isGood),cfg.polyOrder);

if cfg.doPlotAlignCoarse || cfg.doPlotAll
    figure
    clear ax
    ax(1) = subplottight(3,1,1);
    ii = f.rangeCoarse<fit.maxDepth*1.1;
    plotAmp(f,g,ii);
    title('Profile co-registration - amplitude cross-correlation')
    
    ax(2) = subplottight(3,1,2);
    %sh = surf(AC.RANGE',AC.LAGS',AC.AMPCOR','edgecol','none');
    sh = pcolor(AC.RANGE',AC.LAGS',AC.AMPCOR'); % ,'edgecol','none'
    shading interp
    view(0,90)
    caxis([0.9 1])
    ylabel('bin lag')
    ch = colorbar('East');
    ylabel(ch,'amplitude correlation')
    set(gca,'ydir','normal')
    hold on
    plot3(AC.range,AC.lags,ones(size(AC.range)),'w.') % all lags
    plot3(AC.range(AC.isGood),AC.lags(AC.isGood),ones(1,sum(AC.isGood)),'g.') % only good lags
    AC.lagsPoly = polyval(AC.P,f.rangeCoarse(ii)); % generate smoothed lags
    plot3(f.rangeCoarse(ii),AC.lagsPoly,ones(1,sum(ii)),'g') % poly fit to lags
    
    ax(3) = subplottight(3,1,3);
    plot(AC.range,AC.ampCor)
    ylabel('correlation')
    xlabel('depth')
    
    %set(gcf,'pos',[232 554 560 420])
    linkaxes(ax,'x')
    
    %keyboard
end

%% Estimate error
switch cfg.errorMethod
    case 'emperical'
        % Process each shot separately to get phase standard deviation at each depth
        % shot 1 (f.specCor)
        [~,~,F] = fmcw_range(f,cfg.p,cfg.maxRange,cfg.winFun);
        f.phaseStdDev = std(F)./abs(f.specCor); % phase standard deviation of the burst
        f.phaseStdError = f.phaseStdDev/sqrt(size(f.vif,1)); % phase standard error of the mean shot - using sqrt(n)
        % shot 2 (g.specCor)
        [~,~,G] = fmcw_range(g,cfg.p,cfg.maxRange,cfg.winFun);
        g.phaseStdDev = std(G)./abs(g.specCor); % phase standard deviation of the burst
        g.phaseStdError = g.phaseStdDev/sqrt(size(g.vif,1)); % phase standard error of the mean shot - using sqrt(n)
    case 'assumedNoiseFloor'
        % Error estimate by assuming a noise level and calculating phase
        % noise from this
        noiseFloor = 10.^(cfg.noiseFloordB/20);
        f.phaseStdError = noiseFloor./abs(f.specCor); % phase error shot 1
        g.phaseStdError = noiseFloor./abs(g.specCor); % phase error shot 2
end
%f.rangeError = angle(f.specCor)./((4*pi/f.lambdac) - (4*f.rangeCoarse*f.K/f.ci^2));
%g.rangeError = angle(g.specCor)./((4*pi/g.lambdac) - (4*g.rangeCoarse*g.K/g.ci^2));
f.rangeError = fmcw_phase2range(f.phaseStdError,f.lambdac,f.rangeCoarse,f.K,f.ci);
g.rangeError = fmcw_phase2range(g.phaseStdError,g.lambdac,g.rangeCoarse,g.K,g.ci);

%% ALIGN FINE: Estimate phase shift between profile segments (complex xcor)
switch cfg.phaseDiffMethod
    case 'peakDiff'
        % Directly difference phase using value at peak and one bin either
        % side for error
        [peaks,peaki] = findpeaks(abs(f.specCor).*double(f.rangeCoarse<fit.maxDepth));
        for ii = 1:length(peaks);
            % Get coarse offsets at bin centres
            if cfg.doPolySmoothCoarseOffset
                peaklag = round(polyval(AC.P,f.rangeCoarse(peaki(ii)))); % generate smoothed lags
            else
                peaklag = round(interp1(AC.range,AC.lags,f.rangeCoarse(peaki(ii)),'linear','extrap')); % offset between shots at this peak (m)
            end
%            peaklag = round(AC.dhInterp/dr); % bin lag
            jj = peaki(ii)-1:peaki(ii)+1;
            fg = f.specCor(jj).*conj(g.specCor((jj)+peaklag)); % one point either side of peak
            AF.phasedh = -angle(fg)*(f.lambdac/(4*pi)); % note: phase changes in opposite sense to range
            AF.lagdh = peaklag*dr; %
            AF.DH(ii,:) = AF.lagdh + AF.phasedh; % range change between shots % Brennan et al. eq 14 and 15;
            
            if cfg.doPlotAlignFine  || cfg.doPlotAll
                if mod(ii,23) == 0 % only plot some
                
                    figure
                    clear ax
                    ax(1) = subplottight(2,1,1);
                    wi = peaki(ii)-100:peaki(ii)+100;
                    %jj = peaki(ii)-1:peaki(ii)+1;
                    plot(f.rangeCoarse(wi),20*log10(abs(f.specCor(wi))),'r');
                    hold on
                    plot(f.rangeCoarse(jj),20*log10(abs(f.specCor(jj))),'ro');
                    plot(f.rangeCoarse(peaki(ii)),20*log10(abs(f.specCor(peaki(ii)))),'r+');
                    plot(g.rangeCoarse(wi),20*log10(abs(g.specCor(wi+peaklag))),'b');
                    plot(f.rangeCoarse(jj),20*log10(abs(g.specCor(jj+peaklag))),'bo');
                    plot(f.rangeCoarse(peaki(ii)),20*log10(abs(g.specCor(peaki(ii)+peaklag))),'b+');
                    ylabel('Vrms (dB)')
                    %xlabel('range (m)')
                    %legend({name1,name2},'Location','SouthWest','interpreter','none')
                    title('Profile phase difference - peaks only (peak 1)')
                    
                    ax(2) = subplottight(2,1,2);
                    plot(f.rangeCoarse(wi),angle(f.specRaw(wi)),'r');
                    hold on
                    plot(g.rangeCoarse(wi),angle(g.specRaw(wi+peaklag)),'b');
                    ylabel('phase (rad)')
                    linkaxes(ax,'x')
                    
                    %keyboard
                end
            end
        end
        dh  = AF.DH(:,2);
        dhe = (AF.DH(:,3)-AF.DH(:,1))/2; % average of upper and lower error...
        range = transpose(f.rangeCoarse(peaki));
        AF.coherence = ones(size(range));
        
    case 'xcor'
        % complex correlation: xcor in small chunks to get good depth resolution
        % this also gives us the AF.coherence of the segments
        stepSizeM = cfg.chunkWidth; % cfg.chunkWidth/2; cfg.chunkWidth;
        binStart = [cfg.minDepth:stepSizeM:fit.maxDepth-cfg.chunkWidth]; % measure offset over a wider range to plot
        %OffsetRange = AC.maxOffset;
        for ii = 1:numel(binStart);
            depthRange = [binStart(ii) binStart(ii) + cfg.chunkWidth];
            binDepth = mean(depthRange);
            maxlag = ceil(AC.maxOffset/dr); % max bin lags
            fi = find((f.rangeCoarse>=min(depthRange) & f.rangeCoarse<max(depthRange))); % depth bins to use (f)
            [AF.RANGEIND(ii,:),AF.AMPCOR(ii,:),AF.COR(ii,:),AF.LAGS(ii,:),AF.PE(ii,:),AF.PSE(ii,:)] = fmcw_xcorr(f.specCor,g.specCor,fi,maxlag,f.phaseStdError,g.phaseStdError,cfg.p);
            AF.RANGE(ii,:) = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,AF.RANGEIND(ii,:));
            
            if cfg.doUseCoarseOffset % Define the bin lag from the coarse correlation
                % Get coarse offsets at bin centres
                if cfg.doPolySmoothCoarseOffset
                    AC.dhInterp = dr*polyval(AC.P,binDepth); % generate smoothed lags
                else
                    AC.dhInterp = interp1(AC.range,AC.dh,binDepth,'linear','extrap'); %
                end
                [~,AF.mci(ii)] = min(abs(AC.dhInterp/dr-AF.LAGS(ii,:))); % bin lags index
            else
                [~,AF.mci(ii)] = max(AF.AMPCOR(ii,:)); % use best lag from fine cor
            end
            AF.ampCor(ii) = AF.AMPCOR(ii,AF.mci(ii));
            AF.cor(ii) = AF.COR(ii,AF.mci(ii)); % complex correlation at best amp correlation point
            % will be overwritten later...
            range(ii) = AF.RANGE(ii,AF.mci(ii)); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
            % will be overwritten later...
        end
        AF.PHASECOR = abs(AF.COR)./AF.AMPCOR;
        AF.lagvec = AF.LAGS(1,:);
        
        % Unwrap
        if cfg.doSmartUnwrap
            % redefine best lags as those closest to zero phase diff
            dphi = 2*pi*f.fc/(f.B*cfg.p); % phase change over 1 bin (brennan eq 16)
            mci_pm = AF.mci - fix(angle(AF.cor)/dphi); % lag to give smallest phase diff near good correlation point
            mci_pm(mci_pm<1) = mci_pm(mci_pm<1) + round(2*pi/dphi); % deal with edge effects
            mci_pm(mci_pm>size(AF.COR,2)) = mci_pm(mci_pm>size(AF.COR,2)) - round(2*pi/dphi); % deal with edge effects
            %mci_pmu = round(unwrap(mci_pm*dphi)/dphi); % lag with no large steps
            
            % Try following phase difference minimum
            % start from maximum correlation point(s)
%             if ~cfg.doBulkAllignment
                [~,si] = max(AF.ampCor); % index of best correlation point
                mci_pmu(si) = mci_pm(si);
%             else
%                 % the above method is suseptable to starting in the wrong phase
%                 % catchment (integer wavelength error) if binWidth is small. Better to
%                 % use the 0 lag at the centre of the bulk alignment region as
%                 % this has used a much larger xcor window
%                 [~,si] = min(abs(range-mean(cfg.bulkAlignRange))); % range bin closest to centre of bulk alignment range
%                 mci_pmu(si) = find(AF.lagvec==0);
%             end
            % Forwards to end
            for ii = si+1:length(AF.mci)
                [~,pki] = findpeaks(-abs(angle(AF.COR(ii,:)))); % index of phase difference minimum
                [~,bpki] = min(abs(mci_pmu(ii-1)-pki)); % closest minimum to last
                mci_pmu(ii) = pki(bpki);
                %mci_pmu(ii) = mci_pmu(ii-1) + round(angle(c(mci_pm(ii-1),ii))/dphi); % lag to give smallest phase diff
            end
            % Backwards to start
            for ii = si-1:-1:1
                [~,pki] = findpeaks(-abs(angle(AF.COR(ii,:)))); % index of phase difference minimum
                [~,bpki] = min(abs(mci_pmu(ii+1)-pki)); % closest minimum to last
                mci_pmu(ii) = pki(bpki);
                %mci_pmu(ii) = mci_pmu(ii-1) + round(angle(c(mci_pm(ii-1),ii))/dphi); % lag to give smallest phase diff
            end
            AF.mciu = mci_pmu;
            
            % Now check that starting at the bulk alignment centre got the
            % same interger pathlength/wavelength as most of the local amp
            % correlations.
            binLagError = mci_pm - mci_pmu;
            numWavelengthError = sum(abs(binLagError)>0);
        end
        
        % Extract values at chosen offsets
        if cfg.doSmartUnwrap
            igood = AF.mciu;
        else
            igood = AF.mci;
        end
        for ii = 1:length(igood)
            AF.cor(ii) = AF.COR(ii,igood(ii));
            range(ii) = AF.RANGE(ii,igood(ii)); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
            AF.lags(ii) = AF.LAGS(ii,igood(ii));
            AF.ampCor(ii) = AF.AMPCOR(ii,igood(ii));
            AF.coherence(ii) = AF.COR(ii,igood(ii));
            AF.phaseCor(ii) = AF.PHASECOR(ii,igood(ii));
            AF.pe(ii) = AF.PE(ii,igood(ii));
            AF.pse(ii) = AF.PSE(ii,igood(ii));
        end
        
        % Calculate the total depth shift from the integer bin lags and the phase shift
        AF.lagdh = AF.lags*dr; % 
        AF.phasedh = -angle(AF.cor)*(f.lambdac/(4*pi)); % note: phase changes in opposite sense to range
        %AF.phasedh = -angle(AF.cor)./((4*pi/f.lambdac) - (4*AF.RANGE*f.K/f.ci^2)); % this is the full equation including the term generated by the last term in (13)
        dh = AF.lagdh + AF.phasedh; % range change between shots % Brennan et al. eq 14 and 15;
        dhe = AF.pse*(f.lambdac/(4*pi)); % using the standard error of the phase difference estimated across the range bin
end

figure
clear ax
ax(1) = subplot(2,1,1);
plotAmp(f,g);
set(gca,'xlim',[0 cfg.maxDepthConfig])
box on
title('Range and vertical velocities')

ax(2) = subplot(2,1,2);
erbar(range,dh./dt,dhe./dt,-dhe./dt,'k','k'); % errors
ylabel('vertical velocity (m day^{-1})')
box on
% Try rudimentary fit 
c = polyfit(range, dh./dt, 1);
hold on;
v_est = polyval(c,range);
plot(ax(2),range,v_est,'r-');
hold off;
linkaxes(ax,'x')



function Disp(text) % Only diplay output if cfg.verbose on
global cfg
if cfg.verbose
    disp(text)
end


function h = fillXRange(x,varargin)
y = get(gca,'ylim');
h = patch([x(1) x(2) x(2) x(1)],[y(1) y(1) y(2) y(2)],[0.6 0.6 0.6]);
set(h,varargin{:})


function plotAmp(f,g,ii)
if nargin<3
    ii = 1:numel(f.rangeCoarse);
end
%plot the standard amplitude profile
fah = plot(f.rangeCoarse(ii),20*log10(abs(f.specCor(ii))),'r');
hold on
gah = plot(g.rangeCoarse(ii),20*log10(abs(g.specCor(ii))),'b');
ylabel('Vrms (dB)')
xlabel('range (m)')
legend([fah(1) gah(1)],{datestr(f.TimeStamp),datestr(g.TimeStamp)},'Location','SouthWest')
ylim([-140 -20])

