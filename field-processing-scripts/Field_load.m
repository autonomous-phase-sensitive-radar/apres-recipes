function vdat = Field_load(filename,burst)

% vdat = fmcw_load(filename,burst)

% Load FMCW radar burst and metadata into Field_fmcw_plot_v_fft_hist 
vdat = Field_LoadBurstRMB5(filename, burst);

% Check file was found
switch(vdat.Code)
    case 2
        fprintf('Too few measurements in burst\n');
        return
    case 1
        fprintf('Corrupt header\n');
        return
    case -1
        fprintf('Unable to open file: %s\n',filename);
        return
    case -2
        fprintf('Corrupt header in burst - fatal %d\n',vdat.Burst);
        return
    case -4
        %        disp(['Burst ' int2str(burst) ' not found in file ' filename]);
        return
end

% Extract just good chirp data from voltage record and rearrange into
% matrix with one chirp per row

%% Add metadata to structure
AttSet = vdat.Attenuator_1 + 1i*vdat.Attenuator_2; % unique code for attenuator setting

% Sampling parameters
vdat.filename = filename;
vdat.processing = {};
vdat.SamplesPerChirp = vdat.Nsamples; %if vdat.SamplesPerChirp=vdat.Nsamples, why keep both

H = Field_ParametersRMB2(vdat.filename);
    
vdat.K = H.K;
vdat.f0 = H.startFreq;
vdat.fs = H.fs;
vdat.f1 = H.startFreq + H.chirpLength * H.K/2/pi;
vdat.SamplesPerChirp = round(H.chirpLength * H.fs);
vdat.T = H.chirpLength;
vdat.B = H.chirpLength * H.K/2/pi;
vdat.fc = H.startFreq + vdat.B/2;

vdat.dt = 1/H.fs;
vdat.er = 3.18;
vdat.ci = 3e8/sqrt(vdat.er);
vdat.lambdac = vdat.ci/vdat.fc;
%vdat.Nsamples = H.nchirpSamples; %hasn't this already been set by N_ADC_SAMPLES ?

% Load each chirp into a row
vdat.Endind = vdat.Startind + vdat.SamplesPerChirp - 1;
vdat.vif = zeros(vdat.ChirpsInBurst,vdat.SamplesPerChirp); % preallocate array
chirpInterval = 1/(24*3600); % time of a chirp in units of days
for chirp = 1:vdat.ChirpsInBurst
    vdat.vif(chirp,:) = vdat.v(vdat.Startind(chirp):vdat.Endind(chirp));
    vdat.chirpNum(chirp,1) = chirp; % chirp number in burst
    vdat.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); % attenuator setting for chirp
    vdat.chirpTime(chirp,1) = vdat.TimeStamp + chirpInterval*(chirp-1); % time of chirp
end


% Create time and frequency stamp for samples
vdat.t = vdat.dt*(0:size(vdat.vif,2)-1); % sampling times (rel to first)
vdat.f = vdat.f0 + vdat.t.*vdat.K/(2*pi);
