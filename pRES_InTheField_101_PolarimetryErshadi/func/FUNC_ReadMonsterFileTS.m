% Written by Reza Ershadi
% May 2023
% University of Tübingen - Geophysics group
function [DAT] = FUNC_ReadMonsterFileTS(filePath,maxRange,att,Burst)
vdat = fmcw_load(filePath,Burst);
%%
p=1; % pad factor (i.e. level of interpolation to use during fft)
winFun=@blackman; % window function handle
frange=[2e8 , 4e8]; %Normally [2e8,4e8]
rhos=407.0613; % surface density (kg/m^3)
rhoi=907.7165; % ice density (kg/m^3)
Lrho=39.5512;  % half depth decay (m) (assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho))
% Lrho=20;  % half depth decay (m) (assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho))
nI=1.68;
% n=3;
%%
tic
nSubBurst = vdat.SubBurstsInBurst; % number of sub-bursts in the burst
nTx = length(vdat.TxAnt); % number of Tx antenna
nRx = length(vdat.RxAnt); % number of Rx antenna
nAtt = vdat.NAttenuators; % number of attenuator
n = 0;
% Index order: sub-burst --> Tx --> Rx --> att
for i = 1:nSubBurst
    for j = 1:nTx
        for k =1:nRx
            for l = 1:nAtt
                n = n+1;
                Guid(n,1) = i;
                Guid(n,2) = j;
                Guid(n,3) = k;
                Guid(n,4) = l;
                Guid(n,5) = n;
            end
        end
    end
end
G = sortrows(Guid,[2,3,4,1]);
n = 0;
% Creating new parameters
vdat.TIME = []; % Time
vdat.Z = []; % Depth
vdat.ZT = []; % Corrected depth
vdat.Signal = []; % Signal
for i = 1:nSubBurst:vdat.ChirpsInBurst
    n = n+1;
    ii = i:(i+nSubBurst)-1;
    vDats(n) = vdat;
    vDats(n).NAttenuators = 1;
    vDats(n).TxAnt = G(i,2);
    vDats(n).RxAnt = G(i,3);
    vDats(n).ChirpsInBurst = nSubBurst;
    vDats(n).Attenuator_1 = vdat.Attenuator_1(G(i,4));
    vDats(n).Attenuator_2 = vdat.Attenuator_2(G(i,4));
    if size(vdat.vif,1) == 1
        vDats(n).vif = vdat.vif;
        vDats(n).chirpNum = vdat.chirpNum;
        vDats(n).chirpAtt = vdat.chirpAtt;
        vDats(n).chirpTime = vdat.chirpTime;
    else
        vDats(n).vif = vdat.vif(G(ii,5),:);
        vDats(n).chirpNum = vdat.chirpNum(G(ii,5));
        vDats(n).chirpAtt = vdat.chirpAtt(G(ii,5));
        vDats(n).chirpTime = vdat.chirpTime(G(ii,5));
    end

    % --------
%     vDats(n) = fmcw_cull_freq(vDats(n),frange);
    % --------
    vDats(n) = fmcw_burst_mean(vDats(n));
    vDats(n).TIME = datetime(vDats(n).TimeStamp,'ConvertFrom','datenum');
    % -------- Signal and depth
    [Range,Rfine,SpecCor,spec] = fmcw_range(vDats(n),p,maxRange,winFun);
    % -------- Depth correction
    ShouldBeZero=@(d,dI,L,RhoSp) -dI+d+L*(nI-1)/nI*(1-RhoSp)*(exp(-d/L)-1);
    TrueDepthfun=@(RhoSp,L,dI) FUNC_bisection(@(d) ShouldBeZero(d,dI,L,RhoSp),0,max(dI)+20);
    TrueDepth = TrueDepthfun(rhos/rhoi,Lrho,Range);
    % --------
    Range(1) = 1e-20;
    TrueDepth(1) = 1e-20;
    vDats(n).Z = Range';
    vDats(n).ZT = TrueDepth';
    vDats(n).Signal = SpecCor.';
end

sAtt1 = vDats(att).Attenuator_1;
sAtt2 = vDats(att).Attenuator_2;

n = 1;
for i = 1:size(vDats,2)
    if vDats(i).Attenuator_1 == sAtt1 && vDats(i).Attenuator_2 == sAtt2
        DAT(n) = vDats(i);
        n = n+1;
    end
end



