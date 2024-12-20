% This script:
% reads the pRES data (mono & mimo)
% Cut the data according to the bed it finds
% Applys Smoothing and then synthesizes the signal
% Calculates the necessary parameters to analyse fabric
% Plot the most important ones

% Written by Reza Ershadi
% May 2023
% University of Tübingen - Geophysics group
%%
clear
close all
clc
restoredefaultpath
addpath(genpath('fmcw'));
addpath(genpath('func'));
%% Check these values ************************ IMPORTANT
maxRange = 1301;
whichATT = 1;
whichBURST = 1;
SF = [1 1 -1 -1]; % a minus factor to correct the antenna orientation on [RR RU UR UU]
% use [1 1 1 1] for the mono test file
% use [-1 -1 1 1] for the mimo test file
%%
fn = '40_SubZero__195754.40_T1HT2VR1HR2V.dat';
% fn = '10_SubZero__193608.20_T1HT2VR1HR2V.dat';

fd = '/Users/reza/pCloud/IceTub/Devices/ApRES/projects/Antarctica_ReMeltRadar2122/data/radar/Rover/Along_A_26_25_24_23_22/';
% [fn, fd] = uigetfile({'*.dat';'*.DAT'}, 'Select a file','MultiSelect','on');
fp = fullfile(fd,fn);
if ischar(fp)
    FileType = "MIMO";
elseif iscell(fp) && (size(fp,2) == 3 || size(fp,2) == 4)
    FileType = "MONO";
end
%% Read the data
switch FileType
    case "MIMO"
        DAT = FUNC_ReadMonsterFileTS(fp,maxRange,whichATT,whichBURST);
        sRR = DAT(1).Signal; 
        sRU = DAT(2).Signal; 
        sUR = DAT(3).Signal; 
        sUU = DAT(4).Signal;
        z = DAT(1).Z;
        zt = DAT(1).ZT;
    case "MONO"
        iRR = find(contains(fn, ['HH',"RR"]));
        DAT_RR = FUNC_ReadMonsterFileTS(fp{iRR},maxRange,whichATT,whichBURST);
        sRR = DAT_RR.Signal;
        iUU = find(contains(fn, ['VV',"UU"]));
        DAT_UU = FUNC_ReadMonsterFileTS(fp{iUU},maxRange,whichATT,whichBURST);
        sUU = DAT_UU.Signal;
        if size(fn,2) == 3
            iRU = find(contains(fn, ['HV',"VH","RU","UR"]));
            DAT_RU = FUNC_ReadMonsterFileTS(fp{iRU},maxRange,whichATT,whichBURST);
            DAT_UR = DAT_RU;
            sRU = DAT_RU.Signal;
            sUR = DAT_UR.Signal;
        elseif size(fn,2) == 4
            iRU = find(contains(fn, ['HV',"RU"]));
            DAT_RU = FUNC_ReadMonsterFileTS(fp{iRU},maxRange,whichATT,whichBURST);
            sRU = DAT_RU.Signal;
            iUR = contains(fn, ['VH',"UR"]);
            DAT_UR = FUNC_ReadMonsterFileTS(fp{iUR},maxRange,whichATT,whichBURST);
            sUR = DAT_UR.Signal;
        end
        z = DAT_RR.Z;
        zt = DAT_RR.ZT;
end
Z = z;
[sRR,sRU,sUR,sUU,Z,zBED] = cut2bed(sRR,sRU,sUR,sUU,Z);
sRR = SF(1) .* sRR; 
sRU = SF(2) .* sRU; 
sUR = SF(3) .* sUR; 
sUU = SF(4) .* sUU;
ZB = round(mean(zBED));
maxRange = max(Z);
%% Denoise and synthesize
dA = 1;
ao = 0:dA:179; 
f = 3.0000e+08;
C_DepthWin = maxRange * 0.05;
C_ConvWin = maxRange * 0.05;
DenoisingFlag.PA = [  "1", "MovingAverage"  , string(maxRange*0.05) ;
                      "0", "Conv1D"         , string(maxRange*0.1) ;
                      "2", "Conv2D"         , string(maxRange*0.05) ;
                      "0", "DenoisePCA"     , string(1)];
DenoisingFlag.PD = [  "1", "MovingAverage"  , string(maxRange*0.05) ;
                      "0", "Conv1D"         , string(maxRange*0.01) ;
                      "0", "Conv2D"         , string(maxRange*0.01) ;
                      "0", "DenoisePCA"     , string(1)];
[HH,VH,HV,VV] = QuadpoleSynthesizer(sRR,sUR,sUR,sUU,ao,0);
Obs = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
pwrRR = 20.*log10(abs(Obs{1}(:,1)));
pwrRU = 20.*log10(abs(Obs{2}(:,1)));
pwrUR = 20.*log10(abs(Obs{3}(:,1)));
pwrUU = 20.*log10(abs(Obs{4}(:,1)));
pwrLim = round([min([pwrRR;pwrRU;pwrUR;pwrUU])-20 max([pwrRR;pwrRU;pwrUR;pwrUU])+20],0);
ll = -0.01; % girdle strength lower limit
dL = Obs{18};
dL(dL<ll) = nan;
dL(dL>1) = nan;  
%% Plot
load('SeismicColorMap100.mat');
load('ColorMapAbsC.mat');
load('ColorMapDLAMBDA.mat');
 
fntsz = 14;
fig0 = CLASS_FixedPlot.SetFigureSize(0.025,0.025,0.5,0.9);
rw = 2;
cl = 12;
ax{1} = subplot(rw,cl,1);
ii = 1;
hold(ax{ii},'on')
plot(ax{ii},pwrRR,Z,'-b')
plot(ax{ii},pwrLim,[ZB ZB],'--k','LineWidth',2)
grid(ax{ii},'on')
set(ax{ii},'YDir','reverse')
xlim(ax{ii},pwrLim)
ylim(ax{ii},[0 Z(end)])
xticks(ax{ii},pwrLim)
xtickangle(ax{ii},90)
ylabel(ax{ii},'Depth [m]')
title(ax{ii},'RR')
set(ax{ii},'FontSize',fntsz)
%---
ax{2} = subplot(rw,cl,2);
ii = 2;
hold(ax{ii},'on')
plot(ax{ii},pwrRU,Z,'-b')
plot(ax{ii},pwrLim,[ZB ZB],'--k','LineWidth',2)
grid(ax{ii},'on')
set(ax{ii},'YDir','reverse')
yticklabels(ax{ii},[])
xlim(ax{ii},pwrLim)
ylim(ax{ii},[0 Z(end)])
xticks(ax{ii},pwrLim)
xtickangle(ax{ii},90)
title(ax{ii},'RU')
set(ax{ii},'FontSize',fntsz)
%---
ax{3} = subplot(rw,cl,3);
ii = 3;
hold(ax{ii},'on')
plot(ax{ii},pwrUR,Z,'-b')
plot(ax{ii},pwrLim,[ZB ZB],'--k','LineWidth',2)
grid(ax{ii},'on')
set(ax{ii},'YDir','reverse')
yticklabels(ax{ii},[])
xlim(ax{ii},pwrLim)
ylim(ax{ii},[0 Z(end)])
xticks(ax{ii},pwrLim)
xtickangle(ax{ii},90)
title(ax{ii},'UR')
set(ax{ii},'FontSize',fntsz)
%---
ax{4} = subplot(rw,cl,4);
ii = 4;
hold(ax{ii},'on')
plot(ax{ii},pwrUU,Z,'-b')
plot(ax{ii},pwrLim,[ZB ZB],'--k','LineWidth',2)
grid(ax{ii},'on')
set(ax{ii},'YDir','reverse')
yticklabels(ax{ii},[])
xlim(ax{ii},pwrLim)
ylim(ax{ii},[0 Z(end)])
xticks(ax{ii},pwrLim)
xtickangle(ax{ii},90)
title(ax{ii},'UU')
set(ax{ii},'FontSize',fntsz)
%---
ax{5} = subplot(rw,cl,5:8);
ii = 5;
hold(ax{ii},'on')
pc = pcolor(ax{ii},ao,Z,Obs{5});
set(pc, 'EdgeColor', 'none');
plot(ax{ii},[0 179],[ZB ZB],'--k','LineWidth',2)
xticks(ax{ii},[0 90 179])
% xticklabels(ax{ii},[])
yticklabels(ax{ii},[])
set(ax{ii},'YDir','reverse')
xlim(ax{ii},[0 179])
ylim(ax{ii},[0 Z(end)])
colormap(ax{ii},cm_rb)
caxis(ax{ii},[-5 5])
colorbar(ax{ii})
title(ax{ii},'\deltaP_{RR}')
set(ax{ii},'FontSize',fntsz)
%---
ax{6} = subplot(rw,cl,9:12);
ii = 6;
hold(ax{ii},'on')
pc = pcolor(ax{ii},ao,Z,Obs{7});
set(pc, 'EdgeColor', 'none');
plot(ax{ii},[0 179],[ZB ZB],'--k','LineWidth',2)
xticks(ax{ii},[0 90 179])
% xticklabels(ax{ii},[])
yticklabels(ax{ii},[])
set(ax{ii},'YDir','reverse')
xlim(ax{ii},[0 179])
ylim(ax{ii},[0 Z(end)])
colormap(ax{ii},cm_rb)
caxis(ax{ii},[-5 5])
colorbar(ax{ii})
title(ax{ii},'\deltaP_{RU}')
set(ax{ii},'FontSize',fntsz)
%---
ax{7} = subplot(rw,cl,13:16);
ii = 7;
hold(ax{ii},'on')
pc = pcolor(ax{ii},ao,Z,Obs{13});
plot(ax{ii},[0 179],[ZB ZB],'--k','LineWidth',2)
set(pc, 'EdgeColor', 'none');
set(ax{ii},'YDir','reverse')
xlim(ax{ii},[0 179])
xticks(ax{ii},[0 90 179])
ylim(ax{ii},[0 Z(end)])
colormap(ax{ii},cm_CMAG)
caxis(ax{ii},[0 1])
cb = colorbar(ax{ii});
cb.Ticks = [0 0.4 1];
xlabel(ax{ii},'\alpha [deg]')
ylabel(ax{ii},'Depth [m]')
title(ax{ii},'|C_{RR,UU}|')
set(ax{ii},'FontSize',fntsz)
%---
ax{8} = subplot(rw,cl,17:20);
ii = 8;
hold(ax{ii},'on')
pc = pcolor(ax{ii},ao,Z,Obs{14});
plot(ax{ii},[0 179],[ZB ZB],'--k','LineWidth',2)
set(pc, 'EdgeColor', 'none');
yticklabels(ax{ii},[])
set(ax{ii},'YDir','reverse')
xlim(ax{ii},[0 179])
xticks(ax{ii},[0 90 179])
ylim(ax{ii},[0 Z(end)])
colormap(ax{ii},cm_rb)
caxis(ax{ii},[-pi pi])
colorbar(ax{ii})
xlabel(ax{ii},'\alpha [deg]')
title(ax{ii},'\phi_{C_{RR,UU}}')
set(ax{ii},'FontSize',fntsz)
%---
ax{9} = subplot(rw,cl,21:24);
ii = 9;
hold(ax{ii},'on')
pc = pcolor(ax{ii},ao,Z,dL);
set(pc, 'EdgeColor', 'none');
plot(ax{ii},[0 179],[ZB ZB],'--k','LineWidth',2)
yticklabels(ax{ii},[])
set(ax{ii},'YDir','reverse')
xlim(ax{ii},[0 179])
xticks(ax{ii},[0 90 179])
ylim(ax{ii},[0 Z(end)])
colormap(ax{ii},cm_dLambda)
caxis(ax{ii},[0 1])
cb = colorbar(ax{ii});
cb.Ticks = [0:0.1:1];
xlabel(ax{ii},'\alpha [deg]')
title(ax{ii},'\Delta\lambda')
set(ax{ii},'FontSize',fntsz)

function [sRR,sRU,sUR,sUU,z,zBED] = cut2bed(sRR,sRU,sUR,sUU,z)
%     dz = mean(diff(z));
%     [zBED(1),ibed(1)] = FUNC_FindBed_RE(sRR,z);
%     [zBED(2),ibed(2)] = FUNC_FindBed_RE(sRU,z);
%     [zBED(3),ibed(3)] = FUNC_FindBed_RE(sUR,z);
%     [zBED(4),ibed(4)] = FUNC_FindBed_RE(sUU,z);
%     iBed = round(mean(ibed),0);
%     ExtraDepth = iBed + round((25/dz),0);
%     if ExtraDepth > length(z)
%         ExtraDepth = length(z);
%     end
    ExtraDepth = length(z);
    sRR = sRR(1:ExtraDepth);
    sRU = sRU(1:ExtraDepth);
    sUR = sUR(1:ExtraDepth);
    sUU = sUU(1:ExtraDepth);
    z = z(1:ExtraDepth);
    zBED = max(z);
end







