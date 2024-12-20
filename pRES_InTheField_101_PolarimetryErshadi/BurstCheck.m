% This script:
% Is useful to quickly check the quality of your data in the field
% Note that it doesnt blut every single chirp
% You only see the mean of all the chirps in a burst
% It reads a single pRES file (mono & mimo)
% plots the voltage, power and phase of the signal

% Written by Reza Ershadi
% May 2023
% University of Tübingen - Geophysics group
%%
clear
% close all
clc
restoredefaultpath
addpath(genpath('fmcw'));
addpath(genpath('func'));
%%
[fn, fd] = uigetfile({'*.dat';'*.DAT'}, 'Select a file');
fp = fullfile(fd,fn);
%%
maxRange = 2501;
whichATT = 1;
whichBURST = 1;
DAT = FUNC_ReadMonsterFileTS(fp,maxRange,whichATT,whichBURST);
FileFlag = "normal";
if size(DAT,2) >= 3
    FileFlag = "mimo";
end
%%
fprintf('In file: %s \n',fn)
fprintf('Burst: %i \n',whichBURST)
fprintf('Attenuator: %i \n',DAT(1).Attenuator_1)
fprintf('Gain: %i \n',DAT(1).Attenuator_2)
for i = 1:size(DAT,2)
    fprintf('INDEX %i --> Tx:%i & Rx:%i \n',i,DAT(i).TxAnt,DAT(i).RxAnt)
    t(:,i) = DAT(i).t; % time
    v(:,i) = DAT(i).vif; % voltage
    pwr(:,i) = 20.*log10(abs(DAT(i).Signal)); % power
    phs(:,i) = angle(DAT(i).Signal);
    z(:,i) = DAT(i).Z;
    zt(:,i) = DAT(i).ZT;
    zBed(:,i) = FUNC_FindBed_RE(DAT(i).Signal,z(:,i));
end
%%
if FileFlag == "normal"
    data_plt(t,v,pwr,phs,z,zBed,DAT(1).TxAnt,DAT(1).RxAnt)
else
    for i = 1:size(DAT,2)
        data_plt(t(:,i),v(:,i),pwr(:,i),phs(:,i),z(:,i),zBed(:,i),DAT(i).TxAnt,DAT(i).RxAnt)
    end
end
%%
function [] = data_plt(t,v,pwr,phs,z,zBed,Tx,Rx)
    f1 = CLASS_FixedPlot.SetFigureSize(0,0,0.5,0.8);
    ax{1} = subplot(1,10,[1:3]);
    ax{2} = subplot(1,10,[5:7]);
    ax{3} = subplot(1,10,[8:10]);
    fsz = 14;
%     v = v+1;

    ii = 1;
    % 0.5<v<2 green
    % 0<v<0.5   &  2<v<2.5 yellow
    % v<0 & v>2.5 red
    iyellow = (v>2 & v<=2.5) | (v>0 & v<=0.5);
    ired = v>2.5 | v<0;
    vyellow = nan(size(v));
    vred = nan(size(v));
    vyellow(iyellow) = v(iyellow);
    vred(ired) = v(ired);
    plot(ax{ii},v(:,1),t(:,1),'g');
    hold(ax{ii},'on')
    plot(ax{ii},vyellow(:,1),t(:,1),'y');
    plot(ax{ii},vred(:,1),t(:,1),'r');
    plot(ax{ii},[2 2],[0 1],':w')
    plot(ax{ii},[0.5 0.5],[0 1],':w')
    plot(ax{ii},[0 0],[0 1],':w')
    plot(ax{ii},[2.5 2.5],[0 1],':w')
    set(ax{ii},'YDir','reverse')
    xlim(ax{ii},[-0.5 3])
    ylim(ax{ii},[0 1])
    ax{ii}.Color = 'k';
    xticks(ax{ii},[-0.25 0.25 1.25 2.25 2.75])
    xticklabels(ax{ii},["RED","OK","GREEN","OK","RED"])
    xtickangle(ax{ii},90)
    ylabel(ax{ii},'Time [s]')
    title(ax{ii},'Voltage')
    set(ax{ii},'FontSize',fsz)
    
    ii = 2;
    pwrlim = round([min(pwr(:)) max(pwr(:))],0);
    plot(ax{ii},pwr(:,1),z(:,1),'.-w');
    hold(ax{ii},'on')
    plot(ax{ii},[pwrlim(1)-20 pwrlim(2)+20],[zBed zBed],'--c','linewidth',2);
    set(ax{ii},'YDir','reverse')
    xlim(ax{ii},[pwrlim(1)-20 pwrlim(2)+20])
    xticks(ax{ii},[pwrlim(1) mean(pwrlim) pwrlim(2)])
    xtickangle(ax{ii},90)
    ylim(ax{ii},[0 max(z)])
    xlabel(ax{ii},'[dB]')
    ax{ii}.Color = 'k';
    ylabel(ax{ii},'Depth [m]')
    title(ax{ii},'Power')
    set(ax{ii},'FontSize',fsz)
    
    ii = 3;
    plot(ax{ii},phs(:,1),z(:,1),'.-r');
    hold(ax{ii},'on')
    plot(ax{ii},[-4 4],[zBed zBed],'--c','linewidth',2);
    set(ax{ii},'YDir','reverse')
    plot(ax{ii},[0 0],[z(1) z(end)],':w')
    plot(ax{ii},[-pi -pi],[z(1) z(end)],':w')
    plot(ax{ii},[pi pi],[z(1) z(end)],':w')
    xlim(ax{ii},[-4 4])
    xlabel(ax{ii},'[Rad]')
    ylim(ax{ii},[0 max(z)])
    yticks(ax{ii},[])
    ax{ii}.Color = 'k';
    xticks(ax{ii},[-pi 0 pi])
    xticklabels(ax{ii},["-\pi",'0',"\pi"])
    title(ax{ii},'Phase')
    set(ax{ii},'FontSize',fsz)

    linkaxes([ax{2}, ax{3}], 'y');
    annotation(f1,'textbox',...
        [0.00655555555555571 0.909255898366606 0.0934444444444443 0.0834845735027223],...
        'String',{"Tx: "+string(Tx),"Rx: "+string(Rx)},...
        'FontWeight','bold',...
        'FontSize',fsz,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor','none');

    print(f1,"Priestley_2022_DC_UR.png",'-dpng','-r300');

end
