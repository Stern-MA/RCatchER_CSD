%SCRIPT TO GENERATE Analysis and FIGURES for RCatchER Paper

%% Figure 1: Version s56 210831 bl
load('s56_210831_bl/S56base1.mat')
aID=[S56base1,S56base1];
fs_2p=30;
kk=1;

%% Fig1F:Traces
%

figure
plot(aID(kk).speedDataTS1,aID(kk).speedData1/100+1,'color',[0 0 0])
hold on
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt1,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt1,1)-1,'color',[1 0 1]);
ylim([-1.5 3.5])
xlim([0 300])
yticks([])
pbaspect([1 .5 1])
hold off

%scale bars
sb_x=20; %20seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_yv=.25; %25mm/s (=#*100; .25*10=25mm/s) accounting for scaling below divide by 100
sb_omg=[0 -0.5];
sb_omr=[0 -1.5];
sb_ov=[0 .5];

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars Velocity
line([sb_ov(1) sb_ov(1)+sb_x],[sb_ov(2) sb_ov(2)],'LineWidth',1,'color',[0 0 0])
line([sb_ov(1) sb_ov(1)],[sb_ov(2) sb_ov(2)+sb_yv],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

% Rasters
figure
kk=1;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:300])
xticklabels([0:50:300])
yticks([253:50:507])
yticklabels([0:50:500])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([0 300])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:300])
xticklabels([0:50:300])
yticks([253:50:507])
yticklabels([0:50:500])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([0 300])


%% Fig 1: Indiv  cells...could add in
cellList=[62,167,183];

figure
hold on
cc=1;

for ii=cellList
    plot(aID(kk).F1_ts,aID(kk).Fc1Rdff_flt2(ii,:)+cc,'color', [0.75 0 0.75],'linewidth',1);
    plot(aID(kk).F1_ts,aID(kk).Fc1Gdff_flt2(ii,:)+cc+1,'color', [0 0.75 0],'linewidth',1);
    line([aID(kk).RecTime_spks{ii,:};aID(kk).RecTime_spks{ii,:}],repmat([cc+1;cc+1.5],1,size(aID(kk).RecTime_spks{ii,:},2)),'color',[0 0 .5],'LineWidth',1.25)
    cc=cc+3.5;
end

xlim([125 250])
ylim([0 cc-.5])

%scale bars
sb_x=10; %10seconds
sb_ymg=.25; %0.25dff(=#*1)
sb_ymr=.25; %0.25dff
sb_omg=[130 1.25];
sb_omr=[130 .5];

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1.5,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1.5,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1.5,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1.5,'color',[0.75 0 0.75])

set(gca,'Visible','off')


print('../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_IndivTraces.svg','-dsvg')


%% Fig 1: Spike Raster and Example with Calcium Level Change

%Vectorize the RecTimes
tempVec=cell2mat(cellfun(@(x) [x nan(1,10-numel(x))],aID(kk).RecTime_spks,'uni',0));
tempVec=transpose(tempVec(~isnan(tempVec)));


figure('Position',[100 100 500 400])
% subplot(4,1,1)
% plot(aID(kk).speedDataTS1,aID(kk).speedData1,'color',[0 0 0])
% xlim([0 300])
%Spike Frequency
subplot(3,1,1)
histogram(tempVec,[1:1:300],'FaceColor',[.5 .5 .5],'EdgeColor',[.25 .25 .25],'Normalization','count')
xlim([0 300])
xticks([])
set(gca,'box','off')

%Spike raster
subplot(3,1,2:3)
for jj=1:size(aID(kk).RecTime_spks,1)
line(repmat(aID(kk).RecTime_spks{jj},2,1),repmat([jj-0.5;jj+0.5],1,length(aID(kk).RecTime_spks{jj})),'color',[0 0 0],'linewidth',1)
hold on
end
clear jj
ylim([0 size(aID(kk).RecTime_spks,1)])
yticks(0:50:300)
xlim([0 300])
set(gca,'box','on')

% set(1,'renderer','painters')
% print('../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_SpkRaster.eps','-depsc')
print('../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_SpkRaster.svg','-dsvg')

%Spike Average
win_pre=1;
win_post=3;

figure('Position',[100 100 400 400])
tempX=[-win_pre:1/fs_2p:win_post];
tempYlG=(aID(kk).spks_trace_avgWtdG-aID(kk).spks_trace_stdpG);
tempYuG=(aID(kk).spks_trace_avgWtdG+aID(kk).spks_trace_stdpG);
tempYlR=(aID(kk).spks_trace_avgWtdR-aID(kk).spks_trace_stdpR);
tempYuR=(aID(kk).spks_trace_avgWtdR+aID(kk).spks_trace_stdpR);
%Green
plot(tempX,aID(kk).spks_trace_avgWtdG,'color',[0 .7 0],'linewidth',2.5)
hold on
%plot(tempX,tempYlG,'color',[0 1 0],'linestyle',':')
%plot(tempX,tempYuG,'color',[0 1 0],'linestyle',':')
patch([tempX fliplr(tempX)], [tempYlG fliplr(tempYuG)], [.5 1 .5],'FaceAlpha',0.2,'EdgeColor','none')
%Red
plot(tempX,aID(kk).spks_trace_avgWtdR,'color',[.7 0 .7],'linewidth',2.5)
%plot(tempX,tempYlR,'color',[1 0 1],'linestyle',':')
%plot(tempX,tempYuR,'color',[1 0 1],'linestyle',':')
patch([tempX fliplr(tempX)], [tempYlR fliplr(tempYuR)], [1 .5 1],'FaceAlpha',0.2,'EdgeColor','none')
xlim([-win_pre win_post])
ylim([-.5 2.1])

% set(2,'renderer','painters')
% print('../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_AvgSpk.eps','-depsc')
print('../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_AvgSpk.svg','-dsvg')

%Ca Levels Plot (Pooled Mean+Pooled SEM)

[Ca_Spks_means,Ca_Spks_std]=deal(nan(size(aID(kk).RecTime_spks,1),2));
Ca_Spks_means(:,1)=cellfun(@nanmean,aID(kk).spks_traceCaDiffG);
Ca_Spks_means(:,2)=cellfun(@nanmean,aID(kk).spks_traceCaDiffR);
Ca_Spks_std(:,1)=cellfun(@nanstd,aID(kk).spks_traceCaDiffG);
Ca_Spks_std(:,2)=cellfun(@nanstd,aID(kk).spks_traceCaDiffR);
Ca_Spks_df=aID(kk).spks_trace_df;

% Pooled SE
Ca_Spks_stdp=sqrt(nansum(Ca_Spks_df.*Ca_Spks_std.^2,1)./nansum(Ca_Spks_df,1));
%Ca_Spks_sep=Ca_Spks_stdp.*sqrt(nansum(1./(Ca_Spks_df+1)));


figure('Position',[100 100 200 400])
plotYvals=[nanmean(Ca_Spks_means(:,1)); nanmean(Ca_Spks_means(:,2))];
%plotYerr=[Ca_Spks_sep(1); Ca_Spks_sep(2)];
plotYerr=[Ca_Spks_stdp(1); Ca_Spks_stdp(2)];
barPlot=bar(plotYvals,'grouped','FaceColor','flat','FaceAlpha',.2,'EdgeColor','flat','LineWidth',1);
barPlot.CData = [0 .9 0 ; 0.9 0 .9];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
%ylim([-.25 .75])
xticks(1.5)
xticklabels({'Spks'})
yticks(-.25:.25:1.75)
ylim([-.25 1.75])
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_SpkCaLvl.svg','-dsvg')



%% Figure 1: Version s102 -Baseline
load('s102_230530_bl/S102base1.mat');
aID=[S102base1,S102base1];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;

%% Fig1F:Traces
% 
figure
kk=1;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+2.25,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./10+1.75,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt1,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt1,1)-1,'color',[1 0 1]);
ylim([-1.5 3])
xlim([0 300])
yticks([])
pbaspect([1 .5 1])
hold off

%scale bars
sb_x=20; %20seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .25uV*2000=.5mV) accounting for scaling below divide by 2000
sb_yd=.25; %2.5mV (=#*10; .25mV*10=2.5mV) accounting for scaling below divide by 20
sb_omg=[0 -0.5];
sb_omr=[0 -1.5];
sb_oe=[0 1.5];
sb_od=[0 .5]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

% Rasters
figure
kk=1;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:300])
xticklabels([0:50:300])
yticks([350:50:701])
yticklabels([0:50:500])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([0 300])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:300])
xticklabels([0:50:300])
yticks([350:50:701])
yticklabels([0:50:500])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([0 300])


%% Fig 1: Indiv  cells...could add in
cellList=[5,6,7];

figure
hold on
cc=1;

for ii=cellList
%     plot(aID(kk).F1_ts,SternRollAvg(aID(kk).Fc1Rdff_flt2(ii,:),10)+cc,'color', [0.75 0 0.75]);
%     plot(aID(kk).F1_ts,SternRollAvg(aID(kk).Fc1Gdff_flt2(ii,:),10)/2+cc+1,'color', [0 0.75 0]);
    plot(aID(kk).F1_ts,aID(kk).Fc1Rdff_flt2(ii,:)+cc,'color', [0.75 0 0.75]);
    plot(aID(kk).F1_ts,aID(kk).Fc1Gdff_flt2(ii,:)/2+cc+1,'color', [0 0.75 0]);
    cc=cc+2;
end

plot(aID(kk).EEG_ts,aID(kk).EEG./2000+cc+1.5,'color',[0.5 0.5 0.5]);
plot(aID(kk).DC_ts,aID(kk).DC./20+cc+1,'color',[0 0 0]);
xlim([125 225])
ylim([0 cc+3])
pbaspect([1 1 1])

%scale bars
sb_x=5; %20seconds
sb_ymg=.25; %0.5dff(=#*2; .25dff*2=.5dff)
sb_ymr=.25; %0.25dff
sb_ye=.25; %250uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %2.5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[125 1.5];
sb_omr=[125 .5];
sb_oe=[125 8];
sb_od=[125 7.2]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars DC
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')






%% Figure 2: Version s89
load('s89_221027_2_ptz1/S89ptz1.mat');
load('s89_221101_ptz2/S89ptz2.mat');
load('s89_221109_ptz3/S89ptz3.mat');
aID=[S89ptz1,S89ptz2,S89ptz3];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;

%% Fig2A-C:Traces
% Pre-Ictal Only
figure
kk=1;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([900 1400])
yticks([])
pbaspect([1 .5 1])
hold off

%scale bars
sb_x=30; %30seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[905 -0.4];
sb_omr=[905 -1.2];
sb_oe=[905 2.1];
sb_od=[905 .9]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_89TraceSubGen.svg','-dsvg')

%Seizure w/o CSD
figure
kk=2;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([250 750])
yticks([])
pbaspect([1 .5 1])
hold off
set(gca,'Visible','off')

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_89TraceSzOnly.svg','-dsvg')

%Seizure w/ CSD
figure
kk=3;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([500 1000])
yticks([])
pbaspect([1 .5 1])
hold off

% %scale bars
% sb_x=30; %30seconds
% sb_ymg=.25; %0.25dff
% sb_ymr=.25; %0.25dff
% sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
% sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
% sb_omg=[525 -0.5];
% sb_omr=[525 -1.3];
% sb_oe=[525 2];
% sb_od=[525 1.1]; %
% 
% %Scale Bars Mean Green
% line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
% line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])
% 
% %Scale Bars Mean Red
% line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
% line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])
% 
% %Scale Bars EEG
% line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
% line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])
% 
% %Scale Bars EEG
% line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
% line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_89TraceSzCSD.svg','-dsvg')


%% Fig2 A-C: Rasters
%PIS only
figure
kk=1;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1500])
%xticklabels([0:50:1500])
xticklabels([])
yticks([156:30:313])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([900 1400])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1500])
xticklabels([0:50:1500])
yticks([156:30:313])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([900 1400])

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_89RasterSubGen.svg','-dsvg')

%Sz w/o CSD
figure
kk=2;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
%xticklabels([0:50:1050])
xticklabels([])
yticks([56:30:113])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([250 750])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
xticklabels([0:50:1050])
yticks([56:30:113])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([250 750])

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_89RasterSzOnly.svg','-dsvg')

%Sz w/ CSD
figure
kk=3;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
%xticklabels([0:50:1050])
xticklabels([])
yticks([141:30:283])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([500 1000])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
xticklabels([0:50:1050])
yticks([141:30:283])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([500 1000])

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_89RasterSzCSD.svg','-dsvg')


%% FigS3A-C: Bar Plots of Calcium Changes for S89

figure('Position',[100 100 800 300])
subplot(1,6,1)
kk=1;
% plotYvals=[mean(aID(kk).AvgPISCaChangeG(:,1)); -mean(aID(kk).AvgPISCaChangeR(:,1))];
% plotYerr=[std(aID(kk).AvgPISCaChangeG(:,1)); std(aID(kk).AvgPISCaChangeR(:,1))];
plotYvals=[mean(aID(kk).AvgPISCaChangeG(:,2)); -mean(aID(kk).AvgPISCaChangeR(:,2))];
plotYerr=[std(aID(kk).AvgPISCaChangeG(:,2)); std(aID(kk).AvgPISCaChangeR(:,2))];
barPlot=bar(plotYvals,'grouped','FaceColor','flat');
barPlot.CData = [0 .7 0 ; 0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
%ylim([-.25 1.25])
ylim([-.5 1.5])
xticks(1.5)
xticklabels({'PIS'})
set(gca,'box','off')

subplot(1,6,2:3)
kk=2;
% plotYvals=[mean(aID(kk).AvgPISCaChangeG(:,1)), mean(aID(kk).SzMeanDffDiffG); -mean(aID(kk).AvgPISCaChangeR(:,1)), -mean(aID(kk).SzMeanDffDiffR)]';
% plotYerr=[std(aID(kk).AvgPISCaChangeG(:,1)), std(aID(kk).SzMeanDffDiffG); std(aID(kk).AvgPISCaChangeR(:,1)), std(aID(kk).SzMeanDffDiffR)]';
plotYvals=[mean(aID(kk).AvgPISCaChangeG(:,2)), mean(aID(kk).SzMeanDffDiffG(aID(kk).isRecruitedG)); -mean(aID(kk).AvgPISCaChangeR(:,2)), -mean(aID(kk).SzMeanDffDiffR(aID(kk).isRecruitedG))]';
plotYerr=[std(aID(kk).AvgPISCaChangeG(:,2)), std(aID(kk).SzMeanDffDiffG(aID(kk).isRecruitedG)); std(aID(kk).AvgPISCaChangeR(:,2)), std(aID(kk).SzMeanDffDiffR(aID(kk).isRecruitedG))]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
%ylim([-.25 1.25])
ylim([-.5 1.5])
xticklabels({'PIS';'Sz'})
set(gca,'box','off')

subplot(1,6,4:6)
kk=3;
% plotYvals=[mean(aID(kk).AvgPISCaChangeG(:,1)), mean(aID(kk).SzMeanDffDiffG), mean(aID(kk).CSDMeanDffDiffG); -mean(aID(kk).AvgPISCaChangeR(:,1)), -mean(aID(kk).SzMeanDffDiffR), -mean(aID(kk).CSDMeanDffDiffR)]';
% plotYerr=[std(aID(kk).AvgPISCaChangeG(:,1)), std(aID(kk).SzMeanDffDiffG), std(aID(kk).CSDMeanDffDiffG); std(aID(kk).AvgPISCaChangeR(:,1)), std(aID(kk).SzMeanDffDiffR), std(aID(kk).CSDMeanDffDiffR)]';
plotYvals=[mean(aID(kk).AvgPISCaChangeG(:,2)), mean(aID(kk).SzMeanDffDiffG(aID(kk).isRecruitedG)), mean(aID(kk).CSDMeanDffDiffG(aID(kk).isRecruitedG_csd)); -mean(aID(kk).AvgPISCaChangeR(:,2)), -mean(aID(kk).SzMeanDffDiffR(aID(kk).isRecruitedG)), -mean(aID(kk).CSDMeanDffDiffR(aID(kk).isRecruitedR_csd))]';
plotYerr=[std(aID(kk).AvgPISCaChangeG(:,2)), std(aID(kk).SzMeanDffDiffG(aID(kk).isRecruitedG)), std(aID(kk).CSDMeanDffDiffG(aID(kk).isRecruitedG_csd)); std(aID(kk).AvgPISCaChangeR(:,2)), std(aID(kk).SzMeanDffDiffR(aID(kk).isRecruitedG)), std(aID(kk).CSDMeanDffDiffR(aID(kk).isRecruitedR_csd))]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
%ylim([-.25 1.25])
ylim([-.5 1.5])
xticklabels({'PIS';'Sz';'CSD'})

set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/FigureS2/FigureS2_89CaLvl.svg','-dsvg')

%% FigS3D: Bar Plots of Calcium Changes Post Seizure 
%Plots the max (G) and min(R) sustained calcium change in each post ictal
%period following  seizures to demsotnrate the post ictal calcium changes
%during CSD that do not occur without CSD.
%nned to invert R channel as it is the max of the inverted signal

figure('Position',[100 100 400 400])
plotYvals=[mean(aID(2).PostIctalCaGmax), mean(aID(3).PostIctalCaGmax); -mean(aID(2).PostIctalCaRmax), -mean(aID(3).PostIctalCaRmax)]';
plotYerr=[std(aID(2).PostIctalCaGmax), std(aID(3).PostIctalCaGmax); std(aID(2).PostIctalCaRmax), std(aID(3).PostIctalCaRmax)]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.5 1.5])
xticklabels({'Sz w/o CSD';'Sz w/ CSD'})
title('Post-Ictal Ca')
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/FigureS2/FigureS2_89PostIctalCaLvl.svg','-dsvg')







%% Get some values
% loads the calcium data for generating the glm table to determine
% %statistical significance

StructFieldNames={'FOV_dim','death','firstSpkI','fs_2p',...
    'fs_eeg','EEGbias','EEG_ts','F1_ts','EEG','PSD_F','PSD_P','DC','fs_dc','DC_ts',...
    'F1G','Fneu1G','stat','xCoor','yCoor','xCoorI','yCoorI','iscell','iscellI',...
    'Fb1G','Fneub1G','Fc1G','Fb1Gdff','Fneub1Gdff','Fc1Gdff','Fneub1Gdff_flt1',...
    'Fc1Gdff_flt1','Fneub1Gdff_flt1','Fc1Gdff_flt2','F1R','Fneu1R','Fb1R','Fneub1R',...
    'Fc1R','Fb1Rdff','Fneub1Rdff','Fc1Rdff','Fneub1Rdff_flt1','Fc1Rdff_flt1',...
    'Fneub1Rdff_flt1','Fc1Rdff_flt2','EEG_PIS_times2','EEG_PIS_times3',...
    'MeanCaG_PIS_times','MeanCaR_PIS_times','PISisRec','RecTimeG_PIS','RecTimeR_PIS',...
    'MeanSzRecTimeR','StimStartTime','StimEndTime','AvgPISCaChangeG','AvgPISCaChangeR',...
    'SzMeanDffDiffG','SzMeanDffDiffG','MeanCSDRecTimeR','CSDMeanDffDiffG','CSDMeanDffDiffR',...
    'CSDoffTimeDC','MeanCSDendTimeR','RecTimeG_sz2','RecTimeR_sz2',...
    'RecTimeG_csd','RecTimeR_csd','RecTimeG_csdRTN','RecTimeR_csdRTN','isRecG_csdRTN','isRecR_csdRTN'};

%szID={'S60ptz1','S60ptz2','S86ptz1','S86ptz2','S86ptz3','S86ptz4','S86ptz5','S87ptz1','S87ptz2','S87ptz3','S87ptz4','S87ptz5','S87ptz6','S87ptz7','S89ptz1','S89ptz2','S89ptz3','S102ptz1','S102ptz2','S102ptz3','S103ptz1','S104ptz1','S104ptz2','S105ptz1','S105ptz2','S105ptz3'};
%szID={'S60ptz1','S60ptz2','S86ptz2','S86ptz3','S86ptz4','S86ptz5','S87ptz4','S87ptz5','S87ptz7','S89ptz2','S89ptz3','S102ptz2','S102ptz3','S103ptz1','S104ptz1','S104ptz2','S105ptz1','S105ptz2','S105ptz3'};
szID={'S60ptz1','S60ptz2','S89ptz3','S102ptz3','S103ptz1','S105ptz3'};


struct_string='aID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii
clear *ptz*

szID={'S86ptz2','S86ptz3','S86ptz4','S86ptz5','S87ptz4','S87ptz5','S87ptz7','S89ptz2','S102ptz2','S104ptz1','S104ptz2','S105ptz1','S105ptz2'};


struct_string='bID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii

fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
clear *ptz*

cID=[aID,bID];

%Seizure start relative to PTZ
%spiking start time
%seizure length
%CSD length
%CSD latency


%% Fig 2 Additional Relevant Values (DC Shift mag, length)
AA=nan(1,length(cID));%w/ and w/o csd
for kk=1:length(cID)
    AA(kk)=cID(kk).firstSpkT;
end
disp({'time to PIS',[];'mean',nanmean(AA)+60;'std',nanstd(AA);'ste',nanstd(AA)/length(AA)^.5;'min',nanmin(AA)+60;'max',nanmax(AA)+60})

AA=nan(1,length(cID));%w/ and w/o csd
for kk=1:length(cID)
    AA(kk)=cID(kk).MeanSzRecTimeG;
end
disp({'time to Sz',[];'mean',nanmean(AA)+60;'std',nanstd(AA);'ste',nanstd(AA)/length(AA)^.5;'min',nanmin(AA)+60;'max',nanmax(AA)+60})

AA=nan(1,length(cID));%w/ and w/o csd
for kk=1:length(cID)
    AA(kk)=cID(kk).SzEndTime-cID(kk).MeanSzRecTimeG;
end
disp({'length of sz (all)',[];'mean',nanmean(AA);'std',nanstd(AA);'ste',nanstd(AA)/length(AA)^.5;'min',nanmin(AA);'max',nanmax(AA)})

AA=nan(1,length(aID));%w/ and w/o csd
for kk=1:length(aID)
    AA(kk)=aID(kk).SzEndTime-aID(kk).MeanSzRecTimeG;
end
disp({'length of sz (w/)',[];'mean',nanmean(AA);'std',nanstd(AA);'ste',nanstd(AA)/length(AA)^.5;'min',nanmin(AA);'max',nanmax(AA)})

BB=nan(1,length(bID));%w/ and w/o csd
for kk=1:length(bID)
    BB(kk)=bID(kk).SzEndTime-bID(kk).MeanSzRecTimeG;
end
disp({'length of sz (w/o)',[];'mean',nanmean(BB);'std',nanstd(BB);'ste',nanstd(BB)/length(BB)^.5;'min',nanmin(BB);'max',nanmax(BB)})

disp(['ranksum of w/ and w/o length p='])
disp(ranksum(AA,BB))

AA=nan(1,length(aID));%only w/ csd
for kk=1:length(aID)
    AA(kk)=aID(kk).MeanCSDRecTimeG-aID(kk).SzEndTime;
end
disp({'time sz to csd',[];'mean',nanmean(AA);'std',nanstd(AA);'ste',nanstd(AA)/length(AA)^.5;'min',nanmin(AA);'max',nanmax(AA)})


%%


%% Figure 2D(CSD) GLME for Calcium within subject (w/CSD)
% loads the calcium data for generating the glm table to determine
% %statistical significance

StructFieldNames={'FOV_dim','death','deathTime','firstSpkT','firstSpkI','fs_2p',...
    'fs_eeg','EEGbias','EEG_ts','F1_ts','EEG','PSD_F','PSD_P','DC','fs_dc','DC_ts',...
    'F1G','Fneu1G','stat','xCoor','yCoor','xCoorI','yCoorI','iscell','iscellI',...
    'Fb1G','Fneub1G','Fc1G','Fb1Gdff','Fneub1Gdff','Fc1Gdff','Fneub1Gdff_flt1',...
    'Fc1Gdff_flt1','Fneub1Gdff_flt1','Fc1Gdff_flt2','F1R','Fneu1R','Fb1R','Fneub1R',...
    'Fc1R','Fb1Rdff','Fneub1Rdff','Fc1Rdff','Fneub1Rdff_flt1','Fc1Rdff_flt1',...
    'Fneub1Rdff_flt1','Fc1Rdff_flt2','EEG_PIS_times2','EEG_PIS_times3',...
    'MeanCaG_PIS_times','MeanCaR_PIS_times','PISisRec','RecTimeG_PIS','RecTimeR_PIS',...
    'MeanSzRecTimeG','MeanSzRecTimeR','SzEndTime','StimStartTime','StimEndTime',...
    'MeanCSDRecTimeG','MeanCSDRecTimeR','CSDstartTimeDC','CSDendTimeDC',...
    'CSDoffTimeDC','MeanCSDendTimeG','MeanCSDendTimeR','RecTimeG_sz2','RecTimeR_sz2',...
    'RecTimeG_csd','RecTimeR_csd','RecTimeG_csdRTN','RecTimeR_csdRTN','isRecG_csdRTN','isRecR_csdRTN'};

szID={'S60ptz1','S60ptz2','S89ptz3','S102ptz3','S103ptz1','S105ptz3'};

struct_string='aID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii

fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
clear *ptz*


% store relevant data
%invert the r channel as all the values are determiend based upon an
%inverted signal



%%
[Ca_PIS,Ca_Sz,Ca_CSD]=deal(cell(length(aID),2));
% [bothIndex_Sz,bothIndex_CSD]=deal(cell(length(aID),1));
for kk=1:length(aID)
%     Ca_PIS{kk,1}=aID(kk).AvgPISCaChangeG(:,2);
%     Ca_PIS{kk,2}=-aID(kk).AvgPISCaChangeR(:,2);
%     temp=double(aID(kk).isRecruitedG);
%     temp(temp==0)=NaN;
%     Ca_Sz{kk,1}=aID(kk).SzMeanDffDiffG.*temp;
%     Ca_Sz{kk,2}=-aID(kk).SzMeanDffDiffR.*temp;
%     temp=double(aID(kk).isRecruitedG_csd);
%     temp(temp==0)=NaN;
%     Ca_CSD{kk,1}=aID(kk).CSDMeanDffDiffG.*temp;
%     temp=double(aID(kk).isRecruitedR_csd);
%     temp(temp==0)=NaN;
%     Ca_CSD{kk,2}=-aID(kk).CSDMeanDffDiffR.*temp;
    
    Ca_PIS{kk,1}=aID(kk).AvgPISCaChangeG(:,1);
    Ca_PIS{kk,2}=-aID(kk).AvgPISCaChangeR(:,1);
    Ca_Sz{kk,1}=aID(kk).SzMeanDffDiffG;
    Ca_Sz{kk,2}=-aID(kk).SzMeanDffDiffR;
    Ca_CSD{kk,1}=aID(kk).CSDMeanDffDiffG;
    Ca_CSD{kk,2}=-aID(kk).CSDMeanDffDiffR;
end
clear kk

%% GLME for difference between events (all data; treating subject as random effect) w/CSD
%uses dummy variaible coding on effects to model the categorical variable
%relative to the average of the data across categories, enabling pairwise
%comparisons and the estimation of real effects rather than referencial
%dummy encoding which only allows comparision to the first variable set as
%as a reference.
%random effects are assigned to subject with each cell being treated as a
%replicate value within that subject for each condition, with no repeated
%measurment factored in otherwise.
%Notes: Using all data and using cellID as random effect made little
%difference so not using it; also using only recruited cells for the data
%made the outcomes often more significant but perhaps slightly unusual
%(p<1E-50)

GLMtable=cell(1,length(aID));
varNames={'CaLevelG','CaLevelR','Event','CellID','Subject'};%Event 1:PIS 2:Sz 3:CSD
cc=0;
for kk=1:length(aID)
    tempCaG=[Ca_PIS{kk,1};Ca_Sz{kk,1};Ca_CSD{kk,1}];
    tempCaR=[Ca_PIS{kk,2};Ca_Sz{kk,2};Ca_CSD{kk,2}];
    %tempEv=[ones(length(Ca_PIS{kk,1}),1);2*ones(length(Ca_Sz{kk,1}),1);3*ones(length(Ca_CSD{kk,1}),1)];
    tempEv=[repmat("PIS",length(Ca_PIS{kk,1}),1);repmat("Sz",length(Ca_Sz{kk,1}),1);repmat("CSD",length(Ca_CSD{kk,1}),1)];
    tempCID=[cc+1:cc+length(Ca_PIS{kk,1})]';
    GLMtable{kk}=table(tempCaG,tempCaR,tempEv,[tempCID;tempCID;tempCID],kk*ones(length(tempCaG),1),'VariableNames',varNames);
    cc=cc+length(Ca_PIS{kk,1});
end
CaLevel_GLMdata=vertcat(GLMtable{:});

CaLevelG_GLMout=fitglme(CaLevel_GLMdata,'CaLevelG ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');
%CaLevelG_GLMout=fitglme(CaLevel_GLMdata,'CaLevelG ~ Event + (Event|Subject) + (Event|CellID)','Distribution','normal','DummyVarCoding','effects');
CaLevelR_GLMout=fitglme(CaLevel_GLMdata,'CaLevelR ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');
%CaLevelR_GLMout=fitglme(CaLevel_GLMdata,'CaLevelR ~ Event + (Event|Subject) + (Event|CellID)','Distribution','normal','DummyVarCoding','effects');

%Test for multiple comparisions 
%H is the contrast matrix coorepsonding to the coef outputs of model
%if using dummy variables coefs are [intercept PIS Sz]; PIS+Sz+CSD=0 therefore CSD=-PIS-Sz
Hmat=[1 0 0; %Intercept
      0 1 0; %PIS
      0 0 1; %Sz
      0 -1 -1; %CSD=-PIS-Sz
      0 1 -1; %PIS-Sz
      0 1 2; %Sz-CSD Sz-(-PIS-Sz)=2Sz+PIS
      0 2 1]; %PIS-CSD PIS-(-PIS-Sz)=2PIS+Sz
CaLevel_GLMpp=nan(size(Hmat,1),2);
for ii=1:size(Hmat,1)
    CaLevel_GLMpp(ii,1)=coefTest(CaLevelG_GLMout,Hmat(ii,:));
    CaLevel_GLMpp(ii,2)=coefTest(CaLevelR_GLMout,Hmat(ii,:));
end
 clear ii


 %% Figure 2D GLME plot (w/ CSD)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
[Ca_PIS_means,Ca_PIS_std,Ca_PIS_df,Ca_Sz_means,Ca_Sz_std,Ca_Sz_df,Ca_CSD_means,Ca_CSD_std,Ca_CSD_df,]=deal(nan(length(aID),2));
for kk=1:length(aID)
    Ca_PIS_means=cellfun(@nanmean,Ca_PIS);
    Ca_PIS_std=cellfun(@nanstd,Ca_PIS);
    Ca_PIS_df=cellfun(@length,Ca_PIS)-1;
    Ca_Sz_means=cellfun(@nanmean,Ca_Sz);
    Ca_Sz_std=cellfun(@nanstd,Ca_Sz);
    Ca_Sz_df=cellfun(@length,Ca_Sz)-1;
    Ca_CSD_means=cellfun(@nanmean,Ca_CSD);
    Ca_CSD_std=cellfun(@nanstd,Ca_CSD);
    Ca_CSD_df=cellfun(@length,Ca_CSD)-1;
end

% Pooled SE
Ca_PIS_stdp=sqrt(sum(Ca_PIS_df.*Ca_PIS_std.^2,1)./sum(Ca_PIS_df,1));
Ca_PIS_sep=Ca_PIS_stdp.*sqrt(sum(1./(Ca_PIS_df+1)));
Ca_Sz_stdp=sqrt(sum(Ca_Sz_df.*Ca_Sz_std.^2,1)./sum(Ca_Sz_df,1));
Ca_Sz_sep=Ca_Sz_stdp.*sqrt(sum(1./(Ca_Sz_df+1)));
Ca_CSD_stdp=sqrt(sum(Ca_CSD_df.*Ca_CSD_std.^2,1)./sum(Ca_CSD_df,1));
Ca_CSD_sep=Ca_CSD_stdp.*sqrt(sum(1./(Ca_CSD_df+1)));


figure(100)
set(gcf,'Position',[100 100 1000 500])
subplot(1,6,4:6)
plotYvals=[mean(Ca_PIS_means(:,1)), mean(Ca_Sz_means(:,1)), mean(Ca_CSD_means(:,1)); mean(Ca_PIS_means(:,2)), mean(Ca_Sz_means(:,2)), mean(Ca_CSD_means(:,2))]';
plotYerr=[Ca_PIS_sep(1), Ca_Sz_sep(1), Ca_CSD_sep(1); Ca_PIS_sep(2), Ca_Sz_sep(2), Ca_CSD_sep(2)]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.3 .6])
xticklabels({'PIS';'Sz';'CSD'})
set(gca,'box','off')


%% Figure 2D(Sz) GLME for Calcium within subject (w/o CSD)
%loads the calcium data for generating the glm table to determine
%statistical significance 

StructFieldNames={'FOV_dim','death','deathTime','firstSpkT','firstSpkI','fs_2p',...
    'fs_eeg','EEGbias','EEG_ts','F1_ts','EEG','PSD_F','PSD_P','DC','fs_dc','DC_ts',...
    'F1G','Fneu1G','stat','xCoor','yCoor','xCoorI','yCoorI','iscell','iscellI',...
    'Fb1G','Fneub1G','Fc1G','Fb1Gdff','Fneub1Gdff','Fc1Gdff','Fneub1Gdff_flt1',...
    'Fc1Gdff_flt1','Fneub1Gdff_flt1','Fc1Gdff_flt2','F1R','Fneu1R','Fb1R','Fneub1R',...
    'Fc1R','Fb1Rdff','Fneub1Rdff','Fc1Rdff','Fneub1Rdff_flt1','Fc1Rdff_flt1',...
    'Fneub1Rdff_flt1','Fc1Rdff_flt2','EEG_PIS_times2','EEG_PIS_times3',...
    'MeanCaG_PIS_times','MeanCaR_PIS_times','PISisRec','RecTimeG_PIS','RecTimeR_PIS',...
    'MeanSzRecTimeG','MeanSzRecTimeR','SzEndTime','StimStartTime','StimEndTime',...
    'MeanCSDRecTimeG','MeanCSDRecTimeR','CSDstartTimeDC','CSDendTimeDC',...
    'CSDoffTimeDC','MeanCSDendTimeG','MeanCSDendTimeR','RecTimeG_sz2','RecTimeR_sz2',...
    'RecTimeG_csd','RecTimeR_csd','RecTimeG_csdRTN','RecTimeR_csdRTN','isRecG_csdRTN','isRecR_csdRTN'};


%szID={'S89ptz2','S102ptz2','S105ptz1'};%Only the ones we have in the w/ CSD group
%szID={'S86ptz2','S86ptz3','S87ptz4','S87ptz5','S87ptz7','S89ptz2','S102ptz2','S104ptz1','S105ptz1'};%Not including the ones we stimulated
szID={'S86ptz3','S86ptz4','S86ptz5','S87ptz4','S87ptz5','S87ptz7','S89ptz2','S102ptz2','S104ptz1','S104ptz2','S105ptz1','S105ptz2'};
%szID={'S86ptz3','S86ptz4','S86ptz5','S89ptz2','S102ptz2','S104ptz1','S104ptz2','S105ptz1','S105ptz2'};%no s87 (movement artifact during seizure)
%szID={'S86ptz3','S89ptz2','S102ptz2','S104ptz1','S105ptz1'};%no s87 or stim

struct_string='aID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii
clear *ptz*

% store relevant data (w/o CSD)
%invert the r channel as all the values are determiend based upon an
%inverted signal
%%
[Ca_PIS,Ca_Sz]=deal(cell(length(aID),2));
for kk=1:length(aID)
%     Ca_PIS{kk,1}=aID(kk).AvgPISCaChangeG(:,2);
%     Ca_PIS{kk,2}=-aID(kk).AvgPISCaChangeR(:,2);
%     temp=double(aID(kk).isRecruitedG);
%     temp(temp==0)=NaN;
%     Ca_Sz{kk,1}=aID(kk).SzMeanDffDiffG.*temp;
%     Ca_Sz{kk,2}=-aID(kk).SzMeanDffDiffR.*temp;
    
    Ca_PIS{kk,1}=aID(kk).AvgPISCaChangeG(:,1);
    Ca_Sz{kk,1}=aID(kk).SzMeanDffDiffG;
    Ca_PIS{kk,2}=-aID(kk).AvgPISCaChangeR(:,1);
    Ca_Sz{kk,2}=-aID(kk).SzMeanDffDiffR;
end
clear kk

%% GLME for difference between events (all data; treating subject as random effect) (w/o CSD)
%uses dummy variable coding on effects to model the categorical variable
%relative to the average of the data across catagories, enabling pairwise
%comparisons and the estimation of real effects rather than referencial
%dummy recording which only allows comaprision to the first variable set as
%as a reference...however as there are only two options as categorical
%variables this is essentaily the same as doing referential coding
%random effects are assigned to subject with each cell being treated as a
%replicate value within that subject for each condition, with no repeted
%measurment factored in otherwise.

GLMtable=cell(1,length(aID));
varNames={'CaLevelG','CaLevelR','Event','CellID','Subject'};%Event 1:PIS 2:Sz 3:CSD
cc=0;
for kk=1:length(aID)
    tempCaG=[Ca_PIS{kk,1};Ca_Sz{kk,1}];
    tempCaR=[Ca_PIS{kk,2};Ca_Sz{kk,2}];
    tempEv=[repmat("PIS",length(Ca_PIS{kk,1}),1);repmat("Sz",length(Ca_Sz{kk,1}),1)];
    tempCID=[cc+1:cc+length(Ca_PIS{kk,1})]';
    GLMtable{kk}=table(tempCaG,tempCaR,tempEv,[tempCID;tempCID],kk*ones(length(tempCaG),1),'VariableNames',varNames);
    cc=cc+length(Ca_PIS{kk,1});
end
CaLevel_GLMdata=vertcat(GLMtable{:});

CaLevelG_GLMout=fitglme(CaLevel_GLMdata,'CaLevelG ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');
CaLevelR_GLMout=fitglme(CaLevel_GLMdata,'CaLevelR ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');

% Test for multiple comparisions (w/o CSD)
%H is the contrast matrix coorepsonding to the coef outputs of model
%if using dummy variables coefs are [intercept PIS Sz]; PIS+Sz+CSD=0 therefore CSD=-PIS-Sz

Hmat=[1 0 ; %Intercept
      0 1 ; %PIS
      0 -1; %Sz : PTZ + Sz = 0 Sz = -PTZ
      0 2]; %PIS-Sz PTZ-(-PTZ)
CaLevel_GLMpp=nan(size(Hmat,1),2);
for ii=1:size(Hmat,1)
    CaLevel_GLMpp(ii,1)=coefTest(CaLevelG_GLMout,Hmat(ii,:));
    CaLevel_GLMpp(ii,2)=coefTest(CaLevelR_GLMout,Hmat(ii,:));
end
 clear ii


 %% GLME plot (w/o CSD)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
[Ca_PIS_means,Ca_PIS_std,Ca_PIS_df,Ca_Sz_means,Ca_Sz_std,Ca_Sz_df]=deal(nan(length(aID),2));
for kk=1:length(aID)
    Ca_PIS_means=cellfun(@nanmean,Ca_PIS);
    Ca_PIS_std=cellfun(@nanstd,Ca_PIS);
    Ca_PIS_df=cellfun(@length,Ca_PIS)-1;
    Ca_Sz_means=cellfun(@nanmean,Ca_Sz);
    Ca_Sz_std=cellfun(@nanstd,Ca_Sz);
    Ca_Sz_df=cellfun(@length,Ca_Sz)-1;
end

% Pooled SE
Ca_PIS_stdp=sqrt(sum(Ca_PIS_df.*Ca_PIS_std.^2,1)./sum(Ca_PIS_df,1));
Ca_PIS_sep=Ca_PIS_stdp.*sqrt(sum(1./(Ca_PIS_df+1)));
Ca_Sz_stdp=sqrt(sum(Ca_Sz_df.*Ca_Sz_std.^2,1)./sum(Ca_Sz_df,1));
Ca_Sz_sep=Ca_Sz_stdp.*sqrt(sum(1./(Ca_Sz_df+1)));                           

figure(100)
subplot(1,6,2:3)
plotYvals=[mean(Ca_PIS_means(:,1)), mean(Ca_Sz_means(:,1)); mean(Ca_PIS_means(:,2)), mean(Ca_Sz_means(:,2))]';
plotYerr=[Ca_PIS_sep(1), Ca_Sz_sep(1); Ca_PIS_sep(2), Ca_Sz_sep(2)]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.3 .6])
xticklabels({'PIS';'Sz'})
set(gca,'box','off')
 
 
%% Figure 2 GLME for Calcium within subject (subgen)
%loads the calcium data for generating the glm table to determine
%statistical significance 

StructFieldNames={'FOV_dim','death','deathTime','firstSpkT','firstSpkI','fs_2p',...
    'fs_eeg','EEGbias','EEG_ts','F1_ts','EEG','PSD_F','PSD_P','DC','fs_dc','DC_ts',...
    'F1G','Fneu1G','stat','xCoor','yCoor','xCoorI','yCoorI','iscell','iscellI',...
    'Fb1G','Fneub1G','Fc1G','Fb1Gdff','Fneub1Gdff','Fc1Gdff','Fneub1Gdff_flt1',...
    'Fc1Gdff_flt1','Fneub1Gdff_flt1','Fc1Gdff_flt2','F1R','Fneu1R','Fb1R','Fneub1R',...
    'Fc1R','Fb1Rdff','Fneub1Rdff','Fc1Rdff','Fneub1Rdff_flt1','Fc1Rdff_flt1',...
    'Fneub1Rdff_flt1','Fc1Rdff_flt2','EEG_PIS_times2','EEG_PIS_times3',...
    'MeanCaG_PIS_times','MeanCaR_PIS_times','PISisRec','RecTimeG_PIS','RecTimeR_PIS',...
    'MeanSzRecTimeG','MeanSzRecTimeR','SzEndTime','StimStartTime','StimEndTime',...
    'MeanCSDRecTimeG','MeanCSDRecTimeR','CSDstartTimeDC','CSDendTimeDC',...
    'CSDoffTimeDC','MeanCSDendTimeG','MeanCSDendTimeR','RecTimeG_sz2','RecTimeR_sz2',...
    'RecTimeG_csd','RecTimeR_csd','RecTimeG_csdRTN','RecTimeR_csdRTN','isRecG_csdRTN','isRecR_csdRTN'};

%szID={'S89ptz1','S102ptz1'};
szID={'S86ptz1','S87ptz1','S87ptz2','S87ptz3','S87ptz6','S89ptz1','S102ptz1'};

struct_string='aID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii
clear *ptz*

% load('s89_221027_2_ptz1/S89ptz1.mat');
% load('s102_230522_ptz1/S102ptz1.mat');
% aID=[S89ptz1,S102ptz1]; %sz w/o CSD

% store relevant data (subgen)
%invert the r channel as all the values are determiend based upon an
%inverted signal

Ca_PIS=cell(length(aID),2);
for kk=1:length(aID)
    Ca_PIS{kk,1}=aID(kk).AvgPISCaChangeG(:,1);
    Ca_PIS{kk,2}=-aID(kk).AvgPISCaChangeR(:,1);
end
clear kk


%% GLME plot (subgen; no stats)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
[Ca_PIS_means,Ca_PIS_std,Ca_PIS_df]=deal(nan(length(aID),2));
for kk=1:length(aID)
    Ca_PIS_means=cellfun(@nanmean,Ca_PIS);
    Ca_PIS_std=cellfun(@nanstd,Ca_PIS);
    Ca_PIS_df=cellfun(@length,Ca_PIS)-1;
end

% Pooled SE
Ca_PIS_stdp=sqrt(sum(Ca_PIS_df.*Ca_PIS_std.^2,1)./sum(Ca_PIS_df,1));
Ca_PIS_sep=Ca_PIS_stdp.*sqrt(sum(1./(Ca_PIS_df+1)));


figure(100)
subplot(1,6,1)
plotYvals=[mean(Ca_PIS_means(:,1)); mean(Ca_PIS_means(:,2))];
plotYerr=[Ca_PIS_sep(1); Ca_PIS_sep(2)];
barPlot=bar(plotYvals,'grouped','FaceColor','flat');
barPlot.CData = [0 .7 0 ; 0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.3 .6])
xticks(1.5)
xticklabels({'SWD'})
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_CaLvl.svg','-dsvg')

%% Figure 2E: Post Ictal Ca Level
%loads the calcium data for generating the glm table to determine
%statistical significance

StructFieldNames={'FOV_dim','death','deathTime','firstSpkT','firstSpkI','fs_2p',...
    'fs_eeg','EEGbias','EEG_ts','F1_ts','EEG','PSD_F','PSD_P','DC','fs_dc','DC_ts',...
    'F1G','Fneu1G','stat','xCoor','yCoor','xCoorI','yCoorI','iscell','iscellI',...
    'Fb1G','Fneub1G','Fc1G','Fb1Gdff','Fneub1Gdff','Fc1Gdff','Fneub1Gdff_flt1',...
    'Fc1Gdff_flt1','Fneub1Gdff_flt1','Fc1Gdff_flt2','F1R','Fneu1R','Fb1R','Fneub1R',...
    'Fc1R','Fb1Rdff','Fneub1Rdff','Fc1Rdff','Fneub1Rdff_flt1','Fc1Rdff_flt1',...
    'Fneub1Rdff_flt1','Fc1Rdff_flt2','EEG_PIS_times2','EEG_PIS_times3',...
    'MeanCaG_PIS_times','MeanCaR_PIS_times','PISisRec','RecTimeG_PIS','RecTimeR_PIS',...
    'MeanSzRecTimeG','MeanSzRecTimeR','SzEndTime','StimStartTime','StimEndTime',...
    'MeanCSDRecTimeG','MeanCSDRecTimeR','CSDstartTimeDC','CSDendTimeDC',...
    'CSDoffTimeDC','MeanCSDendTimeG','MeanCSDendTimeR','RecTimeG_sz2','RecTimeR_sz2',...
    'RecTimeG_csd','RecTimeR_csd','RecTimeG_csdRTN','RecTimeR_csdRTN','isRecG_csdRTN','isRecR_csdRTN'};

%CSD
%szID={'S89ptz3','S102ptz3','S103ptz1','S105ptz3'};
szID={'S60ptz1','S60ptz2','S89ptz3','S102ptz3','S103ptz1','S105ptz3'};

struct_string='aID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii
clear *ptz*

%
%Sz only
%szID={'S86ptz3','S89ptz2','S104ptz1','S105ptz1'};%-only good recordings
%szID={'S86ptz3','S87ptz4','S87ptz5','S87ptz7','S89ptz2','S104ptz1','S105ptz1'};%Not including the ones we stimulated
szID={'S86ptz3','S86ptz4','S86ptz5','S87ptz4','S87ptz5','S87ptz7','S89ptz2','S104ptz1','S104ptz2','S105ptz1','S105ptz2'};% omit S102ptz2 due to split recording and s86ptz2 due to objective moving during recording

struct_string='bID=[';

for ii=1:length(szID)
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
PrmTbl=removevars(PrmTbl,{'szID'});
file_path=PrmTbl{szID{ii},'file_path'};
load(strcat(file_path{1,1},'/',szID{ii},'.mat'));
eval(strcat(szID{ii},'=rmfield(',szID{ii},',StructFieldNames);'));
struct_string=strcat(struct_string,szID{ii});
if ii<length(szID)
struct_string=strcat(struct_string,',');
end
end
struct_string=strcat(struct_string,']');
eval(struct_string);
clear ii
clear *ptz*

% load('s89_221101_ptz2/S89ptz2.mat');
% load('s105_230526/S105ptz1.mat'); 
%s102ptz2 does not have post ictal phase to be measured
%s86ptz2,s86ptz3,s87ptz4,s87ptz5,s87ptz7,s104ptz1)
% bID=[S89ptz2,S105ptz1]; %sz w/o CSD

fs_2p=30;
fs_eeg=2000;
fs_dc=2000;

% store relevant data (w/ & w/o CSD)
%invert the r channel as all the values are determiend based upon an
%inverted signal
%%
Ca_PostIct_CSD=cell(length(aID),2);
for kk=1:length(aID)
    Ca_PostIct_CSD{kk,1}=aID(kk).PostIctalCaGmax;
    Ca_PostIct_CSD{kk,2}=-aID(kk).PostIctalCaRmax;
    %Ca_PostIct_CSD{kk,1}=aID(kk).PostIctalCaGmax(aID(kk).isRecruitedG);
    %Ca_PostIct_CSD{kk,2}=-aID(kk).PostIctalCaRmax(aID(kk).isRecruitedG);
    %Ca_PostIct_CSD{kk,1}=aID(kk).PostIctalCaGmax(aID(kk).isRecruitedG_csd);
    %Ca_PostIct_CSD{kk,2}=-aID(kk).PostIctalCaRmax(aID(kk).isRecruitedG_csd);
end
clear kk


Ca_PostIct_noCSD=cell(length(bID),2);
for kk=1:length(bID)
    Ca_PostIct_noCSD{kk,1}=bID(kk).PostIctalCaGmax;
    Ca_PostIct_noCSD{kk,2}=-bID(kk).PostIctalCaRmax;
    %Ca_PostIct_noCSD{kk,1}=bID(kk).PostIctalCaGmax(bID(kk).isRecruitedG);
    %Ca_PostIct_noCSD{kk,2}=-bID(kk).PostIctalCaRmax(bID(kk).isRecruitedG);
end
clear kk


%% Fig 2E: GLME for difference between events (all data; treating subject as random effect) post ictal

GLMtable=cell(1,length(bID)+length(aID));
varNames={'CaLevelG','CaLevelR','Event','Subject'};%Event 1:PIS 2:Sz 3:CSD
for kk=1:length(bID)
    tempCaG=[Ca_PostIct_noCSD{kk,1}];
    tempCaR=[Ca_PostIct_noCSD{kk,2}];
    tempEv=repmat("noCSD",length(tempCaG),1);
    GLMtable{kk}=table(tempCaG,tempCaR,tempEv,kk*ones(length(tempCaG),1),'VariableNames',varNames);
end
clear kk
for kk=1:length(aID)
    tempCaG=[Ca_PostIct_CSD{kk,1}];
    tempCaR=[Ca_PostIct_CSD{kk,2}];
    tempEv=repmat("CSD",length(tempCaG),1);
    GLMtable{length(bID)+kk}=table(tempCaG,tempCaR,tempEv,(length(bID)+kk)*ones(length(tempCaG),1),'VariableNames',varNames);
end
CaLevel_GLMdata=vertcat(GLMtable{:});

CaLevelG_GLMout=fitglme(CaLevel_GLMdata,'CaLevelG ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');
CaLevelR_GLMout=fitglme(CaLevel_GLMdata,'CaLevelR ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');


% Fig 2E: Test for multiple comparisions post ictal
%H is the contrast matrix coorepsonding to the coef outputs of model
%if using dummy variables coefs are [intercept PIS Sz]; PIS+Sz+CSD=0 therefore CSD=-PIS-Sz

Hmat=[1 0 ; %Intercept
      0 1 ; %noCSD
      0 -1; %CSD : noCSD + CSD = 0 CSD = -noCSD
      0 2]; %noCSD-CSD noCSD-(-noCSD)
CaLevel_GLMpp=nan(size(Hmat,1),2);
for ii=1:size(Hmat,1)
    CaLevel_GLMpp(ii,1)=coefTest(CaLevelG_GLMout,Hmat(ii,:));
    CaLevel_GLMpp(ii,2)=coefTest(CaLevelR_GLMout,Hmat(ii,:));
end
 clear ii


 %% Fig 2E: GLME plot (w/ & w/o CSD)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
[Ca_PostIct_noCSD_means,Ca_PostIct_noCSD_std,Ca_PostIct_noCSD_df,Ca_PostIct_CSD_means,Ca_PostIct_CSD_std,Ca_PostIct_CSD_df]=deal(nan(length(aID),2));
for kk=1:length(aID)
    Ca_PostIct_noCSD_means=cellfun(@mean,Ca_PostIct_noCSD);
    Ca_PostIct_noCSD_std=cellfun(@std,Ca_PostIct_noCSD);
    Ca_PostIct_noCSD_df=cellfun(@length,Ca_PostIct_noCSD)-1;
    Ca_PostIct_CSD_means=cellfun(@mean,Ca_PostIct_CSD);
    Ca_PostIct_CSD_std=cellfun(@std,Ca_PostIct_CSD);
    Ca_PostIct_CSD_df=cellfun(@length,Ca_PostIct_CSD)-1;
end

% Pooled SE
Ca_PostIct_noCSD_stdp=sqrt(sum(Ca_PostIct_noCSD_df.*Ca_PostIct_noCSD_std.^2,1)./sum(Ca_PostIct_noCSD_df,1));
Ca_PostIct_noCSD_sep=Ca_PostIct_noCSD_stdp.*sqrt(sum(1./(Ca_PostIct_noCSD_df+1)));
Ca_PostIct_CSD_stdp=sqrt(sum(Ca_PostIct_CSD_df.*Ca_PostIct_CSD_std.^2,1)./sum(Ca_PostIct_CSD_df,1));
Ca_PostIct_CSD_sep=Ca_PostIct_CSD_stdp.*sqrt(sum(1./(Ca_PostIct_CSD_df+1)));


%Plots the max (G) and min(R) sustained calcium change in each post ictal
%period following  seizures to demsotnrate the post ictal calcium changes
%during CSD that do not occur without CSD.
%nned to invert R channel as it is the max of the inverted signal

figure ('Position', [100 100 300 500])
plotYvals=[mean(Ca_PostIct_noCSD_means(:,1)), mean(Ca_PostIct_CSD_means(:,1)); mean(Ca_PostIct_noCSD_means(:,2)), mean(Ca_PostIct_CSD_means(:,2))]';
plotYerr=[Ca_PostIct_noCSD_sep(1), Ca_PostIct_CSD_sep(1); Ca_PostIct_noCSD_sep(2), Ca_PostIct_CSD_sep(2)]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.5 .5])
xticklabels({'Sz w/o CSD';'Sz w/ CSD'})
title('Post-Ictal Ca')
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/Figure2/Figure2_PostIctCaLvl.svg','-dsvg')










%% Extra Panel: DC = CSD
load('s105_230609_ptz3/S105ptz3.mat');
aID=[S105ptz3,S105ptz3];
fs_2p=30;
fs_eeg=2000;
kk=1;


%% Mean Traces and EEG

figure
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+2.5,'color',[0.3 0.3 0.3]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1)-1,'color', [0 .75 0]);
ylim([-1.5 3.5])
xlim([105 380])
yticks([])
pbaspect([1 .6 1])
hold off

%scale bars
sb_x=30; %30seconds
sb_ymg=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[125 -1.4];
sb_oe=[125 1.4];
sb_od=[125 .4];

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars DC
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')


print('../../RCatchER_CSD_Paper/Figures/Figure2/S105ptz3_Traces.svg','-dsvg')


% Movie S1: Rasters
figure
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
% c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
% c1.Label.String='dF/F';
caxis([-1 3])
xticks([0:60:1050])
xticklabels([0:60:1050])
yticks([259:60:519])
yticklabels([0:60:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([105 380])

print('../../RCatchER_CSD_Paper/Figures/Figure2/S105ptz3_rasters.svg','-dsvg')










%% FIGURE 3: Vectors of Propogation
% Load Data
load('s60_210915_ptz1/S60ptz1.mat');
load('s60_211008_ptz2/S60ptz2.mat');
load('s89_221109_ptz3/S89ptz3.mat');
load('s102_230530_ptz2/S102ptz3.mat');
load('s103_230531_ptz1/S103ptz1.mat');
load('s105_230609_ptz3/S105ptz3.mat');
aID=[S60ptz1,S60ptz2,S89ptz3,S102ptz3,S103ptz1,S105ptz3];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
cMap1=colormap(lines);


%% Linear Regression
% Propogation Dynamics
[t_csdG,t_csdR,betaG,VG,pG,velocityG,thetaG,betaR,VR,pR,velocityR,thetaR,RelRec_csd,length_csd,latency_csdStart,latency_csdEnd,bothIndex,both2Index]=deal(cell(length(aID),2));
[roi_coordsG,roi_coordsR,length_csdDC,CSD_DCshift]=deal(cell(length(aID),1));

for kk=1:length(aID)

pxlen=aID(kk).FOV_dim/512;%um:px


roi_coordsG{kk} = [aID(kk).xCoor*pxlen, (aID(kk).yCoor)*pxlen];
roi_coordsR{kk} = roi_coordsG{kk}; 

%CSD Start
t_csdG{kk,1} = aID(kk).RecTimeG_csd;
t_csdR{kk,1} = aID(kk).RecTimeR_csd;

bothIndex{kk,1} = and(aID(kk).isRecruitedG_csd,aID(kk).isRecruitedR_csd);

tempt = num2cell(t_csdG{kk,1});
[betaG{kk,1}, VG{kk,1}, pG{kk,1}] = SpatialLinearRegression(tempt(aID(kk).isRecruitedG_csd),roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityG{kk,1}=sqrt(VG{kk,1}(1)^2+VG{kk,1}(2)^2);
if VG{kk,1}(1)>=0 && VG{kk,1}(2)>=0
    thetaG{kk,1}=atand(VG{kk,1}(2)/VG{kk,1}(1));
elseif VG{kk,1}(1)>=0 && VG{kk,1}(2)<0
    thetaG{kk,1}=360-atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
elseif VG{kk,1}(1)<0 && VG{kk,1}(2)<0
    thetaG{kk,1}=180+atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
else
    thetaG{kk,1}=180-atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
end

tempt = num2cell(t_csdR{kk,1});
[betaR{kk,1}, VR{kk,1}, pR{kk,1}] = SpatialLinearRegression(tempt(aID(kk).isRecruitedR_csd),roi_coordsR{kk}(aID(kk).isRecruitedR_csd,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityR{kk,1}=sqrt(VR{kk,1}(1)^2+VR{kk,1}(2)^2);
if VR{kk,1}(1)>=0 && VR{kk,1}(2)>=0
    thetaR{kk,1}=atand(VR{kk,1}(2)/VR{kk,1}(1));
elseif VR{kk,1}(1)>=0 && VR{kk,1}(2)<0
    thetaR{kk,1}=360-atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
elseif VR{kk,1}(1)<0 && VR{kk,1}(2)<0
    thetaR{kk,1}=180+atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
else
    thetaR{kk,1}=180-atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
end

RelRec_csd{kk,1}=t_csdR{kk,1}(bothIndex{kk,1})-t_csdG{kk,1}(bothIndex{kk,1});

%CSD End
t_csdG{kk,2} = aID(kk).RecTimeG_csdRTN;
t_csdR{kk,2} = aID(kk).RecTimeR_csdRTN;

bothIndex{kk,2} = and(aID(kk).isRecG_csdRTN,aID(kk).isRecR_csdRTN);

tempt = num2cell(t_csdG{kk,2});
[betaG{kk,2}, VG{kk,2}, pG{kk,2}] = SpatialLinearRegression(tempt(aID(kk).isRecG_csdRTN),roi_coordsG{kk}(aID(kk).isRecG_csdRTN,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityG{kk,2}=sqrt(VG{kk,2}(1)^2+VG{kk,2}(2)^2);
if VG{kk,2}(1)>=0 && VG{kk,2}(2)>=0
    thetaG{kk,2}=atand(VG{kk,2}(2)/VG{kk,2}(1));
elseif VG{kk,2}(1)>=0 && VG{kk,2}(2)<0
    thetaG{kk,2}=360-atand(abs(VG{kk,2}(2)/VG{kk,2}(1)));
elseif VG{kk,2}(1)<0 && VG{kk,2}(2)<0
    thetaG{kk,2}=180+atand(abs(VG{kk,2}(2)/VG{kk,2}(1)));
else
    thetaG{kk,2}=180-atand(abs(VG{kk,2}(2)/VG{kk,2}(1)));
end

tempt = num2cell(t_csdR{kk,2});
[betaR{kk,2}, VR{kk,2}, pR{kk,2}] = SpatialLinearRegression(tempt(aID(kk).isRecR_csdRTN),roi_coordsR{kk}(aID(kk).isRecR_csdRTN,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityR{kk,2}=sqrt(VR{kk,2}(1)^2+VR{kk,2}(2)^2);
if VR{kk,2}(1)>=0 && VR{kk,2}(2)>=0
    thetaR{kk,2}=atand(VR{kk,2}(2)/VR{kk,2}(1));
elseif VR{kk,2}(1)>=0 && VR{kk,2}(2)<0
    thetaR{kk,2}=360-atand(abs(VR{kk,2}(2)/VR{kk,2}(1)));
elseif VR{kk,2}(1)<0 && VR{kk,2}(2)<0
    thetaR{kk,2}=180+atand(abs(VR{kk,2}(2)/VR{kk,2}(1)));
else
    thetaR{kk,2}=180-atand(abs(VR{kk,2}(2)/VR{kk,2}(1)));
end

RelRec_csd{kk,2}=t_csdR{kk,2}(bothIndex{kk,2})-t_csdG{kk,2}(bothIndex{kk,2});

both2Index{kk,1}=and(aID(kk).isRecruitedG_csd,aID(kk).isRecG_csdRTN);%column 1 is G and 2 is R
length_csd{kk,1}=t_csdG{kk,2}(both2Index{kk,1})-t_csdG{kk,1}(both2Index{kk,1});

both2Index{kk,2}=and(aID(kk).isRecruitedR_csd,aID(kk).isRecR_csdRTN);%column 1 is G and 2 is R
length_csd{kk,2}=t_csdR{kk,2}(both2Index{kk,2})-t_csdR{kk,1}(both2Index{kk,2});

latency_csdStart{kk,1}=t_csdG{kk,1}(aID(kk).isRecruitedG_csd)-aID(kk).CSDstartTimeDC;
latency_csdEnd{kk,1}=t_csdG{kk,2}(aID(kk).isRecG_csdRTN)-aID(kk).CSDendTimeDC;
latency_csdStart{kk,2}=t_csdR{kk,1}(aID(kk).isRecruitedR_csd)-aID(kk).CSDstartTimeDC;
latency_csdEnd{kk,2}=t_csdR{kk,2}(aID(kk).isRecR_csdRTN)-aID(kk).CSDendTimeDC;

length_csdDC{kk}=aID(kk).CSDendTimeDC-aID(kk).CSDstartTimeDC;
CSD_DCshift{kk}=aID(kk).CSD_dcShift(1);

end


%% GLM for difference between groups

%CSD Start
GLMtable=cell(1,length(aID));
varNames={'RecTime','RCatchER','Subject'};
for kk=1:length(aID)
GLMtable{kk}=table([t_csdG{kk,1}(bothIndex{kk,1});t_csdR{kk,1}(bothIndex{kk,1})],[zeros(sum(bothIndex{kk,1}),1);ones(sum(bothIndex{kk,1}),1)],kk*ones(sum(bothIndex{kk,1})*2,1),'VariableNames',varNames);
end
RelRec_GLMdata=vertcat(GLMtable{:});

RelRec_GLMout=fitglme(RelRec_GLMdata,'RecTime ~ 1 + RCatchER + (RCatchER|Subject) + (1|Subject)','Distribution','normal');

%CSD End
GLMtable2=cell(1,length(aID));
varNames={'RecTime','RCatchER','Subject'};
for kk=1:length(aID)
GLMtable2{kk}=table([t_csdG{kk,2}(bothIndex{kk,2});t_csdR{kk,2}(bothIndex{kk,2})],[zeros(sum(bothIndex{kk,2}),1);ones(sum(bothIndex{kk,2}),1)],kk*ones(sum(bothIndex{kk,2})*2,1),'VariableNames',varNames);
end
RelRec_GLMdata2=vertcat(GLMtable2{:});

RelRec_GLMout2=fitglme(RelRec_GLMdata2,'RecTime ~ 1 + RCatchER + (RCatchER|Subject) + (1|Subject)','Distribution','normal');


%% Fig3G: bar chart of velocity
figure('position',[100 100 500 400])
plotYvals=[mean(cell2mat(velocityG)); mean(cell2mat(velocityR))]';
plotYerr=[std(cell2mat(velocityG))/sqrt(length(aID)); std(cell2mat(velocityR))/sqrt(length(aID))]';
barPlot=bar(plotYvals,'grouped','FaceColor','flat','FaceAlpha',0.2,'EdgeColor','flat','LineWidth',1);
barPlot(1).FaceColor=[0 .9 0];
barPlot(2).FaceColor=[0.9 0 .9];
barPlot(1).EdgeColor=[0 .9 0];
barPlot(2).EdgeColor=[0.9 0 .9];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
scatter(ones(length(aID),1)*errXvals(1,1),cell2mat(velocityG(:,1)),40,cMap1(1:length(aID),:),'filled')
scatter(ones(length(aID),1)*errXvals(2,1),cell2mat(velocityR(:,1)),40,cMap1(1:length(aID),:),'filled')
scatter(ones(length(aID),1)*errXvals(1,2),cell2mat(velocityG(:,2)),40,cMap1(1:length(aID),:),'filled')
scatter(ones(length(aID),1)*errXvals(2,2),cell2mat(velocityR(:,2)),40,cMap1(1:length(aID),:),'filled')
hold off
%ylim([-.5 .5])
xticklabels({'CSD Start';'CSD End'})
title('CSD Wavefront Speed')
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_Velocity.svg','-dsvg')

%% Fig3H: Boxplot of theta difference between channels
thetaDiff=abs(cell2mat(thetaG(:,1))-cell2mat(thetaR(:,1)));
thetaDiff(thetaDiff>180)=abs(thetaDiff(thetaDiff>180)-360);
figure('position',[100 100 325 400])
boxchart(ones(length(aID),1),thetaDiff(:,1),'BoxFaceColor',[0.25 0.25 0.25]);
hold on
scatter(ones(length(aID),1),thetaDiff(:,1),[],[1:1:length(aID)],'filled')
boxchart(ones(length(aID),1)*2,thetaDiff(:,2),'BoxFaceColor',[0.75 0.75 0.75]);
scatter(ones(length(aID),1)*2,thetaDiff(:,2),[],[1:1:length(aID)],'filled')
colormap(cMap1)
xlim([0 3])
ylim([0 90])
%yticks([-45:5:45])
%yticklabels([-45:5:45])
xticklabels({'CSD Start';'CSD End'})
title('theta b/n channels')
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_ThetaBnCh.svg','-dsvg')


%% Fig3I: Boxplot of theta difference between start and end by channel
thetaDiffG=abs(diff(cell2mat(thetaG),1,2));
thetaDiffG=abs(thetaDiffG-360*(thetaDiffG>180));
thetaDiffR=abs(diff(cell2mat(thetaR),1,2));
thetaDiffR=abs(thetaDiffR-360*(thetaDiffR>180));
figure('position',[100 100 325 400])
boxchart(ones(length(aID),1),thetaDiffG,'BoxFaceColor',[0 1 0]);
hold on
scatter(ones(length(aID),1),thetaDiffG,[],[1:1:length(aID)],'filled')
boxchart(ones(length(aID),1)*2,thetaDiffR,'BoxFaceColor',[1 0 1]);
scatter(ones(length(aID),1)*2,thetaDiffR,[],[1:1:length(aID)],'filled')
colormap(cMap1)
xlim([0 3])
ylim([0 180])
%yticks([-45:5:45])
%yticklabels([-45:5:45])
title('theta b/n start and end')
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_ThetaBnOnOff.svg','-dsvg')


%% Fig3F: Box plots of relative recruitment (only Start)
%end is not performed as it is arbitrarily based on selection point of curve
figure('position',[100 100 350 400])
hold on
for ii=1:length(aID)
boxchart(ones(size(RelRec_csd{ii,1}))*ii,RelRec_csd{ii,1},'BoxFaceColor',cMap1(ii,:),'MarkerColor','none');
end
line([0,7],[0,0],'color',[.5 0 0],'LineStyle','--')
xlim([0 7])
ylim([-6 6])
title('Relative Recruitment by Seizure')
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_RelRecIndiv.svg','-dsvg')


%% plot of population relative recruitment vertical
clear hh* clear pp*

RelRec_csdMean=cellfun(@nanmean,RelRec_csd);
RelRec_csd_std=cellfun(@nanstd,RelRec_csd);
RelRec_csd_df=cellfun(@length,RelRec_csd)-1;

% Pooled SE
RelRec_csd_stdp=sqrt(sum(RelRec_csd_df.*RelRec_csd_std.^2,1)./sum(RelRec_csd_df,1));
RelRec_csd_sep=RelRec_csd_stdp.*sqrt(sum(1./(RelRec_csd_df+1)));
    
figure('position',[100 100 140 400])
barPlot=bar(1,mean(RelRec_csdMean(:,1)),'FaceColor','flat','FaceAlpha',0.2,'EdgeColor','flat','LineWidth',1);
barPlot.CData = [0.25 .25 0.25];
hold on
scatter(ones(length(aID),1),RelRec_csdMean(:,1),[],[1:1:length(aID)],'filled')
errorbar(1,mean(RelRec_csdMean(:,1)),RelRec_csd_sep(:,1),RelRec_csd_sep(:,1),'color',[0.25 0.25 0.25])
colormap(cMap1)
xlim([0 2])
ylim([-4 4])
title('Relative Recruitment')
set(gca,'box','off','xcolor','none')


print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_RelRecGrp.svg','-dsvg')


%% Fig3E: CSD length (pooled means and indiv sz plots; no stats)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
length_csdMean=cellfun(@nanmean,length_csd);
length_csd_std=cellfun(@nanstd,length_csd);
length_csd_df=cellfun(@length,length_csd)-1;
length_csd_GrandMean=mean(length_csdMean,1);

% Pooled SE
length_csd_stdp=sqrt(sum(length_csd_df.*length_csd_std.^2,1)./sum(length_csd_df,1));
length_csd_sep=length_csd_stdp.*sqrt(sum(1./(length_csd_df+1)));

%figure
figure('position',[100 100 250 400])
plotYvals=length_csd_GrandMean';
plotYerr=length_csd_sep';
barPlot=bar(plotYvals,'grouped','FaceColor','flat','FaceAlpha',0.2,'EdgeColor','flat','LineWidth',1);
barPlot.CData = [0 .7 0 ; .7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
title('CSD Length')

% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25],'linewidth',1.25);
scatter(ones(length(aID),1),length_csdMean(:,1),40,cMap1(1:length(aID),:),'filled')
colormap(cMap1)
scatter(ones(length(aID),1)*2,length_csdMean(:,2),40,cMap1(1:length(aID),:),'filled')
colormap(cMap1)
hold off
%ylim([-.25 .75])
set(gca,'box','off')
print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_CSDLength.svg','-dsvg')

figure('position',[100 100 650 400])
hold on
for ii=1:length(aID)
boxchart(ones(size(length_csd{ii,1}))*ii*2-1,length_csd{ii,1},'BoxFaceColor',[0 .7 0],'MarkerColor','none');
boxchart(ones(size(length_csd{ii,2}))*ii*2,length_csd{ii,2},'BoxFaceColor',[.7 0 .7],'MarkerColor','none');
text(ii*2-.5,85,'seizure '+string(ii),'color',cMap1(ii,:),'HorizontalAlignment','center')
end
xlim([0 13])
ylim([30 90])
title('CSD Length by Seizure')
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_CSDLength_Indiv.svg','-dsvg')

%% Fig3B-D: Rainbow ROI Field (actual times axis for red)
kk=4;

% [RecTimeGsortReG2,RecTimeGReG_I2]=sort(aID(kk).RecTimeG_csd); %sorts cells by rec times according to green
% [RecTimeRsortReR2,RecTimeRReR_I2]=sort(aID(kk).RecTimeR_csd); %sorts cells by rec times according to red
% [RecTimeGsortReG2_end,RecTimeGReG_I2_end]=sort(aID(kk).RecTimeG_csdRTN);
% [RecTimeRsortReR2_end,RecTimeRReR_I2_end]=sort(aID(kk).RecTimeR_csdRTN);

minOnT=min([aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd);aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd)]);
maxOnT=max([aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd);aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd)]);
minOffT=min([aID(kk).RecTimeG_csdRTN(aID(kk).isRecG_csdRTN);aID(kk).RecTimeR_csdRTN(aID(kk).isRecR_csdRTN)]);
maxOffT=max([aID(kk).RecTimeG_csdRTN(aID(kk).isRecG_csdRTN);aID(kk).RecTimeR_csdRTN(aID(kk).isRecR_csdRTN)]);

figure('position',[100 100 600 600])
ax(1)=subplot(2,2,1);
xx = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,1);
yy = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,2);
scatter(xx,yy,[],aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd),'filled')
% xx = roi_coordsG{kk}(RecTimeGReG_I2(aID(kk).isRecruitedG_csd),1); %NOTE this section is the old method where there was a mismatch between coordiantes and times 
% yy = roi_coordsG{kk}(RecTimeGReG_I2(aID(kk).isRecruitedG_csd),2);
% scatter(xx,yy,[],RecTimeGsortReG2(aID(kk).isRecruitedG_csd),'filled')
hold on
line([10,60],[15,15],'color',[0 0 0],'LineStyle','-','LineWidth',2) %scale bar 50um
hold off
colormap(ax(1),parula)
%colorbar('location','southoutside','Ticks',[],'TickLabels',{})
colorbar('location','southoutside')
caxis([floor(minOnT) ceil(maxOnT)])
pbaspect([1 1 1])
%title('Non-vGAT (jYCaMP1s)','fontsize',12)
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

ax(2)=subplot(2,2,2);
xx = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,1);
yy = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,2);
scatter(xx,yy,[],aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd),'filled')
colormap(ax(2),parula)
colorbar('location','southoutside','Ticks',[],'TickLabels',{})
caxis([floor(minOnT) ceil(maxOnT)])
pbaspect([1 1 1])
%title('VGAT (jRGECO1a)','fontsize',12)
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

ax(3)=subplot(2,2,3);
xx = roi_coordsG{kk}(aID(kk).isRecG_csdRTN,1);
yy = roi_coordsG{kk}(aID(kk).isRecG_csdRTN,2);
scatter(xx,yy,[],aID(kk).RecTimeG_csdRTN(aID(kk).isRecG_csdRTN),'filled')
colormap(ax(3),parula)
colorbar('location','southoutside')
caxis([5*floor(minOffT/5) 5*ceil(maxOffT/5)])
pbaspect([1 1 1])
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

ax(4)=subplot(2,2,4);
xx = roi_coordsR{kk}(aID(kk).isRecR_csdRTN,1);
yy = roi_coordsR{kk}(aID(kk).isRecR_csdRTN,2);
scatter(xx,yy,[],aID(kk).RecTimeR_csdRTN(aID(kk).isRecR_csdRTN),'filled')
colormap(ax(4),parula)
colorbar('location','southoutside','Ticks',[],'TickLabels',{})
caxis([5*floor(minOffT/5) 5*ceil(maxOffT/5)])
pbaspect([1 1 1])
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_ROIcMap.svg','-dsvg')

%Vector Plots (using compass)

pthresh=0.05;
max_lim = 100;
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];

figure('position',[100 100 400 400])
h_fake=compass(x_fake,y_fake);
hold on;
c1=compass(VG{kk,1}(1),VG{kk,1}(2));
c2=compass(VR{kk,1}(1),VR{kk,1}(2));
c3=compass(VG{kk,2}(1),VG{kk,2}(2));
c4=compass(VR{kk,2}(1),VR{kk,2}(2));

set(h_fake,'Visible','off');


if pG{kk,1}>=pthresh
   c1.Color=[0.9 1 0.9];
else
   c1.Color=[0 0.75 0];
end
c1.LineWidth=4;
c1.LineJoin='chamfer';
c1.LineStyle='-';
if pR{kk,1}>=pthresh
   c2.Color=[1 0.9 1];
else
   c2.Color=[0.75 0 0.75];
end
c2.LineWidth=4;
c2.LineJoin='chamfer';
c2.LineStyle='-';

if pG{kk,2}>=pthresh
   c3.Color=[0.9 1 0.9];
else
   c3.Color=[0 0.75 0];
end
c3.LineWidth=4;
c3.LineJoin='round';
c3.LineStyle=':';

if pR{kk,2}>=pthresh
   c4.Color=[1 0.9 1];
else
   c4.Color=[0.75 0 0.75];
end
c4.LineWidth=4;
c4.LineJoin='round';
c4.LineStyle=':';

hold off

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_Polar.svg','-dsvg')


%Propogation Vector Recruitment Times

%order of recruitments by subtype projected onto vector of propogation in
%green

figure('position',[100 100 600 400])
subplot(1,2,1)
P_v_axisG = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:)*betaG{kk,1}(1:2)/norm(betaG{kk,1}(1:2));
P_v_axisG = P_v_axisG - min(P_v_axisG);
P_v_axisR = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,:)*betaG{kk,1}(1:2)/norm(betaG{kk,1}(1:2));
P_v_axisR = P_v_axisR - min(P_v_axisR);
plot(t_csdG{kk,1}(aID(kk).isRecruitedG_csd),P_v_axisG,'.','MarkerSize',11,'color',[0 0.75 0]);
hold on
plot(t_csdR{kk,1}(aID(kk).isRecruitedR_csd),P_v_axisR,'.','MarkerSize',11,'color',[0.75 0 0.75]);
plot(roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:)*betaG{kk,1}(1:2)+betaG{kk,1}(end),P_v_axisG,'color',[0 0 .5]);
%ylabel('um');
%xlabel('seconds')
%xlim([750 780])
ax2 = gca;
ax2.FontSize = 12;
pbaspect([1 1.5 1])
box off
set(gca,'Visible','on')
Xtimes=[t_csdG{kk,1}(aID(kk).isRecruitedG_csd)',t_csdR{kk,1}(aID(kk).isRecruitedR_csd)'];
xlim([min(Xtimes)-2 max(Xtimes)+2])
ylim([0 500])
title('CSD Start')
clear Xtimes

subplot(1,2,2)
P_v_axisG = roi_coordsG{kk}(aID(kk).isRecG_csdRTN,:)*betaG{kk,2}(1:2)/norm(betaG{kk,2}(1:2));
P_v_axisG = P_v_axisG - min(P_v_axisG);
P_v_axisR = roi_coordsR{kk}(aID(kk).isRecR_csdRTN,:)*betaG{kk,2}(1:2)/norm(betaG{kk,2}(1:2));
P_v_axisR = P_v_axisR - min(P_v_axisR);
plot(t_csdG{kk,2}(aID(kk).isRecG_csdRTN),P_v_axisG,'.','MarkerSize',11,'color',[0 0.75 0]);
hold on
plot(t_csdR{kk,2}(aID(kk).isRecR_csdRTN),P_v_axisR,'.','MarkerSize',11,'color',[0.75 0 0.75]);
plot(roi_coordsG{kk}(aID(kk).isRecG_csdRTN,:)*betaG{kk,2}(1:2)+betaG{kk,2}(end),P_v_axisG,'color',[0 0 .5]);
%ylabel('um');
%xlabel('seconds')
%xlim([750 780])
ax2 = gca;
ax2.FontSize = 12;
pbaspect([1 1.5 1])
box off
set(gca,'Visible','on')
Xtimes=[t_csdG{kk,2}(aID(kk).isRecG_csdRTN)',t_csdR{kk,2}(aID(kk).isRecR_csdRTN)'];
xlim([min(Xtimes)-2 max(Xtimes)+2])
ylim([0 450])
title('CSD End')
clear Xtimes

print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_RelRecAxis.svg','-dsvg')



%% Figure S3: All CSD Vectors
max_lim = 100;
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];

figure('Position',[100,100,1000,1000])

for kk=1:6
subplot(2,2,1)
h_fake=compass(x_fake,y_fake);
hold on;
c1=compass(VG{kk,1}(1),VG{kk,1}(2));
c1.Color=cMap1(kk,:);
c1.LineWidth=4;
c1.LineJoin='chamfer';
c1.LineStyle='-';
set(h_fake,'Visible','off');

subplot(2,2,2)
h_fake=compass(x_fake,y_fake);
hold on;
c2=compass(VR{kk,1}(1),VR{kk,1}(2));
c2.Color=cMap1(kk,:);
c2.LineWidth=4;
c2.LineJoin='chamfer';
c2.LineStyle='-';
set(h_fake,'Visible','off');

subplot(2,2,3)
h_fake=compass(x_fake,y_fake);
hold on;
c3=compass(VG{kk,2}(1),VG{kk,2}(2));
c3.Color=cMap1(kk,:);
c3.LineWidth=4;
c3.LineJoin='chamfer';
c3.LineStyle='-';
set(h_fake,'Visible','off');

subplot(2,2,4)
h_fake=compass(x_fake,y_fake);
hold on;
c4=compass(VR{kk,2}(1),VR{kk,2}(2));
c4.Color=cMap1(kk,:);
c4.LineWidth=4;
c4.LineJoin='chamfer';
c4.LineStyle='-';
set(h_fake,'Visible','off');

end

print('../../RCatchER_CSD_Paper/Figures/FigureS3/FigureS3_PolarPlots.svg','-dsvg')




%% Fig3A: Individual traces with detection during CSD
%For the code used to figure out which traces to use see below in archive
kk=4;

%green
[Fb1Gdff_flt2,Fb1Rdff_flt2,Fc1Gdff_flt3,Fc1Rdff_flt3]=deal(zeros(size(aID(kk).Fc1Gdff)));
for ii=1:size(aID(kk).Fc1Gdff,1)
    Fb1Gdff_flt2(ii,:)=lofi(aID(kk).Fb1Gdff(ii,:),10^6/fs_2p,1,'verbose',0);
    Fb1Rdff_flt2(ii,:)=lofi(aID(kk).Fb1Rdff(ii,:),10^6/fs_2p,1,'verbose',0);
    %Fc1Gdff_flt3(ii,:)=lofi(aID(kk).Fc1Gdff(ii,:),10^6/fs_2p,.5,'verbose',0);
    %Fc1Rdff_flt3(ii,:)=lofi(aID(kk).Fc1Rdff(ii,:),10^6/fs_2p,.5,'verbose',0);
end

clear lpfiltG ii


%%
cellList=[84,156,290,340];
figure('Position',[100 100 500 800])
hold on
cc=1;

for ii=cellList %NOTE FLIPPED CELL PLOTTING ORDER IN ADOBE
    plot(aID(kk).F1_ts,Fb1Rdff_flt2(ii,:)+cc,'color', [1 0 1]);
    plot(aID(kk).F1_ts,Fb1Gdff_flt2(ii,:)*.5+cc+.85,'color', [0 .75 0]);
    %plot(aID(kk).F1_ts,aID(kk).Fc1Rdff_flt2(ii,:)+cc,'color', [1 0 1]);
    %plot(aID(kk).F1_ts,aID(kk).Fc1Gdff_flt2(ii,:)*.5+cc+.85,'color', [0 .75 0]);
    %plot(aID(kk).F1_ts,Fc1Rdff_flt3(ii,:)/2+cc,'color', [.9 0 .9]);
    %plot(aID(kk).F1_ts,Fc1Gdff_flt3(ii,:)+cc+1,'color', [0 .9 0]);
    
    line([aID(kk).RecTimeG_csd(ii);aID(kk).RecTimeG_csd(ii)],repmat([cc+.85;cc+1.1],1,2),'color',[0.7 .3 0],'LineWidth',2)
    line([aID(kk).RecTimeG_csdRTN(ii);aID(kk).RecTimeG_csdRTN(ii)],repmat([cc+.85;cc+1.1],1,2),'color',[0 0 .8],'LineWidth',2)
    
    line([aID(kk).RecTimeR_csd(ii);aID(kk).RecTimeR_csd(ii)],repmat([cc;cc+.25],1,2),'color',[0.7 .3 0],'LineWidth',2)
    line([aID(kk).RecTimeR_csdRTN(ii);aID(kk).RecTimeR_csdRTN(ii)],repmat([cc;cc+.25],1,2),'color',[0 0 .8],'LineWidth',2)
    
    cc=cc+2;
end

plot(aID(kk).EEG_ts,aID(kk).EEG./2000+cc+2.2,'color',[0.5 0.5 0.5]);
plot(aID(kk).DC_ts,aID(kk).DC./20+cc+1,'color',[0 0 0]);
xlim([10 260])
ylim([0 cc+3])

%scale bars
sb_x=30; %30seconds
sb_ymg=.25*.75; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .25uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[25 1.5];
sb_omr=[25 .5];
sb_oe=[25 9.9];
sb_od=[25 9]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

%print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_IndivTraces.svg','-dsvg')

figure
scatter(aID(kk).xCoor,aID(kk).yCoor,100,[0.9 0.9 0.9])
hold on
scatter(aID(kk).xCoor(cellList),aID(kk).yCoor(cellList),100,cMap1(1:length(cellList),:))
text(aID(kk).xCoor(cellList)+10,aID(kk).yCoor(cellList),string(1:4))
xlim([0 512])
ylim([0 512])
pbaspect([1 1 1])
yticks([])
xticks([])
box on
set(gca,'Visible','on')

%print('../../RCatchER_CSD_Paper/Figures/Figure3/Figure3_ROIs.svg','-dsvg')


%% Fig 3 Additional Relevant Values (DC Shift mag, length)
CSDlengthDC_mean=nanmean(cell2mat(length_csdDC));
CSDlengthDC_std=nanstd(cell2mat(length_csdDC));
CSDlengthDC_ste=CSDlengthDC_std/sqrt(sum(~isnan(cell2mat(length_csdDC))));

DCshift_mean=nanmean(cell2mat(CSD_DCshift));
DCshift_std=nanstd(cell2mat(CSD_DCshift));
DCshift_ste=DCshift_std/sqrt(sum(~isnan(cell2mat(CSD_DCshift))));

latency_csdStart_mean=cellfun(@nanmean,latency_csdStart);
latency_csdStart_std=cellfun(@nanstd,latency_csdStart);
latency_csdStart_df=cellfun(@length,latency_csdStart)-1;
latency_csdEnd_mean=cellfun(@nanmean,latency_csdEnd);
latency_csdEnd_std=cellfun(@nanstd,latency_csdEnd);
latency_csdEnd_df=cellfun(@length,latency_csdEnd)-1;

%Mean of Means
latency_csdStart_GrandMean=nanmean(latency_csdStart_mean,1);
latency_csdEnd_GrandMean=nanmean(latency_csdEnd_mean,1);

% Pooled Ste
latency_csdStart_stdp=sqrt(nansum(latency_csdStart_df.*latency_csdStart_std.^2,1)./nansum(latency_csdStart_df,1));
latency_csdEnd_stdp=sqrt(nansum(latency_csdEnd_df.*latency_csdEnd_std.^2,1)./nansum(latency_csdEnd_df,1));
latency_csdStart_sep=latency_csdStart_stdp.*sqrt(sum(1./(latency_csdStart_df+1)));
latency_csdEnd_sep=latency_csdEnd_stdp.*sqrt(sum(1./(latency_csdEnd_df+1)));









%% Figure 4: Terminal Spreading Depolarizations

%***Note for this figure the data is from no neuropil subtracted traces and filt1:.5hz and filt2:.1hz*****

%s56 no DC or EEG
%s86 looks bad re traces
%s87 has motion artifact
%s104 may be the best bet

load('s56_210831_ptz1/S56ptz1.mat');
load('s86_230104_ptz6/S86ptz6.mat');
load('s87_221219_ptz8/S87ptz8.mat');
load('s104_230608_ptz3/S104ptz3.mat');
aID=[S56ptz1,S86ptz6,S87ptz8,S104ptz3];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
cMap1=colormap(lines);


%% Fig4A: Traces (s104)
%lofi(aID(kk).Fc1Gdff_Mean,10^6/fs_2p,lpfilt,'verbose',0)
%Seizure w/ Terminal CSD
figure
kk=4;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]); %flt1 and flt2 are both filtered at 1hz.
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([200 700])
yticks([])
pbaspect([1 .5 1])
hold off

%scale bars
sb_x=30; %30seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[225 -0.5];
sb_omr=[225 -1.3];
sb_oe=[225 2];
sb_od=[225 1.1]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_104_TracesTSD.svg','-dsvg')

%Sz w/ CSD
figure
kk=4;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
xticklabels([0:50:1050])
yticks([35:10:69])
yticklabels([0:10:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([5 1 1])
xlim([200 700])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
xticklabels([0:50:1050])
yticks([35:10:69])
yticklabels([0:10:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([5 1 1])
xlim([200 700])

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_104RasterTSD.svg','-dsvg')


%% Fig4B: Ca Levels during seizure (summary data)

% store relevant data (CSD)
%invert the r channel as all the values are determined based upon an
%inverted signal

%[Ca_Sz,Ca_CSD]=deal(cell(length(aID),2));
[Ca_Sz,Ca_CSD]=deal(cell(3,2));
%for kk=1:length(aID)
cc=1;
for kk=[1 2 4]%omit s87 movement artifact
    Ca_Sz{cc,1}=aID(kk).SzMeanDffDiffG;
    Ca_Sz{cc,2}=-aID(kk).SzMeanDffDiffR;
    Ca_CSD{cc,1}=aID(kk).CSDMeanDffDiffG;
    Ca_CSD{cc,2}=-aID(kk).CSDMeanDffDiffR;
    cc=cc+1;
end
clear kk


%% GLME for difference between events (all data; treating subject as random effect)

%GLMtable=cell(1,length(aID));
GLMtable=cell(1,3);
varNames={'CaLevelG','CaLevelR','Event','CellID','Subject'};%Event 1:PIS 2:Sz 3:CSD
cc=0;
%for kk=1:length(aID)
for kk = 1:3 %omit s87 movement artifact
    tempCaG=[Ca_Sz{kk,1};Ca_CSD{kk,1}];
    tempCaR=[Ca_Sz{kk,2};Ca_CSD{kk,2}];
    tempEv=[repmat("Sz",length(Ca_Sz{kk,1}),1);repmat("CSD",length(Ca_CSD{kk,1}),1)];
    tempCID=[cc+1:cc+length(Ca_Sz{kk,1})]';
    GLMtable{kk}=table(tempCaG,tempCaR,tempEv,[tempCID;tempCID],kk*ones(length(tempCaG),1),'VariableNames',varNames);
    cc=cc+length(Ca_Sz{kk,1});
end
CaLevel_GLMdata=vertcat(GLMtable{:});

CaLevelG_GLMout=fitglme(CaLevel_GLMdata,'CaLevelG ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');
CaLevelR_GLMout=fitglme(CaLevel_GLMdata,'CaLevelR ~ 1 + Event + (Event|Subject) + (1|Subject)','Distribution','normal','DummyVarCoding','effects');

% Test for multiple comparisions (w/o CSD)
%H is the contrast matrix coorepsonding to the coef outputs of model
%if using dummy variables coefs are [intercept PIS Sz]; PIS+Sz+CSD=0 therefore CSD=-PIS-Sz

Hmat=[1 0 ; %Intercept
      0 1 ; %Sz
      0 -1; %CSD : Sz + CSD = 0 CSD = -Sz
      0 2]; %Sz-CSD Sz-(-Sz)
CaLevel_GLMpp=nan(size(Hmat,1),2);
for ii=1:size(Hmat,1)
    CaLevel_GLMpp(ii,1)=coefTest(CaLevelG_GLMout,Hmat(ii,:));
    CaLevel_GLMpp(ii,2)=coefTest(CaLevelR_GLMout,Hmat(ii,:));
end
 clear ii
 

%% GLME plot (w/o CSD)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate

Ca_Sz_means=cellfun(@nanmean,Ca_Sz);
Ca_Sz_std=cellfun(@nanstd,Ca_Sz);
Ca_Sz_df=cellfun(@length,Ca_Sz)-1;
Ca_CSD_means=cellfun(@nanmean,Ca_CSD);
Ca_CSD_std=cellfun(@nanstd,Ca_CSD);
Ca_CSD_df=cellfun(@length,Ca_CSD)-1;

% Pooled SE
Ca_Sz_stdp=sqrt(sum(Ca_Sz_df.*Ca_Sz_std.^2,1)./sum(Ca_Sz_df,1));
Ca_Sz_sep=Ca_Sz_stdp.*sqrt(sum(1./(Ca_Sz_df+1)));
Ca_CSD_stdp=sqrt(sum(Ca_CSD_df.*Ca_CSD_std.^2,1)./sum(Ca_CSD_df,1));
Ca_CSD_sep=Ca_CSD_stdp.*sqrt(sum(1./(Ca_CSD_df+1)));                           

figure('position',[100 100 450 400])
plotYvals=[mean(Ca_Sz_means(:,1)), mean(Ca_CSD_means(:,1)); mean(Ca_Sz_means(:,2)), mean(Ca_CSD_means(:,2))]';
plotYerr=[Ca_Sz_sep(1), Ca_CSD_sep(1); Ca_Sz_sep(2), Ca_CSD_sep(2)]';
barPlot=bar(plotYvals);
barPlot(1).FaceColor=[0 .7 0];
barPlot(2).FaceColor=[0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.3 .6])
xticklabels({'Sz';'CSD'})
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_CaLvl.svg','-dsvg')


%% Linear Regression for velocity fig 4
% Propogation Dynamics
[t_csdG,t_csdR,betaG,VG,pG,velocityG,thetaG,betaR,VR,pR,velocityR,thetaR,bothIndex,latency_csdStart]=deal(cell(length(aID),2));
[roi_coordsG,roi_coordsR,RelRec_csd,CSD_DCshift]=deal(cell(length(aID),1));


for kk=1:length(aID)

pxlen=aID(kk).FOV_dim/512;%um:px

roi_coordsG{kk} = [aID(kk).xCoor*pxlen, (aID(kk).yCoor)*pxlen];
roi_coordsR{kk} = roi_coordsG{kk}; 

%CSD Start
t_csdG{kk,1} = aID(kk).RecTimeG_csd;
t_csdR{kk,1} = aID(kk).RecTimeR_csd;

bothIndex{kk,1} = and(aID(kk).isRecruitedG_csd,aID(kk).isRecruitedR_csd);

tempt = num2cell(t_csdG{kk,1});
[betaG{kk,1}, VG{kk,1}, pG{kk,1}] = SpatialLinearRegression(tempt(aID(kk).isRecruitedG_csd),roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityG{kk,1}=sqrt(VG{kk,1}(1)^2+VG{kk,1}(2)^2);
if VG{kk,1}(1)>=0 && VG{kk,1}(2)>=0
    thetaG{kk,1}=atand(VG{kk,1}(2)/VG{kk,1}(1));
elseif VG{kk,1}(1)>=0 && VG{kk,1}(2)<0
    thetaG{kk,1}=360-atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
elseif VG{kk,1}(1)<0 && VG{kk,1}(2)<0
    thetaG{kk,1}=180+atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
else
    thetaG{kk,1}=180-atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
end

tempt = num2cell(t_csdR{kk,1});
[betaR{kk,1}, VR{kk,1}, pR{kk,1}] = SpatialLinearRegression(tempt(aID(kk).isRecruitedR_csd),roi_coordsR{kk}(aID(kk).isRecruitedR_csd,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityR{kk,1}=sqrt(VR{kk,1}(1)^2+VR{kk,1}(2)^2);
if VR{kk,1}(1)>=0 && VR{kk,1}(2)>=0
    thetaR{kk,1}=atand(VR{kk,1}(2)/VR{kk,1}(1));
elseif VR{kk,1}(1)>=0 && VR{kk,1}(2)<0
    thetaR{kk,1}=360-atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
elseif VR{kk,1}(1)<0 && VR{kk,1}(2)<0
    thetaR{kk,1}=180+atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
else
    thetaR{kk,1}=180-atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
end

RelRec_csd{kk,1}=t_csdR{kk,1}(bothIndex{kk,1})-t_csdG{kk,1}(bothIndex{kk,1});


latency_csdStart{kk,1}=t_csdG{kk,1}(aID(kk).isRecruitedG_csd)-aID(kk).CSDstartTimeDC;
latency_csdStart{kk,2}=t_csdR{kk,1}(aID(kk).isRecruitedR_csd)-aID(kk).CSDstartTimeDC;

CSD_DCshift{kk}=aID(kk).CSD_dcShift(1);

end

%% Fig4C-E: Vector of propogation plots(actual times axis for red) (s87)
kk=3;

%[RecTimeGsortReG2,RecTimeGReG_I2]=sort(aID(kk).RecTimeG_csd); %sorts cells by rec times according to green
%[RecTimeRsortReR2,RecTimeRReR_I2]=sort(aID(kk).RecTimeR_csd); %sorts cells by rec times according to red

minOnT=min([aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd);aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd)]);
maxOnT=max([aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd);aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd)]);

figure('position',[100 100 600 300])
ax(1)=subplot(1,2,1);
xx = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,1);
yy = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,2);
scatter(xx,yy,[],aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd),'filled')
hold on
line([10,60],[15,15],'color',[0 0 0],'LineStyle','-','LineWidth',2) %scale bar 50um
hold off
colormap(ax(1),parula)
colorbar('location','southoutside','Ticks',[3*floor(minOnT/3):3:3*ceil(maxOnT/3)])
caxis([3*floor(minOnT/3) 3*ceil(maxOnT/3)])
pbaspect([1 1 1])
%title('Non-vGAT (jYCaMP1s)','fontsize',12)
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')


ax(2)=subplot(1,2,2);
xx = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,1);
yy = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,2);
scatter(xx,yy,[],aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd),'filled')
colormap(ax(2),parula)
colorbar('location','southoutside','Ticks',[],'TickLabels',{})
caxis([3*floor(minOnT/3) 3*ceil(maxOnT/3)])
pbaspect([1 1 1])
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_ROImap.svg','-dsvg')



%Vector Plots (using compass)

pthresh=0.05;
max_lim = 100;
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];

figure('position',[100 100 400 400])
h_fake=compass(x_fake,y_fake);
hold on;
c1=compass(VG{kk,1}(1),VG{kk,1}(2));
c2=compass(VR{kk,1}(1),VR{kk,1}(2));

set(h_fake,'Visible','off');


if pG{kk,1}>=pthresh
   c1.Color=[0.9 1 0.9];
else
   c1.Color=[0 0.75 0];
end
c1.LineWidth=4;
c1.LineJoin='chamfer';
c1.LineStyle='-';
if pR{kk,1}>=pthresh
   c2.Color=[1 0.9 1];
else
   c2.Color=[0.75 0 0.75];
end
c2.LineWidth=4;
c2.LineJoin='chamfer';
c2.LineStyle='-';

hold off

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_Polar.svg','-dsvg')


%Propogation Vector Recruitment Times

%order of recruitments by subtype projected onto vector of propogation in
%green

figure('position',[100 100 600 400])
P_v_axisG = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:)*betaG{kk,1}(1:2)/norm(betaG{kk,1}(1:2));
P_v_axisG = P_v_axisG - min(P_v_axisG);
P_v_axisR = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,:)*betaG{kk,1}(1:2)/norm(betaG{kk,1}(1:2));
P_v_axisR = P_v_axisR - min(P_v_axisR);
plot(t_csdG{kk,1}(aID(kk).isRecruitedG_csd),P_v_axisG,'.','MarkerSize',11,'color',[0 0.75 0]);
hold on
plot(t_csdR{kk,1}(aID(kk).isRecruitedR_csd),P_v_axisR,'.','MarkerSize',11,'color',[0.75 0 0.75]);
plot(roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:)*betaG{kk,1}(1:2)+betaG{kk,1}(end),P_v_axisG,'color',[0 0 .5]);
%ylabel('um');
%xlabel('seconds')
%xlim([750 780])
ax2 = gca;
ax2.FontSize = 12;
pbaspect([1 1.5 1])
box off
set(gca,'Visible','on')
Xtimes=[t_csdG{kk,1}(aID(kk).isRecruitedG_csd)',t_csdR{kk,1}(aID(kk).isRecruitedR_csd)'];
xlim([min(Xtimes)-2 max(Xtimes)+2])
ylim([0 550])
clear Xtimes

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_VectorAxis.svg','-dsvg')



%% Fig4F-H: Bar/Box plots of Vector properities (Summary)

% Bar chart Speed by channel
figure('position',[100 100 200 400])
bar(1,mean(cell2mat(velocityG)),'FaceColor',[0 .9 0],'FaceAlpha',0.2,'EdgeColor','flat')
hold on
scatter(ones(length(aID),1),cell2mat(velocityG),40,cMap1(1:length(aID),:),'filled')
errorbar(1,mean(cell2mat(velocityG)),std(cell2mat(velocityG))/sqrt(length(aID)),std(cell2mat(velocityG))/sqrt(length(aID)),'color',[0.25 0.25 0.25])
bar(2,mean(cell2mat(velocityR)),'FaceColor',[.9 0 .9],'FaceAlpha',0.2,'EdgeColor','flat')
scatter(ones(length(aID),1)*2,cell2mat(velocityR),40,cMap1(1:length(aID),:),'filled')
errorbar(2,mean(cell2mat(velocityR)),std(cell2mat(velocityR))/sqrt(length(aID)),std(cell2mat(velocityR))/sqrt(length(aID)),'color',[0.25 0.25 0.25])
xlim([0 3])
ylim([0 100])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_Velocity.svg','-dsvg')


%% Box plot speed by channel
figure
boxchart(ones(size(cell2mat(velocityG))),cell2mat(velocityG),'BoxFaceColor',[0 1 0]);
hold on
scatter(ones(length(aID),1),cell2mat(velocityG),[],cMap1(1:length(aID),:),'filled')
boxchart(ones(size(cell2mat(velocityR)))*2,cell2mat(velocityR),'BoxFaceColor',[1 0 1]);
scatter(ones(length(aID),1)*2,cell2mat(velocityR),[],cMap1(1:length(aID),:),'filled')
xlim([0 3])
ylim([0 110])
pbaspect([1 3 1])
set(gca,'box','off','xcolor','none')

%% Theta difference between channels
thetaDiff=abs(cell2mat(thetaG(:,1))-cell2mat(thetaR(:,1)));
thetaDiff(thetaDiff>180)=abs(thetaDiff(thetaDiff>180)-360);
figure('position',[100 100 200 400])
boxchart(ones(size(thetaDiff)),thetaDiff,'BoxFaceColor',[0.5 0.5 0.5]);
hold on
scatter(ones(length(aID),1),thetaDiff,[],cMap1(1:length(aID),:),'filled')
colormap(cMap1)
xlim([0 2])
ylim([0 90])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_theta.svg','-dsvg')


%% GLM for difference between groups onset

%CSD Start
GLMtable=cell(1,length(aID));
varNames={'RecTime','RCatchER','Subject'};
for kk=1:length(aID)
GLMtable{kk}=table([t_csdG{kk,1}(bothIndex{kk,1});t_csdR{kk,1}(bothIndex{kk,1})],[zeros(sum(bothIndex{kk,1}),1);ones(sum(bothIndex{kk,1}),1)],kk*ones(sum(bothIndex{kk,1})*2,1),'VariableNames',varNames);
end
RelRec_GLMdata=vertcat(GLMtable{:});

RelRec_GLMout=fitglme(RelRec_GLMdata,'RecTime ~ 1 + RCatchER + (RCatchER|Subject) + (1|Subject)','Distribution','normal');


%% Box plots of relative recruitment
figure('position',[100 100 250 400])
hold on
for ii=1:length(RelRec_csd)
boxchart(ones(size(RelRec_csd{ii}))*ii,RelRec_csd{ii},'BoxFaceColor',cMap1(ii,:),'MarkerColor','none');
end
line([0,5],[0,0],'color',[.5 0 0],'LineStyle','--')
xlim([0 length(aID)+1])
ylim([-8 8])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_RelRecIndiv.svg','-dsvg')


%% plot of population relative recruitment vertical
clear hh* clear pp*
RelRec_csdMean=cellfun(@nanmean,RelRec_csd);
RelRec_csd_std=cellfun(@nanstd,RelRec_csd);
RelRec_csd_df=cellfun(@length,RelRec_csd)-1;

% Pooled SE
RelRec_csd_stdp=sqrt(sum(RelRec_csd_df.*RelRec_csd_std.^2,1)./sum(RelRec_csd_df,1));
RelRec_csd_sep=RelRec_csd_stdp.*sqrt(sum(1./(RelRec_csd_df+1)));

figure('position',[100 100 140 400])
barPlot=bar(1,mean(RelRec_csdMean),'FaceColor','flat','FaceAlpha',0.2,'EdgeColor','flat','LineWidth',1);
barPlot.CData = [0.25 .25 0.25];
hold on
scatter(ones(length(aID),1),RelRec_csdMean,[],cMap1(1:length(aID),:),'filled')
errorbar(1,mean(RelRec_csdMean),RelRec_csd_sep,RelRec_csd_sep,'color',[0.25 0.25 0.25])
colormap(cMap1)
xlim([0 2])
ylim([-4 4])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure4/Figure4_TSD_RelRecGrp.svg','-dsvg')


%% Fig 4 Additional Relevant Values (DC Shift mag, length)

DCshift_mean=nanmean(cell2mat(CSD_DCshift));
DCshift_std=nanstd(cell2mat(CSD_DCshift));
DCshift_ste=DCshift_std/sqrt(sum(~isnan(cell2mat(CSD_DCshift))));

latency_csdStart_mean=cellfun(@nanmean,latency_csdStart);
latency_csdStart_std=cellfun(@nanstd,latency_csdStart);
latency_csdStart_df=cellfun(@length,latency_csdStart)-1;

%Mean of Means
latency_csdStart_GrandMean=nanmean(latency_csdStart_mean,1);

% Pooled Ste
latency_csdStart_stdp=sqrt(nansum(latency_csdStart_df.*latency_csdStart_std.^2,1)./nansum(latency_csdStart_df,1));
latency_csdStart_sep=latency_csdStart_stdp.*sqrt(sum(1./(latency_csdStart_df+1)));








%% Figure 5 Stim
load('s86_221024_1_stim1/S86stim1.mat')
load('s87_221005_stim1/S87stim1.mat')
load('s88_221013_stim1/S88stim1.mat')
load('s89_221027_1_stim1/S89stim1.mat')
load('s102_230605_stim1/S102stim1.mat')
load('s103_230607_stim1/S103stim1.mat')
load('s105_230616_stim1/S105stim1.mat')
aID=[S86stim1,S87stim1,S88stim1,S89stim1,S102stim1,S103stim1,S105stim1];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
cMap1=colormap(lines);


%% Fig 5B: Mean Traces, DC and EEG

kk=6;
figure
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+2,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1)-.2,'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-.9,'color',[1 0 1]);
ylim([-2 5])
xlim([10 310])
yticks([])
pbaspect([1 .6 1])
hold off

%scale bars
sb_x=15; %30seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[15 -0.6];
sb_omr=[15 -1.3];
sb_oe=[15 2.4];
sb_od=[15 1.2]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_103_TracesStim.svg','-dsvg')

%% Fig5B: Rasters
figure
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:30:1050])
xticklabels([])
yticks([349:50:699])
yticklabels([0:50:600])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([10 310])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt2,1),aID(kk).Fc1Rdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:30:1050])
xticklabels([0:30:1050])
yticks([349:50:699])
yticklabels([0:50:600])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([10 310])


print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_103RasterStim.svg','-dsvg')

%% Figure 5C: CSD Ca Levels

% store relevant data (CSD)
%invert the r channel as all the values are determined based upon an
%inverted signal

Ca_CSD=cell(length(aID),2);
for kk=1:length(aID)
    Ca_CSD{kk,1}=aID(kk).CSDMeanDffDiffG;
    Ca_CSD{kk,2}=-aID(kk).CSDMeanDffDiffR;
end
clear kk

% GLME plot
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
Ca_CSD_means=cellfun(@mean,Ca_CSD);
Ca_CSD_std=cellfun(@std,Ca_CSD);
Ca_CSD_df=cellfun(@length,Ca_CSD)-1;
Ca_CSD_GrandMean=mean(Ca_CSD_means,1);

% Pooled SE
Ca_CSD_stdp=sqrt(sum(Ca_CSD_df.*Ca_CSD_std.^2,1)./sum(Ca_CSD_df,1));
Ca_CSD_sep=Ca_CSD_stdp.*sqrt(sum(1./(Ca_CSD_df+1)));


figure('Position',[100 100 220 400])
plotYvals=Ca_CSD_GrandMean';
plotYerr=Ca_CSD_sep';
barPlot=bar(plotYvals,'grouped','FaceColor','flat');
barPlot.CData = [0 .7 0 ; 0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.4 .8])
xticks(1.5)
xticklabels({'CSD'})
set(gca,'box','off')


print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_CaLvl.svg','-dsvg')


%% Linear Regression
% Propogation Dynamics
[t_csdG,t_csdR,betaG,VG,pG,velocityG,thetaG,betaR,VR,pR,velocityR,thetaR,RelRec_csd,RelRec_csdMean,length_csd,latency_csdStart,latency_csdEnd,bothIndex,both2Index]=deal(cell(length(aID),2));
[roi_coordsG,roi_coordsR,length_csdDC,CSD_DCshift]=deal(cell(length(aID),1));

for kk=1:length(aID)

pxlen=aID(kk).FOV_dim/512;%um:px


roi_coordsG{kk} = [aID(kk).xCoor*pxlen, (aID(kk).yCoor)*pxlen];
roi_coordsR{kk} = roi_coordsG{kk}; 

%CSD Start
t_csdG{kk,1} = aID(kk).RecTimeG_csd;
t_csdR{kk,1} = aID(kk).RecTimeR_csd;

bothIndex{kk,1} = and(aID(kk).isRecruitedG_csd,aID(kk).isRecruitedR_csd);

tempt = num2cell(t_csdG{kk,1});
[betaG{kk,1}, VG{kk,1}, pG{kk,1}] = SpatialLinearRegression(tempt(aID(kk).isRecruitedG_csd),roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityG{kk,1}=sqrt(VG{kk,1}(1)^2+VG{kk,1}(2)^2);
if VG{kk,1}(1)>=0 && VG{kk,1}(2)>=0
    thetaG{kk,1}=atand(VG{kk,1}(2)/VG{kk,1}(1));
elseif VG{kk,1}(1)>=0 && VG{kk,1}(2)<0
    thetaG{kk,1}=360-atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
elseif VG{kk,1}(1)<0 && VG{kk,1}(2)<0
    thetaG{kk,1}=180+atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
else
    thetaG{kk,1}=180-atand(abs(VG{kk,1}(2)/VG{kk,1}(1)));
end

tempt = num2cell(t_csdR{kk,1});
[betaR{kk,1}, VR{kk,1}, pR{kk,1}] = SpatialLinearRegression(tempt(aID(kk).isRecruitedR_csd),roi_coordsR{kk}(aID(kk).isRecruitedR_csd,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityR{kk,1}=sqrt(VR{kk,1}(1)^2+VR{kk,1}(2)^2);
if VR{kk,1}(1)>=0 && VR{kk,1}(2)>=0
    thetaR{kk,1}=atand(VR{kk,1}(2)/VR{kk,1}(1));
elseif VR{kk,1}(1)>=0 && VR{kk,1}(2)<0
    thetaR{kk,1}=360-atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
elseif VR{kk,1}(1)<0 && VR{kk,1}(2)<0
    thetaR{kk,1}=180+atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
else
    thetaR{kk,1}=180-atand(abs(VR{kk,1}(2)/VR{kk,1}(1)));
end

RelRec_csd{kk,1}=t_csdR{kk,1}(bothIndex{kk,1})-t_csdG{kk,1}(bothIndex{kk,1});

%CSD End
t_csdG{kk,2} = aID(kk).RecTimeG_csdRTN;
t_csdR{kk,2} = aID(kk).RecTimeR_csdRTN;

bothIndex{kk,2} = and(aID(kk).isRecG_csdRTN,aID(kk).isRecR_csdRTN);

tempt = num2cell(t_csdG{kk,2});
[betaG{kk,2}, VG{kk,2}, pG{kk,2}] = SpatialLinearRegression(tempt(aID(kk).isRecG_csdRTN),roi_coordsG{kk}(aID(kk).isRecG_csdRTN,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityG{kk,2}=sqrt(VG{kk,2}(1)^2+VG{kk,2}(2)^2);
if VG{kk,2}(1)>=0 && VG{kk,2}(2)>=0
    thetaG{kk,2}=atand(VG{kk,2}(2)/VG{kk,2}(1));
elseif VG{kk,2}(1)>=0 && VG{kk,2}(2)<0
    thetaG{kk,2}=360-atand(abs(VG{kk,2}(2)/VG{kk,2}(1)));
elseif VG{kk,2}(1)<0 && VG{kk,2}(2)<0
    thetaG{kk,2}=180+atand(abs(VG{kk,2}(2)/VG{kk,2}(1)));
else
    thetaG{kk,2}=180-atand(abs(VG{kk,2}(2)/VG{kk,2}(1)));
end

tempt = num2cell(t_csdR{kk,2});
[betaR{kk,2}, VR{kk,2}, pR{kk,2}] = SpatialLinearRegression(tempt(aID(kk).isRecR_csdRTN),roi_coordsR{kk}(aID(kk).isRecR_csdRTN,:),'switch_plot',0,'Lossfun','L1','n_shuffle',1000); % Least absolute deviation estimator

velocityR{kk,2}=sqrt(VR{kk,2}(1)^2+VR{kk,2}(2)^2);
if VR{kk,2}(1)>=0 && VR{kk,2}(2)>=0
    thetaR{kk,2}=atand(VR{kk,2}(2)/VR{kk,2}(1));
elseif VR{kk,2}(1)>=0 && VR{kk,2}(2)<0
    thetaR{kk,2}=360-atand(abs(VR{kk,2}(2)/VR{kk,2}(1)));
elseif VR{kk,2}(1)<0 && VR{kk,2}(2)<0
    thetaR{kk,2}=180+atand(abs(VR{kk,2}(2)/VR{kk,2}(1)));
else
    thetaR{kk,2}=180-atand(abs(VR{kk,2}(2)/VR{kk,2}(1)));
end

RelRec_csd{kk,2}=t_csdR{kk,2}(bothIndex{kk,2})-t_csdG{kk,2}(bothIndex{kk,2});

both2Index{kk,1}=and(aID(kk).isRecruitedG_csd,aID(kk).isRecG_csdRTN);%column 1 is G and 2 is R
length_csd{kk,1}=t_csdG{kk,2}(both2Index{kk,1})-t_csdG{kk,1}(both2Index{kk,1});

both2Index{kk,2}=and(aID(kk).isRecruitedR_csd,aID(kk).isRecR_csdRTN);%column 1 is G and 2 is R
length_csd{kk,2}=t_csdR{kk,2}(both2Index{kk,2})-t_csdR{kk,1}(both2Index{kk,2});

latency_csdStart{kk,1}=t_csdG{kk,1}(aID(kk).isRecruitedG_csd)-aID(kk).CSDstartTimeDC;
latency_csdEnd{kk,1}=t_csdG{kk,2}(aID(kk).isRecG_csdRTN)-aID(kk).CSDendTimeDC;
latency_csdStart{kk,2}=t_csdR{kk,1}(aID(kk).isRecruitedR_csd)-aID(kk).CSDstartTimeDC;
latency_csdEnd{kk,2}=t_csdR{kk,2}(aID(kk).isRecR_csdRTN)-aID(kk).CSDendTimeDC;

length_csdDC{kk}=aID(kk).CSDendTimeDC-aID(kk).CSDstartTimeDC;
CSD_DCshift{kk}=aID(kk).CSD_dcShift(1);

end


%% GLM for difference between groups
%CSD Start
GLMtable=cell(1,length(aID));
varNames={'RecTime','RCatchER','Subject'};
for kk=1:length(aID)
GLMtable{kk}=table([t_csdG{kk,1}(bothIndex{kk,1});t_csdR{kk,1}(bothIndex{kk,1})],[zeros(sum(bothIndex{kk,1}),1);ones(sum(bothIndex{kk,1}),1)],kk*ones(sum(bothIndex{kk,1})*2,1),'VariableNames',varNames);
end
RelRec_GLMdata=vertcat(GLMtable{:});

RelRec_GLMout=fitglme(RelRec_GLMdata,'RecTime ~ 1 + RCatchER + (RCatchER|Subject) + (1|Subject)','Distribution','normal');


%% Fig5C-E: Rainbow ROI Field (actual times axis for red)

kk=6;

%[RecTimeGsortReG2,RecTimeGReG_I2]=sort(aID(kk).RecTimeG_csd); %sorts cells by rec times according to green
%[RecTimeRsortReR2,RecTimeRReR_I2]=sort(aID(kk).RecTimeR_csd); %sorts cells by rec times according to red

minOnT=min([aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd);aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd)]);
maxOnT=max([aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd);aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd)]);

figure('position',[100 100 600 300])
ax(1)=subplot(1,2,1);
xx = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,1);
yy = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,2);
scatter(xx,yy,[],aID(kk).RecTimeG_csd(aID(kk).isRecruitedG_csd),'filled')
hold on
line([10,60],[15,15],'color',[0 0 0],'LineStyle','-','LineWidth',2) %scale bar
hold off
colormap(ax(1),parula)
colorbar('location','southoutside','Ticks',[2*floor(minOnT/2):2:2*ceil(maxOnT/2)])
caxis([2*floor(minOnT/2) 2*ceil(maxOnT/2)])
pbaspect([1 1 1])
%title('Non-vGAT (jYCaMP1s)','fontsize',12)
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

ax(2)=subplot(1,2,2);
xx = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,1);
yy = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,2);
scatter(xx,yy,[],aID(kk).RecTimeR_csd(aID(kk).isRecruitedR_csd),'filled')
colormap(ax(2),parula)
colorbar('location','southoutside','Ticks',[],'TickLabels',{})
caxis([2*floor(minOnT/2) 2*ceil(maxOnT/2)])
pbaspect([1 1 1])
%title('VGAT (jRGECO1a)','fontsize',12)
yticks([])
xticks([])
xlim([0 450])
ylim([0 450])
box on
set(gca,'Visible','on')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_ROImap.svg','-dsvg')


%Vector Plots (using compass)

pthresh=0.05;
max_lim = 100;
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];

figure('position',[100 100 400 400])
h_fake=compass(x_fake,y_fake);
hold on;
c1=compass(VG{kk,1}(1),VG{kk,1}(2));
c2=compass(VR{kk,1}(1),VR{kk,1}(2));

set(h_fake,'Visible','off');


if pG{kk,1}>=pthresh
   c1.Color=[0.9 1 0.9];
else
   c1.Color=[0 0.75 0];
end
c1.LineWidth=4;
c1.LineJoin='chamfer';
if pR{kk,1}>=pthresh
   c2.Color=[1 0.9 1];
else
   c2.Color=[0.75 0 0.75];
end
c2.LineWidth=4;
c2.LineJoin='chamfer';
hold off

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_Polar.svg','-dsvg')

%Propogation Vector Recruitment Times

%order of recruitments by subtype projected onto vector of propogation in
%green

figure('position',[100 100 600 400]);
P_v_axisG = roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:)*betaG{kk,1}(1:2)/norm(betaG{kk,1}(1:2));
P_v_axisG = P_v_axisG - min(P_v_axisG);
P_v_axisR = roi_coordsR{kk}(aID(kk).isRecruitedR_csd,:)*betaG{kk,1}(1:2)/norm(betaG{kk,1}(1:2));
P_v_axisR = P_v_axisR - min(P_v_axisR);
plot(t_csdG{kk,1}(aID(kk).isRecruitedG_csd),P_v_axisG,'.','MarkerSize',11,'color',[0 0.75 0]);
hold on
plot(t_csdR{kk,1}(aID(kk).isRecruitedR_csd),P_v_axisR,'.','MarkerSize',11,'color',[0.75 0 0.75]);
plot(roi_coordsG{kk}(aID(kk).isRecruitedG_csd,:)*betaG{kk,1}(1:2)+betaG{kk,1}(end),P_v_axisG,'color',[0 0 .5]);
ax2 = gca;
ax2.FontSize = 12;
pbaspect([1 1.5 1])
box off
set(gca,'Visible','on')
Xtimes=[t_csdG{kk,1}(aID(kk).isRecruitedG_csd)',t_csdR{kk,1}(aID(kk).isRecruitedR_csd)'];
xlim([min(Xtimes)-2 max(Xtimes)+2])
ylim([0 550])
clear Xtimes

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_VectorAxis.svg','-dsvg')


%% Fig5F-H: Bar/Box plots of Vector properities (Summary)

% Bar chart Speed by channel
figure('position',[100 100 200 400])
bar(1,mean(cell2mat(velocityG(:,1))),'FaceColor',[0 .9 0],'FaceAlpha',0.2,'EdgeColor','flat')
hold on
scatter(ones(length(aID),1),cell2mat(velocityG(:,1)),[],cMap1(1:length(aID),:),'filled')
errorbar(1,mean(cell2mat(velocityG(:,1))),std(cell2mat(velocityG(:,1)))/sqrt(length(aID)),std(cell2mat(velocityG(:,1)))/sqrt(length(aID)),'color',[0.25 0.25 0.25])
bar(2,mean(cell2mat(velocityR(:,1))),'FaceColor',[.9 0 .9],'FaceAlpha',0.2,'EdgeColor','flat')
scatter(ones(length(aID),1)*2,cell2mat(velocityR(:,1)),[],cMap1(1:length(aID),:),'filled')
errorbar(2,mean(cell2mat(velocityR(:,1))),std(cell2mat(velocityR(:,1)))/sqrt(length(aID)),std(cell2mat(velocityR(:,1)))/sqrt(length(aID)),'color',[0.25 0.25 0.25])
xlim([0 3])
ylim([0 110])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_Velocity.svg','-dsvg')


%% Box plot Speed by channel
figure
boxchart(ones(size(cell2mat(velocityG(:,1)))),cell2mat(velocityG(:,1)),'BoxFaceColor',[0 1 0]);
hold on
scatter(ones(length(aID),1),cell2mat(velocityG(:,1)),[],cMap1(1:length(aID),:),'filled')
boxchart(ones(size(cell2mat(velocityR(:,1))))*2,cell2mat(velocityR(:,1)),'BoxFaceColor',[1 0 1]);
scatter(ones(length(aID),1)*2,cell2mat(velocityR(:,1)),[],cMap1(1:length(aID),:),'filled')
xlim([0 3])
ylim([0 110])
set(gca,'box','off','xcolor','none')

%% Theta difference between channels
thetaDiff=abs(cell2mat(thetaG(:,1))-cell2mat(thetaR(:,1)));
thetaDiff(thetaDiff>180)=abs(thetaDiff(thetaDiff>180)-360);
figure('position',[100 100 200 400])
boxchart(ones(size(thetaDiff)),thetaDiff,'BoxFaceColor',[0.5 0.5 0.5]);
hold on
scatter(ones(length(aID),1),thetaDiff,[],cMap1(1:length(aID),:),'filled');%,'Xjitter','randn', 'XJitterWidth', .5)
colormap(cMap1)
xlim([0 2])
ylim([0 90])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_theta.svg','-dsvg')


%% Box plots of relative recruitment
figure('position',[100 100 400 400])
hold on
for ii=1:size(RelRec_csd,1)
boxchart(ones(size(RelRec_csd{ii,1}))*ii,RelRec_csd{ii,1},'BoxFaceColor',cMap1(ii,:),'MarkerColor','none');
end
line([0,size(RelRec_csd,1)+1],[0,0],'color',[.5 0 0],'LineStyle','--')
xlim([0 size(RelRec_csd,1)+1])
ylim([-6 6])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_RelRecIndiv.svg','-dsvg')


%% plot of population relative recruitment vertical
clear hh* clear pp*

RelRec_csdMean=cellfun(@nanmean,RelRec_csd);
RelRec_csd_std=cellfun(@nanstd,RelRec_csd);
RelRec_csd_df=cellfun(@length,RelRec_csd)-1;

% Pooled SE
RelRec_csd_stdp=sqrt(sum(RelRec_csd_df.*RelRec_csd_std.^2,1)./sum(RelRec_csd_df,1));
RelRec_csd_sep=RelRec_csd_stdp.*sqrt(sum(1./(RelRec_csd_df+1)));


figure('position',[100 100 140 400])
barPlot=bar(1,mean(RelRec_csdMean(:,1)),'FaceColor','flat','FaceAlpha',0.2,'EdgeColor','flat','LineWidth',1);
barPlot.CData = [0.25 .25 0.25];
hold on
scatter(ones(length(aID),1),RelRec_csdMean(:,1),[],cMap1(1:length(aID),:),'filled')
errorbar(1,mean(RelRec_csdMean(:,1)),RelRec_csd_sep(:,1),RelRec_csd_sep(:,1),'color',[0.25 0.25 0.25])
colormap(cMap1)
xlim([0 2])
ylim([-4 4])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_Stim_RelRecGrp.svg','-dsvg')


%% Fig5: CSD length (pooled means and indiv sz plots; no stats)
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate

length_csdMean=cellfun(@nanmean,length_csd);
length_csd_std=cellfun(@std,length_csd);
length_csd_df=cellfun(@length,length_csd)-1;
length_csd_GrandMean=mean(length_csdMean,1);

% Pooled SE
length_csd_stdp=sqrt(sum(length_csd_df.*length_csd_std.^2,1)./sum(length_csd_df,1));
length_csd_sep=length_csd_stdp.*sqrt(sum(1./(length_csd_df+1)));

%figure
figure('Position',[100 100 250,450])
plotYvals=length_csd_GrandMean';
plotYerr=length_csd_sep';
barPlot=bar(plotYvals,'grouped','FaceColor','flat','FaceAlpha',0.2,'EdgeColor','flat','LineWidth',1);
barPlot.CData = [0 .9 0 ; .9 0 .9];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25],'linewidth',1.25);
scatter(ones(length(aID),1),length_csdMean(:,1),40,cMap1(1:length(aID),:),'filled')
colormap(cMap1)
scatter(ones(length(aID),1)*2,length_csdMean(:,2),40,cMap1(1:length(aID),:),'filled')
colormap(cMap1)
hold off
%ylim([-.25 .75])
set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_CSDLength.svg','-dsvg')


figure('position',[100 100 650 400])
hold on
for ii=1:length(aID)
boxchart(ones(size(length_csd{ii,1}))*ii*2-1,length_csd{ii,1},'BoxFaceColor',[0 .9 0],'MarkerColor','none');
boxchart(ones(size(length_csd{ii,2}))*ii*2,length_csd{ii,2},'BoxFaceColor',[.9 0 .9],'MarkerColor','none');
text(ii*2-.5,85,'animal '+string(ii),'color',cMap1(ii,:),'HorizontalAlignment','center')
end
xlim([0 2*length(aID)+1])
ylim([20 90])
set(gca,'box','off','xcolor','none')

print('../../RCatchER_CSD_Paper/Figures/Figure5/Figure5_CSDLength_Indiv.svg','-dsvg')


%% Fig 5 Additional Relevant Values (DC Shift mag, length)

CSDlengthDC_mean=nanmean(cell2mat(length_csdDC));
CSDlengthDC_std=nanstd(cell2mat(length_csdDC));
CSDlengthDC_ste=CSDlengthDC_std/sqrt(sum(~isnan(cell2mat(length_csdDC))));

DCshift_mean=nanmean(cell2mat(CSD_DCshift));
DCshift_std=nanstd(cell2mat(CSD_DCshift));
DCshift_ste=DCshift_std/sqrt(sum(~isnan(cell2mat(CSD_DCshift))));

latency_csdStart_mean=cellfun(@nanmean,latency_csdStart);
latency_csdStart_std=cellfun(@nanstd,latency_csdStart);
latency_csdStart_df=cellfun(@length,latency_csdStart)-1;
latency_csdEnd_mean=cellfun(@nanmean,latency_csdEnd);
latency_csdEnd_std=cellfun(@nanstd,latency_csdEnd);
latency_csdEnd_df=cellfun(@length,latency_csdEnd)-1;

%Mean of Means
latency_csdStart_GrandMean=nanmean(latency_csdStart_mean,1);
latency_csdEnd_GrandMean=nanmean(latency_csdEnd_mean,1);

% Pooled Ste
latency_csdStart_stdp=sqrt(nansum(latency_csdStart_df.*latency_csdStart_std.^2,1)./nansum(latency_csdStart_df,1));
latency_csdEnd_stdp=sqrt(nansum(latency_csdEnd_df.*latency_csdEnd_std.^2,1)./nansum(latency_csdEnd_df,1));
latency_csdStart_sep=latency_csdStart_stdp.*sqrt(sum(1./(latency_csdStart_df+1)));
latency_csdEnd_sep=latency_csdEnd_stdp.*sqrt(sum(1./(latency_csdEnd_df+1)));










%% Figure s4 (RNS)
%Cumulative propability plot

plotVal_X={'0','25','50','100','250','500','750','1000'};
%plotVal_Y=[250,500,750;250,750,750];
plotVal_Y=[4,5,6;4,6,6];
figure('position',[100 100 600 300])
subplot(1,3,1:2)
hh1=cdfplot(plotVal_Y(:));
set(hh1,'LineWidth',1)
xlim([0 7])
xticks(0:7)
xticklabels(plotVal_X)
xlabel('Threshold Current (uA)')
ylabel('Cumlative Probability')
title('Stimulation to Induce CSD (RNS Parameters)')
set(gca,'box','off')
grid off

subplot(1,3,3)
scatter([1,2,1,2,1,2],plotVal_Y(:),40,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold on
line([1,1,1;2,2,2],plotVal_Y,'color',[0.5 0 0],'LineWidth',1)
ylim([0 7])
xlim([.5 2.5])
yticks(0:7)
yticklabels(plotVal_X)
xlabel('Trail')
ylabel('Threshold Current (uA)')

set(gca,'box','off')

print('../../RCatchER_CSD_Paper/Figures/FigureS4/FigureS4_RNSstimCummInc.svg','-dsvg')


%% load neuropacepulse.mat from thomass eggers
%200 hz with 1 pulse over 160us (pulse width) burts for 100ms (so 20 pulses) and 5 seconds between each burst for 5 bursts.
ampY=35000;
figure('position',[100 100 1000 200])
subplot(1,4,1:2)
rectangle('Position',[2452.6,-ampY,2650-2452.6,2*ampY],'FaceColor',[0.9 0.9 0.9],'linestyle','none')
hold on
plot(2*10^-3:2*10^-3:5*2*length(neuropacepulse)/10^3,repmat(neuropacepulse(:,1),5,1),'color', [0 0 0])
set(gca,'Visible','off')
subplot(1,4,3)
rectangle('Position',[2504.18,-ampY,2505.5-2504.18,2*ampY],'FaceColor',[0.9 0.9 0.9],'linestyle','none')
hold on
plot(2*10^-3:2*10^-3:2*length(neuropacepulse)/10^3,neuropacepulse(:,1),'color', [0 0 0])
xlim([2450 2650])
set(gca,'Visible','off')
subplot(1,4,4)
plot(2*10^-3:2*10^-3:2*length(neuropacepulse)/10^3,neuropacepulse(:,1),'color', [0 0 0])
xlim([2504.18 2505.5])
set(gca,'Visible','off')
print('../../RCatchER_CSD_Paper/Figures/FigureS4/FigureS4_RNSstimWaveForm.svg','-dsvg')

%%
figure
tt = 0:1/1e4:10;
yy = square(2*pi*2000*tt);
plot(tt,yy)



%%
%%

%% Figure 6 Stim
load('s102_230628_cicr1/S102cicr1.mat')
load('s105_230629_cicr1/S105cicr1.mat')
load('s102_230605_stim1/S102stim1.mat')
load('s105_230616_stim1/S105stim1.mat')
aID=[S102stim1,S102cicr1,S105stim1,S105cicr1];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
cMap1=colormap(lines);

    
%% cicr block study: Mean Traces, DC and EEG

kk=1;
figure
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+2,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt1,1)-.2,'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt1,1)-.9,'color',[1 0 1]);
ylim([-2 5])
xlim([10 310])
yticks([])
pbaspect([1 .6 1])
hold off

%scale bars
sb_x=10; %30seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[15 -0.6];
sb_omr=[15 -1.3];
sb_oe=[15 2.4];
sb_od=[15 1.2]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

%% cicr block study: Rasters
figure
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt1,1),aID(kk).Fc1Gdff_flt1(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:30:1050])
xticklabels([])
yticks([349:50:699])
yticklabels([0:50:600])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([10 310])

ax(2)=subplot(2,1,2);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Rdff_flt1,1),aID(kk).Fc1Rdff_flt1(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 5])
xticks([0:30:1050])
xticklabels([0:30:1050])
yticks([349:50:699])
yticklabels([0:50:600])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([10 310])


%% cicr block study: CSD Ca Levels


%aID=[S102stim1,S105stim1];
%aID=[S102cicr1,S105cicr1];

% store relevant data (CSD)
%invert the r channel as all the values are determined based upon an
%inverted signal

Ca_CSD=cell(length(aID),2);
for kk=1:length(aID)
    Ca_CSD{kk,1}=aID(kk).CSDMeanDffDiffG;
    Ca_CSD{kk,2}=-aID(kk).CSDMeanDffDiffR;
end
clear kk

% GLME plot
 % do a mean of means bar plot series for the figure E above but now with  all data
 % grouped and do SEM for error bars...could look at dta to see if bar plot
 % is more appropriate
 
Ca_CSD_means=cellfun(@mean,Ca_CSD);
Ca_CSD_std=cellfun(@std,Ca_CSD);
Ca_CSD_df=cellfun(@length,Ca_CSD)-1;
Ca_CSD_GrandMean=mean(Ca_CSD_means,1);

% Pooled SE
Ca_CSD_stdp=sqrt(sum(Ca_CSD_df.*Ca_CSD_std.^2,1)./sum(Ca_CSD_df,1));
Ca_CSD_sep=Ca_CSD_stdp.*sqrt(sum(1./(Ca_CSD_df+1)));


figure('Position',[100 100 200 500])
plotYvals=Ca_CSD_GrandMean';
plotYerr=Ca_CSD_sep';
barPlot=bar(plotYvals,'grouped','FaceColor','flat');
barPlot.CData = [0 .7 0 ; 0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
ylim([-.5 1])
xticks(1.5)
xticklabels({'CSD'})
set(gca,'box','off')

%% cicr block study
aID=[S102stim1,S102cicr1,S105stim1,S105cicr1];
figure('Position',[100 100 800 300])
for kk=1:4
subplot(1,4,kk)
plotYvals=[mean(aID(kk).CSDMeanDffDiffG); -mean(aID(kk).CSDMeanDffDiffR)];
plotYerr=[std(aID(kk).CSDMeanDffDiffG); std(aID(kk).CSDMeanDffDiffR)];
barPlot=bar(plotYvals,'grouped','FaceColor','flat');
barPlot.CData = [0 .7 0 ; 0.7 0 .7];
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(plotYvals);
% Get the x coordinate of the bars
errXvals = nan(nbars, ngroups);
for ii = 1:nbars
    errXvals(ii,:) = barPlot(ii).XEndPoints;
end
% Plot the errorbars
errorbar(errXvals',plotYvals,plotYerr,'k','linestyle','none','color',[.25 .25 .25]);
hold off
%ylim([-.25 1.25])
ylim([-.5 1.5])
xticks(1.5)
xticklabels({'PIS'})
title(string(aID(kk).szID))
set(gca,'box','off')
end

%% Fig 5 Additional Relevant Values (DC Shift mag, length)

CSDlengthDC_mean=nanmean(cell2mat(length_csdDC));
CSDlengthDC_std=nanstd(cell2mat(length_csdDC));
CSDlengthDC_ste=CSDlengthDC_std/sqrt(sum(~isnan(cell2mat(length_csdDC))));

DCshift_mean=nanmean(cell2mat(CSD_DCshift));
DCshift_std=nanstd(cell2mat(CSD_DCshift));
DCshift_ste=DCshift_std/sqrt(sum(~isnan(cell2mat(CSD_DCshift))));

latency_csdStart_mean=cellfun(@nanmean,latency_csdStart);
latency_csdStart_std=cellfun(@nanstd,latency_csdStart);
latency_csdStart_df=cellfun(@length,latency_csdStart)-1;
latency_csdEnd_mean=cellfun(@nanmean,latency_csdEnd);
latency_csdEnd_std=cellfun(@nanstd,latency_csdEnd);
latency_csdEnd_df=cellfun(@length,latency_csdEnd)-1;

%Mean of Means
latency_csdStart_GrandMean=nanmean(latency_csdStart_mean,1);
latency_csdEnd_GrandMean=nanmean(latency_csdEnd_mean,1);

% Pooled Ste
latency_csdStart_stdp=sqrt(nansum(latency_csdStart_df.*latency_csdStart_std.^2,1)./nansum(latency_csdStart_df,1));
latency_csdEnd_stdp=sqrt(nansum(latency_csdEnd_df.*latency_csdEnd_std.^2,1)./nansum(latency_csdEnd_df,1));
latency_csdStart_sep=latency_csdStart_stdp.*sqrt(sum(1./(latency_csdStart_df+1)));
latency_csdEnd_sep=latency_csdEnd_stdp.*sqrt(sum(1./(latency_csdEnd_df+1)));



%%









%% Figure Methods single trace detection
load('s103_230531/S103ptz1.mat');
aID=[S103ptz1,S103ptz1];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
kk=1;

%%
cellList=[9,13];

figure
hold on
cc=1;

for ii=cellList
%     plot(aID(kk).F1_ts,SternRollAvg(aID(kk).Fc1Rdff_flt2(ii,:),10)+cc,'color', [0.75 0 0.75]);
%     plot(aID(kk).F1_ts,SternRollAvg(aID(kk).Fc1Gdff_flt2(ii,:),10)/2+cc+1,'color', [0 0.75 0]);
    plot(aID(kk).F1_ts,aID(kk).Fc1Rdff_flt2(ii,:)+cc,'color', [0.75 0 0.75]);
    plot(aID(kk).F1_ts,aID(kk).Fc1Gdff_flt2(ii,:)/2+cc+1,'color', [0 0.75 0]);
    line([aID(kk).RecTimeG_csd(ii);aID(kk).RecTimeG_csd(ii)],[cc+.75;cc+1.25],'color',[0 .25 0],'linewidth',.75)
    line([aID(kk).RecTimeR_csd(ii);aID(kk).RecTimeR_csd(ii)],[cc-.25;cc+.25],'color',[.25 0 .25],'linewidth',.75)
    cc=cc+2;
end

plot(aID(kk).EEG_ts,aID(kk).EEG./2000+cc+1.5,'color',[0.5 0.5 0.5]);
plot(aID(kk).DC_ts,aID(kk).DC./20+cc+1,'color',[0 0 0]);
xlim([275 350])
ylim([0 cc+3])
pbaspect([1 1.5 1])
set(gca,'Visible','off')










%% ARCHIVE

%% Figure 2: Version s89
load('s89_221027_2_ptz1/S89ptz1.mat');
load('s89_221101_ptz2/S89ptz2.mat');
load('s89_221109_ptz3/S89ptz3.mat');
aID=[S89ptz1,S89ptz2,S89ptz3];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;
%%
% Pre-Ictal Only
figure
kk=1;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([900 1400])
yticks([])
pbaspect([1 .5 1])
hold off

%scale bars
sb_x=30; %30seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[925 -0.4];
sb_omr=[925 -1.2];
sb_oe=[925 2.1];
sb_od=[925 .9]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

%Seizure w/o CSD
figure
kk=2;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([250 750])
yticks([])
pbaspect([1 .5 1])
hold off
set(gca,'Visible','off')


%Seizure w/ CSD
figure
kk=3;
plot(aID(kk).EEG_ts,aID(kk).EEG./2000+3.1,'color',[0.5 0.5 0.5]);
hold on
plot(aID(kk).DC_ts,aID(kk).DC./20+1.7,'color',[0 0 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Gdff_flt2,1),'color', [0 .75 0]);
plot(aID(kk).F1_ts,mean(aID(kk).Fc1Rdff_flt2,1)-0.75,'color',[1 0 1]);
ylim([-1.5 4])
xlim([500 1000])
yticks([])
pbaspect([1 .5 1])
hold off
set(gca,'Visible','off')


%% Figure3A:indiv traces, used fo figuring out the traces i wanted to use
kk=4;
lpfilt=.05;
%green
[Fc1Gdff_flt3,Fc1Rdff_flt3]=deal(zeros(size(aID(kk).Fc1Gdff)));
for ii=1:size(aID(kk).Fc1Gdff,1)
    Fb1Gdff_flt2(ii,:)=lofi(aID(kk).Fb1Gdff(ii,:),10^6/fs_2p,1,'verbose',0);
    Fb1Rdff_flt2(ii,:)=lofi(aID(kk).Fb1Rdff(ii,:),10^6/fs_2p,1,'verbose',0);
    Fb1Gdff_flt3(ii,:)=lofi(aID(kk).Fb1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fb1Rdff_flt3(ii,:)=lofi(aID(kk).Fb1Rdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Gdff_flt3(ii,:)=lofi(aID(kk).Fc1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Rdff_flt3(ii,:)=lofi(aID(kk).Fc1Rdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
end

clear lpfiltG ii

%
for jj=1
%cellList=[jj,jj+1,jj+2];
%cellList=[36,43,61,84,90,140,144,149,154,156,158,212,215,221,223,235,236,240,274,275,276,277,279,280,290,291,292,337,340,347,352];
cellList=[84,156,290,340];
figure('Position',[100 100 400 800])
hold on
cc=1;

for ii=cellList
    plot(aID(kk).F1_ts,Fb1Rdff_flt2(ii,:)+cc,'color', [.9 0.75 .9]);
    plot(aID(kk).F1_ts,Fb1Gdff_flt2(ii,:)+cc+1,'color', [0.75 .9 0.75]);
    plot(aID(kk).F1_ts,Fb1Rdff_flt3(ii,:)+cc,'color', [0.75 0 0.75],'linewidth',1);
    plot(aID(kk).F1_ts,Fb1Gdff_flt3(ii,:)+cc+1,'color', [0 0.75 0],'linewidth',1);
    
    line([aID(kk).RecTimeG_csd(ii);aID(kk).RecTimeG_csd(ii)],repmat([cc+1;cc+1.25],1,2),'color',[0 1 1],'LineWidth',2)
    line([aID(kk).RecTimeG_csdRTN(ii);aID(kk).RecTimeG_csdRTN(ii)],repmat([cc+1;cc+1.25],1,2),'color',[0 0 1],'LineWidth',2)
    if aID(kk).isRecruitedG_csd(ii)==1
    line([aID(kk).RecTimeG_csd(ii);aID(kk).RecTimeG_csd(ii)],repmat([cc+1;cc+1.25],1,2),'color',[0 .5 .5],'LineWidth',2)
    end
    if aID(kk).isRecG_csdRTN(ii)==1
    line([aID(kk).RecTimeG_csdRTN(ii);aID(kk).RecTimeG_csdRTN(ii)],repmat([cc+1;cc+1.25],1,2),'color',[0 0 .5],'LineWidth',2)
    end
    
    
    line([aID(kk).RecTimeR_csd(ii);aID(kk).RecTimeR_csd(ii)],repmat([cc;cc+.25],1,2),'color',[0 1 1],'LineWidth',2)
    line([aID(kk).RecTimeR_csdRTN(ii);aID(kk).RecTimeR_csdRTN(ii)],repmat([cc;cc+.25],1,2),'color',[0 0 1],'LineWidth',2)
    if aID(kk).isRecruitedR_csd(ii)==1
    line([aID(kk).RecTimeR_csd(ii);aID(kk).RecTimeR_csd(ii)],repmat([cc;cc+.25],1,2),'color',[0 .5 .5],'LineWidth',2)
    end
    if aID(kk).isRecR_csdRTN(ii)==1
    line([aID(kk).RecTimeR_csdRTN(ii);aID(kk).RecTimeR_csdRTN(ii)],repmat([cc;cc+.25],1,2),'color',[0 0 .5],'LineWidth',2)
    end
    
    cc=cc+2;
end

plot(aID(kk).EEG_ts,aID(kk).EEG./2000+cc+1.5,'color',[0.5 0.5 0.5]);
plot(aID(kk).DC_ts,aID(kk).DC./20+cc+1,'color',[0 0 0]);
xlim([10 260])
ylim([0 cc+3])
title(string(jj))

%scale bars
sb_x=30; %20seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %250uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %2.5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[25 1.5];
sb_omr=[25 .5];
sb_oe=[25 8];
sb_od=[25 7.2]; %

%Scale Bars Mean Green
line([sb_omg(1) sb_omg(1)+sb_x],[sb_omg(2) sb_omg(2)],'LineWidth',1,'color',[0 0.75 0])
line([sb_omg(1) sb_omg(1)],[sb_omg(2) sb_omg(2)+sb_ymg],'LineWidth',1,'color',[0 0.75 0])

%Scale Bars Mean Red
line([sb_omr(1) sb_omr(1)+sb_x],[sb_omr(2) sb_omr(2)],'LineWidth',1,'color',[0.75 0 0.75])
line([sb_omr(1) sb_omr(1)],[sb_omr(2) sb_omr(2)+sb_ymr],'LineWidth',1,'color',[0.75 0 0.75])

%Scale Bars EEG
line([sb_oe(1) sb_oe(1)+sb_x],[sb_oe(2) sb_oe(2)],'LineWidth',1,'color',[0.5 0.5 0.5])
line([sb_oe(1) sb_oe(1)],[sb_oe(2) sb_oe(2)+sb_ye],'LineWidth',1,'color',[0.5 0.5 0.5])

%Scale Bars EEG
line([sb_od(1) sb_od(1)+sb_x],[sb_od(2) sb_od(2)],'LineWidth',1,'color',[0 0 0])
line([sb_od(1) sb_od(1)],[sb_od(2) sb_od(2)+sb_yd],'LineWidth',1,'color',[0 0 0])

set(gca,'Visible','off')

end

figure
scatter(xCoor(cellList),yCoor(cellList),100,cMap1(1:length(cellList),:),'filled')
hold on
text(xCoor(cellList),yCoor(cellList),string(1:4))
