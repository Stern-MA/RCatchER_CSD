%% Figure 2: Version s89
load('S89ptz3.mat');
aID=[S89ptz3,S89ptz3];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;

%% Fig2A-C:Traces
figure
kk=1;

%Seizure w/ CSD
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

%scale bars
sb_x=30; %30seconds
sb_ymg=.25; %0.25dff
sb_ymr=.25; %0.25dff
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV) accounting for scaling below divide by 2000
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV) accounting for scaling below divide by 20
sb_omg=[525 -0.5];
sb_omr=[525 -1.3];
sb_oe=[525 2];
sb_od=[525 1.1]; %

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

print('89ptz3_TraceSzCSD.svg','-dsvg')


%% Fig2 A-C: Rasters
%Sz w/ CSD
figure
kk=1;
ax(1)=subplot(2,1,1);
%Green Ca Raster
[~,xCoorI]=sort(aID(kk).xCoor);
imagesc(aID(kk).F1_ts,size(aID(kk).Fc1Gdff_flt2,1),aID(kk).Fc1Gdff_flt2(xCoorI,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
%c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
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
%c1=colorbar('Ticks',[-3:1:6],'TickLabels',[-3:1:6]);
%c1.Label.String='dF/F';
caxis([-2 6])
xticks([0:50:1050])
xticklabels([0:50:1050])
yticks([141:30:283])
yticklabels([0:30:300])
ax1=gca;
ax1.FontSize=11;
pbaspect([3 1 1])
xlim([500 1000])

print('89ptz3_RasterSzCSD.svg','-dsvg')



