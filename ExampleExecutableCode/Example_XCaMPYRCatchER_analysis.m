%EXAMPLE SCRIPT FOR ANALYZING XCaMP-Y-RCatchER PTZ DATA
%
% Code to generate figure panels depicting the example recording of a
% seizure with a post-ictal CSD including the changes in calcium and the
% vectors of propegation detrmined from the recruitment times of individual
% cells to the onset and offset of the CSD event.
%
% Load the output from the processing script (must have run processing
% script with the output file 'S89ptz3.mat' saved in the directory where
% this code is being executed.
%
% Make sure all the SpatialLinearRegression directory scripts have been
% added to the MATLAB path


%% Figure 2: Version s89
load('S89ptz3.mat');
aID=[S89ptz3,S89ptz3];
fs_2p=30;
fs_eeg=2000;
fs_dc=2000;

%% Traces
%EEG (grey) DC (black) and mean population calcium traces (XCaMP-Y in green
%and RCatchER in magenta)
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
sb_ye=.25; %500uV (=#*2000; .5uV*2000=1mV)
sb_yd=.25; %5mV (=#*20; .25mV*20=5mV)
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

print('s89ptz3_TraceSzCSD.svg','-dsvg')


%% Rasters
%Colored raster plots of individual cell calcium transient data
%
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

print('s89ptz3_RasterSzCSD.svg','-dsvg')


%% Spatial Linear Regression
%Spatial Linear Regression method for microeelctrode arrays adapted from Liou et al 2016 J Neural Eng
%14:044001 for calcium data. Here vectors of propgation during CSDs are
%determined in each channel for CSD onset and offset


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


%% ROI Field of View Colormap of Recrutiment Times
%2x2 array with onset in the top row and offset in the bottom row XCaMP-Y
%in the left column and RCatchER in the right column

kk=1;

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

print('s89ptz3_ROIcMap.svg','-dsvg')

%% Vector Plots (using compass)
% compass plot depicting the velocity (direction and speed) of the CSD  
% onset (solid line) and offet vectors (dashed line) in both the XCaMP-Y 
% (green) and RCatchER (magenta) channels

pthresh=0.05;
max_lim = 100;
x_fake=[0 max_lim 0 -max_lim];
y_fake=[max_lim 0 -max_lim 0];

figure('position',[100 100 600 600])
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

print('s89ptz3_Polar.svg','-dsvg')


%% Projection of Recruitment Times on Propogation Axis
% plots depicting the recruimtent times of indivudal cells ordered by their 
% position along the CSD vector of propegation axis
% onset (left plot) and offset (right plot) in both the XCaMP-Y in green
% and RCatchER in magenta


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

print('s89ptz3_RelRecAxis.svg','-dsvg')




