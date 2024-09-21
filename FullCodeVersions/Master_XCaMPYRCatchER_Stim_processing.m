%MASTER SCRIPT FOR PROCESSING XCaMP-Y-RCatchER STIM DATA

%*****___STIM DATA___*****

%extracted from interleaved two channel 2PCI using suite2P
%***Fall is ch2 (r-catcher)***
%***F_chan2 is XCaMP-Y***

%use Suite2P to extract caclium transients by channel manually removing
%yellow bleedthrough cells in red channel
%load colors seperatly below


%
%
%Input .mat file from Suite2P
%Output is .mat file structure with fields to work with data for anaysis


%Navigate to parent suite2P directory for seizure processing
%Directory Tree ch1 is in suite2P from directory 1 and ch 2 has its own

%Outline
%  0) Define Variables
%  1) Open EEG from EDF and annotation files
%  2) Open DC data from .abf files
%  3) Open F and Fneu values for cells
%  4) Sort the F values into true cells and generate sorting indicies
%  5) Subtract background signal and neuropil contamination
%  6) Determine dF/F
%  7) Generate population mean traces
%  8) Power Spectral Density Determinations
%  9) Mean CSD Seed Time Determination
% 10) Filtering of Individual Traces to 1 Hz
% 11) Individual trace seizure recruitment determination (f') ---RecTimeG_csd
% 12) Define if a cell is recruited during seizure ---isRecruitedG_csd
%
%
%
%
%

%% START Loop

szList={'S86stim1','S87stim1','S88stim1','S89stim1','S102stim1','S103stim1','S105stim1'};
%szList={'S102cicr1','S105cicr1'};
%szList={'S105stim1'};

for mm=1:length(szList)

    
%szID='S60ptz2';
szID=szList{mm};
disp(szID)
plotON=0;

%PrmTbl=readtable('../RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
file_path=PrmTbl{szID,'file_path'};
PrmTbl=removevars(PrmTbl,{'szID'});
FOV_dim=PrmTbl{szID,'FOV_dim'};
EEG_On=PrmTbl{szID,'EEG_On'};
DC_On=PrmTbl{szID,'DC_On'};
EEGbias=PrmTbl{szID,'EEG_bias'};

cd(file_path{1,1})

%% Variables
fs_eeg=2000;
fs_2p=30;

%% Loader EEG from EDF file
if EEG_On==1
[EEG_all]=ParseEDF('export.edf','annotations.txt',fs_eeg,plotON);
EEG=EEG_all.part{1};
EEG_ts=[1/fs_eeg:1/fs_eeg:length(EEG)/fs_eeg];
else
EEG=NaN;
EEG_ts=NaN;
end

%% Loader DC from .abf file
if DC_On==1
    fileName = dir('*.abf');
    [DC,DC_tinv,~]=abfload(fileName.name);%DC is the vector DC_tinv is the time interval in us
    fs_dc=10^6/DC_tinv;
    if and(~isnan(EEG),length(DC)>length(EEG))
        DC=DC(1:floor(length(EEG)/fs_eeg*fs_dc));
    end
    DC_ts=[1/fs_dc:1/fs_dc:length(DC)/fs_dc];
else
   [DC,DC_ts,fs_dc]=deal(NaN);
end

%% Loader Ca Data (Red (ch 1); Green (ch 2))
load('suite2p/plane0/Fall.mat')
load('suite2p/plane0/F_chan2.mat')
load('suite2p/plane0/Fneu_chan2.mat')
F_R=F; F_G=F2; Fneu_R=Fneu; Fneu_G=F2neu;
clear F F2 Fneu F2neu spks


%% Select Only Selected ROIs
% GREEN CHANNEL
F1G=nan(sum(iscell(:,1)),size(F_G,2)-1); %generate matrix to fill with only cell ROIs F
F1R=nan(sum(iscell(:,1)),size(F_R,2)-1);
Fneu1G=nan(sum(iscell(:,1)),size(Fneu_G,2)-1); %generate matrix to fill with only cell ROIs neu
Fneu1R=nan(sum(iscell(:,1)),size(Fneu_R,2)-1);
xCoor=nan(sum(iscell(:,1)),1); %generate matrix to fill with x coordinate values for only cells green
yCoor=nan(sum(iscell(:,1)),1); %generate matrix to fill with y coordinate values for only cells green
iscellI=nan(size(iscell)); %generate matrix for a key of the iscell ROIs' original ROI ID numbers

jj=1;%counter total will equal number of true cells
for ii=1:size(F_R,1)%number of ROIs (both true and not true cells)
    if iscell(ii)==1
        F1R(jj,:)=F_R(ii,1:end-1);
        F1G(jj,:)=F_G(ii,1:end-1); %populate new cell only matrix with traces
        Fneu1R(jj,:)=Fneu_R(ii,1:end-1);
        Fneu1G(jj,:)=Fneu_G(ii,1:end-1);
        xCoor(jj)=stat{ii}.med(2); %MODIFIED TO FIX X AND Y FLIP
        yCoor(jj)=513-stat{ii}.med(1);%Modified to flip the y axis to match the images without needign to change axis orientation 
        iscellI(ii,:)=[ii,jj];
        jj=jj+1;
    end
end
clear ii
clear jj

[~,xCoorI]=sort(xCoor); %get index values for a sorted dataset with respect x coordinate
%can sort cells by this index using F(xCorrI,:); %reorder F based upon x
%coordinant
[~,yCoorI]=sort(yCoor); %same for y


%Write a timestamp for plotting
F1_ts=[1/fs_2p:1/fs_2p:length(F1G)/fs_2p];


%% Baseline shift and neuropil subtraction
%must subtract the baseline or at least get transients to roughly same
%positive range
%Fb: indicates background subtracted; Fc neuropil corrected
Fb1G=F1G-min(min(Fneu1G));%shifts F trace by the global Fneu minimum (assumed as background)
Fb1R=F1R-min(min(Fneu1R));
Fneub1G=Fneu1G-min(min(Fneu1G));% must shift the neuropil too so it can be properly scaled
Fneub1R=Fneu1R-min(min(Fneu1R));

Fc1G=Fb1G-(0.7.*Fneub1G);%gives neuropil subtracted data (serves as pure F)
Fc1R=Fb1R-(0.7.*Fneub1R);


%% Generate Mean Population Traces
F1G_Mean=mean(F1G,1);%population mean of raw F
F1R_Mean=mean(F1R,1);
Fneu1G_Mean=mean(Fneu1G,1);%population mean of raw F
Fneu1R_Mean=mean(Fneu1R,1);

Fb1G_Mean=mean(Fb1G,1);%cell population mean of background subtracted data 
Fb1R_Mean=mean(Fb1R,1);
Fneub1G_Mean=mean(Fneub1G,1);%neuropil population mean of background subtracted data 
Fneub1R_Mean=mean(Fneub1R,1);

Fc1G_Mean=mean(Fc1G,1);%population mean of background and neuropil subtracted data 
Fc1R_Mean=mean(Fc1R,1);


%% Generate dF/F
%(F-F0)/F0 with F0 being the firat 30 seconds of data
F0sec=30; %time period for F0 in seconds
Fc1Gdff=(Fc1G-mean(Fc1G(:,1:round(F0sec*fs_2p)),2))./mean(Fc1G(:,1:round(F0sec*fs_2p)),2);%dF/F
Fc1Rdff=(Fc1R-mean(Fc1R(:,1:round(F0sec*fs_2p)),2))./mean(Fc1R(:,1:round(F0sec*fs_2p)),2);
Fc1Gdff_Mean=mean(Fc1Gdff,1);%population mean of dF/F
Fc1Rdff_Mean=mean(Fc1Rdff,1);
Fc1Gdff_std=std(Fc1Gdff);%standard deviation of the population mean of dF/F
Fc1Rdff_std=std(Fc1Rdff);

%background subtracted soma F
Fb1Gdff=(Fb1G-mean(Fb1G(:,1:round(F0sec*fs_2p)),2))./mean(Fb1G(:,1:round(F0sec*fs_2p)),2);
Fb1Rdff=(Fb1R-mean(Fb1R(:,1:round(F0sec*fs_2p)),2))./mean(Fb1R(:,1:round(F0sec*fs_2p)),2);
Fb1Gdff_Mean=mean(Fb1Gdff,1);
Fb1Rdff_Mean=mean(Fb1Rdff,1);
Fb1Gdff_std=std(Fb1Gdff);
Fb1Rdff_std=std(Fb1Rdff);

%background subtracted neuropil Fneu
Fneub1Gdff=(Fneub1G-mean(Fneub1G(:,1:round(F0sec*fs_2p)),2))./mean(Fneub1G(:,1:round(F0sec*fs_2p)),2);
Fneub1Rdff=(Fneub1R-mean(Fneub1R(:,1:round(F0sec*fs_2p)),2))./mean(Fneub1R(:,1:round(F0sec*fs_2p)),2);
Fneub1Gdff_Mean=mean(Fneub1Gdff,1);
Fneub1Rdff_Mean=mean(Fneub1Rdff,1);
Fneub1Gdff_std=std(Fneub1Gdff);
Fneub1Rdff_std=std(Fneub1Rdff);


%% Filter (low pass) version 1: on dffs calculated from raw unfiltered data
lpfilt=5;
%hpfilt=.1;
%green
[Fb1Gdff_flt1,Fneub1Gdff_flt1,Fc1Gdff_flt1]=deal(zeros(size(F1G)));
for ii=1:size(F1G,1)
    Fb1Gdff_flt1(ii,:)=lofi(Fb1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fneub1Gdff_flt1(ii,:)=lofi(Fneub1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Gdff_flt1(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    %Fc1Gdff_flt1hp(ii,:)=hifi(Fc1Gdff_flt1(ii,:),10^6/F1_fs,hpfilt);
end

%red
[Fb1Rdff_flt1,Fneub1Rdff_flt1,Fc1Rdff_flt1]=deal(zeros(size(F1R)));
for jj=1:size(F1R,1)
    Fb1Rdff_flt1(jj,:)=lofi(Fb1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fneub1Rdff_flt1(jj,:)=lofi(Fneub1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Rdff_flt1(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    %Fc1Rdff_flt1hp(jj,:)=hifi(Fc1Rdff_flt1(jj,:),10^6/F1_fs,hpfilt);
end

clear lpfiltG hpfilt
clear ii jj

% Generate Means for filtered data version 1
Fb1Gdff_flt1_Mean=mean(Fb1Gdff_flt1,1);
Fb1Rdff_flt1_Mean=mean(Fb1Rdff_flt1,1);
Fneub1Gdff_flt1_Mean=mean(Fneub1Gdff_flt1,1);
Fneub1Rdff_flt1_Mean=mean(Fneub1Rdff_flt1,1);
Fc1Gdff_flt1_Mean=mean(Fc1Gdff_flt1,1);
Fc1Rdff_flt1_Mean=mean(Fc1Rdff_flt1,1);

%Fc1Gdff_flt1hp_Mean=mean(Fc1Gdff_flt1hp,1);
%Fc1Rdff_flt1hp_Mean=mean(Fc1Rdff_flt1hp,1);


%% Mean Trace Seizure and CSD Time Seeds (max first derivative block method)

%Mean Trace Filter
lpfilt=1;

%green
FGdff_MeanSm=lofi(Fc1Gdff_Mean,10^6/fs_2p,lpfilt,'verbose',0); %smoothed

%red
FRdff_MeanSm=lofi(-Fc1Rdff_Mean,10^6/fs_2p,lpfilt,'verbose',0); %smoothed

clear lpfilt

% determine all the points that are greater than the average of the min and
% max values of the trace

temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),max(FGdff_MeanSm)]);
temp1R=FRdff_MeanSm>=mean([min(FRdff_MeanSm),max(FRdff_MeanSm)]);

% find block bounds
temp2G=diff(temp1G,1,2);
temp2R=diff(temp1R,1,2);

temp3G=cell(1,2);
temp3R=cell(1,2);
temp3G{1}=find(temp2G==1);
temp3G{2}=find(temp2G==-1);
temp3R{1}=find(temp2R==1);
temp3R{2}=find(temp2R==-1);

%integrate over blocks
temp4G=cell(1,1);
temp4R=cell(1,1);
if temp3G{1}(1)>temp3G{2}(1)%correct for first slope value index being negative
    temp3G{2}=temp3G{2}(2:end);
end

if temp3G{1}(end)>temp3G{2}(end)%correct for extra final slope value being positive
    temp3G{1}=temp3G{1}(1:end-1);
end

for jj=1:length(temp3G{1}) %find the area under the curve of each region
    temp4G{1}(jj)=trapz(FGdff_MeanSm(1,[temp3G{1}(jj):temp3G{2}(jj)]));
end
clear jj

if temp3R{1}(1)>temp3R{2}(1)%correct for first slope value index being negative
    temp3R{2}=temp3R{2}(2:end);
end

if temp3R{1}(end)>temp3R{2}(end)%correct for extra final slope value being positive
    temp3R{1}=temp3R{1}(1:end-1);
end

for jj=1:length(temp3R{1}) %find the area under the filtered curve of each positive region
    temp4R{1}(jj)=trapz(FRdff_MeanSm(1,[temp3R{1}(jj):temp3R{2}(jj)]));
end

clear ii jj


% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(1,2);%gives the block in each cell that is seizure
temp5R=cell(1,2);%gives the block in each cell that is seizure

[temp5G{1},temp5G{2}] = sort(temp4G{1},1,'descend');%sort blocks by size (block index is cell 2)

%determine the block bounds of the seizure and CSD blocks
MeanCSDRecTimeG = (temp3G{1}(temp5G{2}(1))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks

[temp5R{1},temp5R{2}] = sort(temp4R{1},1,'descend');

%determine the block bounds of the seizure and CSD blocks
MeanCSDRecTimeR = (temp3R{1}(temp5R{2}(1))+1)/fs_2p;


%% Filtering for Individual Trace Recruitment Time Determination
% filtering to smooth for recruitment detection
lpfilt=1;
%hpfilt=0.1;
Fc1Gdff_flt2=zeros(size(Fc1Gdff));
%Fc1Gdff_flt2_hp=zeros(size(Fc1Gdff_flt2));
Fc1Rdff_flt2=zeros(size(Fc1Rdff));
%Fc1Rdff_flt2_hp=zeros(size(Fc1Rdff_flt2));
%green
for ii=1:size(Fc1Gdff_flt1,1)
    Fc1Gdff_flt2(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    %Fc1Gdff_flt2_hp(ii,:)=hifi(Fc1Gdff_flt2(ii,:),10^6/fs_2p,hpfilt);
end

%red
for jj=1:size(Fc1Rdff_flt1,1)
    Fc1Rdff_flt2(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    %Fc1Rdff_flt2_hp(jj,:)=hifi(Fc1Rdff_flt2(jj,:),10^6/fs_2p,hpfilt);
end

clear lpfilt
clear ii jj


%% Individual CSD recruitment detection

MeanCSDRecTimeR=MeanCSDRecTimeG;

RecTimeG_csd_struct=IndivRecTimes3(Fc1Gdff_flt2,MeanCSDRecTimeG,fs_2p,5);%average speed of 40um/s gives means it takes 7s to cross 280um feidl of view, therefore a sigma of 5s gives us a 10s window over which 
RecTimeR_csd_struct=IndivRecTimes3(-Fc1Rdff_flt2,MeanCSDRecTimeR,fs_2p,10);%by seeding R with G mean want to increse window for detection a bit to account for possibliy of variation in time.
RecTimeG_csd=RecTimeG_csd_struct.time;
RecTimeR_csd=RecTimeR_csd_struct.time;


%% Recruited cell or not (CSD)
%trace needs to have atleast a 20% increase in signal during CSD relative to before and the
%CSD detected time needs to be determined to be within 10 seconds of
%the seed
%NOTE:Invert the Red channel signal for CSD

tprd=5;%time period (s)
gthresh=0.2; 
rthresh=0.15;
CSDMeanDffDiffG=nan(size(RecTimeG_csd));
for ii=1:length(RecTimeG_csd)
    preCSDrange=floor([(RecTimeG_csd(ii)-tprd)*fs_2p:RecTimeG_csd(ii)*fs_2p]);
    postCSDrange=floor([RecTimeG_csd(ii)*fs_2p:(RecTimeG_csd(ii)+tprd)*fs_2p]);
    preCSDMeanDff=mean(Fc1Gdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff=mean(Fc1Gdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffG(ii)=postCSDMeanDff-preCSDMeanDff;
end
isRecruitedG_csd=and(CSDMeanDffDiffG>gthresh,abs(RecTimeG_csd-MeanCSDRecTimeG)<FOV_dim/45); 
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff

CSDMeanDffDiffR=nan(size(RecTimeR_csd));
for ii=1:length(RecTimeR_csd)
    preCSDrange=floor([(RecTimeR_csd(ii)-tprd)*fs_2p:RecTimeR_csd(ii)*fs_2p]);
    postCSDrange=floor([RecTimeR_csd(ii)*fs_2p:(RecTimeR_csd(ii)+tprd)*fs_2p]);
    preCSDMeanDff=mean(-Fc1Rdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff=mean(-Fc1Rdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffR(ii)=postCSDMeanDff-preCSDMeanDff;
end
isRecruitedR_csd=and(CSDMeanDffDiffR>rthresh,abs(RecTimeR_csd-MeanCSDRecTimeR)<FOV_dim/45);
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff
clear tprd gthresh rthresh




%% Stim End Time

% PSD initial processing for feature detection
%used for finding pre-ictal spikes and also death following seizure death
if DC_On==1
PSDbin=1024;
[~,PSD_F]=pwelch(DC,[],[],PSDbin,fs_eeg);
PSDwin=1;%in seconds

DCstart=[1:PSDwin*fs_dc/2:length(DC)];%sliding window width 2x PSDwin shifting by 1/2 PSDwin
PSD_P=nan(length(PSD_F),length(DCstart));
for ii = 1:numel(DCstart)-2
    DC1=DC(DCstart(ii):DCstart(ii)+PSDwin*fs_dc-1);
    [PSD_P(:,ii),~]=pwelch(DC1,[],[],PSDbin,fs_dc);
    clear EEG1
end

elseif EEG_On==1
PSDbin=1024;
[~,PSD_F]=pwelch(EEG,[],[],PSDbin,fs_eeg);
PSDwin=1;%in seconds

EEGstart=[1:PSDwin*fs_eeg/2:length(EEG)];%sliding window width 2x PSDwin shifting by 1/2 PSDwin
PSD_P=nan(length(PSD_F),length(EEGstart));
for ii = 1:numel(EEGstart)-2
    EEG1=EEG(EEGstart(ii):EEGstart(ii)+PSDwin*fs_eeg-1);
    [PSD_P(:,ii),~]=pwelch(EEG1,[],[],PSDbin,fs_eeg);
    clear EEG1
end

end

%
PSD_totP=sum(PSD_P,1);
PSD_stimArt=sum(PSD_P(or(and(PSD_F>50,PSD_F<70),and(PSD_F>170,PSD_F<190)),:),1);% +/- 30 second window roll avg (time step here is 0.5s)

%determines stim start and end time using normalized power in the artifact frequency range
StimStartTime=find(PSD_stimArt(1:60/PSDwin*2)./max(PSD_stimArt)>0.5,1)*PSDwin/2;
StimEndTime=find(PSD_stimArt(StimStartTime/PSDwin*2:60/PSDwin*2)./max(PSD_stimArt)<0.5,1)*PSDwin/2+StimStartTime;


%% DC anaysis and Recovery signal from CSD
%grabbing the approx end of the csd return of signal to baseline Ca


% filtering to smooth for recruitment detection
lpfilt=.05;
F1Gdff_csdEnd=zeros(size(Fc1Gdff));
F1Rdff_csdEnd=zeros(size(Fc1Rdff));

FGdff_MeanSm=lofi(Fb1Gdff_Mean,10^6/fs_2p,.1,'verbose',0);
FRdff_MeanSm=lofi(Fb1Rdff_Mean,10^6/fs_2p,.1,'verbose',0);

%green
for ii=1:size(F1Gdff_csdEnd,1)
    F1Gdff_csdEnd(ii,:)=lofi(Fb1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
end
%red
for jj=1:size(F1Rdff_csdEnd,1)
    F1Rdff_csdEnd(jj,:)=lofi(Fb1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
end
clear lpfilt
clear ii jj

if DC_On==1
[~,DCpostCSDmaxI]=max(SternRollAvg(DC(round((MeanCSDRecTimeG+30)*fs_dc):round((MeanCSDRecTimeG+130)*fs_dc)),fs_dc));
CSDoffTimeDC=MeanCSDRecTimeG+30+DCpostCSDmaxI/fs_dc;

tempRange=round(StimEndTime*fs_dc):round((MeanCSDRecTimeG)*fs_dc);
if strcmp(szID,'S89stim1')
    tempRange=round(StimEndTime*fs_dc):round((MeanCSDRecTimeG+10)*fs_dc);
end
tempTrace=SternRollAvg(DC(tempRange),2*fs_dc);
tempTrace_ds=downsample(tempTrace,fs_dc);
tempTrace_sp=spline(1:1:length(tempTrace_ds),tempTrace_ds,1:1/10:length(tempTrace_ds));
[~,DCstartCSDmaxI]=min(diff(tempTrace_sp,2,2));%finds the elbow of the dc trace at start of csd

tempRange2=round((MeanCSDRecTimeG+10)*fs_dc):round((MeanCSDRecTimeG+60)*fs_dc);
if strcmp(szID,'S89stim1')
    tempRange2=round((MeanCSDRecTimeG+30)*fs_dc):round((MeanCSDRecTimeG+60)*fs_dc);
elseif strcmp(szID,'S102stim1')
    tempRange2=round((MeanCSDRecTimeG+20)*fs_dc):round((MeanCSDRecTimeG+60)*fs_dc);
end
tempTrace2=SternRollAvg(DC(tempRange2),2*fs_dc);
tempTrace2_ds=downsample(tempTrace2,fs_dc);
tempTrace2_sp=spline(1:1:length(tempTrace2_ds),tempTrace2_ds,1:1/10:length(tempTrace2_ds));
[~,DCendCSDmaxI]=max(diff(tempTrace2_sp(11:end-10),2,2));%finds the elbow of the dc trace at start of csd
[~,DCpostCSDmaxI]=max(SternRollAvg(DC(round((MeanCSDRecTimeG+30)*fs_dc):round((MeanCSDRecTimeG+130)*fs_dc)),fs_dc));

CSDstartTimeDC=StimEndTime+(DCstartCSDmaxI-1)/fs_dc*200;
CSDendTimeDC=MeanCSDRecTimeG+10+(DCendCSDmaxI+10)/fs_dc*200;
if strcmp(szID,'S89stim1')
    CSDendTimeDC=CSDendTimeDC+20;
elseif strcmp(szID,'S102stim1')
    CSDendTimeDC=CSDendTimeDC+10;
end
CSDoffTimeDC=MeanCSDRecTimeG+30+(DCpostCSDmaxI-1)/fs_dc;

DCShiftPre=mean(DC(round(StimEndTime*fs_dc):round(CSDstartTimeDC*fs_dc)));
DCShiftPost=min(SternRollAvg(DC(round(CSDstartTimeDC*fs_dc):round(CSDendTimeDC*fs_dc)),fs_dc));
CSD_dcShift=[DCShiftPost-DCShiftPre,DCShiftPost];%the difference and the raw average shift
else %manually seed times if missing DC trace
    if strcmp(szID,'S88stim1')
    CSDoffTimeDC=160;
    CSDstartTimeDC=nan;
    CSDendTimeDC=nan;
    CSD_dcShift=nan;
    end
end

%Find all contiguous periods of positive slope within a range of time
%find the largest of these periods (slope integral)
%index this by the max second derivative during this time

SeedTime=CSDoffTimeDC; %seconds
sigma1=40; %seconds for length of standard deviation around seed time for detecting signal

%generate gaussian centered at seed time
trace1=1:length(FGdff_MeanSm);
gauss1=exp(-(trace1-SeedTime*fs_2p).^2/(sigma1*fs_2p)^2);
gauss2=[gauss1(1:ceil(SeedTime*fs_2p)),zeros(1,length(gauss1)-ceil(SeedTime*fs_2p))];%half gaussian as peak dc is after calcium signal

%calculate the derivates of the filtered traces
FGdff_MeanSm_df=diff(FGdff_MeanSm,1,2);
FGdff_MeanSm_d2f=diff(FGdff_MeanSm,2,2);
FGdff_MeanSm_d2fgauss=FGdff_MeanSm_d2f.*gauss2(3:end);
FRdff_MeanSm_df=diff(FRdff_MeanSm,1,2);
FRdff_MeanSm_d2f=diff(FRdff_MeanSm,2,2);
FRdff_MeanSm_d2fgauss=FRdff_MeanSm_d2f.*gauss2(3:end);

% determine all the points that are negative for green and positive for red
temp1G=FGdff_MeanSm_df<=0;
temp1R=FRdff_MeanSm_df>=0;

% find block bounds
temp2G=diff(temp1G,1,2);
temp2R=diff(temp1R,1,2);

temp3G=cell(1,2);
temp3G{1,1}=find(temp2G==1);
temp3G{1,2}=find(temp2G==-1);
temp3R=cell(1,2);
temp3R{1,1}=find(temp2R==1);
temp3R{1,2}=find(temp2R==-1);

%integrate over blocks
[temp4G,minD2FwinG,temp5G,temp6G]=deal(cell(1,1));
[temp4R,maxD2FwinR,temp5R,temp6R]=deal(cell(1,1));

if temp3G{1,1}(1)>temp3G{1,2}(1)%correct for first slope value index being negative
    temp3G{1,2}=temp3G{1,2}(2:end);
end
    
if temp3G{1,1}(end)>temp3G{1,2}(end)%correct for extra final slope value being positive
        temp3G{1,1}=temp3G{1,1}(1:end-1);
end

if temp3R{1,1}(1)>temp3R{1,2}(1)%correct for first slope value index being negative
    temp3R{1,2}=temp3R{1,2}(2:end);
end
    
if temp3R{1,1}(end)>temp3R{1,2}(end)%correct for extra final slope value being positive
        temp3R{1,1}=temp3R{1,1}(1:end-1);
end
    
for jj=1:length(temp3G{1,1}) %find the area under the first derivative curve of each postitive region
    temp4G{1}(jj)=trapz(FGdff_MeanSm_df([temp3G{1,1}(jj):temp3G{1,2}(jj)]));
    [minD2FwinG{1}(jj,1),minD2FwinG{1}(jj,2)]=min(FGdff_MeanSm_d2fgauss([temp3G{1,1}(jj):temp3G{1,2}(jj)]));
end
clear jj

for jj=1:length(temp3R{1,1}) %find the area under the first derivative curve of each postitive region
    temp4R{1}(jj)=trapz(FRdff_MeanSm_df([temp3R{1,1}(jj):temp3R{1,2}(jj)]));
    [maxD2FwinR{1}(jj,1),maxD2FwinR{1}(jj,2)]=max(FRdff_MeanSm_d2fgauss([temp3R{1,1}(jj):temp3R{1,2}(jj)]));
end
clear jj

%recalculate block index times based upon max second derivative time 
temp5G{1,1}=temp3G{1,1}+minD2FwinG{1}(:,2)';
temp5R{1,1}=temp3R{1,1}+maxD2FwinR{1}(:,2)';
%recalculate the features weighted by the gaussian kernel
temp6G{1,1}=temp4G{1}.*gauss2(temp5G{1,1});
temp6R{1,1}=temp4R{1}.*gauss2(temp5R{1,1});

% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp5 to find the index (time point) of the event in the origional trace
[~, feats_minIG] = min(temp6G{1,1});
MeanCSDendTimeG = (temp5G{1,1}(feats_minIG)+1)/fs_2p;

[~, feats_maxIR] = max(temp6R{1,1});
MeanCSDendTimeR = (temp5R{1,1}(feats_maxIR)+1)/fs_2p;


if abs(MeanCSDendTimeR-MeanCSDendTimeG)>7
MeanCSDendTimeR=MeanCSDendTimeG;
end

%MeanCSDendTimeR=MeanCSDendTimeG;

%Closest 2nd deriv peak detection indiv
CSDend_trange=FOV_dim/30;%time range dimension scaled for FOV size (equivalent to half the time needed to cross the hypotenuse of the field if traveling at 20um/s (lower speed limit)

tprdG=CSDend_trange*2;
CSDtempI=round((mean([MeanCSDendTimeG,MeanCSDendTimeR]))*fs_2p);
CSDtempIstart=floor(CSDtempI-tprdG*fs_2p);

RecTimeG_csdRTN=nan(size(F1Gdff_csdEnd,1),1);
tempTraceMat=diff(F1Gdff_csdEnd(:,CSDtempIstart:floor((CSDtempI+tprdG*fs_2p))),2,2);
for ii=1:size(tempTraceMat,1)
    [~,tempPeaksI]=findpeaks(-tempTraceMat(ii,:),'MinPeakHeight',std(tempTraceMat(ii,:)));%negative to grab local minima
    if ~isempty(tempPeaksI)
        [~,peakMinDistI]=sort(abs(tempPeaksI-CSDtempI+CSDtempIstart+2),'ascend');
        RecTimeG_csdRTN(ii)=(CSDtempIstart+tempPeaksI(peakMinDistI(1))+2)/fs_2p;
    else
        RecTimeG_csdRTN(ii)=nan;
    end  
end

tprdR=tprdG;
RecTimeR_csdRTN=nan(size(F1Rdff_csdEnd,1),1);
tempTraceMat=diff(F1Rdff_csdEnd(:,CSDtempIstart:floor((CSDtempI+tprdR*fs_2p))),2,2);
for ii=1:size(tempTraceMat,1)
    [~,tempPeaksI]=findpeaks(tempTraceMat(ii,:),'MinPeakHeight',std(tempTraceMat(ii,:)));
    if ~isempty(tempPeaksI)
        [~,peakMinDistI]=sort(abs(tempPeaksI-CSDtempI+CSDtempIstart+2),'ascend');
        RecTimeR_csdRTN(ii)=(CSDtempIstart+tempPeaksI(peakMinDistI(1))+2)/fs_2p;
    else
        RecTimeR_csdRTN(ii)=nan;
    end  
end

%GREEN min 2nd deriv SEEDED RED DETECTION INDIV
% Indiv RecTimes 
% elbow of signal (2nd deriv max within range around mean seed times)
% tprdG=CSDend_trange*2;%default ptz 30 default stim 20
% CSDtempI=round((MeanCSDendTimeG)*fs_2p);
% [~,RecTime_chng]=min(diff(F1Gdff_csdEnd(:,(CSDtempI-tprdG*fs_2p):(CSDtempI+tprdG*fs_2p)),2,2),[],2);
% RecTimeG_csdRTN=(CSDtempI-tprdG*fs_2p+RecTime_chng+2)/fs_2p;
% 
% 
% tprdR=5;
% RecTimeR_csdRTN=nan(size(F1Rdff_csdEnd,1),1);
% CSDtempI=round(RecTimeG_csdRTN*fs_2p);
% for ii=1:size(F1Rdff_csdEnd,1)
%     [~,RecTime_chng]=max(diff(F1Rdff_csdEnd(ii,(CSDtempI(ii)-tprdR*fs_2p):(CSDtempI(ii)+tprdR*fs_2p)),2));
%     RecTimeR_csdRTN(ii)=(CSDtempI(ii)-tprdR*fs_2p+RecTime_chng+2)/fs_2p;
%     clear  RecTime_chng
% end

%RED DETECTION min 2nd deriv INDIV
% tprdR=CSDend_trange*2;%default ptz 30 default stim 20
% CSDtempI=round((MeanCSDendTimeR)*fs_2p);
% [~,RecTime_chng]=max(diff(F1Rdff_csdEnd(:,(CSDtempI-tprdR*fs_2p):(CSDtempI+tprdR*fs_2p)),2,2),[],2);
% RecTimeR_csdRTN=(CSDtempI-tprdR*fs_2p+RecTime_chng+2)/fs_2p;

%GAUSSIAN WEIGHTED VERSION INDIV

% gauss3=exp(-(trace1-MeanCSDendTimeG*fs_2p).^2/(15*fs_2p)^2);
% [~,CSDtempI]=min(diff(F1Gdff_csdEnd,2,2).*gauss3(3:end),[],2);
% [~,RecTime_chng]=min(diff(F1Gdff_csdEnd(:,(CSDtempI-5*fs_2p):(CSDtempI+5*fs_2p)),2,2),[],2);
% RecTimeG_csdRTN=(CSDtempI-5*fs_2p+RecTime_chng+2)./fs_2p;
% 
% gauss3=exp(-(trace1-MeanCSDendTimeR*fs_2p).^2/(15*fs_2p)^2);
% [~,CSDtempI]=max(diff(F1Rdff_csdEnd,2,2).*gauss3(3:end),[],2);
% [~,RecTime_chng]=max(diff(F1Rdff_csdEnd(:,(CSDtempI-5*fs_2p):(CSDtempI+5*fs_2p)),2,2),[],2);
% RecTimeR_csdRTN=(CSDtempI-5*fs_2p+RecTime_chng+2)./fs_2p;


% gauss3=exp(-(trace1-MeanCSDendTimeR*fs_2p).^2/(10*fs_2p)^2);
% figure
% plot(diff(F1Rdff_csdEnd(1,:),2,2)*10000,'r')
% hold on
% plot(F1Rdff_csdEnd(1,:),'b')
% plot(gauss3(3:end),'g')
% plot(diff(F1Rdff_csdEnd(1,:),2,2).*gauss3(3:end)*100000,'m')

%
% % IsRecCSD RTN
% % difference before and after within time range of mean exceeds threshold
% if strcmp(szID,'S105stim1') || strcmp(szID,'S88stim1')%to restrict to the most robust cells in this noisy recording
%     gthresh=0.05; %.05
%     rthresh=0.125; %.125 
% elseif strcmp(szID,'S89stim1') || strcmp(szID,'S103stim1')%s89,s103 ;;;;;decide on s102 parameters
%     gthresh=0.1;
%     rthresh=0.05;
% elseif strcmp(szID,'S87stim1')
%     gthresh=0.1;
%     rthresh=0.07;
% else
%     gthresh=0.05;
%     rthresh=0.05;
% end
% 
%     gthresh=0.01;
%     rthresh=0.08;

%gthresh=(mean(FGdff_MeanSm2(floor([(MeanCSDendTimeG-CaLvl_tprd)*fs_2p:MeanCSDendTimeG*fs_2p])),2)-mean(FGdff_MeanSm2(floor([MeanCSDendTimeG*fs_2p:(MeanCSDendTimeG+CaLvl_tprd)*fs_2p])),2))*1;
%rthresh=(mean(FRdff_MeanSm2(floor([MeanCSDendTimeR*fs_2p:(MeanCSDendTimeR+CaLvl_tprd)*fs_2p])),2)-mean(FRdff_MeanSm2(floor([(MeanCSDendTimeR-CaLvl_tprd)*fs_2p:MeanCSDendTimeR*fs_2p])),2))*1;

CaLvl_tprd=30;
CSDrtnDiffG=nan(size(RecTimeG_csdRTN));
for ii=1:length(RecTimeG_csdRTN)
    if ~isnan(RecTimeG_csdRTN(ii))
    preCSDrange=floor([(RecTimeG_csdRTN(ii)-CaLvl_tprd)*fs_2p:RecTimeG_csdRTN(ii)*fs_2p]);
    postCSDrange=floor([RecTimeG_csdRTN(ii)*fs_2p:(RecTimeG_csdRTN(ii)+CaLvl_tprd)*fs_2p]);
    preCSDMeanDff=mean(F1Gdff_csdEnd(ii,preCSDrange),2);
    postCSDMeanDff=mean(F1Gdff_csdEnd(ii,postCSDrange),2);
    CSDrtnDiffG(ii)=preCSDMeanDff-postCSDMeanDff; %Pre>Post
    else
    CSDrtnDiffG(ii)=nan;
    end
end
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff

CSDrtnDiffR=nan(size(RecTimeR_csdRTN));
for ii=1:length(RecTimeR_csdRTN)
    if ~isnan(RecTimeR_csdRTN(ii))
    preCSDrange=floor([(RecTimeR_csdRTN(ii)-CaLvl_tprd)*fs_2p:RecTimeR_csdRTN(ii)*fs_2p]);
    postCSDrange=floor([RecTimeR_csdRTN(ii)*fs_2p:(RecTimeR_csdRTN(ii)+CaLvl_tprd)*fs_2p]);
    preCSDMeanDff=mean(F1Rdff_csdEnd(ii,preCSDrange),2);
    postCSDMeanDff=mean(F1Rdff_csdEnd(ii,postCSDrange),2);
    CSDrtnDiffR(ii)=postCSDMeanDff-preCSDMeanDff; %Pre<Post
    else
    CSDrtnDiffR(ii)=nan;
    end
end
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff

gthresh=std(CSDrtnDiffG(and(CSDrtnDiffG<mean(CSDrtnDiffG,'omitnan'),CSDrtnDiffG>0)));
rthresh=std(CSDrtnDiffR(and(CSDrtnDiffR<mean(CSDrtnDiffR,'omitnan'),CSDrtnDiffR>0)));

gthresh=.01;
rthresh=.02;

isRecG_csdRTN=and(CSDrtnDiffG>gthresh,abs(RecTimeG_csdRTN-mean(RecTimeG_csdRTN,'omitnan'))<CSDend_trange);%15
isRecR_csdRTN=and(CSDrtnDiffR>rthresh,abs(RecTimeR_csdRTN-mean(RecTimeR_csdRTN,'omitnan'))<CSDend_trange);%15

%isRecG_csdRTN=and(CSDrtnDiffG>gthresh,abs(RecTimeG_csdRTN-MeanCSDendTimeG)<CSDend_trange);%15
%isRecR_csdRTN=and(CSDrtnDiffR>rthresh,abs(RecTimeR_csdRTN-MeanCSDendTimeR)<CSDend_trange);%15
%isRecR_csdRTN=and(CSDrtnDiffR>rthresh,and(RecTimeR_csdRTN-MeanCSDendTimeR<5,RecTimeR_csdRTN-MeanCSDendTimeR>-10));

%% RecTime PLOTS
if plotON==1
tempG=Fc1Gdff_flt2; tempR=Fc1Rdff_flt2;

%Mean PIS, Sz and CSD traces and times 
figure
suptitle('Means')
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
plot(DC_ts,DC./10+2,'color',[0.2 0.2 0.2]);
plot(F1_ts,Fc1Gdff_flt1_Mean,'color', [0 1 0]);
plot(F1_ts,Fc1Rdff_flt1_Mean-1,'color',[1 0 1]);
line([MeanSzRecTimeG;MeanSzRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-')
line([MeanSzRecTimeR;MeanSzRecTimeR],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','-')
line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-.')
line([MeanCSDRecTimeR;MeanCSDRecTimeR],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','-.')

% Plot CSD indiv rec times

% Individual cell csd/death recruitment times
figure
suptitle('CSD')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    if isRecruitedG_csd(ii)==1 %incldG(ii)==1
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    else
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 1]);
    end
    line([RecTimeG_csd(ii);RecTimeG_csd(ii)],repmat([ii*10;ii*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end

subplot(1,2,2)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for jj=1:size(Fc1Rdff_flt2,1)
    if isRecruitedR_csd(jj)==1
        plot(F1_ts,tempR(jj,:)+jj*10,'color',[1 0.25 1]);
    else
        plot(F1_ts,tempR(jj,:)+jj*10,'color',[1 0.75 1]);
    end
    line([RecTimeR_csd(jj);RecTimeR_csd(jj)],repmat([jj*10;jj*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanCSDRecTimeR;MeanCSDRecTimeR],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj

end
% Generate Data Structure
%General Variables
StructFieldNames={'szID'
    'FOV_dim'
    %Time/Sampling
    'fs_2p'
    'fs_eeg'
    'EEG_ts'
    'F1_ts'
    %EEG
    'EEG'
    %DC
    'DC'
    'fs_dc'
    'DC_ts'
    %F1 processing
    'stat'
    'xCoor'
    'yCoor'
    'xCoorI'
    'yCoorI'
    'iscell'
    'iscellI'
    'F1G'
    'Fneu1G'
    'Fb1G'
    'Fneub1G'
    'Fc1G'
    'Fb1Gdff'
    'Fneub1Gdff'
    'Fc1Gdff'
    'Fneub1Gdff_flt1'
    'Fc1Gdff_flt1'
    'Fneub1Gdff_flt1'
    'Fc1Gdff_flt2'
    'F1R'
    'Fneu1R'
    'Fb1R'
    'Fneub1R'
    'Fc1R'
    'Fb1Rdff'
    'Fneub1Rdff'
    'Fc1Rdff'
    'Fneub1Rdff_flt1'
    'Fc1Rdff_flt1'
    'Fneub1Rdff_flt1'
    'Fc1Rdff_flt2'
    %Features
    'StimStartTime'
    'StimEndTime'
    'MeanCSDRecTimeG'
    'MeanCSDRecTimeR'
    'CSDoffTimeDC'
    'CSDstartTimeDC'
    'CSDendTimeDC'
    'CSD_dcShift'
    'MeanCSDendTimeG'
    'MeanCSDendTimeR'
    'RecTimeG_csd'
    'RecTimeR_csd'
    'RecTimeG_csdRTN'
    'RecTimeR_csdRTN'
    'isRecruitedG_csd'
    'isRecruitedR_csd'
    'isRecG_csdRTN'
    'isRecR_csdRTN'
    'CSDMeanDffDiffG'
    'CSDMeanDffDiffR'};    
for kk=1:length(StructFieldNames)
    eval(strcat('szStruct.',StructFieldNames{kk},'=',StructFieldNames{kk},';'));
end

% SAVE Structure Varaible as the SeizureID (szID)
assignin('base', szID, szStruct);
save(strcat(szID,'.mat'),szID,'-v7.3')
disp('saved')

cd ..
clearvars -except mm szList

end
