%MASTER SCRIPT FOR PROCESSING XCaMP-Y-RCatchER PTZ DATA

%*****___PTZ DATA___*****

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
%  2) Open F and Fneu values for cells
%  3) Sort the F values into true cells and generate sorting indicies
%  4) Subtract background signal and neuropil contamination
%  5) Determine dF/F
%  6) Generate population mean traces
%  7) Correct for EEG Bias (time shift) between 2p and eeg
%  8) Power Spectral Density Determinations
%  9) Death Determination
% 10) Mean Seziure and CSD Seed Time Determination
% 11) Mean Pre-ictal Spike Seed Time Determiantion ---EEG_PIS_times
% 12) EEG Kernal Generation and Refinement of EEG PIS detection ---EEG_PIS_times2
% 13) First PIS determination and Final PIS SEED adjustment ---EEG_PIS_times2
% 14) Filtering of Individual Traces to 1 Hz
% 15) Individual trace seizure recruitment determination (f') ---RecTimeG_sz2
% 16) Define if a cell is recruited during seizure ---isRecruitedG
% 17) Individual trace PIS recruitment determination (f') ---RecTimeG_PIS.v3
% 18) Define if spike included in analysis (>10% cell rec) ---EEG_PIS_times3, recTimeG_PIS.true
%
%
%
%
%


%% START Loop

szList={'S56ptz1','S86ptz6','S87ptz8','S104ptz3'};
%szList={'S105stim1'};

for mm=1:length(szList)

    
%szID='S60ptz2';
szID=szList{mm};
disp(szID)
plotON=0;

PrmTbl=readtable('RCatchER_loadParam.csv');%loads the parameters for the recording analysis
PrmTbl.Properties.RowNames = string(PrmTbl{:,'szID'});
file_path=PrmTbl{szID,'file_path'};
PrmTbl=removevars(PrmTbl,{'szID'});
FOV_dim=PrmTbl{szID,'FOV_dim'};
EEG_On=PrmTbl{szID,'EEG_On'};
DC_On=PrmTbl{szID,'DC_On'};
EEGbias=PrmTbl{szID,'EEG_bias'};
PIS_On=PrmTbl{szID,'PIS_On'};
Sz_On=PrmTbl{szID,'Sz_On'};
TSW_On=PrmTbl{szID,'TSW_On'};
split_On=PrmTbl{szID,'split_On'}; %Indicate if the recording needs to be split
split_time=PrmTbl{szID,'split_time'}; %Indicate the time (in seconds) for the split to occur
split_sel=PrmTbl{szID,'split_sel'}; %Select the first or second segment for processing
%s102_230530 split time is 616s between first and second seizure so the
%interictal spieks are included in the second seizure recording

cd(file_path{1,1})

%% Variables
fs_eeg=2000;
fs_2p=30;


%% Loader EEG from EDF file
if EEG_On==1
    [EEG_all]=ParseEDF('export.edf','annotations.txt',fs_eeg,plotON);
    EEG=EEG_all.part{2};
    EEG_ts=[1/fs_eeg:1/fs_eeg:length(EEG)/fs_eeg];
else 
   EEG=nan;
   EEG_ts=nan;
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


%% Split module
% Modifies the loaded varaibles to be a specific segment of the recording
if split_On==1
    if split_sel==1
        DC=DC(1:ceil(split_time*fs_dc));
        EEG=EEG(1:ceil(split_time*fs_eeg));
        F_G=F_G(:,1:ceil(split_time*fs_2p));
        F_R=F_R(:,1:ceil(split_time*fs_2p));
        Fneu_G=Fneu_G(:,1:ceil(split_time*fs_2p));
        Fneu_R=Fneu_R(:,1:ceil(split_time*fs_2p));
    elseif split_sel==2
        DC=DC(ceil(split_time*fs_dc):end);
        EEG=EEG(ceil(split_time*fs_eeg):end);
        F_G=F_G(:,ceil(split_time*fs_2p):end);
        F_R=F_R(:,ceil(split_time*fs_2p):end);
        Fneu_G=Fneu_G(:,ceil(split_time*fs_2p):end);
        Fneu_R=Fneu_R(:,ceil(split_time*fs_2p):end);
    end
    DC_ts=[1/fs_dc:1/fs_dc:length(DC)/fs_dc];
    EEG_ts=[1/fs_eeg:1/fs_eeg:length(EEG)/fs_eeg];
end


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
    if iscell(ii,1)==1
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


%% Filter (low pass): on dffs calculated from raw unfiltered data
lpfilt=1;
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

%Mean Trace Filter
lpfilt=.5;

%green
FGdff_MeanSm=lofi(Fc1Gdff_Mean,10^6/fs_2p,lpfilt,'verbose',0); %smoothed

%red
FRdff_MeanSm=lofi(-Fc1Rdff_Mean,10^6/fs_2p,lpfilt,'verbose',0); %smoothed

clear lpfilt


%% Correction of EEG bias
if EEGbias>0
    EEGa=[EEG(fs_eeg*EEGbias+1:end),zeros(1,fs_eeg*EEGbias)];
    EEG=EEGa;
elseif EEGbias<0
    EEGa=[zeros(1,fs_eeg*EEGbias),EEG(1:end-fs_eeg*EEGbias)];
    EEG=EEGa;
end
EEGb=circshift(EEG,-fs_eeg*EEGbias);

clear EEGa EEGb


%% PSD initial processing for feature detection
%used for finding pre-ictal spikes and also death following seizure death
if EEG_On==1
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

PSD_totP=sum(PSD_P,1);
PSD_theta=sum(PSD_P(and(PSD_F>3,PSD_F<15),:),1);
PSD_gammaL=sum(PSD_P(and(PSD_F>20,PSD_F<55),:),1);
PSD_100plus=sum(PSD_P(PSD_F>=100 & PSD_F<250 ,:),1);
PSD_100less=sum(PSD_P(PSD_F<100,:),1);
Ratio100=PSD_100plus./PSD_100less;
RA_Ratio100=SternRollAvg(Ratio100,20/PSDwin); % +/- 10 second window roll avg (time step here is 0.5s)
else
    PSD_F=nan;
    PSD_P=nan;
end

%% Deriving recording Features (Death) %need to tune to WOD seizures

if EEG_On==1 && TSW_On==1 && max(RA_Ratio100) > 1
    death=1;
    deathTime=find(Ratio100>1.25,1)*PSDwin/2;
    deathTime_eeg=deathTime*fs_eeg;
    deathTime_2p=deathTime*fs_2p;
elseif strcmp(szID,'S56ptz1')
    death=1;
    deathTime=93;%manually seeded due to lack of EEG recording
    deathTime_2p=deathTime*fs_2p;
else
    death=0;
    [deathTime,deathTime_2p,deathTime_eeg]=deal(NaN); 
end


%% Mean Trace Seizure and CSD Time Seeds (max first derivative block method)
if Sz_On==1
% determine all the points that are greater than the average of the min and
% max values of the trace

if death==1 %to account for an animal dying the data 3 min beyond death is excluded (3 min is to account for hypoxic wave to propogate uncertatinty in exact time of death)
    postDeathT=deathTime_2p+fs_2p*60*3;
    temp1G=FGdff_MeanSm(1:postDeathT)>=mean([min(FGdff_MeanSm(1:postDeathT)),max(FGdff_MeanSm(1:postDeathT))]);
    temp1R=FRdff_MeanSm(1:postDeathT)>=mean([min(FRdff_MeanSm(1:postDeathT)),max(FRdff_MeanSm(1:postDeathT))]);
else
    temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),max(FGdff_MeanSm)]);
    temp1R=FRdff_MeanSm>=mean([min(FRdff_MeanSm),max(FRdff_MeanSm)]);   
end

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

%determine the block bounds of the seizure and CSD blocks
if TSW_On==1
    [temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)
    MeanSzRecTimeG = min(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p;%defines the Sz block as the first block of the top 2 blocks; get the first bound of that block which is the max recruitment time (add 1 due to use of diff fn)
    MeanCSDRecTimeG = max(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
    
    [temp5R{1},temp5R{2}] = sort(temp4R{1},2,'descend');
    %determine the block bounds of the seizure and CSD blocks
    MeanSzRecTimeR = min(temp3R{1}(temp5R{2}(1:2))+1)/fs_2p;
    MeanCSDRecTimeR = max(temp3R{1}(temp5R{2}(1:2))+1)/fs_2p;

else
    [temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)
    %determine the block bounds of the seizure blocks
    MeanSzRecTimeG = (temp3G{1}(temp5G{2}(1))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
    MeanCSDRecTimeG = NaN;
    
    [temp5R{1},temp5R{2}] = sort(temp4R{1},2,'descend');
    %determine the block bounds of the seizure blocks
    MeanSzRecTimeR = (temp3R{1}(temp5R{2}(1))+1)/fs_2p;
    MeanCSDRecTimeR = NaN;
end

newMaxG = max(FGdff_MeanSm(round((MeanSzRecTimeG-2)*fs_2p):round((MeanSzRecTimeG+5)*fs_2p)));
newMaxR = max(FRdff_MeanSm(round((MeanSzRecTimeR-2)*fs_2p):round((MeanSzRecTimeR+5)*fs_2p)));

if strcmp(szID,'S56ptz1')
newMaxR=max(FRdff_MeanSm);
end

% Iteration 2

if death==1 %to account for an animal dying the data 3 min beyond death is excluded (3 min is to account for hypoxic wave to propogate uncertatinty in exact time of death)
    postDeathT=round(deathTime_2p+fs_2p*60*3);
    temp1G=FGdff_MeanSm(1:postDeathT)>=mean([min(FGdff_MeanSm(1:postDeathT)),newMaxG]);
    temp1R=FRdff_MeanSm(1:postDeathT)>=mean([min(FRdff_MeanSm(1:postDeathT)),newMaxR]);
else
    temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),newMaxG]);
    temp1R=FRdff_MeanSm>=mean([min(FRdff_MeanSm),newMaxR]);   
end

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

%determine the block bounds of the seizure blocks
if TSW_On==1
    [temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)
    MeanSzRecTimeG = min(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p;%defines the Sz block as the first block of the top 2 blocks; get the first bound of that block which is the max recruitment time (add 1 due to use of diff fn)
       
    [temp5R{1},temp5R{2}] = sort(temp4R{1},2,'descend');
    MeanSzRecTimeR = min(temp3R{1}(temp5R{2}(1:2))+1)/fs_2p;

else
    [temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)
    MeanSzRecTimeG = (temp3G{1}(temp5G{2}(1))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
    
    [temp5R{1},temp5R{2}] = sort(temp4R{1},2,'descend');
    MeanSzRecTimeR = (temp3R{1}(temp5R{2}(1))+1)/fs_2p;
end

else
MeanSzRecTimeG=NaN;
MeanSzRecTimeR=NaN;
MeanCSDRecTimeG=NaN;
MeanCSDRecTimeR=NaN;
end


%% End of Seizure
%determines seizure end time using normalized total power below 100Hz
if Sz_On==1 && death==0 && EEG_On==1
    SzEndTime=find(PSD_100less(floor(MeanSzRecTimeG/PSDwin*2):end)./max(PSD_100less)<0.05,1)*PSDwin/2+floor(MeanSzRecTimeG);
else
    SzEndTime=deal(NaN);
end

%% Preictal Mean RecTime Generation
% use the green channel times for calculations

if PIS_On==1
PreSpkRecTimeG=MeanRecTimes3pre(FGdff_MeanSm,fs_2p,15);
PreSpkRecTimeR=MeanRecTimes3pre(FRdff_MeanSm,fs_2p,15);

if strcmp(szID,'S102ptz3')%need to modify for split recording issues
PreSpkRecTimeG=MeanRecTimes3pre(FGdff_MeanSm,fs_2p,35);
PreSpkRecTimeR=MeanRecTimes3pre(FRdff_MeanSm,fs_2p,35);
end

%calculate an integral of the time 1 second following recruitment to filter
%out false positives
for ii=1:length(PreSpkRecTimeG{1,1})
    tempTrace=FGdff_MeanSm(round(PreSpkRecTimeG{1,1}(ii)*fs_2p):(round(PreSpkRecTimeG{1,1}(ii)+1)*fs_2p));
    valPIS1sG(ii)=trapz(tempTrace-min(tempTrace));
    clear tempTrace
end
clear ii

for ii=1:length(PreSpkRecTimeR{1,1})
    tempTrace=FRdff_MeanSm(round(PreSpkRecTimeR{1,1}(ii)*fs_2p):round((PreSpkRecTimeR{1,1}(ii)+1)*fs_2p));
    valPIS1sR(ii)=trapz(tempTrace-min(tempTrace));
    clear tempTrace
end
clear ii

if isnan(MeanSzRecTimeG)
    CaisLrgG=valPIS1sG>std(FGdff_MeanSm)*fs_2p/3/0.3989; %gaussian, height=std of preictal period, sigma=1/3 second in frames (1s width at 1.5 sigma)
    CaisLrgR=valPIS1sR>std(FRdff_MeanSm)*fs_2p/3/0.3989;
else
    CaisLrgG=valPIS1sG>std(FGdff_MeanSm(1:round(MeanSzRecTimeG*fs_2p)))*fs_2p/3/0.3989; %gaussian, height=std of preictal period, sigma=1/3 second in frames (1s width at 1.5 sigma)
    CaisLrgR=valPIS1sR>std(FRdff_MeanSm(1:round(MeanSzRecTimeG*fs_2p)))*fs_2p/3/0.3989;
end

for ii=1:4
     PreSpkRecTimeG{1,ii}=PreSpkRecTimeG{1,ii}(CaisLrgG);
     PreSpkRecTimeR{1,ii}=PreSpkRecTimeR{1,ii}(CaisLrgR);
end

S_th = 4; %power threshold based on median for PIS detection
winW = .7; %seconds for spike detection window width

CaDetectWin=.25;

EEG_spk_prelim = BO_IIS_Detect_Power_Rolston_v5(EEG',S_th,fs_eeg,winW,0); %use v4 to not use theta gamma power ratio component
EEG_spkTimes_prelim_old=EEG_spk_prelim.time_Sp_final{1};%grab eeg spikes indexed by largest change in feature point
EEG_spkTimes_prelim=EEG_spk_prelim.time_Sp_minPeak{1};%grab all possible EEG spike times indexed by the min peak (spike) time point

%find Ca spikes in mean traces that coorepsond to EEG spikes
DiffPreEEG_G=abs(PreSpkRecTimeG{1,1}-(EEG_spkTimes_prelim'));%generates a Ca spk by EEG spk matrix of time diff
DiffPreEEG_R=abs(PreSpkRecTimeR{1,1}-(EEG_spkTimes_prelim'));
if isnan(MeanSzRecTimeG)
    MeanPreIsRecG=min(DiffPreEEG_G,[],1)<CaDetectWin;%logical indexing across Ca spk min to indicate if recruited
    MeanPreIsRecR=min(DiffPreEEG_R,[],1)<CaDetectWin;
else
    MeanPreIsRecG=min(DiffPreEEG_G,[],1)<CaDetectWin & PreSpkRecTimeG{1,1}<MeanSzRecTimeG;%logical indexing across Ca spk min to indicate if recruited
    MeanPreIsRecR=min(DiffPreEEG_R,[],1)<CaDetectWin & PreSpkRecTimeR{1,1}<MeanSzRecTimeG;
end
MeanCaG_PIS_times=PreSpkRecTimeG{1,1}(MeanPreIsRecG);
MeanCaR_PIS_times=PreSpkRecTimeR{1,1}(MeanPreIsRecR);

%refine EEG spike times based on Ca recruitment and the seizure recruitment
%time (smaller of the red or green channel)...may need to make it the
%larger of the two instead, REVISIT
if isnan(MeanSzRecTimeG)
    EEGisTrue=min(DiffPreEEG_G,[],2)<CaDetectWin;%final EEG spikes times i.e. correspond to at least one channel of Ca recruitment
else
    EEGisTrue=and(min(DiffPreEEG_G,[],2)<CaDetectWin, EEG_spkTimes_prelim'<MeanSzRecTimeG);%final EEG spikes times i.e. correspond to at least one channel of Ca recruitment
end
EEG_PIS_times=EEG_spkTimes_prelim(EEGisTrue);%pre-ictal spike times

if strcmp(szID,'S102ptz3')%need to manually seed times for s103 recording 3 because of split recording issues)
EEG_PIS_times=[6.1575,10.8495,15.9085,22,0405,28.757,36.421,45.116,52.36];
end

%% EEG spike wave kernel and PIS refinement
EEGswd_KernSet=nan(numel(EEG_PIS_times),1*fs_eeg);
for ii=1:numel(EEG_PIS_times)
    EEGtempI=round((EEG_PIS_times(ii)-.5)*fs_eeg);
    EEGswd_KernSet(ii,:)=EEG(EEGtempI+1:EEGtempI+fs_eeg);
end
clear EEGtempI
EEGswd_kern=mean(EEGswd_KernSet,1);

% Convolve template with eeg
EEGswdConv=conv(flip(EEGswd_kern), EEG);
if isnan(MeanSzRecTimeG)
    [~,EEGswd_spikeTimes]=findpeaks(EEGswdConv,'MinPeakHeight',5*std(EEGswdConv),'MinPeakDistance',fs_eeg);
else
    [~,EEGswd_spikeTimes]=findpeaks(EEGswdConv(1:round(MeanSzRecTimeG*fs_eeg)),'MinPeakHeight',5*std(EEGswdConv(1:round(MeanSzRecTimeG*fs_eeg))),'MinPeakDistance',fs_eeg);
end
EEGswd_spikeTimes=EEGswd_spikeTimes/fs_eeg-0.5;

%combine the new detected spikes into the calcium contstrained spikes
EEGswdPISdiff=EEG_PIS_times-EEGswd_spikeTimes';
EEGswdPISmin=min(abs(EEGswdPISdiff),[],2);
EEG_PIS_new=EEGswd_spikeTimes(EEGswdPISmin>1);
EEG_PIS_times2=sort([EEG_PIS_times,EEG_PIS_new]);
if ~isnan(MeanSzRecTimeG)
    EEG_PIS_times2=EEG_PIS_times2(EEG_PIS_times2<MeanSzRecTimeG); %make sure to grab sentinal spike but not seizure
end

if plotON==1
figure
for ii=1:size(EEGswd_KernSet,1)
    subplot(ceil(sqrt(size(EEGswd_KernSet,1))),ceil(sqrt(size(EEGswd_KernSet,1))),ii)
    plot(EEGswd_KernSet(ii,:))
    title(ii)
end
suptitle('Ca EEG Concordant Spks')
clear ii

figure
plot(EEGswd_kern)
suptitle('EEG PIS Kernel')
end

%% interspike interval and first spike determination
%looks for a rolling average of 10 to designate consistent spiking but then
%looks for the first spike in this train that has an ISI within 15 seconds
%as an upper limit on the allowable ISI
ISI=diff(EEG_PIS_times2);
ISI_RAwin=5;
RAtolerance=10;
IndivTolerance=15;
if numel(ISI)<ISI_RAwin
    firstSpkI=find(ISI<IndivTolerance,1);%if less than 5 spikes allow a slightly larger tolerance for ISI to grab spieks to maximize spikes collected
else
    ISI_RA=SternRollAvg(ISI,ISI_RAwin);
    
    ISI_RAthresh=find(ISI_RA<RAtolerance,1);%finds the spike at the beginning of a train where the roll avg is <10
    if ISI_RAthresh>ISI_RAwin
        firstSpkI=find(ISI(ISI_RAthresh-ISI_RAwin:ISI_RAthresh)<IndivTolerance,1);%finds the first spike with an ISI <10 in the first window where the roll avg is <10.
        firstSpkI=firstSpkI+ISI_RAthresh-ISI_RAwin-1;%neg 1 is needed for fence post (index 1 needs to coorespond to first spk in range considered)
    else
        firstSpkI=find(ISI(1:ISI_RAthresh+ISI_RAwin)<IndivTolerance,1);
    end
end

firstSpkT=EEG_PIS_times2(firstSpkI);

EEG_PIS_times2=EEG_PIS_times2(firstSpkI:end);

end


%% Filtering for Individual Trace Recruitment Time Determination
% filtering to smooth for recruitment detection
lpfilt=.1; %**************USING B TRACE***************

Fb1Gdff_flt2=zeros(size(Fb1Gdff));
Fc1Gdff_flt2=zeros(size(Fc1Gdff));
Fb1Rdff_flt2=zeros(size(Fb1Rdff));
Fc1Rdff_flt2=zeros(size(Fc1Rdff));

%green
for ii=1:size(Fc1Gdff_flt1,1)
    Fb1Gdff_flt2(ii,:)=lofi(Fb1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Gdff_flt2(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,1,'verbose',0);
end

%red
for jj=1:size(Fc1Rdff_flt1,1)
    Fb1Rdff_flt2(jj,:)=lofi(Fb1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Rdff_flt2(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,1,'verbose',0);
end

clear lpfilt
clear ii jj

Fb1Gdff_flt2_Mean=mean(Fb1Gdff_flt2,1);
Fb1Rdff_flt2_Mean=mean(Fb1Rdff_flt2,1);

%% Individual seizure recruitment detection (Version 3)
if Sz_On==1

if strcmp(szID,'S86ptz6')
    %MeanSzRecTimeG=145;
elseif strcmp(szID,'S87ptz8')
    MeanSzRecTimeG=406;
elseif strcmp(szID,'S104ptz3')
    %MeanSzRecTimeG=466;
end    

MeanSzRecTimeR=MeanSzRecTimeG;

RecTimeG_sz3_struct=IndivRecTimes3(Fb1Gdff_flt2,MeanSzRecTimeG,fs_2p,1); %selecting a 1s sigma for gaussian to prevent sentinal spike from being grabbed and also if the average speed of the seizure is 200um/s then it should not take more than 1.5 seconds to cross the field so a 1s sigma gives us a 3 second window at at least 33% (1.5sigma) amplification of feature occurs
RecTimeR_sz3_struct=IndivRecTimes3(-Fb1Rdff_flt2,MeanSzRecTimeR,fs_2p,1);
RecTimeG_sz2=RecTimeG_sz3_struct.time;
RecTimeR_sz2=RecTimeR_sz3_struct.time;

% Recruited cell or not
%trace needs to have atleast a 20% increase in signal during seizure and the
%seizure detected time needs to be determined to be within 3 seconds of
%the seed time (if average speed is 100um/s and field is <300um all values
%should be within 3 seconds of mean giving 6 second window)

tprd=10;%time period (s)
gthresh=0.2; %green threshold typical G threshold .3
rthresh=0.2; %red threshold typical R threshold .3
preSzMeanDffG_b=mean(Fb1Gdff_flt2(:,floor((MeanSzRecTimeG-tprd)*fs_2p):floor(MeanSzRecTimeG*fs_2p)),2);
postSzMeanDffG_b=mean(Fb1Gdff_flt2(:,floor(MeanSzRecTimeG*fs_2p):floor((MeanSzRecTimeG+tprd)*fs_2p)),2);
SzMeanDffDiffG_b=postSzMeanDffG_b-preSzMeanDffG_b;
isRecruitedG=and(SzMeanDffDiffG_b>gthresh,abs(RecTimeG_sz2-MeanSzRecTimeG)<3); 

preSzMeanDffG=mean(Fc1Gdff_flt2(:,floor((MeanSzRecTimeG-tprd)*fs_2p):floor(MeanSzRecTimeG*fs_2p)),2);
postSzMeanDffG=mean(Fc1Gdff_flt2(:,floor(MeanSzRecTimeG*fs_2p):floor((MeanSzRecTimeG+tprd)*fs_2p)),2);
SzMeanDffDiffG=postSzMeanDffG-preSzMeanDffG;
%disp(sum(isRecruitedG))

preSzMeanDffR_b=mean(-Fb1Rdff_flt2(:,floor((MeanSzRecTimeR-tprd)*fs_2p):floor(MeanSzRecTimeR*fs_2p)),2);
postSzMeanDffR_b=mean(-Fb1Rdff_flt2(:,floor(MeanSzRecTimeR*fs_2p):floor((MeanSzRecTimeR+tprd)*fs_2p)),2);
SzMeanDffDiffR_b=postSzMeanDffR_b-preSzMeanDffR_b;
isRecruitedR=and(SzMeanDffDiffR_b>rthresh,abs(RecTimeR_sz2-MeanSzRecTimeR)<3);

preSzMeanDffR=mean(-Fc1Rdff_flt2(:,floor((MeanSzRecTimeR-tprd)*fs_2p):floor(MeanSzRecTimeR*fs_2p)),2);
postSzMeanDffR=mean(-Fc1Rdff_flt2(:,floor(MeanSzRecTimeR*fs_2p):floor((MeanSzRecTimeR+tprd)*fs_2p)),2);
SzMeanDffDiffR=postSzMeanDffR-preSzMeanDffR;
%disp(sum(isRecruitedR))
clear tprd gthresh rthresh

else
[RecTimeG_sz2, RecTimeR_sz2, isRecruitedG, isRecruitedR, SzMeanDffDiffG, SzMeanDffDiffR] = deal(NaN);
end

%% Individual CSD / Death recruitment detection
if TSW_On==1

if strcmp(szID,'S86ptz6')
    MeanCSDRecTimeG=152.5;%had to manual seed due to low signal
elseif strcmp(szID,'S87ptz8')
    MeanCSDRecTimeG=538;%had to manual seed due to motion artifact
elseif strcmp(szID,'S104ptz3')
    %MeanCSDRecTimeG=462;
elseif strcmp(szID,'S56ptz1')
    MeanCSDRecTimeG=96;%had to manual seed due to motion artifact
end


% Individual Trace Recruitment Time Determination csd
% filtering to smooth for recruitment detection


%Temporarily need to use the Yellow CSD time seed the Red until red is fixed
MeanCSDRecTimeR=MeanCSDRecTimeG;

%NOTE: Invert the Red Channel for CSD detection

RecTimeG_csd_struct=IndivRecTimes3(Fb1Gdff_flt2,MeanCSDRecTimeG,fs_2p,5);%average speed of 40um/s gives means it takes 7s to cross 280um feidl of view, therefore a sigma of 5s gives us a 10s window over which 
RecTimeR_csd_struct=IndivRecTimes3(-Fb1Rdff_flt2,MeanCSDRecTimeR,fs_2p,FOV_dim/45);
RecTimeG_csd=RecTimeG_csd_struct.time;
RecTimeR_csd=RecTimeR_csd_struct.time;


% Recruited cell or not (CSD)
%trace needs to have atleast a 20% increase in signal during CSD relative to before and the
%CSD detected time needs to be determined to be within 10 seconds of
%the seed
%NOTE:Invert the Red channel signal for CSD

tprd=5;%time period (s) 5s before in case sz was close before
if strcmp(szID,'S86ptz6') || strcmp(szID,'S56ptz1')
    gthresh=0.05; 
    rthresh=0.1;
else
    gthresh=0.1; 
    rthresh=0.02;
end

CSDMeanDffDiffG=nan(size(RecTimeG_csd));
CSDMeanDffDiffG_b=nan(size(RecTimeG_csd));
for ii=1:length(RecTimeG_csd)
    preCSDrange=floor([(RecTimeG_csd(ii)-tprd)*fs_2p:RecTimeG_csd(ii)*fs_2p]);
    postCSDrange=floor([RecTimeG_csd(ii)*fs_2p:(RecTimeG_csd(ii)+tprd)*fs_2p]);
    preCSDMeanDff_b=mean(Fb1Gdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff_b=mean(Fb1Gdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffG_b(ii)=postCSDMeanDff_b-preCSDMeanDff_b;
    
    preCSDMeanDff=mean(Fc1Gdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff=mean(Fc1Gdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffG(ii)=postCSDMeanDff-preCSDMeanDff;
end
isRecruitedG_csd=and(CSDMeanDffDiffG_b>gthresh,abs(RecTimeG_csd-MeanCSDRecTimeG)<FOV_dim/45); 
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff

CSDMeanDffDiffR=nan(size(RecTimeR_csd));
CSDMeanDffDiffR_b=nan(size(RecTimeG_csd));
for ii=1:length(RecTimeR_csd)
    preCSDrange=floor([(RecTimeR_csd(ii)-tprd)*fs_2p:RecTimeR_csd(ii)*fs_2p]);
    postCSDrange=floor([RecTimeR_csd(ii)*fs_2p:(RecTimeR_csd(ii)+tprd)*fs_2p]);
    preCSDMeanDff_b=mean(-Fb1Rdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff_b=mean(-Fb1Rdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffR_b(ii)=postCSDMeanDff_b-preCSDMeanDff_b;
    
    preCSDMeanDff=mean(-Fc1Rdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff=mean(-Fc1Rdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffR(ii)=postCSDMeanDff-preCSDMeanDff;
end
isRecruitedR_csd=and(CSDMeanDffDiffR_b>rthresh,abs(RecTimeR_csd-MeanCSDRecTimeR)<FOV_dim/45);
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff
clear tprd gthresh rthresh

% isRecruitedG_csd = ~isnan(isRecruitedG_csd);
% isRecruitedR_csd = ~isnan(isRecruitedR_csd);


else
[RecTimeG_csd, RecTimeR_csd, isRecruitedG_csd, isRecruitedR_csd, CSDMeanDffDiffG, CSDMeanDffDiffR] = deal(NaN);
end


%% DC Analysis

if DC_On==1
tempRange=round(deathTime*fs_dc):round((MeanCSDRecTimeG+10)*fs_dc);
tempTrace=SternRollAvg(DC(tempRange),2*fs_dc);
tempTrace_ds=downsample(tempTrace,fs_dc);
tempTrace_sp=spline(1:1:length(tempTrace_ds),tempTrace_ds,1:1/10:length(tempTrace_ds));
[~,DCstartCSDmaxI]=min(diff(tempTrace_sp(1:end-100),2,2));%finds the elbow of the dc trace at start of csd
CSDstartTimeDC=deathTime+(DCstartCSDmaxI-1)/fs_dc*200;

DCShiftPre=mean(DC(round(deathTime*fs_dc):round(CSDstartTimeDC*fs_dc)));
DCShiftPost=min(SternRollAvg(DC(round(CSDstartTimeDC*fs_dc):round((CSDstartTimeDC+60)*fs_dc)),fs_dc));
CSD_dcShift=[DCShiftPost-DCShiftPre,DCShiftPost];%the difference and the raw average shift

else
%manually seed times if missing DC trace
CSDstartTimeDC=nan;
CSD_dcShift=nan;
end


%% Post Ictal calcium level LEGACY NOT USED FOR WOD
%Determines the maximum calcium change during a moving window post seizure
%for 5 minutes as a comparator for the post ictal CSD period
if Sz_On==1 && ~strcmp(szID,'S102ptz2') && death==0
PostIctalWinStart=[SzEndTime:SzEndTime+300]; %selected as the expected post ictal period we would expect to observed a csd in
PostIctalWin=20; %selected to be 1/3 length of the average csd depressed period
PostIctalCaG=nan(size(Fc1Gdff_flt2,1),length(PostIctalWinStart));
PostIctalCaR=nan(size(Fc1Rdff_flt2,1),length(PostIctalWinStart));
for jj=1:length(PostIctalWinStart)
    PostIctalCaGbase=mean(Fc1Gdff_flt2(:,(SzEndTime:SzEndTime+5)*fs_2p),2);%5 second baseline period selected post seizure to reduce change of including csd
    PostIctalCaRbase=mean(Fc1Rdff_flt2(:,(SzEndTime:SzEndTime+5)*fs_2p),2);
    PostIctalCaG(:,jj)=mean(Fc1Gdff_flt2(:,(PostIctalWinStart(jj):PostIctalWinStart(jj)+PostIctalWin)*fs_2p),2)-PostIctalCaGbase;
    PostIctalCaR(:,jj)=mean(Fc1Rdff_flt2(:,(PostIctalWinStart(jj):PostIctalWinStart(jj)+PostIctalWin)*fs_2p),2)-PostIctalCaRbase;
end
PostIctalCaGmax=max(PostIctalCaG,[],2); %max increase in G calcium activity post seizure
PostIctalCaRmax=max(-PostIctalCaR,[],2); %max decrease in R calcium activity post seizure
else
[PostIctalCaGmax,PostIctalCaRmax]=deal(nan);
end


%% Individual Cell PIS Recruitment Time Detection
%taking all PIS including sentintal spike, and using f'-steepest slope for
%time indexing
if PIS_On==1
    if isnan(MeanSzRecTimeG)
        RecTimeG_PIS=IndivRecTimes3pre(Fc1Gdff_flt2,EEG_PIS_times2,fs_2p,20,1,4);%fed the EEG PIS times
        RecTimeR_PIS=IndivRecTimes3pre(-Fc1Rdff_flt2,EEG_PIS_times2,fs_2p,20,1,4);
    else
        RecTimeG_PIS=IndivRecTimes3pre(Fc1Gdff_flt2(:,1:floor(MeanSzRecTimeG*fs_2p)),EEG_PIS_times2,fs_2p,20,1,4);%fed the EEG PIS times
        RecTimeR_PIS=IndivRecTimes3pre(-Fc1Rdff_flt2(:,1:floor(MeanSzRecTimeG*fs_2p)),EEG_PIS_times2,fs_2p,20,1,4);
    end


%% Refine Spikes with min 10% recruitment of cells in green channel
   %And remove spikes where the mean is clearly before or after the EEG
   %event window of 1 seconds (window used for spike seperation)
PISisRecG=nansum(~isnan(RecTimeG_PIS.v3),1)/size(RecTimeG_PIS.v3,1)>.1;
PISisRecG2=abs(nanmean(RecTimeG_PIS.v3,1)-EEG_PIS_times2)<1;
PISisRec=(PISisRecG+PISisRecG2)==2;
RecTimeG_PIS.true=RecTimeG_PIS.v3(:,PISisRec);
RecTimeR_PIS.true=RecTimeR_PIS.v3(:,PISisRec);
EEG_PIS_times3=EEG_PIS_times2(PISisRec);


%% Determine Ca change during spike
[PISCaChangeG,PISCaChangeR]=deal(nan(size(Fc1Gdff_flt2,1),length(EEG_PIS_times2)));
for ii = 1:size(Fc1Gdff_flt2,1)
    tempTraceG=lofi(Fc1Gdff(ii,:),10^6/fs_2p,.5,'verbose',0);
    tempTraceR=lofi(-Fc1Rdff(ii,:),10^6/fs_2p,.5,'verbose',0);
    for jj = 1:length(EEG_PIS_times2)
        spkIval=ceil(EEG_PIS_times2(jj)*fs_2p);
        %green
        beforePISmeanG=max(tempTraceG(spkIval-2*fs_2p:spkIval-fs_2p));
        PISmaxCaG=max(tempTraceG(spkIval-fs_2p:spkIval+fs_2p));
        PISCaChangeG(ii,jj)=PISmaxCaG-beforePISmeanG;
        %red
        beforePISmeanR=max(tempTraceR(spkIval-2*fs_2p:spkIval-fs_2p));
        PISmaxCaR=max(tempTraceR(spkIval-fs_2p:spkIval+fs_2p));
        PISCaChangeR(ii,jj)=PISmaxCaR-beforePISmeanR;
    end
end
AvgPISCaChangeG=mean(PISCaChangeG,2);
AvgPISCaChangeR=mean(PISCaChangeR,2);


else
    [MeanCaG_PIS_times,MeanCaR_PIS_times,PISisRec,firstSpkT,firstSpkI,RecTimeG_PIS,RecTimeR_PIS,EEG_PIS_times2,EEG_PIS_times3,AvgPISCaChangeG,AvgPISCaChangeR]=deal(nan);
end


%% RecTime PLOTS
if plotON==1
tempG=Fc1Gdff_flt2; tempR=Fc1Rdff_flt2;

%Mean PIS, Sz and CSD traces and times 
figure
suptitle('Means')
plot(EEG_ts,EEG./4000+1,'color',[0.7 0.7 0.7]);
hold on
plot(F1_ts,Fc1Gdff_flt1_Mean,'color', [0 1 0]);
plot(F1_ts,Fc1Rdff_flt1_Mean-.5,'color',[1 0 1]);
plot(DC_ts,DC./20+2,'color',[0 0 0]);
%plot(F1_ts,mean(tempG),'color', [0 1 0]);
%plot(F1_ts,mean(tempR),'color',[1 0 1]);
line([MeanSzRecTimeG;MeanSzRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','--')
line([MeanSzRecTimeR;MeanSzRecTimeR],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','--')
line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-.')
line([MeanCSDRecTimeR;MeanCSDRecTimeR],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','-.')
%line([EEG_PIS_times3;EEG_PIS_times3],repmat([-4;5],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
if death==1
    line([deathTime;deathTime],[-4;5],'color',[.5 0 0],'linewidth',2)
end


%% Individual cell sz recruitment times
if Sz_On==1
figure
suptitle('Sz')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    if isRecruitedG(ii)==1 %incldG(ii)==1
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    else
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 1]);
    end
    line([RecTimeG_sz2(ii);RecTimeG_sz2(ii)],repmat([ii*10;ii*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanSzRecTimeG;MeanSzRecTimeG],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end

subplot(1,2,2)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for jj=1:size(Fc1Rdff_flt2,1)
    if isRecruitedR(jj)==1
        plot(F1_ts,tempR(jj,:)+jj*6,'color',[1 0.25 1]);
    else
        plot(F1_ts,tempR(jj,:)+jj*6,'color',[1 0.75 1]);
    end
    line([RecTimeR_sz2(jj);RecTimeR_sz2(jj)],repmat([jj*6;jj*6+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanSzRecTimeR;MeanSzRecTimeR],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj
end

%% Plot CSD indiv rec times
if TSW_On==1
% Individual cell csd/death recruitment times
figure
suptitle('CSD')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    if isRecruitedG_csd(ii)==1 %incldG(ii)==1
        plot(F1_ts,tempG(ii,:)*10+ii*10,'color', [0 1 0]);
    else
        plot(F1_ts,tempG(ii,:)*10+ii*10,'color', [0 1 1]);
    end
    line([RecTimeG_csd(ii);RecTimeG_csd(ii)],repmat([ii*10;ii*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end

subplot(1,2,2)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for jj=1:size(Fc1Rdff_flt2,1)
    if isRecruitedR_csd(jj)==1
        plot(F1_ts,tempR(jj,:)*10+jj*10,'color',[1 0.25 1]);
    else
        plot(F1_ts,tempR(jj,:)*10+jj*10,'color',[1 0.75 1]);
    end
    line([RecTimeR_csd(jj);RecTimeR_csd(jj)],repmat([jj*10;jj*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanCSDRecTimeR;MeanCSDRecTimeR],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj
end

%% pre ictal spikes
if PIS_On==1
figure
suptitle('PIS')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    line([RecTimeG_PIS.v1(ii,:);RecTimeG_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeG_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v3(ii,:);RecTimeG_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)

end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)

subplot(1,2,2)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Rdff_flt2,1)
    plot(F1_ts,tempR(ii,:)+ii*10,'color', [1 0 1]);
    line([RecTimeR_PIS.v1(ii,:);RecTimeR_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeR_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeR_PIS.v3(ii,:);RecTimeR_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeR_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)

end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)

clear ii tempG tempR
end
end

%% Generate Data Structure
%General Variables
StructFieldNames={'szID'
    'FOV_dim'
    %Time/Sampling
    'death'
    'deathTime'
    'firstSpkT'
    'firstSpkI'
    'fs_2p'
    'fs_eeg'
    'EEGbias'
    'EEG_ts'
    'F1_ts'
    %EEG
    'EEG'
    'PSD_F'
    'PSD_P'
    %DC
    'DC'
    'fs_dc'
    'DC_ts'
    %F1 processing
    'F1G'
    'Fneu1G'
    'stat'
    'xCoor'
    'yCoor'
    'xCoorI'
    'yCoorI'
    'iscell'
    'iscellI'
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
    'EEG_PIS_times2'
    'EEG_PIS_times3'
    'MeanCaG_PIS_times'
    'MeanCaR_PIS_times'
    'PISisRec'
    'RecTimeG_PIS'
    'RecTimeR_PIS'
    'MeanSzRecTimeG'
    'MeanSzRecTimeR'
    'SzEndTime'
    'MeanCSDRecTimeG'
    'MeanCSDRecTimeR'
    'RecTimeG_sz2'
    'RecTimeR_sz2'
    'RecTimeG_csd'
    'RecTimeR_csd'
    'isRecruitedG'
    'isRecruitedR'
    'isRecruitedG_csd'
    'isRecruitedR_csd'
    'AvgPISCaChangeG'
    'AvgPISCaChangeR'
    'SzMeanDffDiffG'
    'SzMeanDffDiffR'
    'CSDMeanDffDiffG'
    'CSDMeanDffDiffR'
    'PostIctalCaGmax'
    'PostIctalCaRmax'
    'CSDstartTimeDC'
    'CSD_dcShift'};    
for kk=1:length(StructFieldNames)
    eval(strcat('szStruct.',StructFieldNames{kk},'=',StructFieldNames{kk},';'));
end

%% SAVE Structure Varaible as the SeizureID (szID)
assignin('base', szID, szStruct);
save(strcat(szID,'.mat'),szID,'-v7.3')
disp('saved')
%cannot save file as normal .mat file, needs the compression of v7.3

cd ..
clearvars -except mm szList

end





