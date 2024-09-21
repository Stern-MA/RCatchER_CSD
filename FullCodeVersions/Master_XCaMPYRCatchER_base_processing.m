%MASTER SCRIPT FOR PROCESSING XCaMP-Y-RCatchER BASELINE DATA

%*****___Baseline DATA___*****

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
% 10) Filtering of Individual Traces to 1 Hz
%
%
%
%
%

%%
szID='S56galvo1';
plotON=0;
fs_2p=30;


if strcmp(szID,'S56galvo1')
fs_2p=1.07;%for galvo recording)
end

%% Loader EEG from EDF file
if strcmp(szID,'S102base1')
    fs_eeg=2000;
    [EEG_all]=ParseEDF('export.edf','annotations.txt',fs_eeg,plotON);
    EEG=EEG_all.part{1};
    EEG_ts=[1/fs_eeg:1/fs_eeg:length(EEG)/fs_eeg];
else
    fs_eeg=NaN;
    EEG=NaN;
    EEG_ts=NaN; 
end

%% Loader DC from .abf file
if strcmp(szID,'S102base1')
    
fileName = dir('*.abf');
[DC,DC_tinv,~]=abfload(fileName.name);%DC is the vector DC_tinv is the time interval in us
fs_dc=10^6/DC_tinv;
if and(~isnan(EEG),length(DC)>length(EEG))
    DC=DC(1:floor(length(EEG)/fs_eeg*fs_dc));
end
DC_ts=[1/fs_dc:1/fs_dc:length(DC)/fs_dc];

else
    fs_dc=NaN;
    DC=NaN;
    DC_ts=NaN;
end


%% Loader Motion Tracking 

if ~strcmp(szID,'S102base1')
lpfilt=1;
cd('MotionTracking')
fileName = dir('*.tdms');
my_tdms_struct1=TDMS_getStruct(fileName.name);
speedData1=my_tdms_struct1.Data.Speed_mm_per_s.data; %speed data
speedDataTSstr=my_tdms_struct1.Data.Relative_timestamp.data; %time stamp string for speed data
TSMatrix=datevec(speedDataTSstr); %convert time stamp string vector into numerical matrix, rows are each time stamp entry and columns are y, m, d, h, m, s  
speedDataTS1=nan(size(TSMatrix,1),1);
for ii=1:size(TSMatrix,1) %convert time stamp matrix into time stamp vector in seconds
    speedDataTS1(ii)=TSMatrix(ii,4)*3600+TSMatrix(ii,5)*60+TSMatrix(ii,6);
end
fs_mt=length(speedDataTS1)/speedDataTS1(end);
speedDataTS1=speedDataTS1';
speedData1_flt1=lofi(speedData1,10^6/fs_mt,lpfilt,'verbose',0);

clear ii
clear speedDataTSstr
clear TSMatrix
cd ..

else
    [speedData1,speedDataTS1]=deal(nan);
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
Fc1Gdff_std=std(Fc1Gdff);%standard deviation of the population mean of dF/F overall timesteps 
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
if ~strcmp(szID,'S56galvo1')
lpfilt=5;
lpfilt2=1;
%green
[Fb1Gdff_flt1,Fneub1Gdff_flt1,Fc1Gdff_flt1,Fc1Gdff_flt2]=deal(zeros(size(F1G)));
for ii=1:size(F1G,1)
    Fb1Gdff_flt1(ii,:)=lofi(Fb1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fneub1Gdff_flt1(ii,:)=lofi(Fneub1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Gdff_flt1(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Gdff_flt2(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt2,'verbose',0);
end

%red
[Fb1Rdff_flt1,Fneub1Rdff_flt1,Fc1Rdff_flt1,Fc1Rdff_flt2]=deal(zeros(size(F1R)));
for jj=1:size(F1R,1)
    Fb1Rdff_flt1(jj,:)=lofi(Fb1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fneub1Rdff_flt1(jj,:)=lofi(Fneub1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Rdff_flt1(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Rdff_flt2(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,lpfilt2,'verbose',0);
end

clear lpfiltG hpfilt
clear ii jj
else
    Fb1Gdff_flt1=Fb1Gdff;
    Fneub1Gdff_flt1=Fneub1Gdff;
    [Fc1Gdff_flt1, Fc1Gdff_flt2]=deal(Fc1Gdff);
    Fb1Rdff_flt1=Fb1Rdff;
    Fneub1Rdff_flt1=Fneub1Rdff;
    [Fc1Rdff_flt1, Fc1Rdff_flt2]=deal(Fc1Rdff);
end


%% Spike Detector
%detect individual calcium spikes and grabs there rec times
%grab the traces
%calculate the average peak heights
%calculate change in the red during the same period
win_post=3;%in seconds
win_pre=1;%in seconds

[spksI,RecTime_spks,spks_traceG,spks_traceR,spks_traceCaDiffG,spks_traceCaDiffR]=deal(cell(size(Fc1Gdff_flt2,1),1));
Fc1Gdff_flt2_df=diff(Fc1Gdff_flt2,1,2);
thresh=mean(Fc1Gdff_flt2,2)+4*std(Fc1Gdff_flt2,0,2);
[spks_trace_avgG,spks_trace_stdG,spks_trace_avgR,spks_trace_stdR]=deal(nan(size(Fc1Gdff_flt2,1),floor((win_pre+win_post)*fs_2p)+1));
for ii=1:size(Fc1Gdff_flt2,1)
    [~,spksI{ii}]=findpeaks(Fc1Gdff_flt2(ii,:),'MinPeakHeight',thresh(ii),'MinPeakDistance',1);
    TempTrace=nan(length(spksI{ii}),floor(win_pre*fs_2p)+1);
    for jj=1:length(spksI{ii})
        if spksI{ii}(jj)-win_pre*fs_2p<=0
            TempTrace(jj,abs(spksI{ii}(jj)-floor(win_pre*fs_2p))+2:end)=Fc1Gdff_flt2_df(ii,1:spksI{ii}(jj));
        else
            TempTrace(jj,:)=Fc1Gdff_flt2_df(ii,spksI{ii}(jj)-floor(win_pre*fs_2p):spksI{ii}(jj));
        end
    end
    [~,maxI]=max(TempTrace,[],2);
    RecTime_spks{ii}=maxI'+spksI{ii}-floor(win_pre*fs_2p);
    [spks_traceG{ii},spks_traceR{ii}]=deal(nan(length(spksI{ii}),floor((win_pre+win_post)*fs_2p)+1));
    [spks_traceCaDiffG{ii},spks_traceCaDiffR{ii}]=deal(nan(length(spksI{ii}),1));
    for jj=1:length(spksI{ii})
        if RecTime_spks{ii}(jj)-win_pre*fs_2p<=0
            spks_traceG{ii}(jj,abs(RecTime_spks{ii}(jj)-floor(win_pre*fs_2p))+2:end)=Fc1Gdff_flt2(ii,1:RecTime_spks{ii}(jj)+floor(win_post*fs_2p));
            spks_traceR{ii}(jj,abs(RecTime_spks{ii}(jj)-floor(win_pre*fs_2p))+2:end)=Fc1Rdff_flt2(ii,1:RecTime_spks{ii}(jj)+floor(win_post*fs_2p));
        elseif RecTime_spks{ii}(jj)+win_post*fs_2p>=size(Fc1Gdff_flt2,2)
            spks_traceG{ii}(jj,1:end+size(Fc1Gdff_flt2,2)-(RecTime_spks{ii}(jj)+floor(win_post*fs_2p)))=Fc1Gdff_flt2(ii,RecTime_spks{ii}(jj)-floor(win_pre*fs_2p):end);
            spks_traceR{ii}(jj,1:end+size(Fc1Rdff_flt2,2)-(RecTime_spks{ii}(jj)+floor(win_post*fs_2p)))=Fc1Rdff_flt2(ii,RecTime_spks{ii}(jj)-floor(win_pre*fs_2p):end);
        else
            spks_traceG{ii}(jj,:)=Fc1Gdff_flt2(ii,RecTime_spks{ii}(jj)-floor(win_pre*fs_2p):RecTime_spks{ii}(jj)+floor(win_post*fs_2p));
            spks_traceR{ii}(jj,:)=Fc1Rdff_flt2(ii,RecTime_spks{ii}(jj)-floor(win_pre*fs_2p):RecTime_spks{ii}(jj)+floor(win_post*fs_2p));
        end
        %Calcium Levels
        if RecTime_spks{ii}(jj)-.5*fs_2p<=0
            spks_traceCaDiffG{ii}(jj)=mean(Fc1Gdff_flt2(ii,RecTime_spks{ii}(jj):RecTime_spks{ii}(jj)+ceil(.5*fs_2p)))-mean(Fc1Gdff_flt2(ii,1:RecTime_spks{ii}(jj)));
            spks_traceCaDiffR{ii}(jj)=mean(Fc1Rdff_flt2(ii,RecTime_spks{ii}(jj):RecTime_spks{ii}(jj)+ceil(.5*fs_2p)))-mean(Fc1Rdff_flt2(ii,1:RecTime_spks{ii}(jj)));
        else
            spks_traceCaDiffG{ii}(jj)=mean(Fc1Gdff_flt2(ii,RecTime_spks{ii}(jj):RecTime_spks{ii}(jj)+ceil(.5*fs_2p)))-mean(Fc1Gdff_flt2(ii,RecTime_spks{ii}(jj)-ceil(.5*fs_2p):RecTime_spks{ii}(jj)));
            spks_traceCaDiffR{ii}(jj)=mean(Fc1Rdff_flt2(ii,RecTime_spks{ii}(jj):RecTime_spks{ii}(jj)+ceil(.5*fs_2p)))-mean(Fc1Rdff_flt2(ii,RecTime_spks{ii}(jj)-ceil(.5*fs_2p):RecTime_spks{ii}(jj)));
        end
    end
    RecTime_spks{ii}=RecTime_spks{ii}./fs_2p;
    spks_trace_avgG(ii,:)=nanmean(spks_traceG{ii},1);
    spks_trace_avgR(ii,:)=nanmean(spks_traceR{ii},1);
    spks_trace_stdG(ii,:)=nanstd(spks_traceG{ii},0,1);
    spks_trace_stdR(ii,:)=nanstd(spks_traceR{ii},0,1);
end
spks_trace_avgWtdG=nansum((spks_trace_avgG.*cellfun(@length,RecTime_spks)),1)./sum(cellfun(@length,RecTime_spks));
spks_trace_avgWtdR=nansum((spks_trace_avgR.*cellfun(@length,RecTime_spks)),1)./sum(cellfun(@length,RecTime_spks));
spks_trace_df=cellfun(@length,RecTime_spks)-1;
spks_trace_empty=cellfun(@isempty,RecTime_spks);
spks_trace_df(spks_trace_df==-1)=NaN;%correct for empty cells

% Pooled Standard devitation (not sure if can calculate a pooled standard error)
spks_trace_stdpG=sqrt(nansum(spks_trace_df.*spks_trace_stdG.^2,1)./nansum(spks_trace_df,1));
spks_trace_stdpR=sqrt(nansum(spks_trace_df.*spks_trace_stdR.^2,1)./nansum(spks_trace_df,1));


%%
if plotON==1

%% indiv spiks
figure('Position',[50 50 1600 800])
cc=1;
for ii=1:size(spks_traceG,1)
    for jj=1:length(RecTime_spks{ii})
        if ~isnan(nanmean(spks_traceG{ii}(jj,:))) && cc<=200
            subplot(10,20,cc)
            plot(spks_traceG{ii}(jj,:))
            cc=cc+1;
        end
    end
end


%% indiv cell avg spikes
figure('Position',[50 50 1600 800])
cc=1;
for ii=1:size(spks_traceG,1)
    if ~isnan(nanmean(nanmean(spks_traceG{ii}))) && cc<=200
        subplot(10,20,cc)
        plot(nanmean(spks_traceG{ii},1))
        cc=cc+1;
    end
end
%% Average Spike weighted vs not
figure
plot(spks_trace_avgWtdG,'b')
hold on
plot(nanmean(spks_trace_avgG,1),'r')

%% Average Spike with pooled std
figure
tempX=[-win_pre:1/fs_2p:win_post];
tempYlG=(spks_trace_avgWtdG-spks_trace_stdpG);
tempYuG=(spks_trace_avgWtdG+spks_trace_stdpG);
tempYlR=(spks_trace_avgWtdR-spks_trace_stdpR);
tempYuR=(spks_trace_avgWtdR+spks_trace_stdpR);
%Green
plot(tempX,spks_trace_avgWtdG,'color',[0 .7 0],'linewidth',2)
hold on
plot(tempX,spks_trace_avgWtdG+spks_trace_stdpG,'color',[0 1 0],'linestyle',':')
plot(tempX,spks_trace_avgWtdG-spks_trace_stdpG,'color',[0 1 0],'linestyle',':')
patch([tempX fliplr(tempX)], [tempYlG fliplr(tempYuG)], [.5 1 .5],'FaceAlpha',0.2,'EdgeColor','none')
%Red
plot(tempX,spks_trace_avgWtdR,'color',[.7 0 .7],'linewidth',2)
plot(tempX,spks_trace_avgWtdR+spks_trace_stdpR,'color',[1 0 1],'linestyle',':')
plot(tempX,spks_trace_avgWtdR-spks_trace_stdpR,'color',[1 0 1],'linestyle',':')
patch([tempX fliplr(tempX)], [tempYlR fliplr(tempYuR)], [1 .5 1],'FaceAlpha',0.2,'EdgeColor','none')
xlim([-win_pre win_post])
ylim([-.5 2.1])
pbaspect([1 1 1])

    
%% Indiv rec
figure
plot(speedDataTS1,speedData1,'color',[0 0 0])
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    plot(F1_ts,Fc1Gdff_flt2(ii,:)+ii*10,'color', [0 1 0]);
    line([RecTime_spks{ii,:};RecTime_spks{ii,:}],repmat([ii*10+1;ii*10+1.75],1,size(RecTime_spks{ii,:},2)),'color',[0 0 0],'LineWidth',1)
end

%% Raster

%Vectorize the RecTimes
tempVec=cell2mat(cellfun(@(x) [x nan(1,10-numel(x))],RecTime_spks,'uni',0));
tempVec=transpose(tempVec(~isnan(tempVec)));

figure('Position',[100 100 400 600])
subplot(5,1,1)
plot(speedDataTS1,speedData1,'color',[0 0 0])
xlim([0 300])

subplot(5,1,2)
histogram(tempVec,[1:1:300],'FaceColor',[.5 .5 .5],'EdgeColor',[.25 .25 .25])
xlim([0 300])

subplot(5,1,3:5)
%Spike raster
for jj=1:size(RecTime_spks,1)
line(repmat(RecTime_spks{jj},2,1),repmat([jj-0.5;jj+0.5],1,length(RecTime_spks{jj})),'color',[0 0 0],'linewidth',1)
hold on
end
clear jj
ylim([0 size(RecTime_spks,1)])
yticks(0:50:300)
xlim([0 300])


%% Ca Levels Plot (Pooled Mean+Pooled SEM)

[Ca_Spks_means,Ca_Spks_std]=deal(nan(size(RecTime_spks,1),2));
Ca_Spks_means(:,1)=cellfun(@nanmean,spks_traceCaDiffG);
Ca_Spks_means(:,2)=cellfun(@nanmean,spks_traceCaDiffR);
Ca_Spks_std(:,1)=cellfun(@nanstd,spks_traceCaDiffG);
Ca_Spks_std(:,2)=cellfun(@nanstd,spks_traceCaDiffR);
Ca_Spks_df=spks_trace_df;

% Pooled SE
Ca_Spks_stdp=sqrt(nansum(Ca_Spks_df.*Ca_Spks_std.^2,1)./nansum(Ca_Spks_df,1));
%Ca_Spks_sep=Ca_Spks_stdp.*sqrt(nansum(1./(Ca_Spks_df+1)));


figure('Position',[100 100 200 400])
plotYvals=[nanmean(Ca_Spks_means(:,1)); nanmean(Ca_Spks_means(:,2))];
%plotYerr=[Ca_Spks_sep(1); Ca_Spks_sep(2)];
plotYerr=[Ca_Spks_stdp(1); Ca_Spks_stdp(2)];
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
%ylim([-.25 .75])
xticks(1.5)
xticklabels({'Spks'})
yticks(-.25:.25:1.75)
ylim([-.25 1.75])
set(gca,'box','off')


end


%% Generate Data Structure
%General Variables
StructFieldNames={'szID'
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
    %Motion Tracking
    'speedDataTS1'
    'speedData1'
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
    %Spike Analysis
    'RecTime_spks'
    'spks_traceG'
    'spks_traceR'
    'spks_trace_avgG'
    'spks_trace_avgR'
    'spks_trace_stdG'
    'spks_trace_stdR'
    'spks_trace_avgWtdG'
    'spks_trace_avgWtdR'
    'spks_trace_empty'
    'spks_trace_df'
    'spks_trace_stdpG'
    'spks_trace_stdpR'
    'spks_traceCaDiffG'
    'spks_traceCaDiffR'};
    
for kk=1:length(StructFieldNames)
    eval(strcat('szStruct.',StructFieldNames{kk},'=',StructFieldNames{kk},';'));
end


%% SAVE Structure Varaible as the SeizureID (szID)
assignin('base', szID, szStruct);
save(strcat(szID,'.mat'),szID,'-v7.3')
disp('saved')




