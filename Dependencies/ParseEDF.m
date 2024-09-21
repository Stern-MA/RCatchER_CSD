%This function reads a Pinnacle EEG file in EDF format along with an 
%annotation file in txt format and partitions the EEG file into segments
%Paritioning occurs by taking the EEG segment starting at every odd 
%numbered pulse and ending with the subsequent event pulse .
%Output:
%Cell array where each cell is a row vector of the EEG segment between TTL
%pulses.

%Created 210918-Matthew A. Stern
%Current Update: 210919-MAS

function [EEG] = ParseEDF(EDFfile,AnnotFile,fs,plotFlag)
%% Set Defaults

if ~exist('plotFlag')
    plotFlag = 1;
end

%% Loader
EEG.fs=fs;
[~,EEG.full]=edfread(EDFfile);
EEG.full=EEG.full(1,:);
%load annotations as cell array
Annot = table2cell(readtable(AnnotFile,'Format','%s%s%s%f%s%s'));

%% Extract Annotation Time Stamps
EEG.TTL=Annot(1:5:end-1,4:2:6);
if floor(size(EEG.TTL,1)/2)<size(EEG.TTL,1)/2
    EEG.TTL(size(EEG.TTL,1)+1,:)=Annot(end,4:2:6);
end
%Limit cell array to timestamp and type (rise or fall) columns only

%% Partition EEG
%Seperate EEG files into segments dictated by EEG annotations
EEG.part=cell(size(EEG.TTL,1)/2,1);
for ii=1:2:size(EEG.TTL,1)
    EEG.part((ii+1)/2)={EEG.full(EEG.TTL{ii,1}*fs:EEG.TTL{ii+1,1}*fs)};
end
clear ii

%% Plot EEG with Time Stamp annotations
if plotFlag==1
figure
subplot(size(EEG.part,1),2,[1:2:size(EEG.part,1)*2])
plot([1/fs:1/fs:length(EEG.full)/fs],EEG.full, 'color', [0 0 0])
hold on
line(repmat(cell2mat(EEG.TTL(:,1))',2,1),repmat([min(EEG.full); max(EEG.full)],1,size(EEG.TTL,1)),'color',[1 0 0])
for ii=1:size(EEG.part,1)
    subplot(size(EEG.part,1),2,ii*2)
    plot([1/fs:1/fs:length(EEG.part{ii})/fs],EEG.part{ii}, 'color', [0 0 0])
end
end


