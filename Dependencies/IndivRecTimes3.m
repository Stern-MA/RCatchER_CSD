% Determine Seizure related wavefront recruitment times in individual
% traces of calcium data
% *** version using MAX SLOPE time indexing

%this script take a matrix of calcium traces (cells by row, time by column)
%and a seed time for an event (e.g. TSW or Seizure) and finds the maximum
%period where the integral of first deriviative (positive portions only) 
%of the trace is largest (selecting for the longest sustatined and steepist
%slope). It weighs these features by a guassian kernel (sigmal1) centered
%at the seed time.

%Inputs:
%CaTraceMatrix: matrix of calcium traces (cells by row, time by column)
%SeedTime: seed time for an event (e.g. TSW or Seizure)
%fs_2p: sampling frequency of 2p imaging
%sigma1: Standard deviation for width of gaussian weighting
%
%Output Structure Fields:
%time: recruitment times
%value: slope integral value of trace at recruitment 

% Version 231116 Matthew A. Stern and Eric R. Cole, Emory University
% Contact: matthew.a.stern@emory.edu or matt@matthewastern.com

function RecTimes=IndivRecTimes3(CaTraceMatrix,SeedTime,fs_2p,sigma1)

%generate gaussian centered at seed time
trace1=1:size(CaTraceMatrix,2);
gauss1=exp(-(trace1-SeedTime*fs_2p).^2/(sigma1*fs_2p)^2);

%define trace matrix
Fc1Gdff_flt2=CaTraceMatrix;

%calculate the derivates of the filtered traces
Fc1Gdff_flt2_df=diff(Fc1Gdff_flt2,1,2);
%Fc1Gdff_flt2_d2f=diff(Fc1Gdff_flt2,2,2);

% determine all the points that are positive
temp1G=Fc1Gdff_flt2_df>=0;

% find block bounds
temp2G=diff(temp1G,1,2);

temp3G=cell(size(temp2G,1),2);
for ii=1:size(temp3G,1)%green
    temp3G{ii,1}=find(temp2G(ii,:)==1);
    temp3G{ii,2}=find(temp2G(ii,:)==-1);
end
clear ii

%integrate over blocks
[temp4G,maxDFwin,temp5G,temp6G]=deal(cell(size(temp2G,1),1));

for ii=1:length(temp4G) %GREEN
    if temp3G{ii,1}(1)>temp3G{ii,2}(1)%correct for first slope value index being negative
        temp3G{ii,2}=temp3G{ii,2}(2:end);
    end
    
    if temp3G{ii,1}(end)>temp3G{ii,2}(end)%correct for extra final slope value being positive
        temp3G{ii,1}=temp3G{ii,1}(1:end-1);
    end
    
    for jj=1:length(temp3G{ii,1}) %find the area under the first derivative curve of each postitive region
        temp4G{ii}(jj)=trapz(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
        [maxDFwin{ii}(jj,1),maxDFwin{ii}(jj,2)]=max(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
    end
    %recalculate block index times based upon max slope time 
    temp5G{ii,1}=temp3G{ii,1}+maxDFwin{ii}(:,2)';
    %recalculate the features weighted by the gaussian kernel
    temp6G{ii,1}=temp4G{ii}.*gauss1(temp5G{ii,1});
end
clear ii jj

% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp5 to find the index (time point) of the seizure in the origional trace

RecTimeG_sz=nan(size(temp2G,1),1);
RecTimeG_val=nan(size(temp2G,1),1);

for ii=1:size(temp2G,1)
    % sorting to find top value
    [~, feats_maxI] = max(temp6G{ii,1});
    RecTimeG_sz(ii) = (temp5G{ii,1}(feats_maxI)+1)/fs_2p;
    RecTimeG_val(ii) = temp4G{ii}(feats_maxI);
end
clear ii

RecTimes.time=RecTimeG_sz;
RecTimes.value=RecTimeG_val;
end




