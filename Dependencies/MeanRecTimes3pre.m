% Determine PREICTAL Recruitment Times - Mean traces

%this script take a matrix of calcium traces (cells by row, time by column)
%and finds the points (top percentage: TopP) of steepest recruitment of the 
%trace using the slope integral feature (selecting for the longest 
%sustatined and steepist slope). It then indexes these regions by the point
%of steepest slope (column 1), or the elbow of recruitment (the max 
%concavity; 2nd derviative; column 2). 
%
%Outputs cell array columns"
% 1: Recruitment time by max slope
% 2: Recruitment time by max concavity
% 3: feature value of the slope integral 
% 4: normalized slope integral feature to the max value of that feature

% Version 231116 Matthew A. Stern and Eric R. Cole, Emory University
% Contact: matthew.a.stern@emory.edu or matt@matthewastern.com
%
% If using this code please cite our work.

function RecTimes=MeanRecTimes3pre(CaTraceMatrix,fs_2p,topP)

%define trace matrix
Fc1Gdff_flt2=CaTraceMatrix;

%calculate the derivates of the filtered traces
Fc1Gdff_flt2_df=diff(Fc1Gdff_flt2,1,2);
Fc1Gdff_flt2_d2f=diff(Fc1Gdff_flt2,2,2);

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
temp4G=cell(size(temp2G,1),1);
maxDFwin=cell(size(temp2G,1),1);
for ii=1:length(temp4G)
    if temp3G{ii,1}(1)>temp3G{ii,2}(1)%correct for first slope value index being negative
        temp3G{ii,2}=temp3G{ii,2}(2:end);
    end
    
    if temp3G{ii,1}(end)>temp3G{ii,2}(end)%correct for extra final slope value being positive
        temp3G{ii,1}=temp3G{ii,1}(1:end-1);
    end
    
    for jj=1:length(temp3G{ii,1}) %find the area under the first derivative curve of each positive region
        temp4Ga{ii}(jj)=trapz(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
        temp4Gb{ii}(jj)=mean(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
        temp4G{ii}(jj)=temp4Ga{ii}(jj).*temp4Gb{ii}(jj);
        %Also find the max slope in this region
        [maxDFwin{ii}(jj,1),maxDFwin{ii}(jj,2)]=max(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
        %Also find the max second derivative in this region
        [maxD2Fwin{ii}(jj,1),maxD2Fwin{ii}(jj,2)]=max(Fc1Gdff_flt2_d2f(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
    end
end
clear ii jj

% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(size(temp2G,1),1);%gives the block in each cell that is seizure
RecTimeG=nan(size(temp2G,1),1);
RecTimes=cell(size(temp2G,1),4);

RecTimeG_sz=nan(size(temp2G,1),1);
RecTimeG_sz2=nan(size(temp2G,1),1);

for ii=1:size(temp2G,1)
    % sorting to find top % of values
    [feats_sorted, feats_sorted_inds] = sort(temp4G{ii},2,'descend');
    topN=ceil(length(temp4G{ii})*topP/100);
    feats_sorted = feats_sorted(1:topN);
    feats_sorted_inds = feats_sorted_inds(1:topN);
    
    %recalculate block index times based upon max slope time 
    temp5G{ii,1}=temp3G{ii,1}+maxDFwin{ii}(:,2)';
    temp6G{ii,1}=temp3G{ii,1}+maxD2Fwin{ii}(:,2)';
    RecTimes{ii,1}=(temp5G{ii,1}(feats_sorted_inds)+1)/fs_2p;
    RecTimes{ii,2}=(temp6G{ii,1}(feats_sorted_inds)+2)/fs_2p;%use +2 to correct for diff 2nd deriv
    RecTimes{ii,3}=feats_sorted;
    %normalize feature values to max value
    temp7=log(feats_sorted);
    RecTimes{ii,4}=(temp7-min(temp7))./(max(temp7)-min(temp7));
end

clear ii

end