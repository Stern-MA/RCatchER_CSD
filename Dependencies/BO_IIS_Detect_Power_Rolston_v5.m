
function out = BO_IIS_Detect_Power_Rolston_v5(data,S_thresh,Fs,winW,plot_flag)
%adjusted version of John Rolston's (Rolston 2010 Epilepsia) interictal spike detection algorithm:
%now incorporates MATLAB's findpeaks function with a heuristic about spike
%frequency to improve robustness and avoid over-detecting tiny spikes right
%next to bigger ones

%specifically looks for peaks in theta frequency and then filters these based upon  the
%ratio relative theta and gamma power around them (currently .5 seconds
%ahead and 1 second after then determsin the minimum value
%(negative peak) within a window of the spike detected, and has that as the
%peak point

% Version 231116 Matthew A. Stern and Eric R. Cole, Emory University
% Contact: matthew.a.stern@emory.edu or matt@matthewastern.com
%
% If using this code please cite our work.

if ~exist('W')
    W = 0.5;
end
if ~exist('MW')
    MW = 0.1;
end
if ~exist('Fs')
    Fs = 2000;
end

paramssp.Fs                           = Fs;
paramssp.tapers                       = [3 5]; %delta
%paramssp.fpass                        = [4 100];
paramssp.fpass                        = [3 15]; %extended theta
paramssp.fpass1                        = [20 55];%extended gamma
paramssp.notch                        = [4 14]; %remove theta and alpha

%Find all spikes with high theta power
[b,a] = butter(4,paramssp.fpass/(Fs/2),'bandpass');
[b1,a1] = butter(4,paramssp.fpass1/(Fs/2),'bandpass');
%[b2,a2] = butter(4,paramssp.notch/(Fs/2),'stop');
V_filtered = filter(b,a,data');
%V_filtered_temp = filter(b,a,data');
%V_filtered = filter(b2,a2,V_filtered_temp);

for C=1:1:size(V_filtered,1)
    
    ind_up = [];
    ind_down = [];
    ind_Sp_final = [];
    time_Sp_final = [];
    indraw_Sp_final = [];
    
    S_th = S_thresh*median(abs(V_filtered(C,:)))/0.6745; %percentile threshold for IIS detection
    % ind_Sp = logical(V_filtered(C,:)>S_th); %find the index above threshold
    [~,ind_locs] = findpeaks(V_filtered(C,:),'MinPeakHeight',S_th,'MinPeakDistance',winW*Fs);

% Section filters the above spikes by a theta/subsequent gamma ratio and by the low
% gamma power needing to be above the mean
PSDbin=512;
[~,PSD_F]=pwelch(data,[],[],PSDbin,Fs);

EEGstart=ind_locs-.25*Fs;
EEGstop=ind_locs+.5*Fs;
GammaShift=.25*Fs;
PSD_P=nan(length(PSD_F),length(ind_locs));
PSD_P2=nan(length(PSD_F),length(ind_locs));
for ii = 2:numel(ind_locs)-1
    EEG1=data(EEGstart(ii):EEGstop(ii));
    EEG2=data(EEGstart(ii)+GammaShift:EEGstop(ii)+GammaShift);
    [PSD_P(:,ii),~]=pwelch(EEG1,[],[],PSDbin,Fs);
    [PSD_P2(:,ii),~]=pwelch(EEG2,[],[],PSDbin,Fs);
    clear EEG1
end

PSD_totP=sum(PSD_P,1);
PSD_theta=sum(PSD_P(and(PSD_F>3,PSD_F<15),:),1);
PSD_gammaL=sum(PSD_P2(and(PSD_F>20,PSD_F<55),:),1);%Note shifted PSD for GAMMA

TGratio=PSD_theta./PSD_gammaL; %theta to gamma ratio

ind_Sp = ind_locs(TGratio>20 & PSD_gammaL<(nanmean(PSD_gammaL)+0.25*nanstd(PSD_gammaL)));
    
    ind_Sp_final = ind_Sp;
    time_Sp_final = (1/Fs)*ind_Sp;
    indraw_Sp_final = [];
    
% Index based upon the spike of SWD
    minPeak_win = .25*Fs; %window to find the min peak value around the detected peak in the filtered trace
    sp_minPeak = nan(1,length(ind_Sp));
    minPeak = nan(1,length(ind_Sp));
    for jj = 1:length(ind_Sp)
        [minPeak(jj),minPeakI] = min(data(ind_Sp(jj)-minPeak_win:ind_Sp(jj)+minPeak_win,C));
        sp_minPeak(jj)=(ind_Sp(jj)-minPeak_win+minPeakI)/Fs;
    end
    
%Outputs
    out.ind_Sp_final{C} = ind_Sp_final;
    out.time_Sp_final{C} = time_Sp_final;
    out.indraw_Sp_final{C} = indraw_Sp_final;
    out.time_Sp_minPeak{C} = sp_minPeak;
    out.peak_Sp_minPeak{C} = minPeak;
end


%% plotting
if plot_flag == 1
figure(99)
fig_size = [4,8];
dT = t(2)-t(1);
for i=1:1:size(data,1)
    subplot(fig_size(1),fig_size(2),2*(i-1)+1)
    plot([1:1:size(data,2)]/2000,data(i,:))
    hold on
    subplot(fig_size(1),fig_size(2),2*(i-1)+2)
    plot(t,S_sum(:,i),'b')
    hold on
    for j=1:1:size(out.time_Sp_final{i},1)
        [~,ind1] = min(abs(t-out.time_Sp_final{i}(j,1)));
        [~,ind2] = min(abs(t-out.time_Sp_final{i}(j,2)));
        subplot(fig_size(1),fig_size(2),2*(i-1)+2)
        plot(t(ind1:ind2),S_sum(ind1:ind2,i),'r')
        subplot(fig_size(1),fig_size(2),2*(i-1)+1)
        plot([out.indraw_Sp_final{i}(j,1):out.indraw_Sp_final{i}(j,2)]/2000,data(i,round([out.indraw_Sp_final{i}(j,1):out.indraw_Sp_final{i}(j,2)])),'r')
    end
    subplot(fig_size(1),fig_size(2),2*(i-1)+1)
    hold off
    subplot(fig_size(1),fig_size(2),2*(i-1)+2)
    hold off

end
end