%%
%analysis of EEG in two-photon data to investigate changes in interictal
%activity between pre-seizure and post-CSD
%
%run script to process data and generate figure 6

load('eeg_rcatcher_erc.mat')
    %structure array containing raw EEG, DC recording and metadata describing
    % CSD/not CSD during multiple seizures

%%
code_nums = [1,2,3];
exp_code = {'Seizure only','Stim CSD','Seizure + CSD'};

fs = 2000;
sz_ind = 18;
sd_thresh = 8;

mark_iis = true;

tic
iis_results = BO_IIS_Detect_Power_Rolston_v5(eeg_str(sz_ind).EEG',sd_thresh,fs,0.7,0);
toc

spk_samps = iis_results.ind_Sp_final{1};

t_spk = zeros(size(eeg_str(sz_ind).EEG));
t_spk(spk_samps) = 1;

t_plot = linspace(0,length(eeg_str(sz_ind).EEG')/fs,length(eeg_str(sz_ind).EEG'));
eeg_z = (eeg_str(sz_ind).EEG - mean(eeg_str(sz_ind).EEG))./std(eeg_str(sz_ind).EEG);

[counts,tt] = count_movwind(spk_samps, 5, 0.1, length(eeg_str(sz_ind).EEG),fs);

EEG = eeg_str(sz_ind).EEG;

fs_eeg=2000;
PSDbin=1024;

m_avg = 5;
sd_thresh = 1;
stim_thresh = 4;
seg_length = 10;

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
%RA_Ratio100=SternRollAvg(Ratio100,60/PSDwin); % +/- 30 second window roll avg (time step here is 0.5s)
% SzStartTime = 410;
% SzEndTime=find(PSD_100less(floor(SzStartTime/PSDwin*2):end)./max(PSD_100less)<0.05,1)*PSDwin/2+floor(SzStartTime);

p_mavg = movmean(PSD_totP,m_avg,'omitnan');
p_mavg = (p_mavg - median(p_mavg))./(std(p_mavg));
p_mavg_thr = p_mavg > sd_thresh;

endinds =  reshape(find(diff([0,p_mavg_thr,0])~=0),2,[]); 
for kk = 1:size(endinds,2)
    if (endinds(2,kk) - endinds(1,kk)) > seg_length
        SzStartW = endinds(1,kk);
        SzEndW = endinds(2,kk);
        break
    end
end

SzStartTime = SzStartW/PSDwin/2;
SzEndTime = SzEndW/PSDwin/2;

%% Process seizure data for all files
seg_length = 15;
sd_thresh = 4;
sz_thresh = 1.5;
stim_thresh = 4;
iis_win = 300; %num. seconds before/after seizure 
use_manual_times = true;
use_first_spike = true;
stim_offset = 330; 

sub_spec_thresh = true;

fs = 2000;

sz_times = nan(length(eeg_str),2);
stim_times = nan(length(eeg_str),2);
iis_times = cell(size(eeg_str));
iis_log = cell(size(eeg_str));

iis_sz = nan(length(eeg_str),2);
iis_stim = nan(length(eeg_str),2);
%iis_post_sz = nan(size(eeg_str));

sz_start_og = cell2mat({eeg_str.sz_start})';
sz_end_og = cell2mat({eeg_str.sz_end})';
if use_manual_times
   sz_end_og(10) = 450; 
end

for kk = 1:length(eeg_str)
    tic
    fprintf('Starting file %d\n',kk)
    if (sub_spec_thresh) && (kk == 16 || kk == 18)
        tempthr = 7;
        iis_results = BO_IIS_Detect_Power_Rolston_v5(eeg_str(kk).EEG',tempthr,fs,0.7,0);
    else
        iis_results = BO_IIS_Detect_Power_Rolston_v5(eeg_str(kk).EEG',sd_thresh,fs,0.7,0);
    end
    iis_times{kk} = iis_results.time_Sp_final{1};
    
    EEG = eeg_str(kk).EEG;
    
    temp_log = zeros(size(EEG));
    temp_log(iis_results.ind_Sp_final{1}) = 1;
    iis_log{kk} = logical(temp_log);
    
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
    %RA_Ratio100=SternRollAvg(Ratio100,60/PSDwin); % +/- 30 second window roll avg (time step here is 0.5s)
    % SzStartTime = 410;
    % SzEndTime=find(PSD_100less(floor(SzStartTime/PSDwin*2):end)./max(PSD_100less)<0.05,1)*PSDwin/2+floor(SzStartTime);
    
    p_mavg = movmean(PSD_totP,m_avg,'omitnan');
    p_mavg = (p_mavg - median(p_mavg))./(std(p_mavg));
    p_mavg_thr = p_mavg > sz_thresh;
    
    endinds =  reshape(find(diff([0,p_mavg_thr,0])~=0),2,[]);
    
%     for kl = 1:size(endinds,2)
%         if (endinds(2,kl) - endinds(1,kl)) > seg_length
%             SzStartW = endinds(1,kl);
%             SzEndW = endinds(2,kl);
%             break
%         end
%     end

    [~,closest_sz_ind] = min(abs(mean(endinds)/PSDwin/2-sz_start_og(kk)));
    SzStartW = endinds(1,closest_sz_ind);
    
    [~,closest_sz_ind] = min(abs(mean(endinds)/PSDwin/2-sz_end_og(kk)));
	SzEndW = endinds(2,closest_sz_ind);
    
    sz_times(kk,1)= SzStartW/PSDwin/2;
    sz_times(kk,2) = SzEndW/PSDwin/2;
    
    if use_manual_times
        sz_times(1,:) = [950, 965];
        sz_times(2,:) = [255, 265];
        sz_end_og(10,2) = 450;
        sz_times(12,:) = [602, 707];
    end
    
    if use_first_spike
        %pre_win = sz_times(kk,1) - eeg_str(kk).first_iis;
        if kk == 18
            pre_dur = (sz_times(kk,1) - 0)/60;
            iis_sz(kk,1) = sum(and(iis_times{kk}<sz_times(kk,1),iis_times{kk}>0))/pre_dur;
        else
            pre_dur = (sz_times(kk,1) - eeg_str(kk).first_iis)/60;
            iis_sz(kk,1) = sum(and(iis_times{kk}<sz_times(kk,1),iis_times{kk}>eeg_str(kk).first_iis))/pre_dur;
        end
        
    else
        iis_sz(kk,1) = sum(and(iis_times{kk}<sz_times(kk,1),iis_times{kk}>sz_times(kk,1)-iis_win));
    end
    iis_sz(kk,2) = sum(and(iis_times{kk}>sz_times(kk,2),iis_times{kk}<sz_times(kk,2)+iis_win))/(iis_win/60);
    
    if eeg_str(kk).recType == 2
        p_mavg_thr = p_mavg > stim_thresh;
        
        endinds =  reshape(find(diff([0,p_mavg_thr,0])~=0),2,[]);
        for kl = 1:size(endinds,2)
            if (endinds(2,kl) - endinds(1,kl)) > seg_length
                StimStartW = endinds(1,kl);
                StimEndW = endinds(2,kl);
                break
            end
        end
        
        stim_times(kk,1)= StimStartW/PSDwin/2;
        stim_times(kk,2) = StimEndW/PSDwin/2;
        
%         if use_first_spike
%         %pre_win = sz_times(kk,1) - eeg_str(kk).first_iis;
%             pre_dur = (sz_times(kk,1) - eeg_str(kk).first_iis)/60;
%             iis_stim(kk,1) = sum(and(iis_times{kk}<sz_times(kk,1),iis_times{kk}>eeg_str(kk).first_iis))/pre_dur;
%         else
%             iis_stim(kk,1) = sum(and(iis_times{kk}<stim_times(kk,1),iis_times{kk}>stim_times(kk,1)-iis_win));
%         end
        iis_stim(kk,1) = sum(and(iis_times{kk}<stim_times(kk,1),iis_times{kk}>stim_times(kk,1)-iis_win))/(iis_win/60);
        iis_stim(kk,2) = sum(and(iis_times{kk}>stim_times(kk,2),iis_times{kk}<stim_times(kk,2)+iis_win))/(iis_win/60);
    
    elseif eeg_str(kk).recType == 1
        iis_stim(kk,2) = sum(and(iis_times{kk}>sz_times(kk,1)+stim_offset,iis_times{kk}<sz_times(kk,1)+stim_offset+iis_win))/(iis_win/60);
    
    end
    
    toc
    
end

%% final figure 6 layout:
showmarkers = true;

rectype = [eeg_str.recType];

sp_y1 = [0.74, 0.535, 0.375]; sph = 0.2;
sp_x1 = [.025, 0.3, 0.6]; spw1 = 0.25; spw2 = 0.35;

sp_pos = {[sp_x1(1)    sp_y1(1)    spw1    sph]
          [sp_x1(1)    sp_y1(2)    spw1    sph]
          [sp_x1(1)    sp_y1(3)    spw1    sph*.75]
          [sp_x1(2)    sp_y1(1)    spw1    sph]
          [sp_x1(2)    sp_y1(2)    spw1    sph]
          [sp_x1(2)    sp_y1(3)    spw1    sph*.75]
          [sp_x1(3)    sp_y1(1)    spw2    sph]
          [sp_x1(3)    sp_y1(2)    spw2    sph]
          [sp_x1(3)    sp_y1(3)    spw2    sph*.75]
          [.0575    .075      0.45    0.225]
          [.63    .075      0.32    0.225]};
      
figure('Position',[100 100 800 700])
for kk = 1:length(sp_pos)
subplot('Position',sp_pos{kk})
end

sz_exs = [10, 13, 6]; %seizure, stim, and CSD examples

gk = gausswin(fs*30, 2);
gk_ratio = 2*sum(ones(size(gk)))/sum(gk);
    %adjust weight of gaussian kernel to reflect spikes per minute

pre_lim = 300;
post_lim = 300;
ylim_eeg = [-2000 2000];
ylim_dc = [-20 10];

titles = {'Seizure Only','Seizure with CSD', 'Stimulation'};


for kk = 1:3
    sz_ind = sz_exs(kk);
    %gk = gk/sum(gk)
    t_plot = linspace(0,length(eeg_str(sz_ind).EEG')/fs,length(eeg_str(sz_ind).EEG'));
    
    
    if kk == 3
        %%% to plot just pre- and post- stim block:
%         t_plot = t_plot - stim_times(sz_ind,1);
%         ind1 = find(t_plot>(-pre_lim),1);
%         ind2 = find(t_plot>post_lim,1);
% 
%         sz_len = abs(diff(sz_times(sz_ind,:)));
%         %xlims = [-pre_lim, stim_times(sz_ind,2) - sz_times(sz_ind,1)+post_lim];
%         xlims = [-pre_lim, stim_times(sz_ind,2) - sz_times(sz_ind,1)+post_lim];
%         
%         stim_block = stim_times(sz_ind,:) - sz_times(sz_ind,1);
        
        %%% to plot both seizure and stim:
        t_plot = t_plot - sz_times(sz_ind,1);
        pre_new = -pre_lim;
        post_new = stim_times(sz_ind,2) - sz_times(sz_ind,1)+post_lim;
        ind1 = find(t_plot>(pre_new),1);
        ind2 = find(t_plot>post_new,1);
        
        sz_len = abs(diff(sz_times(sz_ind,:)));
        
        xlims = [pre_new, post_new];
        
        stim_block = stim_times(sz_ind,:) - sz_times(sz_ind,1);
        
    else
        if kk == 1 %plot 1 and 2 on same time scale
            sz_len = abs(diff(sz_times(sz_ind,:)));
        end
        
        xlims = [-pre_lim, sz_len+post_lim];
        
        t_plot = t_plot - sz_times(sz_ind,1);
        ind1 = find(t_plot>(-pre_lim),1);
        ind2 = find(t_plot>sz_len+post_lim,1);
    end
    %sz_block = sz_times(sz_ind,:) - sz_times(sz_ind,1);
    sz_block = [sz_times(sz_ind,1) - sz_times(sz_ind,1), sz_times(sz_ind,2) - sz_times(sz_ind,1)];
 
    fr_mavg = conv(iis_log{sz_ind}, gk, 'same')*gk_ratio;
    
    
    subplot('Position',sp_pos{(kk-1)*3+1})
    
    md = median(eeg_str(sz_ind).EEG(ind1:ind2));
    plot(t_plot(ind1:ind2),eeg_str(sz_ind).EEG(ind1:ind2)-md,'Color',[0.5 0.5 0.5])
    hold on
%     ylim_temp = ylim;
%     ylim([-max(abs(ylim_temp)),max(abs(ylim_temp))])
    ylim(ylim_eeg)
    %patch([0 0 sz_len sz_len], [ylim_temp(1) ylim_temp(2) ylim_temp(2) ylim_temp(1)], [0.8 0.8 0.8])
%     if kk == 3
%         patch([stim_block(1) stim_block(1) stim_block(2) stim_block(2)], [ylim_temp(1) ylim_temp(2) ylim_temp(2) ylim_temp(1)], [0.8 0 0])
%     end
%    plot(t_plot,eeg_str(sz_ind).EEG,'k')
%    xlim(xlims)
    
    hold on
%     xline(sz_times(sz_ind,1),'k-')
%     xline(sz_times(sz_ind,2),'k-')
%     if show_iis
%         for kk = 1:length(iis_times{sz_ind})
%             xline(iis_times{sz_ind}(kk),'r--')
%         end
%     end

    %title(titles{kk},'FontSize',16)
    xlim(xlims)
    
    hold on 
    
    if kk == 1
        y_corner = max(ylim)*0.75; x_corner = 180;
        plot([x_corner,x_corner],[y_corner,y_corner+500],'k','LineWidth',1.5)
        plot([x_corner,x_corner + 60],[y_corner,y_corner],'k','LineWidth',1.5)
        text(x_corner+40,y_corner-.04*diff(ylim),'1 min','HorizontalAlignment', 'center')
        text_y = text(x_corner-30,y_corner-0*diff(ylim),'500 \muV');
        text_y.Rotation = 90;
    end
    
    if kk == 3
        y_corner = max(ylim)*0.75; x_corner = 500;
        plot([x_corner,x_corner],[y_corner,y_corner+500],'k','LineWidth',1.5)
        plot([x_corner,x_corner + 60],[y_corner,y_corner],'k','LineWidth',1.5)
        text(x_corner+40,y_corner-.04*diff(ylim),'1 min','HorizontalAlignment', 'center')
        text_y = text(x_corner-30,y_corner-0*diff(ylim),'500 \muV');
        text_y.Rotation = 90;
    end

    set(gca, 'Visible', 'off')
    
    subplot('Position',sp_pos{(kk-1)*3+2})
    
    md = median(eeg_str(sz_ind).DC(ind1:ind2));
    
    plot(t_plot(ind1:ind2),eeg_str(sz_ind).DC(ind1:ind2)-md,'k')
    xlim(xlims)
    disp(xlims)
    if kk == 3
        %ylim([-10 10])
        ylim(ylim_dc+5)
    else
        ylim(ylim_dc)
    end
    if kk == 1
        hold on
        y_corner = min(ylim)*0.8; x_corner = 180;
        plot([x_corner,x_corner],[y_corner,y_corner+5],'k','LineWidth',1.5)
        plot([x_corner,x_corner + 60],[y_corner,y_corner],'k','LineWidth',1.5)
        text(x_corner+40,y_corner-0.04*diff(ylim),'1 min','HorizontalAlignment', 'center')
        text_y = text(x_corner-30,y_corner+.02*diff(ylim),'5 mV');
        text_y.Rotation = 90;
 
    end
    
    if kk == 3
        hold on
        y_corner = min(ylim)*0.8; x_corner = 500;
        plot([x_corner,x_corner],[y_corner,y_corner+5],'k','LineWidth',1.5)
        plot([x_corner,x_corner + 60],[y_corner,y_corner],'k','LineWidth',1.5)
        text(x_corner+40,y_corner-0.04*diff(ylim),'1 min','HorizontalAlignment', 'center')
        text_y = text(x_corner-30,y_corner+.02*diff(ylim),'5 mV');
        text_y.Rotation = 90;
 
    end
    
    set(gca, 'Visible', 'off')
    
    
    subplot('Position',sp_pos{(kk-1)*3+3})
    
    hold on
%     patch([0 0 sz_len sz_len], [ylim_temp(1) ylim_temp(2) ylim_temp(2) ylim_temp(1)], [0.8 0.8 0.8])
%     if kk == 3
%         patch([stim_block(1) stim_block(1) stim_block(2) stim_block(2)], [ylim_temp(1) ylim_temp(2) ylim_temp(2) ylim_temp(1)], [0.8 0 0])
%     end
    plot(t_plot(ind1:ind2),fr_mavg(ind1:ind2),'b')
    hold on
    
%     xline(sz_times(sz_ind,1),'k-')
%     xline(sz_times(sz_ind,2),'k-')
%     ylabel('Spike rate (count/min)')
%     xlabel('Time (s)')
    xlabel('')
    
    xticks(xlims(1):150:xlims(2));
    xticklabels((xlims(1):150:xlims(2))/60);
    set(gca,'FontSize',14)
    
    if kk > 1
        yticks([])
    end
    ylim([-0.5,50])
    ylimtemp = ylim;
    patch([sz_block(1), sz_block(2),sz_block(2),sz_block(1)], [ylimtemp(1),ylimtemp(1),ylimtemp(2),ylimtemp(2),], [0.4 0.4 0.4])
    if kk==3
        patch([stim_block(1), stim_block(2),stim_block(2),stim_block(1)], [ylimtemp(1),ylimtemp(1),ylimtemp(2),ylimtemp(2),], [0.7 0 0.1])
    end
    
    xlim(xlims)
    box off
    
%     if kk == 1
%         y_corner = max(ylim)*0.55; x_corner = 180;
%         plot([x_corner,x_corner],[y_corner,y_corner+10],'k','LineWidth',1.5)
%         plot([x_corner,x_corner + 60],[y_corner,y_corner],'k','LineWidth',1.5)
%         text(x_corner+30,y_corner-3,'60 s','HorizontalAlignment', 'center')
%         text_y = text(x_corner-30,y_corner-0.04*diff(ylim),'10 SWD/min');
%         text_y.Rotation = 90;
%     end
%     if kk == 3
%         y_corner = max(ylim)*0.55; x_corner = 500;
%         plot([x_corner,x_corner],[y_corner,y_corner+10],'k','LineWidth',1.5)
%         plot([x_corner,x_corner + 60],[y_corner,y_corner],'k','LineWidth',1.5)
%         text(x_corner+45,y_corner-.04*diff(ylim),'60 s','HorizontalAlignment', 'center')
%         text_y = text(x_corner-45,y_corner-.04*diff(ylim),'10 SWD/min');
%         text_y.Rotation = 90;
%     end
end

iis_cat = padcat(iis_sz(rectype==1,1),iis_sz(rectype==1,2), iis_sz(rectype==3,1),iis_sz(rectype==3,2));

xvals = [1,2,4,5];
cols = {[0.4,0,0.5],[0.4,0,0.5], [0,0.5,0.5],[0,0.5,0.5]};

%violinplot(stim_cat)
subplot('Position', sp_pos{10})
hold on
for kk = 1:4
    if rem(kk,2) == 0
        boxchart(xvals(kk)*ones(size(iis_cat,1),1),iis_cat(:,kk),'BoxFaceColor',cols{kk},'MarkerColor',cols{kk},'MarkerStyle','none','HandleVisibility','off')
    else
        boxchart(xvals(kk)*ones(size(iis_cat,1),1),iis_cat(:,kk),'BoxFaceColor',cols{kk},'MarkerColor',cols{kk},'MarkerStyle','none')
    end
    if showmarkers
        hold on
        xdat = xvals(kk)*ones(size(iis_cat,1),1) + 0.075*randn(size(iis_cat,1),1);
        scatter(xdat, iis_cat(:,kk),[],cols{kk},'filled','HandleVisibility','off')
    end
end
%legend({'Without CSD','With CSD'})
xticks(xvals)
ylim([0,35])
ylabel('SWD rate (/min.)')
xticklabels({'Pre','Post','Pre','Post'})
set(gca,'FontSize',14)

stim_cat = [iis_sz(rectype==2,1),iis_stim(rectype==2,1),iis_stim(rectype==2,2)];
bvals = [0.8,0.5,0.2];

subplot('Position', sp_pos{11})
hold on
for kk = 1:3
    boxchart(kk*ones(size(stim_cat,1),1),stim_cat(:,kk),'BoxFaceColor',[0.6,0.2,0]*bvals(kk))
    if showmarkers
        hold on
        xdat = kk*ones(size(stim_cat,1),1) + 0.075*randn(size(stim_cat,1),1);
        scatter(xdat, stim_cat(:,kk),[],[0.6,0.2,0]*bvals(kk),'filled','HandleVisibility','off')
    end
end
ylabel('SWD rate (/min.)')
xticks([1,2,3])
xticklabels({'Pre-seizure','Pre-stim','Post-stim'})
%xtickangle(30)
ylim([0,35])
set(gca,'FontSize',14)
xlim([0.5 3.5])

% annotation('textbox',[.0075 .96 .05 .05], 'String','A','FontSize',36,'EdgeColor','None')
% annotation('textbox',[.0075 .175 .05 .05], 'String','B','FontSize',36,'EdgeColor','None')
% annotation('textbox',[.54 .96 .05 .05], 'String','C','FontSize',36,'EdgeColor','None')
% annotation('textbox',[.54 .175 .05 .05], 'String','D','FontSize',36,'EdgeColor','None')
annotation('textbox',[.11    .001      0.25    0.04], 'String','Without CSD','FontSize',16,'EdgeColor','None')
annotation('textbox',[.345    .001      0.25    0.04], 'String','With CSD','FontSize',16,'EdgeColor','None')
annotation('textbox',[.17    .335      0.35    0.01], 'String','Time from seizure start (min.)','FontSize',16,'EdgeColor','None')
annotation('textbox',[.645    .335      0.35    0.01], 'String','Time from seizure start (min.)','FontSize',16,'EdgeColor','None')

%xlabel('Time from seizure start (min.)');

rectype = [eeg_str.recType]';

sz_diff = iis_sz(rectype==1,2)-iis_sz(rectype==1,1);
stim_diff = iis_stim(rectype==2,2) - iis_stim(rectype==2,1);
csd_diff = iis_sz(rectype==3,2)-iis_sz(rectype==3,1);
iis_diffs = padcat(sz_diff,stim_diff,csd_diff);


p = signrank(sz_diff);
%fprintf('Seizure-only group vs. 0: %1.4f\n',p)

p = signrank(stim_diff);
%fprintf('Stim group vs. 0: %1.4f\n',p)

p = signrank(csd_diff);
%fprintf('CSD group vs. 0: %1.4f\n',p)

[p, table, stats] = friedman(stim_cat);
q = multcompare(stats);

%fprintf('Stim w/ 3 groups: %1.4f\n',p)


function [counts, tt] = count_movwind(events, win_size, overlap, dur, fs)
    %events: vector of samples for detected events
    %win_size: length of moving window in secs
    %overlap: fraction of wsize to advance
    n_winds = floor(dur/(overlap*win_size*fs));
    counts = nan(n_winds,1);
    tt = nan(n_winds,1);
    wintemp = round([0, win_size]*fs);
    win_adv = overlap*win_size*fs;
    for kk = 1:length(counts)
        counts(kk) = sum(and(events>wintemp(1), events < wintemp(2)));
        tt(kk) = (kk-1)*(overlap*win_size) + win_size/2;
        wintemp = wintemp + win_adv;
    end
    
end

