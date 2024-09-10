%Generate spectral plot of XCaMP-Y and RCatchER
%extracted from interleaved two channel imaging in vivo

%use Suite2P to extract caclium transients by channel manually removing
%yellow bleedthrough cells in red channel
%load colors seperatly below

%% Loader Ca Data (Red (ch 1); Green (ch 2))
% load('suite2p/plane0/Fall.mat')
% load('suite2p/plane0/F_chan2.mat')
% load('suite2p/plane0/Fneu_chan2.mat')
load('suite2p2/plane0/Fall.mat')
load('suite2p2/plane0/F_chan2.mat')
load('suite2p2/plane0/Fneu_chan2.mat')
F_R=F; F_G=F2; Fneu_R=Fneu; Fneu_G=F2neu;
clear F F2 Fneu F2neu spks

%% Filter out the manual labels
%assumes the manual cells(puncta contamination) are appended at the begining of the array and they
%have a probability of being a cell (iscell column 2) of 1.
mancellsI=sum(iscell(:,2)==1);
% 6 manual labeled cells are the neuropil for subtraction
F1manG=F_G(1:mancellsI,:);
F1manR=F_R(1:mancellsI,:);
Fneu1manG=Fneu_G(1:mancellsI,:);
Fneu1manR=Fneu_R(1:mancellsI,:);
%rewrite without manual fields
F_G=F_G(mancellsI+1:end,:);
F_R=F_R(mancellsI+1:end,:);
Fneu_G=Fneu_G(mancellsI+1:end,:);
Fneu_R=Fneu_R(mancellsI+1:end,:);
iscell=iscell(mancellsI+1:end,:);
redcell=redcell(mancellsI+1:end,:);
stat=stat(mancellsI+1:end);

%% Select Only Selected ROIs

F1G=nan(sum(iscell(:,1)),size(F_G,2)); %generate matrix to fill with only cell ROIs F
F1R=nan(sum(iscell(:,1)),size(F_R,2));
Fneu1G=nan(sum(iscell(:,1)),size(Fneu_G,2)); %generate matrix to fill with only cell ROIs neu
Fneu1R=nan(sum(iscell(:,1)),size(Fneu_R,2));
xCoor=nan(sum(iscell(:,1)),1); %generate matrix to fill with x coordinate values for only cells green
yCoor=nan(sum(iscell(:,1)),1); %generate matrix to fill with y coordinate values for only cells green
iscellI=nan(size(iscell)); %generate matrix for a key of the iscell ROIs' original ROI ID numbers

jj=1;%counter total will equal number of true cells
for ii=1:size(F_R,1)%number of ROIs (both true and not true cells)
    if iscell(ii,1)==1
        F1R(jj,:)=F_R(ii,:);
        F1G(jj,:)=F_G(ii,:); %populate new cell only matrix with traces
        Fneu1R(jj,:)=Fneu_R(ii,:);
        Fneu1G(jj,:)=Fneu_G(ii,:);
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


%% Normalization of signal and removal of contaminant 

% Baseline shift
Fb1G=F1G-min(min(Fneu1G));%shifts F trace by the global Fneu minimum (assumed as background)
Fb1R=F1R-min(min(Fneu1R));
Fb1manG=F1manG-min(min(Fneu1manG));
Fb1manR=F1manR-min(min(Fneu1manR));

% Normalize traces
F1Gmean_N=(mean(Fb1G,1)-min(mean(Fb1G,1)))./(max(mean(Fb1G,1))-min(mean(Fb1G,1)));
F1Rmean_N=(mean(Fb1R,1)-min(mean(Fb1R,1)))./(max(mean(Fb1R,1))-min(mean(Fb1R,1)));
F1manGmean_N=(mean(Fb1manG,1)-min(mean(Fb1manG,1)))./(max(mean(Fb1manG,1))-min(mean(Fb1manG,1)));
F1manRmean_N=(mean(Fb1manR,1)-min(mean(Fb1manR,1)))./(max(mean(Fb1manR,1))-min(mean(Fb1manR,1)));

%filter out the puncta gfp contaminant
F1Gmean_N=F1Gmean_N.*(.9*-F1manGmean_N+1);%chose 90% of values as it is not 100% contaminant signal
F1Gmean_N=F1Gmean_N./max(F1Gmean_N);
% F1Rmean_N=F1Rmean_N.*(.7*-F1manRmean_N+1);
% F1Rmean_N=F1Rmean_N./max(F1Rmean_N);

% [~,GmaxI]=max(F1G,[],2);
% [~,RmaxI]=max(F1R,[],2);
% incCell=and(GmaxI<=23,GmaxI>=11);
% tempG=mean(F1G(incCell,:),1);
% tempR=mean(F1R(incCell,:),1);
% F1Gmean_N=(tempG-min(tempG))./max(tempG-min(tempG));
% F1Rmean_N=(tempR-min(tempR))./max(tempR-min(tempR));

% F1G_N=(F1G-min(F1G,[],2))./max(F1G-min(F1G,[],2));
% F1R_N=(F1R-min(F1R,[],2))./max(F1R-min(F1R,[],2));

%Individual traces
F1G_N=(Fb1G)./max(Fb1G);
F1R_N=(Fb1R)./max(Fb1R);

% [Gmax,GmaxI]=max(F1G_N,[],2);
% [Rmax,RmaxI]=max(F1R_N,[],2);
% % incCell=and(GmaxI<=23,GmaxI>15) & Gmax>.5;
% incCell=Gmax>.9;
% tempG=mean(F1G(incCell,:),1);
% tempR=mean(F1R(incCell,:),1);
% F1Gmean_N=(tempG-min(tempG))./max(tempG-min(tempG));
% F1Rmean_N=(tempR-min(tempR))./max(tempR-min(tempR));


% spline smoothing (does a cubic interpolation between points
F1Gmean_Nsp=spline([800:10:1250],F1Gmean_N,800:1250);%choose 10-fold increase in sampling frequency
F1Rmean_Nsp=spline([800:10:1250],F1Rmean_N,800:1250);
F1G_Nsp=spline([800:10:1250],F1G_N,800:1250);%chose 10-fold increase in sampling frequency
F1R_Nsp=spline([800:10:1250],F1R_N,800:1250);

%% figure of the plots
figure
rectangle('Position',[1000,0,50,1.1],'FaceColor',[0.9 0.9 0.9],'linestyle','none')
hold on
plot([800:1250],F1Gmean_Nsp, 'color', [1 0.75 0])
plot([800:1250],F1Rmean_Nsp, 'color', [0.6 0 0])
%xlim([800 1250])
ylim([0 1.05])
title('XCaMP-Y and RCatchER 2P Excitation Spectrum')
xlabel('Wavelength (nm)')
ylabel('Normalized F (AU)')
pbaspect([2 1 1])
set(gca, 'box', 'off')
% print('untitled.png', '-dpng', '-r400');
% A = imread('untitled.png');
% imwrite (A, 'XYRCERSpectra.png','png', 'transparency',[1 1 1]);
% clear A
print('../../../RCatchER_CSD_Paper/Figures/Figure1/Figure1_Spectra.svg','-dsvg')

% %% figure of the plots
% cc=0;
% figure
% hold on
% for ii=1:size(F1G,1)
%     if incCell(ii)
%         plot([800:1250],F1G_Nsp(ii,:)+cc, 'color', [.5 0.5 0])
%     else
%         plot([800:1250],F1G_Nsp(ii,:)+cc, 'color', [1 1 0])
%     end
% cc=cc+1;
% end
% xlim([800 1250])
% hold off

