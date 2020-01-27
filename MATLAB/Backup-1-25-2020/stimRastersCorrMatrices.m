clear all
close all
clc

resultspath = 'C:\Users\vince\Documents\'
cd(resultspath)


% PatSepFileListAntoine
% PatSepFilesFalconHawk15b_Kainate
% PatSepFiles_CA3_gzine_30Hz
%PatSepFilesFalconHawk_SFrange_Sa
%NameFile_SH
NameFile_EE


% ok = [1,3:12];
% ok = [4,5,7,8,9];
% ok = 1:2;
ok = 1 %1:4;
dbin = 0.010;   % sec
sttime1 = cputime; 

for catnum = 1:length(Files.Stim.File)
clear stimfilespec
stimfilespec = Files.Stim.File{catnum};

% Load Stim & Data files
[stimtime, stimgroup, S] = parse_axograph( stimfilespec, 0) ;
clear S
time = stimtime;

% Analysis Params
bins = [0:dbin:time(end)];

signaltype  = 'raw';
threshtype  = 'direct';
thresh      = -20e-3; 
peakflag    = 0;
displayflag = 0;
colord = lines(5);
   
%stim raster 
data = cat( 2, stimgroup{:} );
for n = 1:size(data, 2)
    stimtimes{n} = stimtime( find(data(:,n)) );
    FiringRate(n,1)=length(stimtimes{n})./2; %number of impulses divided by 2 sec
end
StimRaster = spikecounts2matrix(bins, stimtimes);
MeanFiringRate = mean(FiringRate);

% STIM Corr
StimCorrCoeffMat = corrcoef(StimRaster);  
nanmask = ones(size(StimCorrCoeffMat));
nanmask = triu(nanmask,1)./triu(nanmask,1);
x = nanmask.*triu(StimCorrCoeffMat, 1);
x = x(~isnan(x));
StimCorrMean = mean( x );
StimCorrSEM(catnum,1)     = std( x )./sqrt(length(x));
StimCorrSd(catnum,1) = std( x );

% variability in percent of the mean (akin to CV or Fano Factor)
StimCorrSEMpercent(catnum,1) = (StimCorrSEM(catnum,1)./StimCorrMean).*100;
StimCorrSDpercent(catnum,1) = (StimCorrSd(catnum,1)./StimCorrMean).*100;


StimTimes{catnum} = stimtimes;
StimRasters{catnum} = StimRaster;
StimCorrCoeffMatrix{catnum} = StimCorrCoeffMat;
StimCorrCoefVect(:, catnum) = x;
StimCorrMeanMeasured(catnum,1) = StimCorrMean;
MeanFiringRateAll(catnum,1) = MeanFiringRate;

% NDP and SF analysis
    for ni = 1:size(StimRaster,2);
        for np = 1:size(StimRaster,2);

                StimDotMat(ni, np) = dot(StimRaster(:,ni),StimRaster(:,np)); %Dot product of column vectors A and B = ||A||*||B||*cos(theta), with theta = angle between vect A and B
                StimNormMat(ni, np) = norm(StimRaster(:,ni)).*norm(StimRaster(:,np)); % ||A||x||B|| = magnitude of the dot product
                StimCosMat(ni, np) = StimDotMat(ni, np)./StimNormMat(ni, np); % cos(theta)

                % compute ||A||./||B|| = scaling factor between the 2 vectors. It's a function of how different the firing rate in each bin is.
                % To get a symetric matrix, always compute ||A||/||B|| with ||A|| < ||B||.
                % This way, 0 < scaling factor < 1
                % When A or B = 0 (i.e when there is no spike in one of the rasters), make the ratio = 0 (and then make it NaN and exclude it from the average)

                if norm(StimRaster(:,ni)) && norm(StimRaster(:,np)) %if both A and B have a norm > 0, i.e at least 1 spike in both of the rasters. 
                    if norm(StimRaster(:,ni)) <= norm(StimRaster(:,np)) 
                       StimScalMat(ni,np) = norm(StimRaster(:,ni))./norm(StimRaster(:,np)); 
                    else  
                       StimScalMat(ni,np) = norm(StimRaster(:,np))./norm(StimRaster(:,ni));
                    end
                else
                       StimScalMat(ni,np) = 0;
                end
        end
    end
    
    nanmask = ones(size(StimDotMat));
    nanmask = triu(nanmask,1)./triu(nanmask,1);
    
    xcos = nanmask.*triu(StimCosMat, 1);
    xcos = xcos(~isnan(xcos));
    StimNDPvect(:, catnum) = xcos;
    StimCosMean(catnum,1) = mean( xcos );
    StimCosSEM(catnum,1) = std( xcos )./sqrt(length(xcos));
    StimCosSD(catnum,1) = std( xcos );
    
    StimCosSEMpercent(catnum,1) = (StimCosSEM(catnum,1)./StimCosMean(catnum,1)).*100;
    StimCosSDpercent(catnum,1) = (StimCosSD(catnum,1)./StimCosMean(catnum,1)).*100;
    
    xscal = nanmask.*triu(StimScalMat, 1);
    xscal(xscal==0) = nan; % to remove the values 0, corresponding to an ill-defined scaling factor (because one of the rasters had no spikes)
    xscal = xscal(~isnan(xscal));
    StimSFvect(:, catnum) = xscal;
    StimScalMean(catnum,1) = mean( xscal );
    StimScalSEM(catnum,1)  = std( xscal )./sqrt(length(xscal));
    StimScalSD(catnum,1)   = std( xscal );
    
    StimScalSEMpercent(catnum,1) = (StimScalSEM(catnum,1)./StimScalMean(catnum,1)).*100
    StimScalSDpercent(catnum,1) = (StimScalSD(catnum,1)./StimScalMean(catnum,1)).*100
  
    % pairwise Burstiness and FR metrics
    for i = 1:size(StimRaster,2) % for each input spiketrain  
    [BcapIn(i), BinfIn(i)] = burstinessbin(StimRaster(:,i));     % burstinessbin
    FR(i) = sum(StimRaster(:,i))./2; %FR of a given sweep, in Hz
        for j = 1:size(StimRaster,2) % for each input spiketrain
%              [NDPmatIn(i,j), SFmatIn(i,j)] = NDP_SF( StimRaster(:,i), StimRaster(:,j)); % compute NDP and SF
%              Rdummy = corrcoef(StimRaster(:,i), StimRaster(:,j)); % compute Pearson's coef of correlation
%              RmatIn(i,j) = Rdummy(1,2); clear Rdummy;     
             
             [BcapIn2(j), BinfIn2(j)] = burstinessbin(StimRaster(:,j));
             BurstMetricCapI(i,j) = abs(BcapIn(i)-BcapIn2(j)); % a burst metric that compares pairs of spiketrains
             BurstMetricInfI(i,j) = abs(BinfIn(i)-BinfIn2(j)); %
             
             FR2(j) = sum(StimRaster(:,j))./2;
             FRmetric(i,j) = abs(FR(i)-FR2(j));
        end  
    end
    
    [ XbcIn, bc_meanIn, bc_semIn ] = SymMat2List(BurstMetricCapI);
    [ XbiIn, bi_meanIn, bi_semIn ] = SymMat2List(BurstMetricInfI);
    [ XfrIn, fr_meanIn, fr_semIn ] = SymMat2List(FRmetric);

    % Burstiness and FRmetric per input set: mean, SD and mean absolute difference (aka Gini mean difference)
    BcapInSD = std(BcapIn); BcapInMean = mean(BcapIn); 
    BinfInSD = std(BinfIn); BinfInMean = mean(BinfIn);
    FRsd = std(FR); FRmean = mean(FR);

end

AverageFRstim = mean(MeanFiringRateAll(ok,1));
SEMfrStim = std(MeanFiringRateAll(ok,1))./sqrt(length(MeanFiringRateAll(ok,1)));

%variability as percent of the mean
StimCorrSEMpercentAverage = mean(StimCorrSEMpercent);
StimCorrSDpercentAverage = mean(StimCorrSDpercent);

StimCosSEMpercentAverage = mean(StimCosSEMpercent);
StimCosSDpercentAverage = mean(StimCosSDpercent);

StimScalSEMpercentAverage = mean(StimScalSEMpercent);
StimScalSDpercentAverage = mean(StimScalSDpercent);
%% Display

% N = length(Files.Stim.File);
N = length(ok);
Black=zeros(5,3);

figure(1)

U = N;
for n = 1:N
subplot(2,N,n)
v = ok(n);
%[spikex, spikey] = makedisplayrasters(StimTimes{v}, 0, Black);
%displaydisplayrasters(spikex, spikey, Black); hold on
testTimes{1, 1} = [0.1, .249, .4, 1, 1.294, 1.8924, 2]
[spikex, spikey] = plotSpikeRaster(testTimes)
grid off
axis square
xlabel('Time, s')
title({'R_{input} = ' num2str(StimCorrMeanMeasured(v,1)) ' +/- SD = ' num2str(StimCorrSd(v,1))})
  
U = U+1;
 
subplot(2,N,U)     
dummy = StimScalMat %StimCorrCoeffMatrix{v};
imagesc(dummy); hold on
axis square
colormap(hot) 
caxis([0 1]);
colorbar;
clear dummy
end

% print('SFinput-2.eps', '-depsc', '-painters')
% print('Rinput.pdf', '-dpdf', '-painters')
% print('Rinput.svg', '-dsvg', '-painters')

figure(2)

    subplot(3,1,1)
    boxplot(StimCorrCoefVect(:, ok), 'labels', StimCorrMeanMeasured(ok) ); hold on
    ylabel('pairwise R');
    xlabel('input set R (mean)')
    title('R, 10 ms')
    set(gca, 'ylim', [0 1])
    
    subplot(3,1,2)
    boxplot(StimSFvect(:, ok), 'labels', StimCorrMeanMeasured(ok) ); hold on 
    set(gca, 'ylim', [0 1])
    ylabel('pairwise SF');
    xlabel('input set R (mean)')
    title('SF, 10 ms')
    
    subplot(3,1,3)
    boxplot(StimNDPvect(:, ok), 'labels', StimCorrMeanMeasured(ok) ); hold on 
    set(gca, 'ylim', [0 1])
    ylabel('pairwise NDP');
    xlabel('input set R (mean)')
    title('NDP, 10 ms')
    
    

figure(3)
N = length(Files.Stim.File);
    subplot(1,3,1)
    for n = 1:N
    Xstim{n} = StimCorrMeanMeasured(n).*ones(size(StimCorrCoefVect(:,n)));
    scatter(Xstim{n}, StimCorrCoefVect(:,n), 'k');hold on
    end
    set(gca, 'ylim', [0 1])
    ylabel('pairwise R');
    xlabel('input set R (mean)')
    title('10 ms')

    subplot(1,3,2)
    for n = 1:N
    SFstim{n} = StimScalMean(n).*ones(size(StimSFvect(:,n)));
    scatter(SFstim{n}, StimSFvect(:,n), 'k');hold on
    end
    set(gca, 'ylim', [0 1])
    ylabel('pairwise SF');
    xlabel('input set SF (mean)')
    title('10 ms')
    
    subplot(1,3,3)
    for n = 1:N
    NDPstim{n} = StimCosMean(n).*ones(size(StimNDPvect(:,n)));
    scatter(NDPstim{n}, StimNDPvect(:,n), 'k');hold on
    end
    set(gca, 'ylim', [0 1])
    ylabel('pairwise NDP');
    xlabel('input set NDP (mean)')
    title('10 ms')

    figure(4) % pairwise R, NDP and SF as a function of pairwise FR, Compactness and occupancy 
    
    subplot(3,3,1)
    scatter(XfrIn,StimCorrCoefVect, 'ok')
    box off
    axis square
%     xlabel('FR difference, Hz')
    ylabel('R')
    set(gca, 'ylim', [-1 1]);
    
    subplot(3,3,2)
    scatter(XbcIn,StimCorrCoefVect, 'ok')
    box off
    axis square
%     xlabel('Compactness difference')
%     ylabel('R')
    set(gca, 'ylim', [-1 1]);
%     title([num2str(dbin*1000) ' ms'])
    
    subplot(3,3,3)
    scatter(XbiIn,StimCorrCoefVect, 'ok')
    box off
    axis square
%     xlabel('Spikes/bin difference')
%     ylabel('R')
    set(gca, 'ylim', [-1 1]);
%     title([num2str(dbin*1000) ' ms'])
    
    subplot(3,3,4)
    scatter(XfrIn,StimNDPvect, 'ok')
    box off
    axis square
%     xlabel('FR difference, Hz')
    ylabel('NDP')
    set(gca, 'ylim', [0 1]);
    
    subplot(3,3,5)
    scatter(XbcIn,StimNDPvect, 'ok')
    box off
    axis square
%     xlabel('Compactness difference')
%     ylabel('NDP')
    set(gca, 'ylim', [0 1]);
%     title([num2str(dbin*1000) ' ms'])
    
    subplot(3,3,6)
    scatter(XbiIn,StimNDPvect, 'ok')
    box off
    axis square
%     xlabel('Spikes/bin difference')
%     ylabel('NDP')
    set(gca, 'ylim', [0 1]);
%     title([num2str(dbin*1000) ' ms'])
    
    subplot(3,3,7)
    scatter(XfrIn,StimSFvect, 'ok')
    box off
    axis square
    xlabel('FR difference, Hz')
    ylabel('SF')
    set(gca, 'ylim', [0 1]);
    
    subplot(3,3,8)
    scatter(XbcIn,StimSFvect, 'ok')
    box off
    axis square
    xlabel('Compactness difference')
%     ylabel('SF')
    set(gca, 'ylim', [0 1]);
%     title([num2str(dbin*1000) ' ms'])
    
    subplot(3,3,9)
    scatter(XbiIn,StimSFvect, 'ok')
    box off
    axis square
    xlabel('Spikes/bin difference')
%     ylabel('SF')
    set(gca, 'ylim', [0 1]);
%     title([num2str(dbin*1000) ' ms'])
    
    


























