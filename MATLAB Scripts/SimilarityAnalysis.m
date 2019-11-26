function CDATA = SimilarityAnalysis(PARAMS);
%Takes input of PARAMS after it has been set by SimilarityAnalysis_Test.m
%Outputs to CDATA variable
% struct including stimfilespec, datafilespec, binsize
%Computes R NDP SF, computes similarity between input spiketrains and output spiketarins

stimfilespec    = PARAMS.stimfilespec;
datafilespec    = PARAMS.datafilespec;
dbin            = PARAMS.binsize;

%Load stim file and data file
%Deviates from old fucntion by using excel input
[stimtime, stimgroup, S] = parse_axograph(stimfilespec, 0);
[datatime, datagroup, S] = parse_axograph(datafilespec, 0);
clear S
time = datatime;

%Analysis
bins = [0:dbin:time(end)]; %creates bin based on data timed passed through parse_excel
signaltype  = 'raw';
threshtype  = 'direct';
thresh      = -20e-3;
peakflag    = 0;
displayflag = 0;

%stim raster
data = cat(2, stimgroup{:} ); %concatenates arrays horizontally (dimension = 2)
nStims = size(data, 2); %returns vector of # data columns
for n = 1:nStims
    stimtimes{n} = stimtime(data(:,n) > 0);
end
%Converts spiketimes to rasters using spikecoutns2matrix function
StimRaster = spikecounts2matrix(bins, stimtimes);
clear data

%data spike times
data = datagroup{:};
[spiketimes, spikelocs, peaktimes, peaklocs] = detectspikes(time, data, signaltype, threshtype, thresh, peakflag, displayflag);

nsweeps = length(spiketimes);
nrepeats = nsweeps./nStims;

for ns = 1:nsweeps
    Nout(ns, 1) = length(spiketimes{ns}); %Matrix Nout of row ns, col 1 = Length of spiketimes{ns}
end
FiringRateSweep = Nout./2; %Vector column from Nout; Divides NOut by 2 (duration of sweep)
FiringRate = mean(FiringRateSweep);

clear n
for n = 1:nStims
    groupedtimes{n} = spiketimes(n:nStims:end);
    GroupRaster{n} = spikecounts2matrix(bins, groupedtimes{n});
end
spiketimesGP = cat(2, groupedtimes{:});
RespRaster = cat(2, GroupRaster{:});

for i = 1:size(StimRaster,2) % for each input spiketrain  
[BcapIn(i), BinfIn(i)] = burstinessbin(StimRaster(:,i));     % burstinessbin
FRi(i) = sum(StimRaster(:,i))./2; %FR of a given sweep, in Hz
XBcapIn{i} = ones(nrepeats, 1).*BcapIn(i); % to have as many values as outputs
XBinfIn{i} = ones(nrepeats, 1).*BinfIn(i);
XFRi{i} = ones(nrepeats, 1).*FRi(i); 
        for j = 1:size(StimRaster,2) % for each input spiketrain
             [NDPmatIn(i,j), SFmatIn(i,j)] = NDP_SF( StimRaster(:,i), StimRaster(:,j)); % compute NDP and SF
             Rdummy = corrcoef(StimRaster(:,i), StimRaster(:,j)); % compute Pearson's coef of correlation
             RmatIn(i,j) = Rdummy(1,2); clear Rdummy;     
             
             [BcapIn2(j), BinfIn2(j)] = burstinessbin(StimRaster(:,j));
             BurstMetricCapI(i,j) = abs(BcapIn(i)-BcapIn2(j)); % a burst metric that compares pairs of spiketrains
             BurstMetricInfI(i,j) = abs(BinfIn(i)-BinfIn2(j)); %
             
             FRi2(j) = sum(StimRaster(:,j))./2;
             FRimetric(i,j) = abs(FRi(i)-FRi2(j));
        end  
    
end
XBcapInAll = cat(1,XBcapIn{:});
XBinfInAll = cat(1,XBinfIn{:});
XFRiAll = cat(1,XFRi{:});

GMDcenterIc = (max(BcapIn)+min(BcapIn))./2;
GMDcenterIi = (max(BinfIn)+min(BinfIn))./2;
GMDcenterIfr = (max(FRi)+min(FRi))./2;

[ XcosIn, ndp_meanIn, ndp_semIn ] = SymMat2List(NDPmatIn);
[ XsfIn, sf_meanIn, sf_semIn ] = SymMat2List(SFmatIn);
[ XrIn, r_meanIn, r_semIn ] = SymMat2List(RmatIn);
[ XbcIn, bc_meanIn, bc_semIn ] = SymMat2List(BurstMetricCapI);
[ XbiIn, bi_meanIn, bi_semIn ] = SymMat2List(BurstMetricInfI);
[ XfrIn, fr_meanIn, fr_semIn ] = SymMat2List(FRimetric);

% Burstiness and FRmetric per input set: mean, SD and mean absolute difference (aka Gini mean difference)
BcapInSD = std(BcapIn); BcapInMean = mean(BcapIn); 
BinfInSD = std(BinfIn); BinfInMean = mean(BinfIn);
FRiSD = std(FRi); FRiMean = mean(FRi);

% Output similarity (NDP, SF, R, STTC + raw Affinitiy): similarity matrices
clear i j
% [STTCmatOut, AffinOut, Pa_Out, Pb_Out, TilesProp_Out] = STTC( time(end), dbin./2, spiketimesGP );
for i = 1:size(RespRaster,2) % for each output spiketrain  
     % burstinessbin
    [BcapOut(i), BinfOut(i)] = burstinessbin(RespRaster(:,i));
    FRo(i) = sum(RespRaster(:,i))./2; %FR of a given sweep, in Hz
    
        for j = 1:size(RespRaster,2) % for each output spiketrain
             [NDPmatOut(i,j), SFmatOut(i,j)] = NDP_SF(RespRaster(:,i), RespRaster(:,j)); % compute NDP and SF
             Rdummy = corrcoef(RespRaster(:,i), RespRaster(:,j)); % compute Pearson's coef of correlation
             RmatOut(i,j) = Rdummy(1,2); clear Rdummy; 
             
             [BcapOut2(j), BinfOut2(j)] = burstinessbin(RespRaster(:,j));
             BurstMetricCapO(i,j) = abs(BcapOut(i)-BcapOut2(j)); % a burst metric that compares pairs of spiketrains
             BurstMetricInfO(i,j) = abs(BinfOut(i)-BinfOut2(j)); %
             
             FRo2(j) = sum(RespRaster(:,j))./2;
             FRometric(i,j) = abs(FRo(i)-FRo2(j));
        end
end
GMDcenterOc = (max(BcapOut)+min(BcapOut))./2;
GMDcenterOi = (max(BinfOut)+min(BinfOut))./2;
GMDcenterOfr = (max(FRo)+min(FRo))./2;

% Output similarity : get Similarity Within (i.e. the spiketrain reliability Sw) and Between (i.e. Soutput)
[Rw, Rb] = SimMat_WandB(nStims, nrepeats, RmatOut);
[SFw, SFb] = SimMat_WandB(nStims, nrepeats, SFmatOut);
[NDPw, NDPb] = SimMat_WandB(nStims, nrepeats, NDPmatOut);
[BuCw, BuCb] = SimMat_WandB(nStims, nrepeats, BurstMetricCapO);
[BuIw, BuIb] = SimMat_WandB(nStims, nrepeats, BurstMetricInfO);
[FRow, FRob] = SimMat_WandB(nStims, nrepeats, FRometric);

BcapOutSD = std(BcapOut, 'omitnan'); BcapOutMean = mean(BcapOut, 'omitnan'); 
BinfOutSD = std(BinfOut, 'omitnan'); BinfOutMean = mean(BinfOut, 'omitnan');
FRoSD = std(FRo); FRoMean = mean(FRo);

%% display

figure
 
subplot(1,2,1)
scatter(0.9 + 0.2.*rand(size(BcapIn)), BcapIn, 50, 'ok'); hold on
scatter(1.9 + 0.2.*rand(size(BcapOut)), BcapOut, 50, 'o','MarkerFaceColor', [0,0.8, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
errorbar([0.8, 1.8], [BcapInMean BcapOutMean], [BcapInSD BcapOutSD], 'x', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1.5) ; hold on 
errorbar([1.2, 2.2], [GMDcenterIc GMDcenterOc], [bc_meanIn BuCb.ListMean], '+', 'color', 'r', 'markerfacecolor', 'r', 'LineWidth', 1.5) ; hold on
legend('mean +/- SD','mean(min,max) +/- ~GMD', 'Location', 'Best');
plot( [0 10], [0 0], 'k')
plot( [0 2.5], [1 1], 'k--');
      axis square
      box off
      ylabel('Compactness');
      xlim([0.5 2.5]);
      ylim([0 2]);
      set(gca,'ysc', 'lin', 'XTick', [1:2], 'XTickLabel', {'Inputs', 'Outputs'}) ;% xtickangle(45); , 'YTick', [1:0.2:2]
      title('nb potential bins / nb of occupied bins')        
      
 subplot(1,2,2)
  scatter(0.9 + 0.2.*rand(size(FRi)),FRi, 50, 'ok'); hold on
  scatter(1.9 + 0.2.*rand(size(FRo)), FRo, 50, 'o','MarkerFaceColor', [0,0.8, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
  errorbar([0.8, 1.8], [FRiMean FRoMean], [FRiSD FRoSD], 'x', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1.5) ; hold on 
  errorbar([1.2, 2.2], [GMDcenterIfr GMDcenterOfr], [fr_meanIn FRob.ListMean], '+', 'color', 'r', 'markerfacecolor', 'r', 'LineWidth', 1.5) ; hold on
 legend('mean +/- SD','mean(min,max) +/- ~GMD', 'Location', 'Best');
         axis square
         box off
        ylabel('FR, Hz');
       xlim([0.5 2.5]);
        ylim([0 32]);
        set(gca,'ysc', 'lin', 'XTick', [1:2], 'XTickLabel', {'Inputs', 'Outputs'}, 'YTick', [0:5:30]) ;% xtickangle(45); ,          title('FR')  

%CData Struct Output
CDATA.binsize = dbin;

CDATA.StimRaster = StimRaster;
CDATA.RespRaster = RespRaster;

CDATA.FiringRateSweep = FiringRateSweep;
CDATA.FiringRate = FiringRate;

CDATA.XcosIn = XcosIn;
CDATA.ndp_meanIn = ndp_meanIn;
CDATA.ndp_semIn =  ndp_semIn;

CDATA.XsfIn = XsfIn;
CDATA.sf_meanIn = sf_meanIn;
CDATA.sf_semIn = sf_semIn;

CDATA.XrIn = XrIn;
CDATA.r_meanIn = r_meanIn;
CDATA.r_semIn = r_semIn;


CDATA.Rw = Rw;
CDATA.Rb = Rb;

CDATA.SFw = SFw;
CDATA.SFb = SFb;

CDATA.NDPw = NDPw;
CDATA.NDPb = NDPb;

CDATA.BcapIn = BcapIn;
CDATA.XBcapInAll = XBcapInAll; % to have as many values as in BcapOut
CDATA.BcapInMean = BcapInMean;
CDATA.BcapInSD = BcapInSD;
CDATA.BcapOut = BcapOut; 
CDATA.BcapOutMean = BcapOutMean;
CDATA.BcapOutSD = BcapOutSD;
CDATA.BurstMetricCapI = BurstMetricCapI;
CDATA.XbcIn = XbcIn;
CDATA.bc_meanIn = bc_meanIn;
CDATA.bc_semIn = bc_semIn;
CDATA.BurstMetricCapO = BurstMetricCapO;
CDATA.BuCw = BuCw;
CDATA.BuCb = BuCb; 

CDATA.BinfIn = BinfIn; 
CDATA.XBinfInAll = XBinfInAll;
CDATA.BinfInMean = BinfInMean;
CDATA.BinfInSD = BinfInSD;
CDATA.BinfOut = BinfOut; 
CDATA.BinfOutMean = BinfOutMean;
CDATA.BinfOutSD = BinfOutSD;
CDATA.BurstMetricInfI = BurstMetricInfI;
CDATA.XbiIn = XbiIn;
CDATA.bi_meanIn = bi_meanIn;
CDATA.bi_semIn = bi_semIn;
CDATA.BurstMetricInfO = BurstMetricInfO;
CDATA.BuIw = BuIw;
CDATA.BuIb = BuIb;

CDATA.Fri = FRi; 
CDATA.XFRiAll = XFRiAll;
CDATA.FRi = FRiMean;
CDATA.FRiSD = FRiSD;
CDATA.FRo = FRo; 
CDATA.FRoMean = FRoMean;
CDATA.FRoSD = FRoSD;
CDATA.FRimetric = FRimetric;
CDATA.XfrIn = XfrIn;
CDATA.fr_meanIn = fr_meanIn;
CDATA.fr_semIn = fr_semIn;
CDATA.FRometric = FRometric;
CDATA.FRow = FRow;
CDATA.Frob = FRob;

figure