clear
close all

% ***** user defined parameters to change depending on the dataset *****

Nmet = 3; % number of metrics to visualize
Npol = 1; % degree of the polynomial to fit to data pat sep distributions. 1 for linear regression, 2 for parabolic regression

% Xbinsize = [10, 20, 100, 250]; % in msec
%Xbinsize = [10, 50, 100, 250]; % in msec -> GZN
% Xbinsize = [10, 100, 500, 2000];
 Xbinsize = [5, 10, 20, 40, 50, 60, 80, 100, 250, 500,1000, 2000];

% XbinsizeCell = {'10', '20', '100', '250'};
XbinsizeCell = {'5', '10', '20', '40', '50', '60', '80', '100', '250', '500','1000', '2000'};
% XbinsizeCell = {'10', '100', '500', '2000'};

% **********************************************************************

% ***** Datasets to select *****
resultspath = 'C:\Users\vince\lab\pierautlab\OriginalTest\result\';

%CTRL datasets
valid = 1:4; %1:4; % the catnum (i.e. input set) we want to keep for plots
% load('/Users/antoine/Google Drive/LabReasearch/Projects/FalconHawk/SFrange/FRvariation_prt/Sa/SimilarityAnalysis_GC_SFrangeFR_BcapIsCompactness_10Hz_5to2000ms.mat');
% load('/Users/antoine/Google Drive/LabReasearch/Projects/FalconHawk/SFrange/BurstVariation_prt/Sa/SimilarityAnalysis_GC_SFrangeBurst_BcapIsCompactness_10Hz_5to2000ms.mat');
% load('/Users/antoine/Google Drive/LabReasearch/Projects/FalconHawk/Poisson10Hzoriginal/PatSep/Pairwise/Sa/SimilarityAnalysis_Sa15b_BcapIsCompactness_10Hz_10_20_100_250ms.mat');
load('C:\Users\vince\lab\pierautlab\OriginalTest\data\SimilarityAnalysis_SH.mat');
Collate{1} = Cdata; clear Cdata
% Datasets to compare to control
valid = 1:4; %1:4;
% load('/Users/antoine/Google Drive/LabReasearch/Projects/FalconHawk/SFrange/FRvariation_prt/KA/SimilarityAnalysis_KA15b_SFrangePoisson_BcapIsCompactness_5to2000ms.mat');
% load('/Users/antoine/Google Drive/LabReasearch/Projects/FalconHawk/SFrange/FRvariation_prt/KA/SimilarityAnalysis_KA_SFrangeBurst_BcapIsCompactness_5to2000ms.mat');
% load('/Users/antoine/Google Drive/LabReasearch/Projects/FalconHawk/Poisson10Hzoriginal/PatSep/Pairwise/KA/SimilarityAnalysis_KA15b_BcapIsCompactness_10Hz_10_20_100_250ms.mat');
load('C:\Users\vince\lab\pierautlab\OriginalTest\data\SimilarityAnalysis_EE.mat');
Collate{2} = Cdata; clear Cdata

Bend = length(Collate{1}); % number of binsize. Should be equal to length(Xbinsize)
Cend = length(Collate); % number of comparison groups (the stats in this script support only 2 for now)
   
for b = 1:Bend % for all binsizes in Ctrl AND Kainate datasets (Collate.Ctrl and Collate.KA should be the same length)
        
   for c = 1:Cend % for Ctrl and KA groups, or whatever the comparison groups are. 
        
       for n = 1:size(Collate{c}{b},1) % n correspondonds to catnum, i.e. for each stim protocol ("valid" will select later the ones we want to keep for plots)
            
            x1 = cat( 2, Collate{c}{b}{n,:} ); %regroup all recordings from same stim family (contained in CELL)in a single struct. Each params will now have multiple values distributed along a row
            % inputs
            x2.NDP{n} = cat(2, x1.XcosIn); % x1.S are column-vectors of pairwise Sinput (for each pair of input spiketrains). So x2{n} is a matrix N x M (N = number of stim traces, M = number of recordings).
            x2.SF{n} = cat(2, x1.XsfIn);
            x2.R{n} = cat(2, x1.XrIn);
%             x2.STTC{n} = cat(2, x1.XsttcIn);
%             x2.rawSTTC{n} = cat(2, x1.XrawsttcIn);

            Stim.NDP{n} = x2.NDP{n}(:,1); % column-vector of all the pairwise Sinput values for each pair of input spiketrains in a given stim protocol. because all stims are the same for each recording, I just took a column from x2{n}
            Stim.SF{n} = x2.SF{n}(:,1);
            Stim.R{n} = x2.R{n}(:,1);
%             Stim.STTC{n} = x2.STTC{n}(:,1);
%             Stim.rawSTTC{n} = x2.rawSTTC{n}(:,1);

            StimAllPoints_NDP{n} = reshape(x2.NDP{n}, length(Stim.NDP{n})*length(x1),1); % column vector containing Sinput values for all output pairwise comparisons. To plot against "RespGPallPoints"   
            StimAllPoints_SF{n} = reshape(x2.SF{n}, length(Stim.SF{n})*length(x1),1);
            StimAllPoints_R{n} = reshape(x2.R{n}, length(Stim.R{n})*length(x1),1);
%             StimAllPoints_STTC{n} = reshape(x2.STTC{n}, length(Stim.STTC{n})*length(x1),1);
%             StimAllPoints_rawSTTC{n} = reshape(x2.rawSTTC{n}, length(Stim.rawSTTC{n})*length(x1),1);

            % outputs 
            NDPout = cat(2, x1.NDPb); % Cell of size N x M = (# of input pairs)x(# of recordings), each element containing a Struct with the mean pairwise Soutput values
            SFout = cat(2, x1.SFb);
            Rout = cat(2, x1.Rb);
%             STTCout = cat(2, x1.STTCb);
%             rawSTTCout = cat(2, x1.rawSTTCb);

            NDPout_GPmat{n} = cat(2, NDPout.ListGPmean); % matrix of size N x M = (# of input pairs)x(# of recordings), containing the mean pairwise Soutput values
            SFout_GPmat{n} = cat(2, SFout.ListGPmean);
            Rout_GPmat{n} = cat(2, Rout.ListGPmean);
%             STTCout_GPmat{n} = cat(2, STTCout.ListGPmean);
%             rawSTTCout_GPmat{n} = cat(2, rawSTTCout.ListGPmean);

            NDPout_GPmean{n} = mean(NDPout_GPmat{n},2); % cell containing column vectors of Soutput for each group, averaged across all M recordings
            SFout_GPmean{n} = mean(SFout_GPmat{n},2);
            Rout_GPmean{n} = mean(Rout_GPmat{n},2);
%             STTCout_GPmean{n} = mean(STTCout_GPmat{n},2);
%             rawSTTCout_GPmean{n} = mean(rawSTTCout_GPmat{n},2);

            NDPout_GPsem{n} = std( NDPout_GPmat{n},0,2 ) ./ sqrt(size(NDPout_GPmat{n},2)); %columns of SEM for each pair of stim, across recordings 
            SFout_GPsem{n} = std( SFout_GPmat{n},0,2 ) ./ sqrt(size(SFout_GPmat{n},2)); 
            Rout_GPsem{n} = std( Rout_GPmat{n},0,2 ) ./ sqrt(size(Rout_GPmat{n},2)); 
%             STTCout_GPsem{n} = std( STTCout_GPmat{n},0,2 ) ./ sqrt(size(STTCout_GPmat{n},2));
%             rawSTTCout_GPsem{n} = std( rawSTTCout_GPmat{n},0,2 ) ./ sqrt(size(rawSTTCout_GPmat{n},2));

            RespGPallPoints_NDP{n} = reshape(NDPout_GPmat{n}, length(Stim.NDP{n})*length(x1),1); % column vector containing pairwise Soutput values for all recordings. To plot against "StimallPoints"
            RespGPallPoints_SF{n} = reshape(SFout_GPmat{n}, length(Stim.SF{n})*length(x1),1);
            RespGPallPoints_R{n} = reshape(Rout_GPmat{n}, length(Stim.R{n})*length(x1),1);
%             RespGPallPoints_STTC{n} = reshape(STTCout_GPmat{n}, length(Stim.STTC{n})*length(x1),1);
%             RespGPallPoints_rawSTTC{n} = reshape(rawSTTCout_GPmat{n}, length(Stim.rawSTTC{n})*length(x1),1);
            
            clear NDPout SFout Rout STTCout rawSTTCout x1 x2 
       end
    
        StimPair_NDP{b}{c} = cat(1,Stim.NDP{valid}); %column-vector with the pairwise Rinput for all stim protocols, to plot with OutPairMean and OutPairSEM
        StimPair_SF{b}{c} = cat(1,Stim.SF{valid});
        StimPair_R{b}{c} = cat(1,Stim.R{valid});
%         StimPair_STTC{b}{c} = cat(1,Stim.STTC{valid});
%         StimPair_rawSTTC{b}{c} = cat(1,Stim.rawSTTC{valid});

        OutPairMean_NDP{b}{c} = cat(1,NDPout_GPmean{valid}); %column with all the average pairwise Rinput for all stim families
        OutPairMean_SF{b}{c} = cat(1,SFout_GPmean{valid}); 
        OutPairMean_R{b}{c} = cat(1,Rout_GPmean{valid}); 
%         OutPairMean_STTC{b}{c} = cat(1,STTCout_GPmean{valid});
%         OutPairMean_rawSTTC{b}{c} = cat(1,rawSTTCout_GPmean{valid});

        OutPairSEM_NDP{b}{c} = cat(1,NDPout_GPsem{valid}); % idem with SEM
        OutPairSEM_SF{b}{c} = cat(1,SFout_GPsem{valid});
        OutPairSEM_R{b}{c} = cat(1,Rout_GPsem{valid});
%         OutPairSEM_STTC{b}{c} = cat(1,STTCout_GPsem{valid});
%         OutPairSEM_rawSTTC{b}{c} = cat(1,rawSTTCout_GPsem{valid});

        StimPairAll_NDP{b}{c} = cat(1, StimAllPoints_NDP{valid});
        StimPairAll_SF{b}{c} = cat(1, StimAllPoints_SF{valid});
        StimPairAll_R{b}{c} = cat(1, StimAllPoints_R{valid});
%         StimPairAll_STTC{b}{c} = cat(1, StimAllPoints_STTC{valid});
%         StimPairAll_rawSTTC{b}{c} = cat(1, StimAllPoints_rawSTTC{valid});

        OutPairAll_NDP{b}{c} = cat(1, RespGPallPoints_NDP{valid});
        OutPairAll_SF{b}{c} = cat(1, RespGPallPoints_SF{valid});
        OutPairAll_R{b}{c} = cat(1, RespGPallPoints_R{valid});
%         OutPairAll_STTC{b}{c} = cat(1, RespGPallPoints_STTC{valid});
%         OutPairAll_rawSTTC{b}{c} = cat(1, RespGPallPoints_rawSTTC{valid});
              
        % compute effective pat sep (Sin - Sout)         
        PATSEP.R{b}{c} = StimPairAll_R{b}{c} - OutPairAll_R{b}{c}; 
        PATSEP.NDP{b}{c} = StimPairAll_NDP{b}{c} - OutPairAll_NDP{b}{c}; 
        PATSEP.SF{b}{c} = StimPairAll_SF{b}{c} - OutPairAll_SF{b}{c}; 
%         PATSEP.rawSTTC{b}{c} = StimPairAll_rawSTTC{b}{c} - OutPairAll_rawSTTC{b}{c}; 
        
        %normalized pat sep ( effective PatSep ./ Sin )
        nPATSEP.R{b}{c} = PATSEP.R{b}{c} ./StimPairAll_R{b}{c};
        nPATSEP.NDP{b}{c} = PATSEP.NDP{b}{c} ./StimPairAll_NDP{b}{c};
        nPATSEP.SF{b}{c} = PATSEP.SF{b}{c} ./StimPairAll_SF{b}{c};
%         nPATSEP.rawSTTC{b}{c} = PATSEP.rawSTTC{b}{c} ./StimPairAll_rawSTTC{b}{c};
       
        clear n Stim NDPout_GPmean SFout_GPmean Rout_GPmean STTCout_GPmean rawSTTCout_GPmean NDPout_GPsem SFout_GPsem Rout_GPsem STTCout_GPsem rawSTTCout_GPsem StimAllPoints_NDP StimAllPoints_SF StimAllPoints_R StimAllPoints_STTC StimAllPoints_rawSTTC RespGPallPoints_NDP RespGPallPoints_SF RespGPallPoints_R RespGPallPoints_STTC RespGPallPoints_rawSTTC
      end  % end of looping through c = datasets to compare to each other (e.g. Ctrl vs KA)     
    
    %% STATs
    
    EffectSize.NDP{b} = OutPairMean_NDP{b}{1} - OutPairMean_NDP{b}{2};
    EffectSize.SF{b} = OutPairMean_SF{b}{1} - OutPairMean_SF{b}{2};
    EffectSize.R{b} = OutPairMean_R{b}{1} - OutPairMean_R{b}{2};
%     EffectSize.rawSTTC{b} = OutPairMean_rawSTTC{b}{1} - OutPairMean_rawSTTC{b}{2};
    
    %polynomial Fitting & Ftest. Ctrl: c = 1 ; Other: c = 2
    [EST_NDP{b}, pF_NDP(b), Fstat_NDP(b)] = PolFitandFtest (Npol, StimPairAll_NDP{b}{2}, OutPairAll_NDP{b}{2}, StimPairAll_NDP{b}{1}, OutPairAll_NDP{b}{1}, 0);
    [EST_SF{b}, pF_SF(b), Fstat_SF(b)] = PolFitandFtest (Npol, StimPairAll_SF{b}{2}, OutPairAll_SF{b}{2}, StimPairAll_SF{b}{1}, OutPairAll_SF{b}{1},  0);
    [EST_R{b}, pF_R(b), Fstat_R(b)] = PolFitandFtest (Npol, StimPairAll_R{b}{2}, OutPairAll_R{b}{2}, StimPairAll_R{b}{1}, OutPairAll_R{b}{1}, 0);
%     [EST_STTC{b}, pF_STTC(b), Fstat_STTC(b)] = PolFitandFtest (Npol, StimPairAll_STTC{b}{2}, OutPairAll_STTC{b}{2}, StimPairAll_STTC{b}{1}, OutPairAll_STTC{b}{1}, 0);
%     [EST_rawSTTC{b}, pF_rawSTTC(b), Fstat_rawSTTC(b)] = PolFitandFtest (Npol, StimPairAll_rawSTTC{b}{2}, OutPairAll_rawSTTC{b}{2}, StimPairAll_rawSTTC{b}{1}, OutPairAll_rawSTTC{b}{1}, 0);
        
    clear DataGroup1Y_NDP DataGroup1Y_SF DataGroup1Y_R G_NDP G_SF G_R DataGroup1X_NDP DataGroup1X_SF DataGroup1X_R TreatmentGroup TreatmentGroup2_Ctrl TreatmentGroup2_KA
end % end of looping through b = binsizes       
        
        
%% display collate

% colormaps
% cmap = hot(Bend + 3);
% SaMap = gray(Bend+2);
% KAmap = autumn(Bend+1);

SaMap = brewermap(Bend+1, '*Greys');
KAmap = brewermap(Bend+1, '*Reds');
cmap = KAmap;

figure(1) % effect size
 clear p
    p = panel(); % panel toolbox (more flexible than subplot). See http://www.mathworks.com/matlabcentral/fileexchange/20003-panel
	p.pack(3, 2); % creates a grid of descendant panels (i.e. subplots)
%     set margins of the descendants panels (i.e. the internal margins) 
    p.de.margin = 17;
%     set margins of the root panel (i.e. the edge margins) p to the bones
    p.margin = [12 12 5 5];
    
p(1, 1).select(); % R    
    for b = 1:Bend
    plot(StimPair_R{b}{1}, EffectSize.R{b}, 'o', 'Color', cmap(b,:)); hold on
    end
    hleg = legend(XbinsizeCell, 'Location', 'BestOutside');
    legend BOXOFF
    htitle = get(hleg,'Title');
    set(htitle,'String','timescale (ms)')
    plot([0 1], [0 0], 'k--');
    box off; axis square; 
    xlabel('Rin')
    ylabel('Rout{Sa}-Rout{KA}')
    set(gca, 'xlim', [0 1]) 
    title('Effect Size: R')

p(2, 1).select(); % NDP    
    for b = 1:Bend
    plot(StimPair_NDP{b}{1}, EffectSize.NDP{b}, 'o', 'Color', cmap(b,:)); hold on
    end
    plot([0 1], [0 0], 'k--');
    box off; axis square; 
    xlabel('NDPin')
    ylabel('NDPout{Sa}-NDPout{KA}')
    set(gca, 'xlim', [0 1]) 
    title('Effect Size: NDP')

p(3, 1).select(); % SF    
    for b = 1:Bend
    plot(StimPair_SF{b}{1}, EffectSize.SF{b}, 'o', 'Color', cmap(b,:)); hold on
    end
    plot([0 1], [0 0], 'k--');
    box off; axis square; 
    xlabel('SFin')
    ylabel('SFout{Sa}-SFout{KA}')
    set(gca, 'xlim', [0 1]) 
    title('Effect Size: SF')

% p(4, 1).select(); % SF    
%     for b = 1:Bend
%     plot(StimPair_rawSTTC{b}{1}, EffectSize.rawSTTC{b}, 'o', 'Color', cmap(b,:)); hold on
%     end
%     plot([0 1], [0 0], 'k--');
%     box off; axis square; 
%     xlabel('rawSTTCin')
%     ylabel('out{Sa}-out{KA}')
%     set(gca, 'xlim', [0 1]) 
%     title('Effect Size: raw STTC (Affinity)')
    
figure(2) % cum distrib of effective pattern separation

if Bend == 4
    clear p
        p = panel(); % panel toolbox (more flexible than subplot). See http://www.mathworks.com/matlabcentral/fileexchange/20003-panel
        p.pack(1, 3); % creates a grid of descendant panels (i.e. subplots)
    %     set margins of the descendants panels (i.e. the internal margins) 
        p.de.margin = 17;
    %     set margins of the root panel (i.e. the edge margins) p to the bones
        p.margin = [12 12 5 5];

    p(1, 1).select();
    for b = 1:Bend
    [cdfSa{b}, StatCdfSa_R{b} ]= cdfplot(PATSEP.R{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    [cdfKA{b}, StatCdfKA_R{b} ] = cdfplot(PATSEP.R{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    end
    % hleg2 = legend(XbinsizeCell, 'Location', 'southeast');
    % legend BOXOFF
    % htitle = get(hleg2,'Title');
    % set(htitle,'String','timescale (ms)')
    set(gca, 'xlim', [-0.5, 1]);
    plot([0 0], [0 1], 'k--');
    box off; grid off; axis square; 
    xlabel('Decorrelation')
    ylabel('cum freq') 
    title('R')

    % p(1, 2).select();
    % clear cdfSa cdfKA
    % for b = 1:Bend
    % cdfSa{b} = cdfplot(nPATSEP.R{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    % cdfKA{b} = cdfplot(nPATSEP.R{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    % end
    % plot([0 0], [0 1], 'k--');
    % box off; grid off; axis square; 
    % xlabel('norm Decorr')
    % ylabel('cum freq') 
    % title('R')

    p(1, 2).select();
    clear cdfSa cdfKA
    for b = 1:Bend
    [cdfSa{b}, StatCdfSa_NDP{b} ] = cdfplot(PATSEP.NDP{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    [cdfKA{b}, StatCdfKA_NDP{b} ] = cdfplot(PATSEP.NDP{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    end
    plot([0 0], [0 1], 'k--');
    set(gca, 'xlim', [-0.2, 1]);
    box off; grid off; axis square; 
    xlabel('Orthogonalization')
    ylabel('cum freq') 
    title('NDP')

    % p(2, 2).select();
    % clear cdfSa cdfKA
    % for b = 1:Bend
    % cdfSa{b} = cdfplot(nPATSEP.NDP{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    % cdfKA{b} = cdfplot(nPATSEP.NDP{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    % end
    % plot([0 0], [0 1], 'k--');
    % box off; grid off; axis square; 
    % xlabel('norm Ortho')
    % ylabel('cum freq') 
    % title('NDP')

    p(1, 3).select();
    clear cdfSa cdfKA
    for b = 1:Bend
    [cdfSa{b}, StatCdfSa_SF{b} ] = cdfplot(PATSEP.SF{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    [cdfKA{b}, StatCdfKA_SF{b} ] = cdfplot(PATSEP.SF{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    end
    plot([0 0], [0 1], 'k--');
    set(gca, 'xlim', [-0.1, 0.5]);
    box off; grid off; axis square; 
    xlabel('Scaling')
    ylabel('cum freq') 
    title('SF')

elseif Bend == 12 % when 12 binsizes have been computed

    clear SamMap KAmap
    tau = [2 3 8 9]; %10, 20, 100, 250 ms
    SaMap = brewermap(length(tau)+1, '*Greys');
    KAmap = brewermap(length(tau)+1, '*Reds');
    cmap = KAmap;
    subplot(1,3,1)
    for b = 2
    ok = tau(b);
    [cdfSa{b}, StatCdfSa_R{b} ] = cdfplot(PATSEP.R{ok}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    [cdfKA{b}, StatCdfKA_R{b} ] = cdfplot(PATSEP.R{ok}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    end
    set(gca, 'xlim', [-0.5, 0.5]);
    plot([0 0], [0 1], 'k--');
    box off; grid off; axis square; 
    xlabel('Decorrelation')
    ylabel('cum freq') 
    title(['R, ms = ' num2str(Xbinsize(b)) ' ms'])

    subplot(1,3,2)
    clear cdfSa cdfKA
    for b = 2
    ok = tau(b);
    [cdfSa{b}, StatCdfSa_NDP{b} ] = cdfplot(PATSEP.NDP{ok}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    [cdfKA{b}, StatCdfKA_NDP{b} ] = cdfplot(PATSEP.NDP{ok}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    end
    plot([0 0], [0 1], 'k--');
    set(gca, 'xlim', [-0.5, 0.5]);
    box off; grid off; axis square; 
    xlabel('Orthogonalization')
    ylabel('cum freq') 
    title(['NDP, ms = ' num2str(Xbinsize(b))])

    subplot(1,3,3)
    clear cdfSa cdfKA
    for b = 2
    ok = tau(b);
    [cdfSa{b}, StatCdfSa_SF{b} ] = cdfplot(PATSEP.SF{ok}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
    [cdfKA{b}, StatCdfKA_SF{b} ] = cdfplot(PATSEP.SF{ok}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
    end
    plot([0 0], [0 1], 'k--');
    set(gca, 'xlim', [-0.5, 0.5]);
    box off; grid off; axis square; 
    xlabel('Scaling')
    ylabel('cum freq') 
    title(['SF, ms = ' num2str(Xbinsize(b))])
    
else 
    disp('pb in number of binsizes - check user defined params and datasets')
end

% p(3, 2).select(); % normalized pat sep
% clear cdfSa cdfKA
% for b = 1:Bend
% cdfSa{b} = cdfplot(nPATSEP.SF{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
% cdfKA{b} = cdfplot(nPATSEP.SF{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
% end
% plot([0 0], [0 1], 'k--');
% box off; grid off; axis square; 
% xlabel('norm Scaling')
% ylabel('cum freq') 
% title('SF')

% p(4, 1).select();
% clear cdfSa cdfKA
% for b = 1:Bend
% cdfSa{b} = cdfplot(PATSEP.rawSTTC{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
% cdfKA{b} = cdfplot(PATSEP.rawSTTC{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
% end
% plot([0 0], [0 1], 'k--');
% box off; grid off; axis square; 
% xlabel('Affinity decrease')
% ylabel('cum freq') 
% title('rawSTTC')
% 
% p(4, 2).select();
% clear cdfSa cdfKA
% for b = 1:Bend
% cdfSa{b} = cdfplot(nPATSEP.rawSTTC{b}{1}); cdfSa{b}.Color = SaMap(b,:); cdfSa{b}.LineWidth = 1.5; hold on
% cdfKA{b} = cdfplot(nPATSEP.rawSTTC{b}{2}); cdfKA{b}.Color = KAmap(b,:); cdfKA{b}.LineWidth = 1.5; hold on
% end
% plot([0 0], [0 1], 'k--');
% box off; grid off; axis square; 
% xlabel('norm Affinity decrease')
% ylabel('cum freq') 
% title('raw STTC')

figure(4)

  subplot(1,Nmet,1)
    for b = 1:4
    plot([1 2],[StatCdfSa_R{b}.median StatCdfKA_R{b}.median], 'k-'); hold on;
    scatter([1 2], [StatCdfSa_R{b}.median, StatCdfKA_R{b}.median], 100, [SaMap(b,:); KAmap(b,:)], 'filled' ); hold on;
    end
%     set(gca, 'xlim', [0 3], 'ylim', [-0.1 0]);
    box off; grid off; axis square; 
    xlabel('Treatment')
    ylabel('PatSep median') 
    title('R')
  subplot(1,Nmet,2)
    for b = 1:4
    plot([1 2],[StatCdfSa_NDP{b}.median StatCdfKA_NDP{b}.median], 'k-'); hold on;
    scatter([1 2], [StatCdfSa_NDP{b}.median, StatCdfKA_NDP{b}.median], 100, [SaMap(b,:); KAmap(b,:)], 'filled' ); hold on;
    end
%     set(gca, 'xlim', [0 3], 'ylim', [-0.1 0.1]);
    box off; grid off; axis square; 
    xlabel('Treatment')
    ylabel('PatSep median') 
    title('NDP')
  subplot(1,Nmet,3)
    for b = 1:4
    plot([1 2],[StatCdfSa_SF{b}.median StatCdfKA_SF{b}.median], 'k-'); hold on;
    scatter([1 2], [StatCdfSa_SF{b}.median, StatCdfKA_SF{b}.median], 100, [SaMap(b,:); KAmap(b,:)], 'filled' ); hold on;
    end
%     set(gca, 'xlim', [0 3], 'ylim', [-0.2 0.1]);
    box off; grid off; axis square; 
    xlabel('Treatment')
    ylabel('PatSep median') 
    title('SF')
    
figure(5) % pval of F-test as a function of binning window, for each metric

    subplot(1, Nmet, 1) % R
            plot( Xbinsize, pF_R, 'ok-' ); hold on
            plot([0 Xbinsize(end)], [0.05 0.05], 'g--'); % significance line
            box off
            axis square
            xlabel('binsize, ms')
            ylabel('F-test pval')
            set(gca, 'xlim', [0 Xbinsize(end)]) 
            title('R')

    subplot(1, Nmet, 2) % NDP
            plot( Xbinsize, pF_NDP, 'ok-' ); hold on
            plot([0 Xbinsize(end)], [0.05 0.05], 'g--'); % significance line
            box off
            axis square
            xlabel('binsize, ms')
            ylabel('F-test pval')
            set(gca, 'xlim', [0 Xbinsize(end)]) 
            title('NDP')

    subplot(1, Nmet, 3) % SF
            plot( Xbinsize, pF_SF, 'ok-' ); hold on
            plot([0 Xbinsize(end)], [0.05 0.05], 'g--'); % significance line
            box off
            axis square
            xlabel('binsize, ms')
            ylabel('F-test pval')
            set(gca, 'xlim', [0 Xbinsize(end)])  
            title('SF')
            
%     subplot(1, Nmet, 4) % rawSTTC
%             plot( Xbinsize, pF_rawSTTC, 'ok-' ); hold on
%             plot([0 Xbinsize(end)], [0.05 0.05], 'g--'); % significance line
%             box off
%             axis square
%             xlabel('binsize, ms')
%             ylabel('F-test pval')
%             set(gca, 'xlim', [0 Xbinsize(end)])  
%             title('raw STTC')
            
