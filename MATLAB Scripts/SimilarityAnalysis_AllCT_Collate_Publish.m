clear
close all

% PatSep 
load('\Users\vince\lab\pierautlab\OriginalTest\data\SimilarityAnalysis_TestOutput.mat');
Collate{1} = Cdata; % GC
Valid{1} = [1, 3:12]; %[6, 10, 11]; %[4,5,7,8,9] %
clearvars -except Collate Valid


Npol = 1; % degree of the polynomial to fit to data pat sep distributions. 1 for linear regression, 2 for parabolic regression
% Xbinsize = [5, 10, 20, 40, 60, 80, 100, 250, 500,1000, 2000]; % in msec
% Xbinsize = [10, 50, 100, 250];
Xbinsize = [5, 10, 20, 40, 50, 60, 80, 100, 250, 500,1000, 2000]; % in msec
Rinterp = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];
RinterpSF = [1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5];

Bend = length(Collate{1}); % number of binsize. Should be equal to length(Xbinsize)
Cend = length(Collate); 

for b = 1:Bend % for all binsizes in Ctrl AND Kainate datasets (Collate.Ctrl and Collate.KA should be the same length)
        
   for c = 1:Cend %
       valid = Valid{c};
       for n = 1:size(Collate{c}{b},1) % n correspondonds to catnum, i.e. for each stim protocol ("valid" will select later the ones we want to keep for plots)
            x1 = cat( 2, Collate{c}{b}{n,:} ); %regroup all recordings from same stim family (contained in CELL)in a single struct. Each params will now have multiple values distributed along a row inputs
            x2.NDP{n} = cat(2, x1.XcosIn); % x1.S are column-vectors of pairwise Sinput (for each pair of input spiketrains). So x2{n} is a matrix N x M (N = number of stim traces, M = number of recordings).
            x2.SF{n} = cat(2, x1.XsfIn);
            x2.R{n} = cat(2, x1.XrIn);
            x2.BuC{n} = cat(2, x1.XbcIn); % compactness 
            x2.BuO{n} = cat(2, x1.XbiIn); % occupancy
            x2.FR{n} = cat(2, x1.XfrIn);


            Stim.NDP{n} = x2.NDP{n}(:,1); % column-vector of all the pairwise Sinput values for each pair of input spiketrains in a given stim protocol. because all stims are the same for each recording, I just took a column from x2{n}
            Stim.SF{n} = x2.SF{n}(:,1);
            Stim.R{n} = x2.R{n}(:,1);
            Stim.BuC{n} = x2.BuC{n}(:,1);
            Stim.BuO{n} = x2.BuO{n}(:,1);
            Stim.FR{n} = x2.FR{n}(:,1);

            StimAllBcapInSD{n} = cat(1, x1.BcapInSD);
            StimAllBinfInSD{n} = cat(1, x1.BinfInSD);
            Stim.BcapInSD(n,1) = StimAllBcapInSD{n}(1);
            Stim.BinfInSD(n,1) = StimAllBinfInSD{n}(1);
            StimAllBuCgmd{n} = cat(1, x1.bc_meanIn);% compactness
            StimAllBuOgmd{n} = cat(1, x1.bi_meanIn);% occupancy
            StimAllFRgmd{n} = cat(1, x1.fr_meanIn);
            
            Stim.BuCgmd(n,1) = StimAllBuCgmd{n}(1); % compactness
            Stim.BuOgmd(n,1) = StimAllBuOgmd{n}(1); % occupancy
            Stim.FRgmd(n,1) = StimAllFRgmd{n}(1);
            
            StimAllPoints_NDP{n} = reshape(x2.NDP{n}, length(Stim.NDP{n})*length(x1),1); % column vector containing Sinput values for all output pairwise comparisons. To plot against "RespGPallPoints"   
            StimAllPoints_SF{n} = reshape(x2.SF{n}, length(Stim.SF{n})*length(x1),1);
            StimAllPoints_R{n} = reshape(x2.R{n}, length(Stim.R{n})*length(x1),1);
            StimAllPoints_BuC{n} = reshape(x2.BuC{n}, length(Stim.BuC{n})*length(x1),1); %compactness
            StimAllPoints_BuO{n} = reshape(x2.BuO{n}, length(Stim.BuO{n})*length(x1),1); %occupancy
            StimAllPoints_FR{n} = reshape(x2.FR{n}, length(Stim.FR{n})*length(x1),1);

            AllSpkBcapIn{n} = cat(1, x1.XBcapInAll); %compactness of all input spiketrains (column) repeated for all recordings (rows)
            AllSpkBinfIn{n} = cat(1, x1.XBinfInAll); %Binf = occupancy
            AllSpkFRi{n} = cat(1, x1.XFRiAll); % FR of all input spiketrains
            AllSpkBcapOut{n} = cat(2, x1.BcapOut); % compactness
            AllSpkBinfOut{n} = cat(2, x1.BinfOut); %Occupancy
            AllSpkFRo{n} = cat(2, x1.FRo);
            

            % outputs 
            NDPout = cat(2, x1.NDPb); % Cell of size N x M = (# of input pairs)x(# of recordings), each element containing a Struct with the mean pairwise Soutput values
            SFout = cat(2, x1.SFb);
            Rout = cat(2, x1.Rb);
            BuCout = cat(2, x1.BuCb); %compactness
            BuOout = cat(2, x1.BuIb); % occupancy: Binf instead of Bcap
            FRout = cat(2, x1.Frob);            
            
            %non pairwise analysis
            NDPout_Listmean{n} = cat(1, NDPout.ListMean); % vector column of average NDP per recording, of length = number of recordings for input set n
            SFout_Listmean{n} = cat(1, SFout.ListMean);
            Rout_Listmean{n} = cat(1, Rout.ListMean);
            BuCout_Listmean{n} = cat(1, BuCout.ListMean); % vector column of average burstiness dispersion per recording, of length = number of recordings for input set n
            BuOout_Listmean{n} = cat(1, BuOout.ListMean);
            FRout_Listmean{n} = cat(1, FRout.ListMean);
            
            NDPout_mean{n} = mean(NDPout_Listmean{n}); 
            SFout_mean{n} = mean(SFout_Listmean{n});
            Rout_mean{n} = mean(Rout_Listmean{n});
            BuCout_mean{n} = mean(BuCout_Listmean{n});
            BuOout_mean{n} = mean(BuOout_Listmean{n});
            FR_mean{n} = mean(FRout_Listmean{n});
            
            NDPout_sem{n} = std(NDPout_Listmean{n}) ./ sqrt(size(NDPout_Listmean{n},2));
            SFout_sem{n} = std(SFout_Listmean{n})./ sqrt(size(SFout_Listmean{n},2));
            Rout_sem{n} = std(Rout_Listmean{n})./ sqrt(size(Rout_Listmean{n},2));
            BuCout_sem{n} = std(BuCout_Listmean{n})./ sqrt(size(BuCout_Listmean{n},2));
            BuOout_sem{n} = std(BuOout_Listmean{n})./ sqrt(size(BuOout_Listmean{n},2));
            FR_sem{n} = std(FRout_Listmean{n})./ sqrt(size(FRout_Listmean{n},2));
            
            BcapOutSD{n} = cat(1, x1.BcapOutSD);
            BcapOutSDmean{n} = mean(BcapOutSD{n});
            BcapOutSDsem{n} = std(BcapOutSD{n}) ./ sqrt(size(BcapOutSD{n},2));
            
            BinfOutSD{n} = cat(1, x1.BinfOutSD);
            BinfOutSDmean{n} = mean(BinfOutSD{n});
            BinfOutSDsem{n} = std(BinfOutSD{n}) ./ sqrt(size(BinfOutSD{n},2));

            %pairwise analysis
            NDPout_GPmat{n} = cat(2, NDPout.ListGPmean); % matrix of size N x M = (# of input pairs)x(# of recordings), containing the mean pairwise Soutput values
            SFout_GPmat{n} = cat(2, SFout.ListGPmean);
            Rout_GPmat{n} = cat(2, Rout.ListGPmean);
            BuCout_GPmat{n} = cat(2, BuCout.ListGPmean);
            BuOout_GPmat{n} = cat(2, BuOout.ListGPmean);
            FRout_GPmat{n} = cat(2, FRout.ListGPmean);
            
            for i = 1:length(Stim.SF{n}) % stat: different from ID line? (i.e. significant pat sep or conv)
            
            [H, Psf2{n}(i,1)] = ttest(SFout_GPmat{n}(i,:)-Stim.SF{n}(i)); clear H
            [H, Pndp2{n}(i,1)] = ttest(NDPout_GPmat{n}(i,:)-Stim.NDP{n}(i)); clear H
            [H, Pr2{n}(i,1)] = ttest(Rout_GPmat{n}(i,:)-Stim.R{n}(i)); clear H
            [H, Pbuc2{n}(i,1)] = ttest(BuCout_GPmat{n}(i,:)-Stim.BuC{n}(i)); clear H
            [H, Pbuo2{n}(i,1)] = ttest(BuOout_GPmat{n}(i,:)-Stim.BuO{n}(i)); clear H
            [H, Pfr2{n}(i,1)] = ttest(FRout_GPmat{n}(i,:)-Stim.FR{n}(i)); clear H
            
            end
    
            NDPout_GPmean{n} = mean(NDPout_GPmat{n},2); % cell containing column vectors of Soutput for each group, averaged across all M recordings
            SFout_GPmean{n} = mean(SFout_GPmat{n},2);
            Rout_GPmean{n} = mean(Rout_GPmat{n},2);
            FRout_GPmean{n} = mean(FRout_GPmat{n},2);
            BuCout_GPmean{n} = mean(BuCout_GPmat{n},2);
            BuOout_GPmean{n} = mean(BuOout_GPmat{n},2);

            NDPout_GPsem{n} = std( NDPout_GPmat{n},0,2 ) ./ sqrt(size(NDPout_GPmat{n},2)); %columns of SEM for each pair of stim, across recordings 
            SFout_GPsem{n} = std( SFout_GPmat{n},0,2 ) ./ sqrt(size(SFout_GPmat{n},2)); 
            Rout_GPsem{n} = std( Rout_GPmat{n},0,2 ) ./ sqrt(size(Rout_GPmat{n},2)); 
            FRout_GPsem{n} = std( FRout_GPmat{n},0,2 ) ./ sqrt(size(FRout_GPmat{n},2));
            BuCout_GPsem{n} = std( BuCout_GPmat{n},0,2 ) ./ sqrt(size(BuCout_GPmat{n},2));
            BuOout_GPsem{n} = std( BuOout_GPmat{n},0,2 ) ./ sqrt(size(BuOout_GPmat{n},2));

            RespGPallPoints_NDP{n} = reshape(NDPout_GPmat{n}, length(Stim.NDP{n})*length(x1),1); % column vector containing pairwise Soutput values for all recordings. To plot against "StimallPoints"
            RespGPallPoints_SF{n} = reshape(SFout_GPmat{n}, length(Stim.SF{n})*length(x1),1);
            RespGPallPoints_R{n} = reshape(Rout_GPmat{n}, length(Stim.R{n})*length(x1),1);
            RespGPallPoints_FR{n} = reshape(FRout_GPmat{n}, length(Stim.FR{n})*length(x1),1);
            RespGPallPoints_BuC{n} = reshape(BuCout_GPmat{n}, length(Stim.BuC{n})*length(x1),1);
            RespGPallPoints_BuO{n} = reshape(BuOout_GPmat{n}, length(Stim.BuO{n})*length(x1),1);
            
            clear NDPout SFout Rout STTCout rawSTTCout x1 x2
       end
        
        StimPair_NDP{c}{b} = cat(1,Stim.NDP{valid}); %column-vector with the pairwise Rinput for all stim protocols, to plot with OutPairMean and OutPairSEM
        StimPair_SF{c}{b} = cat(1,Stim.SF{valid});
        StimPair_R{c}{b} = cat(1,Stim.R{valid});
        StimPair_BuC{c}{b} = cat(1, Stim.BuC{valid});
        StimPair_BuO{c}{b} = cat(1, Stim.BuO{valid});
        StimPair_FR{c}{b} = cat(1,Stim.FR{valid});
        Stim_burstSD{c}{b} = Stim.BcapInSD(valid,1);
        Stim_burstSD2{c}{b} = Stim.BinfInSD(valid,1);
        
        OutPairMean_NDP{c}{b} = cat(1,NDPout_GPmean{valid}); %column with all the average pairwise Rinput for all stim families
        OutPairMean_SF{c}{b} = cat(1,SFout_GPmean{valid}); 
        OutPairMean_R{c}{b} = cat(1,Rout_GPmean{valid}); 
        OutPairMean_FR{c}{b} = cat(1,FRout_GPmean{valid});
        OutPairMean_BuC{c}{b} = cat(1,BuCout_GPmean{valid});
        OutPairMean_BuO{c}{b} = cat(1,BuOout_GPmean{valid});
        OutMean_burstSD{c}{b} = cat(1,BcapOutSDmean{valid});
        OutMean_burstSD2{c}{b} = cat(1,BinfOutSDmean{valid});

        OutPairSEM_NDP{c}{b} = cat(1,NDPout_GPsem{valid}); % idem with SEM
        OutPairSEM_SF{c}{b} = cat(1,SFout_GPsem{valid});
        OutPairSEM_R{c}{b} = cat(1,Rout_GPsem{valid});
        OutPairSEM_FR{c}{b} = cat(1,FRout_GPsem{valid});
        OutPairSEM_BuC{c}{b} = cat(1,BuCout_GPsem{valid});
        OutPairSEM_BuO{c}{b} = cat(1,BuOout_GPsem{valid});
        OutSEM_burstSD{c}{b} = cat(1,BcapOutSDsem{valid});
        OutSEM_burstSD2{c}{b} = cat(1,BinfOutSDsem{valid});
        
        P_sf{c}{b} = cat(1, Psf2{valid});
        P_ndp{c}{b} = cat(1, Pndp2{valid});
        P_r{c}{b} = cat(1, Pr2{valid});
        P_buc{c}{b} = cat(1, Pbuc2{valid});
        P_buo{c}{b} = cat(1, Pbuo2{valid});
        P_fr{c}{b} = cat(1, Pfr2{valid});

        StimPairAll_NDP{c}{b} = cat(1, StimAllPoints_NDP{valid});
        StimPairAll_SF{c}{b} = cat(1, StimAllPoints_SF{valid});
        StimPairAll_R{c}{b} = cat(1, StimAllPoints_R{valid});
        StimPairAll_FR{c}{b} = cat(1, StimAllPoints_FR{valid});
        StimPairAll_BuC{c}{b} = cat(1, StimAllPoints_BuC{valid});
        StimPairAll_BuO{c}{b} = cat(1, StimAllPoints_BuO{valid});

        OutPairAll_NDP{c}{b} = cat(1, RespGPallPoints_NDP{valid});
        OutPairAll_SF{c}{b} = cat(1, RespGPallPoints_SF{valid});
        OutPairAll_R{c}{b} = cat(1, RespGPallPoints_R{valid});
        OutPairAll_FR{c}{b} = cat(1, RespGPallPoints_FR{valid});
        OutPairAll_BuC{c}{b} = cat(1, RespGPallPoints_BuC{valid});
        OutPairAll_BuO{c}{b} = cat(1, RespGPallPoints_BuO{valid});

        StimAll_burstSD{c}{b} = cat(1, StimAllBcapInSD{valid});
        StimAll_BuCgmd{c}{b} = cat(1, StimAllBuCgmd{valid});
        StimAll_burstSD2{c}{b} = cat(1, StimAllBinfInSD{valid});
        StimAll_BuOgmd{c}{b} = cat(1, StimAllBuOgmd{valid});
        StimAll_FRgmd{c}{b} = cat(1, StimAllFRgmd{valid});
                
        OutBcapSD{c}{b} = cat(1, BcapOutSD{valid});
        OutBuCgmd{c}{b} = cat(1, BuCout_Listmean{valid});
        OutBinfSD{c}{b} = cat(1, BinfOutSD{valid});
        OutBuOgmd{c}{b} = cat(1, BuOout_Listmean{valid});
        OutFRgmd{c}{b} = cat(1, FRout_Listmean{valid});
        
        DiffBurstSD{c}{b} = OutBcapSD{c}{b} - StimAll_burstSD{c}{b}; % positive values show pat sep via burstiness variations, negative values crrespond to pattern covergence
        DiffBurstSD2{c}{b} = OutBinfSD{c}{b} - StimAll_burstSD2{c}{b};
        PatSepBuC{c}{b} = OutBuCgmd{c}{b} - StimAll_BuCgmd{c}{b};
        PatSepBuO{c}{b} = OutBuOgmd{c}{b} - StimAll_BuOgmd{c}{b};
        PatSepFR{c}{b} = OutFRgmd{c}{b} - StimAll_FRgmd{c}{b};
        
        PatSepBuC_mean{c}{b} = mean(PatSepBuC{c}{b});
        PatSepBuO_mean{c}{b} = mean(PatSepBuO{c}{b});
        PatSepFR_mean{c}{b} = mean(PatSepFR{c}{b});
        
        PatSepBuC_sem{c}{b} = std(PatSepBuC{c}{b})./ sqrt(size(PatSepBuC{c}{b},2));
        PatSepBuO_sem{c}{b} = std(PatSepBuO{c}{b})./ sqrt(size(PatSepBuO{c}{b},2));
        PatSepFR_sem{c}{b} = std(PatSepFR{c}{b})./ sqrt(size(PatSepFR{c}{b},2));
            
        Stim_AllSpkBcap{c}{b} = cat(1, AllSpkBcapIn{valid});
        Out_AllSpkBcap{c}{b} = cat(2, AllSpkBcapOut{valid});
        Stim_AllSpkBinf{c}{b} = cat(1, AllSpkBinfIn{valid});
        Out_AllSpkBinf{c}{b} = cat(2, AllSpkBinfOut{valid});
        Stim_AllSpkFR{c}{b} = cat(1, AllSpkFRi{valid});
        Out_AllSpkFR{c}{b} = cat(2, AllSpkFRo{valid});
        
        Celltype{c}{b} = c.*ones(size(PatSepBuC{c}{b})); %for boxplots and stats
        
                %% STATs and Fits
       
    %polynomial Fitting 
    EST_NDP{c}{b} = polyfit(StimPairAll_NDP{c}{b}, OutPairAll_NDP{c}{b},Npol);
    EST_SF{c}{b} = polyfit(StimPairAll_SF{c}{b}, OutPairAll_SF{c}{b},Npol);
    EST_R{c}{b} = polyfit(StimPairAll_R{c}{b}, OutPairAll_R{c}{b},Npol);
    EST_FR{c}{b} = polyfit(StimPairAll_FR{c}{b}, OutPairAll_FR{c}{b},Npol);
    EST_BuC{c}{b} = polyfit(StimPairAll_BuC{c}{b}, OutPairAll_BuC{c}{b},Npol);
    EST_BuO{c}{b} = polyfit(StimPairAll_BuO{c}{b}, OutPairAll_BuO{c}{b},Npol);
     
    DTw_NDP{c}(b,:) = Rinterp - polyval(EST_NDP{c}{b},Rinterp);    
    DTw_SF{c}(b,:) = RinterpSF - polyval(EST_SF{c}{b},RinterpSF);
    DTw_R{c}(b,:) = Rinterp - polyval(EST_R{c}{b},Rinterp);    
%     DTw_FR{c}(b,:) = polyval(EST_FR{c}{b},Rinterp) - Rinterp;
%     DTw_BuC{c}(b,:) = polyval(EST_BuC{c}{b},Rinterp) - Rinterp;

    % linear regression (equivalent to polynomial fitting of order 1, but gives 95% CI)
    % Interpolated Pat Sep
    [EST2_NDP{c}{b},BINTndp{c}{b}]  = regress(OutPairAll_NDP{c}{b}, [ ones(size(StimPairAll_NDP{c}{b})) , StimPairAll_NDP{c}{b} ] );
    [EST2_SF{c}{b},BINTsf{c}{b}] = regress(OutPairAll_SF{c}{b}, [ ones(size(StimPairAll_SF{c}{b})) , StimPairAll_SF{c}{b} ] );
    [EST2_R{c}{b},BINTr{c}{b}] = regress(OutPairAll_R{c}{b}, [ ones(size(StimPairAll_R{c}{b})) , StimPairAll_R{c}{b} ] );
    [EST2_FR{c}{b},BINTfr{c}{b}] = regress(OutPairAll_FR{c}{b}, [ ones(size(StimPairAll_FR{c}{b})) , StimPairAll_FR{c}{b} ] );
    [EST2_BuC{c}{b},BINTbuc{c}{b}] = regress(OutPairAll_BuC{c}{b}, [ ones(size(StimPairAll_BuC{c}{b})) , StimPairAll_BuC{c}{b} ] );
    [EST2_BuO{c}{b},BINTbuo{c}{b}] = regress(OutPairAll_BuO{c}{b}, [ ones(size(StimPairAll_BuO{c}{b})) , StimPairAll_BuO{c}{b} ] );

    % different from 0? (ttest or signrank)
    [H, Pbuc{c}{b}] = ttest(PatSepBuC{c}{b}); clear H
    [H, Pbuo{c}{b}] = ttest(PatSepBuO{c}{b}); clear H
    [H, Pfr{c}{b}] = ttest(PatSepFR{c}{b}); clear H
        PSRbuc{c}{b} = signrank(PatSepBuC{c}{b});
        PSRbuo{c}{b} = signrank(PatSepBuO{c}{b});
        PSRfr{c}{b} = signrank(PatSepFR{c}{b});
    
        clear Pbuc2 Pbuo2 Pfr2 Psf2 Pndp2 Pr2 n Stim NDPout_GPmean SFout_GPmean Rout_GPmean STTCout_GPmean rawSTTCout_GPmean NDPout_GPsem SFout_GPsem Rout_GPsem STTCout_GPsem rawSTTCout_GPsem StimAllPoints_NDP StimAllPoints_SF StimAllPoints_R StimAllPoints_STTC StimAllPoints_rawSTTC RespGPallPoints_NDP RespGPallPoints_SF RespGPallPoints_R RespGPallPoints_STTC RespGPallPoints_rawSTTC
    end  % end of looping through c = celltypes 
end

%% display
Absc = -1:0.01:1;
grey = [0.6 0.6 0.6];

figure(1) % just 1 dataset at 1 timescale and 1 metric
b = 2; c= 1; % GCyo 10Hz (10ms) dataset

    subplot(1,3,1) % R
        plot(Absc, polyval(EST_R{c}{b}, Absc),'Color', grey , 'LineWidth',1.5 ); hold on
        XinOK = StimPair_R{c}{b}(P_r{c}{b} <= 0.05); XinNON = StimPair_R{c}{b}(P_r{c}{b} > 0.05);
        YmeanOK = OutPairMean_R{c}{b}(P_r{c}{b} <= 0.05); YmeanNON = OutPairMean_R{c}{b}(P_r{c}{b} > 0.05);
        YsemOK = OutPairSEM_R{c}{b}(P_r{c}{b} <= 0.05); YsemNON = OutPairSEM_R{c}{b}(P_r{c}{b} > 0.05);
        eOK = errorbar(XinOK, YmeanOK, YsemOK, 'x') ; hold on
        eOK.Color = [0 0.5 0.8];
        eNON = errorbar(XinNON, YmeanNON, YsemNON, 'x') ; hold on
        eNON.Color = [0 0 0];
        plot( [-1 1], [-1, 1], 'k--'); hold on
        plot( [-1 1], [0 0], 'Color', grey);
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [-0.1 1], 'ylim', [-0.1 1])
        title(['R' num2str(Xbinsize(b)) 'ms'])
        clear XinOK YmeanOK YsemOK XinNON YmeanNON YsemNON
    subplot(1,3,2)
        plot(Absc, polyval(EST_NDP{c}{b}, Absc),'Color', grey , 'LineWidth',1.5 ); hold on
        XinOK = StimPair_NDP{c}{b}(P_ndp{c}{b} <= 0.05); XinNON = StimPair_NDP{c}{b}(P_ndp{c}{b} > 0.05);
        YmeanOK = OutPairMean_NDP{c}{b}(P_ndp{c}{b} <= 0.05); YmeanNON = OutPairMean_NDP{c}{b}(P_ndp{c}{b} > 0.05);
        YsemOK = OutPairSEM_NDP{c}{b}(P_ndp{c}{b} <= 0.05); YsemNON = OutPairSEM_NDP{c}{b}(P_ndp{c}{b} > 0.05);
        eOK = errorbar(XinOK, YmeanOK, YsemOK, 'x') ; hold on
        eOK.Color = [0 0.5 0.8];
        eNON = errorbar(XinNON, YmeanNON, YsemNON, 'x') ; hold on
        eNON.Color = [0 0 0];
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['NDP' num2str(Xbinsize(b)) 'ms'])
        clear XinOK YmeanOK YsemOK XinNON YmeanNON YsemNON
    subplot(1,3,3)
        plot(Absc, polyval(EST_SF{c}{b}, Absc),'Color', grey , 'LineWidth',1.5 ); hold on
        XinOK = StimPair_SF{c}{b}(P_sf{c}{b} <= 0.05); XinNON = StimPair_SF{c}{b}(P_sf{c}{b} > 0.05);
        YmeanOK = OutPairMean_SF{c}{b}(P_sf{c}{b} <= 0.05); YmeanNON = OutPairMean_SF{c}{b}(P_sf{c}{b} > 0.05);
        YsemOK = OutPairSEM_SF{c}{b}(P_sf{c}{b} <= 0.05); YsemNON = OutPairSEM_SF{c}{b}(P_sf{c}{b} > 0.05);
        eOK = errorbar(XinOK, YmeanOK, YsemOK, 'x') ; hold on
        eOK.Color = [0 0.5 0.8];
        eNON = errorbar(XinNON, YmeanNON, YsemNON, 'x') ; hold on
        eNON.Color = [0 0 0];
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0.7 1], 'ylim', [0 1])
        title(['SF' num2str(Xbinsize(b)) 'ms'])
        
figure(2) % all timescale
colord2 = brewermap(9, 'YlOrRd'); %hot(15); %hot(12) when I was not plotting tw = 250ms
colord = [colord2(3,:); colord2(5,:); colord2(7,:); colord2(9,:)];
ColMap = brewermap(length(Rinterp), 'Spectral');
ColMapSF = brewermap(length(Rinterp), '*YlOrRd');

    subplot(2,3,1) % R
        plot(Absc, polyval(EST_R{c}{1}, Absc),'Color', colord(1,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_R{c}{2}, Absc),'Color', colord(2,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_R{c}{3}, Absc),'Color', colord(3,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_R{c}{8}, Absc),'Color', colord(4,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_R{c}{10}, Absc),'Color', [0.2 0 0.2], 'LineWidth',1.5 ); hold on
%         legend('5ms', '10ms', '20ms', '100ms', '500ms',  'Location', 'northwest')
        scatter(StimPair_R{c}{1}, OutPairMean_R{c}{1},25, 'x','MarkerEdgeColor', colord(1,:));hold on
        scatter(StimPair_R{c}{2}, OutPairMean_R{c}{2},25, 'x','MarkerEdgeColor', colord(2,:));hold on
        scatter(StimPair_R{c}{3}, OutPairMean_R{c}{3},25, 'x','MarkerEdgeColor', colord(3,:));hold on
        scatter(StimPair_R{c}{8}, OutPairMean_R{c}{8},25, 'x','MarkerEdgeColor', colord(4,:));hold on
        scatter(StimPair_R{c}{10}, OutPairMean_R{c}{10},25, 'x','MarkerEdgeColor', [0.2 0 0.2]);hold on
%         errorbar(StimPair_R{c}{1}, OutPairMean_R{c}{1}, OutPairSEM_R{c}{1}, 'x', 'color', colord(1,:), 'markerfacecolor', colord(1,:), 'MarkerSize',5) ; hold on 
%         errorbar(StimPair_R{c}{2}, OutPairMean_R{c}{2}, OutPairSEM_R{c}{2}, 'x', 'color', colord(2,:), 'markerfacecolor', colord(2,:), 'MarkerSize',5) ; hold on 
%         errorbar(StimPair_R{c}{3}, OutPairMean_R{c}{3}, OutPairSEM_R{c}{3}, 'x', 'color', colord(3,:), 'markerfacecolor', colord(3,:), 'MarkerSize',5) ; hold on 
%         errorbar(StimPair_R{c}{4}, OutPairMean_R{c}{4}, OutPairSEM_R{c}{4}, 'x', 'color', colord(4,:), 'markerfacecolor', colord(4,:), 'MarkerSize',5) ; hold on 
        plot( [-1 1], [-1, 1], 'k--'); hold on
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        title('R')
        
    subplot(2,3,4)
        for n=1:length(Rinterp)
        plot(Xbinsize, DTw_R{c}(:,n), 'Color', ColMap(n,:), 'LineWidth',2); hold on
        end
%         legend('1', '0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1', '0', 'Location', 'BestOutside')
        box off
        axis square
        xlabel('bin size, ms', 'fontsize', 14)
        ylabel('Pat Conv <--> Pat Sep', 'fontsize', 14)
        set(gca, 'xlim', [0 1000], 'ylim', [-1 1])
        
    subplot(2,3,2) % NDP
        plot(Absc, polyval(EST_NDP{c}{1}, Absc),'Color', colord(1,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_NDP{c}{2}, Absc),'Color', colord(2,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_NDP{c}{3}, Absc),'Color', colord(3,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_NDP{c}{8}, Absc),'Color', colord(4,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_NDP{c}{10}, Absc),'Color', [0.2 0 0.2], 'LineWidth',1.5 ); hold on
%         legend('5ms', '10ms', '20ms', '100ms', '500ms',  'Location', 'northwest')
        scatter(StimPair_NDP{c}{1}, OutPairMean_NDP{c}{1},25, 'x','MarkerEdgeColor', colord(1,:));hold on
        scatter(StimPair_NDP{c}{2}, OutPairMean_NDP{c}{2},25, 'x','MarkerEdgeColor', colord(2,:));hold on
        scatter(StimPair_NDP{c}{3}, OutPairMean_NDP{c}{3},25, 'x','MarkerEdgeColor', colord(3,:));hold on
        scatter(StimPair_NDP{c}{8}, OutPairMean_NDP{c}{8},25, 'x','MarkerEdgeColor', colord(4,:));hold on
        scatter(StimPair_NDP{c}{10}, OutPairMean_NDP{c}{10},25, 'x','MarkerEdgeColor', [0.2 0 0.2]);hold on
        plot( [-1 1], [-1, 1], 'k--'); hold on
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        title('NDP')
        
    subplot(2,3,5)
        for n=1:length(Rinterp)
        plot(Xbinsize, DTw_NDP{c}(:,n), 'Color', ColMap(n,:), 'LineWidth',2); hold on
        end
%         legend('1', '0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1', '0', 'Location', 'BestOutside')
        box off
        axis square
        xlabel('bin size, ms', 'fontsize', 14)
        ylabel('Pat Conv <--> Pat Sep', 'fontsize', 14)
        set(gca, 'xlim', [0 1000], 'ylim', [-1 1])
        
    subplot(2,3,3)
        %         plot(Absc, polyval(EST_SF{c}{1}, Absc),'Color', colord(1,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_SF{c}{2}, Absc),'Color', colord(2,:), 'LineWidth',1.5 ); hold on
        %         plot(Absc, polyval(EST_SF{c}{3}, Absc),'Color', colord(3,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_SF{c}{8}, Absc),'Color', colord(4,:), 'LineWidth',1.5 ); hold on
        plot(Absc, polyval(EST_SF{c}{10}, Absc),'Color', [0.2 0 0.2], 'LineWidth',1.5 ); hold on

        %         legend('10ms', '100ms', '500ms', 'Location', 'northwest')
        %         scatter(StimPair_SF{c}{1}, OutPairMean_SF{c}{1},25, 'x','MarkerEdgeColor', colord(1,:));hold on
        scatter(StimPair_SF{c}{2}, OutPairMean_SF{c}{2},25, 'x','MarkerEdgeColor', colord(2,:));hold on
        %         scatter(StimPair_SF{c}{3}, OutPairMean_SF{c}{3},25, 'x','MarkerEdgeColor', colord(3,:));hold on
        scatter(StimPair_SF{c}{8}, OutPairMean_SF{c}{8},25, 'x','MarkerEdgeColor', colord(4,:));hold on
        scatter(StimPair_SF{c}{10}, OutPairMean_SF{c}{10},25, 'x','MarkerEdgeColor', [0.2 0 0.2]);hold on
        plot( [-1 1], [-1, 1], 'k--'); hold on
        set(gca, 'xlim', [0.5 1], 'ylim', [0.5 1])
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        title('SF')
        
    subplot(2,3,6)
        for n=1:length(Rinterp)
        plot(Xbinsize,DTw_SF{c}(:,n), 'Color', ColMapSF(n,:), 'LineWidth',2); hold on
        end
        box off
        axis square
        xlabel('bin size, ms', 'fontsize', 14)
        ylabel('Pat Conv <--> Pat Sep', 'fontsize', 14)
        set(gca, 'xlim', [0 1000], 'ylim', [-0.5 0.5])
    
figure(3) %SF input sets A (6) and B (7) i.e. PdeltaFR and B10Hz -> figure 4
b = 2; %10 ms
    subplot(1,3,1)
        plot(Absc, polyval(EST_SF{6}{b}, Absc),'Color', grey , 'LineWidth',1.5 ); hold on
        XinOK = StimPair_SF{6}{b}(P_sf{6}{b} <= 0.05); XinNON = StimPair_SF{6}{b}(P_sf{6}{b} > 0.05);
        YmeanOK = OutPairMean_SF{6}{b}(P_sf{6}{b} <= 0.05); YmeanNON = OutPairMean_SF{6}{b}(P_sf{6}{b} > 0.05);
        YsemOK = OutPairSEM_SF{6}{b}(P_sf{6}{b} <= 0.05); YsemNON = OutPairSEM_SF{6}{b}(P_sf{6}{b} > 0.05);
        eOK = errorbar(XinOK, YmeanOK, YsemOK, 'x') ; hold on
        eOK.Color = [0 0.5 0.8];
        eNON = errorbar(XinNON, YmeanNON, YsemNON, 'x') ; hold on
        eNON.Color = [0 0 0];
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['SF: input set A' num2str(Xbinsize(b)) 'ms'])
    subplot(1,3,2)
        plot(Absc, polyval(EST_SF{7}{b}, Absc),'Color', grey , 'LineWidth',1.5 ); hold on
        XinOK = StimPair_SF{7}{b}(P_sf{7}{b} <= 0.05); XinNON = StimPair_SF{7}{b}(P_sf{7}{b} > 0.05);
        YmeanOK = OutPairMean_SF{7}{b}(P_sf{7}{b} <= 0.05); YmeanNON = OutPairMean_SF{7}{b}(P_sf{7}{b} > 0.05);
        YsemOK = OutPairSEM_SF{7}{b}(P_sf{7}{b} <= 0.05); YsemNON = OutPairSEM_SF{7}{b}(P_sf{7}{b} > 0.05);
        eOK = errorbar(XinOK, YmeanOK, YsemOK, 'x') ; hold on
        eOK.Color = [0 0.5 0.8];
        eNON = errorbar(XinNON, YmeanNON, YsemNON, 'x') ; hold on
        eNON.Color = [0 0 0];
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['SF: input set B' num2str(Xbinsize(b)) 'ms'])
    subplot(1,3,3)
        clear b
        tau = [2, 8, 9];
        okColor = [colord(2,:); colord(4,:); [0.3 0 0.4] ]; 
        for b = 1:length(tau)
            ok = tau(b); 
            BothSF_stim{ok} = [StimPairAll_SF{6}{ok}; StimPairAll_SF{7}{ok}];
            BothSF_out{ok} = [OutPairAll_SF{6}{ok}; OutPairAll_SF{7}{ok}];
            EST_BothSF{ok} = polyfit(BothSF_stim{ok}, BothSF_out{ok},Npol);
            scatter(StimPair_SF{6}{ok}, OutPairMean_SF{6}{ok}, 'x', 'MarkerEdgeColor', okColor(b,:)); hold on % GC adult FR variation
            scatter(StimPair_SF{7}{ok}, OutPairMean_SF{7}{ok}, 'x', 'MarkerEdgeColor', okColor(b,:)); hold on % GC adult Burst variation
            plot(Absc, polyval(EST_BothSF{ok}, Absc),'Color', okColor(b,:), 'LineWidth',1 ); hold on
        end
        plot( [0 4], [0, 4], 'k--'); hold on %identity line
        box off
        axis square
        xlabel('pw BuC in')
        ylabel('pw BuC out')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title('SF: input set A and B')
        
figure(4) % for each spiketrain, plot input compactness VS output burstiness
for b = 1:length(tau)
    ok = tau(b);
    subplot(1,length(tau),b)
        scatter(Stim_AllSpkBcap{6}{ok}, Out_AllSpkBcap{6}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(Stim_AllSpkBcap{7}{ok}, Out_AllSpkBcap{7}{ok}, 15, '<','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5); hold on
        plot( [0 30], [0 30], 'k--'); hold on %identity line
        axis square
        box off
        xlabel('input compactness');
        ylabel('output compactness');
        title([num2str(Xbinsize(ok)) 'ms'])
        xlim([0 1]);
        ylim([0 1]);
end
        
figure(5) % for each spiketrain, plot input FR and occupancy VS output FR and occupancy

for b = 1:length(tau)
    ok = tau(b);
    subplot(1,length(tau)+1,b)
        scatter(Stim_AllSpkBinf{6}{ok}, Out_AllSpkBinf{6}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(Stim_AllSpkBinf{7}{ok}, Out_AllSpkBinf{7}{ok}, 15, '<','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5); hold on
        plot( [0 30], [0 30], 'k--'); hold on %identity line
        axis square
        box off
        xlabel('input occupancy');
        ylabel('output occupancy');
        title([num2str(Xbinsize(ok)) 'ms'])
        xlim([0 10]);
        ylim([0 10]);
end
    subplot(1,length(tau)+1,length(tau)+1)
        scatter(Stim_AllSpkFR{6}{ok}, Out_AllSpkFR{6}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(Stim_AllSpkFR{7}{ok}, Out_AllSpkFR{7}{ok}, 15, '<','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5); hold on
        plot( [0 30], [0 30], 'k--'); hold on %identity line
        axis square
        box off
        ylabel('FRout');
        xlabel('FRin');
        title('FR (2 sec) ')
        xlim([0 33]);
        ylim([0 33]);
        
figure(6)
Absc2 = 0:0.01:30;

subplot(1,3,1) %compactness
for b = 1:length(tau)
    ok = tau(b); 
    BothBuC_stim{ok} = [StimPairAll_BuC{6}{ok}; StimPairAll_BuC{7}{ok}]; %combines datasets from input set A and B
    BothBuC_out{ok} = [OutPairAll_BuC{6}{ok}; OutPairAll_BuC{7}{ok}];
    EST_BothBuC{ok} = polyfit(BothBuC_stim{ok}, BothBuC_out{ok},Npol);
%          errorbar(StimPair_BuC{1}{ok}, OutPairMean_BuC{1}{ok}, OutPairSEM_BuC{1}{ok}, 'x', 'color', [0,1, 0]); hold on % GC
%          scatter(StimPair_BuC{6}{ok}, OutPairMean_BuC{6}{ok}, 'x', 'MarkerEdgeColor', okColor(b,:)); hold on % GC adult FR variation
%          scatter(StimPair_BuC{7}{ok}, OutPairMean_BuC{7}{ok}, 'x', 'MarkerEdgeColor', okColor(b,:)); hold on % GC adult Burst variation
    errorbar(StimPair_BuC{6}{ok}, OutPairMean_BuC{6}{ok}, OutPairSEM_BuC{6}{ok}, 'x', 'color', okColor(b,:)); hold on % GC adult FR variation
    errorbar(StimPair_BuC{7}{ok}, OutPairMean_BuC{7}{ok}, OutPairSEM_BuC{7}{ok}, 'x', 'color', okColor(b,:)); hold on % GC adult Burst variation
    plot(Absc2, polyval(EST_BothBuC{ok}, Absc2),'Color', okColor(b,:), 'LineWidth',1 ); hold on
    plot( [0 4], [0, 4], 'k--'); hold on %identity line
    box off
    axis square
    xlabel('in: pairwise diff')
    ylabel('out: pairwise diff')
    set(gca, 'xlim', [0 1], 'ylim', [0 1])
    title('Compactness')
    
end

subplot(1,3,2) %occupancy
for b = 1:length(tau)
    ok = tau(b); 
    BothBuO_stim{ok} = [StimPairAll_BuO{6}{ok}; StimPairAll_BuO{7}{ok}]; %combines datasets from input set A and B
    BothBuO_out{ok} = [OutPairAll_BuO{6}{ok}; OutPairAll_BuO{7}{ok}];
    EST_BothBuO{ok} = polyfit(BothBuO_stim{ok}, BothBuO_out{ok},Npol);
    scatter(StimPair_BuO{6}{ok}, OutPairMean_BuO{6}{ok}, 'x', 'MarkerEdgeColor', okColor(b,:)); hold on % GC adult FR variation
    scatter(StimPair_BuO{7}{ok}, OutPairMean_BuO{7}{ok}, 'x', 'MarkerEdgeColor', okColor(b,:)); hold on % GC adult Burst variation
    plot(Absc2, polyval(EST_BothBuO{ok}, Absc2),'Color', okColor(b,:), 'LineWidth',1 ); hold on
    plot( [0 10], [0, 10], 'k--'); hold on %identity line
    box off
    axis square
    xlabel('in: pairwise diff')
    ylabel('out: pairwise diff')
    set(gca, 'xlim', [0 6.5], 'ylim', [0 6.5])
    title('Occupancy')
    
end

subplot(1,3,3)
    errorbar(StimPair_FR{6}{b}, OutPairMean_FR{6}{b}, OutPairSEM_FR{6}{b}, 'xk'); hold on % GC adult FR variation
    errorbar(StimPair_FR{7}{b}, OutPairMean_FR{7}{b}, OutPairSEM_FR{7}{b}, '<', 'color', [0,0.5, 0]); hold on % GC adult Burst variation
    plot(Absc2, polyval(EST_FR{6}{b}, Absc2),'Color', grey, 'LineWidth',1 ); hold on
    plot( [0 30], [0 30], 'k--'); hold on %identity line
    box off
    axis square
    xlabel('in: pairwise diff')
    ylabel('out: pairwise diff')
%     set(gca, 'xlim', [0 27], 'ylim', [0 5])
    title('FR')

figure(7) %average pat sep per recording, for Compactness
cet = [1, 6, 7];
for b = 1:length(tau)
    ok = tau(b);
    CelltypesAll{b} = [Celltype{1}{ok}; Celltype{6}{ok}; Celltype{7}{ok}];
    PatSepBuCAll{b} = [PatSepBuC{1}{ok}; PatSepBuC{6}{ok}; PatSepBuC{7}{ok}];
    [Pkl2(b),ANOVATABkl2{b},STATSkl2{b}] = anova1(PatSepBuCAll{b}, CelltypesAll{b}, 'off');
    [pk2(b),tblk2{b},statsk2{b}] = kruskalwallis(PatSepBuCAll{b}, CelltypesAll{b}, 'off');
    Comp{b} = multcompare(statsk2{b}, 'display', 'off');
    
    subplot(1,length(tau),b)
        scatter(0.8 + 0.4.*rand(size(PatSepBuC{1}{ok})), PatSepBuC{1}{ok}, 15, 'o','MarkerFaceColor', [0,1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(1.8 + 0.4.*rand(size(PatSepBuC{6}{ok})), PatSepBuC{6}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(2.8 + 0.4.*rand(size(PatSepBuC{7}{ok})), PatSepBuC{7}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
            for ct = 1:3 %1:length(burstCT)
        okc = cet(ct);        
        errorbar(ct, PatSepBuC_mean{okc}{ok}, PatSepBuC_sem{okc}{ok}, '+', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1) ; hold on
        clear okc
            end
        plot( [0 10], [0 0], 'k')
        axis square
        box off
        ylabel('Compactness disp out - in');
        set(gca,'ysc', 'lin', 'xlim',[0 4],'XTick', [1:3], 'XTickLabel', {'GC', 'GC ad FR', 'GC ad Burst'}); xtickangle(45)
        title([num2str(Xbinsize(ok)) 'ms' ' p = ' num2str(pk2(b))]) 
end

figure(8) %average pat sep per recording, for Occupancy and FR
for b = 1:length(tau)
    ok = tau(b); 
    PatSepBuOAll{b} = [PatSepBuO{1}{ok}; PatSepBuO{6}{ok}; PatSepBuO{7}{ok}];
    [Pkl3(b),ANOVATABkl3{b},STATSkl3{b}] = anova1(PatSepBuOAll{b}, CelltypesAll{b}, 'off');
    [pk3(b),tblk3{b},statsk3{b}] = kruskalwallis(PatSepBuOAll{b}, CelltypesAll{b}, 'off');
    Comp{b} = multcompare(statsk2{b}, 'display', 'off');
    
    subplot(1,length(tau)+1,b)
        scatter(0.8 + 0.4.*rand(size(PatSepBuO{1}{ok})), PatSepBuO{1}{ok}, 15, 'o','MarkerFaceColor', [0,1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(1.8 + 0.4.*rand(size(PatSepBuO{6}{ok})), PatSepBuO{6}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(2.8 + 0.4.*rand(size(PatSepBuO{7}{ok})), PatSepBuO{7}{ok}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
            for ct = 1:3 %1:length(burstCT)
        okc = cet(ct);        
        errorbar(ct, PatSepBuO_mean{okc}{ok}, PatSepBuO_sem{okc}{ok}, '+', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1) ; hold on
        clear okc
            end
        plot( [0 10], [0 0], 'k')
        axis square
        box off
        ylabel('Occupancy disp out - in');
        set(gca,'ysc', 'lin', 'xlim',[0 4], 'XTick', [1:3], 'XTickLabel', {'GC', 'GC ad FR', 'GC ad Burst'}); xtickangle(45)
        title([num2str(Xbinsize(ok)) 'ms' ' p = ' num2str(pk3(b))]) 
end

    PatSepFRAll = [PatSepFR{1}{ok}; PatSepFR{6}{ok}; PatSepFR{7}{ok}];
    [Pkl2FR,ANOVATABkl2FR,STATSkl2FR] = anova1(PatSepFRAll, CelltypesAll{b}, 'off');
    [pk2FR,tblk2FR,statsk2FR] = kruskalwallis(PatSepFRAll, CelltypesAll{b}, 'off');
    CompFR = multcompare(statsk2FR, 'display', 'off');
    
    subplot(1,length(tau)+1,length(tau)+1)
        scatter(0.8 + 0.4.*rand(size(PatSepFR{1}{ok})), PatSepFR{1}{ok}, 15, 'o','MarkerFaceColor', [0,1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on%         scatter(1.8 + 0.4.*rand(size(DiffBurstSD{2}{b})), DiffBurstSD{2}{b}, 50, 'd','MarkerFaceColor', [1,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(1.8 + 0.4.*rand(size(PatSepFR{6}{b})), PatSepFR{6}{b}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(2.8 + 0.4.*rand(size(PatSepFR{7}{b})), PatSepFR{7}{b}, 15, 'o','MarkerFaceColor', [0,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        for ct = 1:3 %1:length(burstCT)
        okc = cet(ct);        
        errorbar(ct, PatSepFR_mean{okc}{b}, PatSepFR_sem{okc}{b}, '+', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1) ; hold on 
        end
        plot( [0 10], [0 0], 'k')
        axis square
        box off
        ylabel('FR disp out - in');
        set(gca,'ysc', 'lin', 'xlim',[0 4], 'XTick', [1:3], 'XTickLabel', {'GC','GC ad FR', 'GC ad Burst'}); xtickangle(45)
        title(['ANOVA1 p = ' num2str(pk2FR)]) 
        
figure(9) % Sin VS Sout: GC before gzn vs GC after gzn
clear b c cet tau
cet = 8:9;
tau = [2, 8, 9];
cetCmap = [0, 0, 0; 1, 0, 0];
ctypes = {'GC before','GC after'};

% R
for b = 1:length(tau)
%ANCOVA separate lines (home-made). Ctrl: c = 1 ; Other: c = 2
    [EST_NDP_gzn{b}, pF_NDP(b), Fstat_NDP(b)] = PolFitandFtest (Npol, StimPairAll_NDP{9}{b}, OutPairAll_NDP{9}{b}, StimPairAll_NDP{8}{b}, OutPairAll_NDP{8}{b}, 0);
    [EST_SF_gzn{b}, pF_SF(b), Fstat_SF(b)] = PolFitandFtest (Npol, StimPairAll_SF{9}{b}, OutPairAll_SF{9}{b}, StimPairAll_SF{8}{b}, OutPairAll_SF{8}{b},  0);
    [EST_R_gzn{b}, pF_R(b), Fstat_R(b)] = PolFitandFtest (Npol, StimPairAll_R{9}{b}, OutPairAll_R{9}{b}, StimPairAll_R{8}{b}, OutPairAll_R{8}{b}, 0);
    
    subplot(3,length(tau),b) 
         okb = tau(b);
        for c = 1:length(cet)
         okc = cet(c);
         plot(Absc, polyval(EST_R{okc}{okb}, Absc),'Color', cetCmap(c,:), 'LineWidth',1 ); hold on
         errorbar(StimPair_R{okc}{okb}, OutPairMean_R{okc}{okb}, OutPairSEM_R{okc}{okb}, 'Marker', 'x', 'MarkerSize', 5, 'Color', cetCmap(c,:),'MarkerFaceColor', cetCmap(c,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none') ; hold on
        end
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['R ' num2str(Xbinsize(okb)) 'ms'])

%NDP
    subplot(3,length(tau),b+3) 
         okb = tau(b);
        for c = 1:length(cet)
         okc = cet(c);
         plot(Absc, polyval(EST_NDP{okc}{okb}, Absc),'Color', cetCmap(c,:), 'LineWidth',1 ); hold on
         errorbar(StimPair_NDP{okc}{okb}, OutPairMean_NDP{okc}{okb}, OutPairSEM_NDP{okc}{okb}, 'Marker', 'x', 'MarkerSize', 5, 'Color', cetCmap(c,:),'MarkerFaceColor', cetCmap(c,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none') ; hold on
        end
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['NDP ' num2str(Xbinsize(okb)) 'ms'])

% SF
    subplot(3,length(tau),b+6) 
         okb = tau(b);
        for c = 1:length(cet)
         okc = cet(c);
         plot(Absc, polyval(EST_SF{okc}{okb}, Absc),'Color', cetCmap(c,:), 'LineWidth',1 ); hold on
         errorbar(StimPair_SF{okc}{okb}, OutPairMean_SF{okc}{okb}, OutPairSEM_SF{okc}{okb}, 'Marker', 'x', 'MarkerSize', 5, 'Color', cetCmap(c,:),'MarkerFaceColor', cetCmap(c,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none') ; hold on
        end
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0.4 1], 'ylim', [0.4 1])
        title(['SF ' num2str(Xbinsize(okb)) 'ms'])
end

figure(10)
    clear ok
    % paired t-tests
    for b = 1:length(tau)
        ok = tau(b);
        [H, Pbuc89(b)] = ttest(PatSepBuC{8}{ok}-PatSepBuC{9}{ok}); clear H
        [H, Pbuo89(b)] = ttest(PatSepBuO{8}{ok}-PatSepBuO{9}{ok}); clear H
        [H, Pfr89(b)] = ttest(PatSepFR{8}{ok}-PatSepFR{9}{ok}); clear H
    end
    clear ok
    ok = 2; % 10ms
    subplot(1,3,1)
      for n = 1:length(PatSepBuC{8}{ok})
      plot(1:2, [PatSepBuC{8}{ok}(n), PatSepBuC{9}{ok}(n)], 'Color', grey); hold on
      end
      errorbar(1, PatSepBuC_mean{8}{ok}, PatSepBuC_sem{8}{ok}, 'ok', 'markerfacecolor', 'k', 'LineWidth', 1)
      errorbar(2, PatSepBuC_mean{9}{ok}, PatSepBuC_sem{9}{ok}, 'or', 'markerfacecolor', 'r', 'LineWidth', 1)
      plot( [0 4], [0 0], 'k')
      axis square
      grid off
      box off
      xlim([0.5,2.5]);
      ylim([-0.05 0.05]);
      ylabel('<Compactness dispersion> {out - in}');
      title({[num2str(Xbinsize(ok)) 'ms'] ; [' p = ' num2str(Pbuc89(1))] } ) 
      set(gca,'XTick', [1:2], 'XTickLabel', {'GC before gzn', 'GC after gzn'}); xtickangle(45);

     subplot(1,3,2)
      for n = 1:length(PatSepBuO{8}{ok})
      plot(1:2, [PatSepBuO{8}{ok}(n), PatSepBuO{9}{ok}(n)], 'Color', grey); hold on
      end
      errorbar(1, PatSepBuO_mean{8}{ok}, PatSepBuO_sem{8}{ok}, 'ok', 'markerfacecolor', 'k', 'LineWidth', 1)
      errorbar(2, PatSepBuO_mean{9}{ok}, PatSepBuO_sem{9}{ok}, 'or', 'markerfacecolor', 'r', 'LineWidth', 1)
      plot( [0 4], [0 0], 'k')
      axis square
      grid off
      box off
      xlim([0.5,2.5]);
      ylim([-0.05 0.05]);
      ylabel('<Occupancy dispersion> {out - in}');
      title({[num2str(Xbinsize(ok)) 'ms'] ; [' p = ' num2str(Pbuo89(1))] } ) 
      set(gca,'XTick', [1:2], 'XTickLabel', {'GC before gzn', 'GC after gzn'}); xtickangle(45);

    subplot(1,3,3)
      for n = 1:length(PatSepBuC{8}{ok})
      plot(1:2, [PatSepFR{8}{ok}(n), PatSepFR{9}{ok}(n)], 'Color', grey); hold on
      end
      errorbar(1, PatSepFR_mean{8}{ok}, PatSepFR_sem{8}{ok}, 'ok', 'markerfacecolor', 'k', 'LineWidth', 1)
      errorbar(2, PatSepFR_mean{9}{ok}, PatSepFR_sem{9}{ok}, 'or', 'markerfacecolor', 'r', 'LineWidth', 1)
    %   plot(1:2, [meanFR, meanFR10ms], 'ko-', 'MarkerFaceColor', 'k');
      plot( [0 4], [0 0], 'k')
      axis square
      grid off
      box off
      xlim([0.5,2.5]);
    %   ylim([-0.05 0.03]);
      ylabel('<FR dispersion> {out - in}');
      title(['FR per recs: '' p = ' num2str(Pfr89(1))] ) 
      set(gca,'XTick', [1:2], 'XTickLabel', {'GC before gzn', 'GC after gzn'}); xtickangle(45);
  
figure(11) % Sin VS Sout: GC before gzn vs GC after gzn
clear b c cet tau
cet = 1:5;
tau = [5, 8, 9];
cetCmap = [0, 1, 0; 1, 0, 0; 1, 0.5, 0; 0, 0.8, 0.8; 0, 0, 0];
markers = {'o','d','^','s','s'};
ctypes = {'GC','FS','MC','GCgzn30','CA3gzn30'};

% R
for b = 1:length(tau)  
    subplot(3,length(tau),b) 
         okb = tau(b);
        for c = 1:length(cet)
         okc = cet(c);
         plot(Absc, polyval(EST_R{okc}{okb}, Absc),'Color', cetCmap(c,:), 'LineWidth',1 ); hold on
         errorbar(StimPair_R{okc}{okb}, OutPairMean_R{okc}{okb}, OutPairSEM_R{okc}{okb}, 'Marker', markers{c}, 'MarkerSize', 5, 'Color', cetCmap(c,:),'MarkerFaceColor', cetCmap(c,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none') ; hold on
        end
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['R ' num2str(Xbinsize(okb)) 'ms'])

%NDP
    subplot(3,length(tau),b+3) 
         okb = tau(b);
        for c = 1:length(cet)
         okc = cet(c);
         plot(Absc, polyval(EST_NDP{okc}{okb}, Absc),'Color', cetCmap(c,:), 'LineWidth',1 ); hold on
         errorbar(StimPair_NDP{okc}{okb}, OutPairMean_NDP{okc}{okb}, OutPairSEM_NDP{okc}{okb}, 'Marker', markers{c}, 'MarkerSize', 5, 'Color', cetCmap(c,:),'MarkerFaceColor', cetCmap(c,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none') ; hold on
        end
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0 1], 'ylim', [0 1])
        title(['NDP ' num2str(Xbinsize(okb)) 'ms'])

% SF
    subplot(3,length(tau),b+6) 
         okb = tau(b);
        for c = 1:length(cet)
         okc = cet(c);
         plot(Absc, polyval(EST_SF{okc}{okb}, Absc),'Color', cetCmap(c,:), 'LineWidth',1 ); hold on
         errorbar(StimPair_SF{okc}{okb}, OutPairMean_SF{okc}{okb}, OutPairSEM_SF{okc}{okb}, 'Marker', markers{c}, 'MarkerSize', 5, 'Color', cetCmap(c,:),'MarkerFaceColor', cetCmap(c,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none') ; hold on
        end
        plot( [-1 1], [-1, 1], 'k--'); hold on
        box off
        axis square
        xlabel('pw Sin')
        ylabel('pw Sout')
        set(gca, 'xlim', [0.4 1], 'ylim', [0.4 1])
        title(['SF ' num2str(Xbinsize(okb)) 'ms'])
end

figure(12) %Ancova pval for interaction model Sin*Celltypes (if significant, that means slopes are different)

% ANCOVA on all celltypes (same as my home-made F test with a linear model, when using a "separate lines' model (i.e. different slopes). 
    % But provides more flexibility and more info for posthoc stats. 
    
    %R
    subplot(1,3,1)
    for b = 1:length(Xbinsize)-1
        okb = b; %tau(b);
        Dumout = cat(1,OutPairAll_R{cet}); % length(cet) * length(Xbinsize) cell
        Sout{b} = cat(1, Dumout{:, okb});
        Dumin = cat(1,StimPairAll_R{cet});
        Sin{b} = cat(1,Dumin{:, okb}); %group1 (list of input corr for each data point)
        
        for c = 1:length(cet)
            okc = cet(c);
            for k=1:length(OutPairAll_R{okc}{okb})
            BlockGroup{c}{k,1}= ctypes{c};
            end    
        end
        BlockGroupCat{b} = cat(1,BlockGroup{:}); clear BlockGroup 

        [Hanc{b},TABanc{b},COEFanc{b},Sanc{b}] = aoctool(Sin{b},Sout{b},BlockGroupCat{b},0.05,'Sin','Sout','Celltype','off','separate lines');
        COMPancS{b} = multcompare(Sanc{b}, 'Estimate', 'slope', 'display','off');
        COMPancInt{b} = multcompare(Sanc{b}, 'Estimate', 'intercept', 'display','off');
        [MAT_pvalS{b}, MAT4_pS{b} ] = pvalmat(COMPancS{b});
        [MAT_pvalI{b}, MAT4_pI{b} ] = pvalmat(COMPancInt{b});
    end

        for b = 1:11
        Panc(b) = TABanc{b}{4,6}; 
        end
        
        plot(Xbinsize(1:11),Panc,'xr-', 'LineWidth',1 );
        box off
        axis square
        xlabel('tau, ms')
        ylabel('ANCOVA pval')
        title('R')
    clear Dumout Sout Sin Dumin
    
    %NDP
    subplot(1,3,2)
    for b = 1:length(Xbinsize)-1
        okb = b; %tau(b);
        Dumout = cat(1,OutPairAll_NDP{cet}); % length(cet) * length(Xbinsize) cell
        Sout{b} = cat(1, Dumout{:, okb});
        Dumin = cat(1,StimPairAll_NDP{cet});
        Sin{b} = cat(1,Dumin{:, okb}); %group1 (list of input corr for each data point)

        [Hanc2{b},TABanc2{b},COEFanc2{b},Sanc2{b}] = aoctool(Sin{b},Sout{b},BlockGroupCat{b},0.05,'Sin','Sout','Celltype','off','separate lines');
        COMPancS2{b} = multcompare(Sanc2{b}, 'Estimate', 'slope', 'display','off');
        COMPancInt2{b} = multcompare(Sanc2{b}, 'Estimate', 'intercept', 'display','off');
        [MAT_pvalS2{b}, MAT4_pS2{b} ] = pvalmat(COMPancS2{b});
        [MAT_pvalI2{b}, MAT4_pI2{b} ] = pvalmat(COMPancInt2{b});
    end

        for b = 1:11
        Panc2(b) = TABanc2{b}{4,6}; 
        end
        
        plot(Xbinsize(1:11),Panc2,'xr-', 'LineWidth',1 );
        box off
        axis square
        xlabel('tau, ms')
        ylabel('ANCOVA pval')
        title('NDP')
        clear Dumout Sout Sin Dumin
        
    %SF
    subplot(1,3,3)
    for b = 1:length(Xbinsize)-1
        okb = b; %tau(b);
        Dumout = cat(1,OutPairAll_SF{cet}); % length(cet) * length(Xbinsize) cell
        Sout{b} = cat(1, Dumout{:, okb});
        Dumin = cat(1,StimPairAll_SF{cet});
        Sin{b} = cat(1,Dumin{:, okb}); %group1 (list of input corr for each data point)

        [Hanc3{b},TABanc3{b},COEFanc3{b},Sanc3{b}] = aoctool(Sin{b},Sout{b},BlockGroupCat{b},0.05,'Sin','Sout','Celltype','off','separate lines');
        COMPancS3{b} = multcompare(Sanc3{b}, 'Estimate', 'slope', 'display','off');
        COMPancInt3{b} = multcompare(Sanc3{b}, 'Estimate', 'intercept', 'display','off');
        [MAT_pvalS3{b}, MAT4_pS3{b} ] = pvalmat(COMPancS3{b});
        [MAT_pvalI3{b}, MAT4_pI3{b} ] = pvalmat(COMPancInt3{b});
    end

        for b = 1:11
        Panc3(b) = TABanc3{b}{4,6}; 
        end
        
        plot(Xbinsize(1:11),Panc3,'xr-', 'LineWidth',1 );
        box off
        axis square
        xlabel('tau, ms')
        ylabel('ANCOVA pval')
        title('SF')

figure(13) % NDP & SF: slope and intercept = f(binsize) for all celltypes on same graph
% clear cet cetCmap; cet = [1, 3:5]; cetCmap = [0, 1, 0; 1, 0.5, 0; 0, 0.8, 0.8; 0, 0, 0];

    subplot(2,2,1) % NDP slope
           for c = 1:length(cet)
            okc = cet(c);
                for b = 1:11
                Slope2{okc}(b, 1) = EST2_NDP{okc}{b}(2);
                CIslopeL{okc}(b, 1) = BINTndp{okc}{b}(2,1);
                CIslopeH{okc}(b, 1) = BINTndp{okc}{b}(2,2);
                end
            okc = cet(c);
            plot_ci(Xbinsize(1:11)',[Slope2{okc},CIslopeL{okc},CIslopeH{okc}],'PatchColor', cetCmap(c,:), 'PatchAlpha', 0.2, ...
            'MainLineWidth', 1, 'MainLineStyle', '-', 'MainLineColor', cetCmap(c,:), ...
            'LineWidth', 1, 'LineStyle','none', 'LineColor', [1 1 1]); hold on
           end
    axis square
    xlabel('tau, ms')
    ylabel('Slope')
    title('NDP')

    subplot(2,2,2)
        for c = 1:length(cet)
        okc = cet(c);
            for b = 1:11
            Slope2{okc}(b, 1) = EST2_SF{okc}{b}(2);
            CIslopeL{okc}(b, 1) = BINTsf{okc}{b}(2,1);
            CIslopeH{okc}(b, 1) = BINTsf{okc}{b}(2,2);
            end
        okc = cet(c);
        plot_ci(Xbinsize(1:11)',[Slope2{okc},CIslopeL{okc},CIslopeH{okc}],'PatchColor', cetCmap(c,:), 'PatchAlpha', 0.2, ...
              'MainLineWidth', 1, 'MainLineStyle', '-', 'MainLineColor', cetCmap(c,:), ...
              'LineWidth', 1, 'LineStyle','none', 'LineColor', [1 1 1]); hold on

        end
    axis square
    xlabel('tau, ms')
    ylabel('Slope')
    xlim([40 250]);
    title('SF zoom')

subplot(2,2,3)
       for c = 1:length(cet)
    okc = cet(c);
        for b = 1:11
        int2{okc}(b, 1) = EST2_NDP{okc}{b}(1); %intercept
        CIintL{okc}(b, 1) = BINTndp{okc}{b}(1,1);%intercept
        CIintH{okc}(b, 1) = BINTndp{okc}{b}(1,2);%intercept
        end
     okc = cet(c);
     plot_ci(Xbinsize(1:11)',[int2{okc},CIintL{okc},CIintH{okc}],'PatchColor', cetCmap(c,:), 'PatchAlpha', 0.2, ...
              'MainLineWidth', 1, 'MainLineStyle', '-', 'MainLineColor', cetCmap(c,:), ...
              'LineWidth', 1, 'LineStyle','none', 'LineColor', [1 1 1]); hold on

       end
axis square
xlabel('tau, ms')
ylabel('Intercept')
title('NDP')

subplot(2,2,4)
       for c = 1:length(cet)
    okc = cet(c);
        for b = 1:11
        int2{okc}(b, 1) = EST2_SF{okc}{b}(1); %intercept
        CIintL{okc}(b, 1) = BINTsf{okc}{b}(1,1);%intercept
        CIintH{okc}(b, 1) = BINTsf{okc}{b}(1,2);%intercept
        end
     okc = cet(c);
     plot_ci(Xbinsize(1:11)',[int2{okc},CIintL{okc},CIintH{okc}],'PatchColor', cetCmap(c,:), 'PatchAlpha', 0.2, ...
              'MainLineWidth', 1, 'MainLineStyle', '-', 'MainLineColor', cetCmap(c,:), ...
              'LineWidth', 1, 'LineStyle','none', 'LineColor', [1 1 1]); hold on

       end
axis square
xlabel('tau, ms')
ylabel('Intercept')
xlim([40 250]);
title('SF zoom')

figure(14) %average pat sep per recording, for Compactness
cet = [1:5];
clear tau;
tau = [2, 8, 9];

for b = 1:length(tau)
    ok = tau(b);
    CelltypesAll{b} = [Celltype{1}{ok}; Celltype{2}{ok}; Celltype{3}{ok}; Celltype{4}{ok}; Celltype{5}{ok}];
    PatSepBuCAll{b} = [PatSepBuC{1}{ok}; PatSepBuC{2}{ok}; PatSepBuC{3}{ok}; PatSepBuC{4}{ok}; PatSepBuC{5}{ok}];
    [Pkl4(b),ANOVATABkl4{b},STATSkl4{b}] = anova1(PatSepBuCAll{b}, CelltypesAll{b}, 'off');
    [pk4(b),tblk4{b},statsk4{b}] = kruskalwallis(PatSepBuCAll{b}, CelltypesAll{b}, 'off');
    Comp4{b} = multcompare(statsk4{b}, 'display', 'off');
    
    subplot(1,length(tau)+1,b)
        scatter(0.8 + 0.4.*rand(size(PatSepBuC{1}{ok})), PatSepBuC{1}{ok}, 15, 'o','MarkerFaceColor', [0,1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(1.8 + 0.4.*rand(size(PatSepBuC{2}{ok})), PatSepBuC{2}{ok}, 15, 'd','MarkerFaceColor', [1,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(2.8 + 0.4.*rand(size(PatSepBuC{3}{ok})), PatSepBuC{3}{ok}, 15, '^','MarkerFaceColor', [1,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(3.8 + 0.4.*rand(size(PatSepBuC{4}{ok})), PatSepBuC{4}{ok}, 15, 's','MarkerFaceColor', [0,0.8, 0.8], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(4.8 + 0.4.*rand(size(PatSepBuC{5}{ok})), PatSepBuC{5}{ok}, 15, 's','MarkerFaceColor', [0,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
            for ct = 1:5
        okc = cet(ct);        
        errorbar(ct, PatSepBuC_mean{okc}{ok}, PatSepBuC_sem{okc}{ok}, '+', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1) ; hold on
        clear okc
            end
        plot( [0 10], [0 0], 'k')
        axis square
        box off
        ylabel('Compactness disp out - in');
        set(gca,'ysc', 'lin', 'xlim', [0 6], 'XTick', [1:5], 'XTickLabel', {'GC', 'FS', 'HMC', 'GC+gzn', 'CA3+gzn'}) ; xtickangle(45);
        title([num2str(Xbinsize(ok)) 'ms' ' p = ' num2str(pk4(b))]) 
end

figure(15) %average pat sep per recording, for Occupancy and FR
for b = 1:length(tau)
    ok = tau(b); 
    PatSepBuOAll{b} = [PatSepBuO{1}{ok}; PatSepBuO{2}{ok}; PatSepBuO{3}{ok}; PatSepBuO{4}{ok}; PatSepBuO{5}{ok}];
    [Pkl5(b),ANOVATABkl5{b},STATSkl4{b}] = anova1(PatSepBuOAll{b}, CelltypesAll{b}, 'off');
    [pk5(b),tblk5{b},statsk5{b}] = kruskalwallis(PatSepBuOAll{b}, CelltypesAll{b}, 'off');
    Comp5{b} = multcompare(statsk2{b}, 'display', 'off');
    
    subplot(1,length(tau)+1,b)
        scatter(0.8 + 0.4.*rand(size(PatSepBuO{1}{ok})), PatSepBuO{1}{ok}, 15, 'o','MarkerFaceColor', [0,1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(1.8 + 0.4.*rand(size(PatSepBuO{2}{ok})), PatSepBuO{2}{ok}, 15, 'o','MarkerFaceColor', [1,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(2.8 + 0.4.*rand(size(PatSepBuO{3}{ok})), PatSepBuO{3}{ok}, 15, 'o','MarkerFaceColor', [1,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(3.8 + 0.4.*rand(size(PatSepBuO{4}{ok})), PatSepBuO{4}{ok}, 15, 's','MarkerFaceColor', [0,0.8, 0.8], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(4.8 + 0.4.*rand(size(PatSepBuO{5}{ok})), PatSepBuO{5}{ok}, 15, 's','MarkerFaceColor', [0,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
            for ct = 1:5
        okc = cet(ct);        
        errorbar(ct, PatSepBuO_mean{okc}{ok}, PatSepBuO_sem{okc}{ok}, '+', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1) ; hold on
        clear okc
            end
        plot( [0 10], [0 0], 'k')
        axis square
        box off
        ylabel('Occupancy disp out - in');
        set(gca,'ysc', 'lin', 'xlim', [0 6], 'XTick', [1:5], 'XTickLabel', {'GC', 'FS', 'HMC', 'GC+gzn', 'CA3+gzn'}) ; xtickangle(45);
        title([num2str(Xbinsize(ok)) 'ms' ' p = ' num2str(pk5(b))]) 
end

    PatSepFRAll = [PatSepFR{1}{ok}; PatSepFR{2}{ok}; PatSepFR{3}{ok}; PatSepFR{4}{ok}; PatSepFR{5}{ok}];
    [Pkl3FR,ANOVATABkl3FR,STATSkl3FR] = anova1(PatSepFRAll, CelltypesAll{b}, 'off');
    [pk3FR,tblk3FR,statsk3FR] = kruskalwallis(PatSepFRAll, CelltypesAll{b}, 'off');
    Comp3FR = multcompare(statsk3FR, 'display', 'off');
    
    subplot(1,length(tau)+1,length(tau)+1)
        scatter(0.8 + 0.4.*rand(size(PatSepFR{1}{ok})), PatSepFR{1}{ok}, 15, 'o','MarkerFaceColor', [0,1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on%         scatter(1.8 + 0.4.*rand(size(DiffBurstSD{2}{b})), DiffBurstSD{2}{b}, 50, 'd','MarkerFaceColor', [1,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(1.8 + 0.4.*rand(size(PatSepFR{2}{ok})), PatSepFR{2}{ok}, 15, 'd','MarkerFaceColor', [1,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(2.8 + 0.4.*rand(size(PatSepFR{3}{ok})), PatSepFR{3}{ok}, 15, '^','MarkerFaceColor', [1,0.5, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(3.8 + 0.4.*rand(size(PatSepFR{4}{ok})), PatSepFR{4}{ok}, 15, 's','MarkerFaceColor', [0,0.8, 0.8], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        scatter(4.8 + 0.4.*rand(size(PatSepFR{5}{ok})), PatSepFR{5}{ok}, 15, 's','MarkerFaceColor', [0,0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
        for ct = 1:5
        okc = cet(ct);        
        errorbar(ct, PatSepFR_mean{okc}{b}, PatSepFR_sem{okc}{b}, '+', 'color', 'k', 'markerfacecolor', 'k', 'LineWidth', 1) ; hold on 
        end
        plot( [0 10], [0 0], 'k')
        axis square
        box off
        ylabel('FR disp out - in');
        set(gca,'ysc', 'lin', 'xlim', [0 6], 'XTick', [1:5], 'XTickLabel', {'GC', 'FS', 'HMC', 'GC+gzn', 'CA3+gzn'}) ; xtickangle(45);
        title(['ANOVA1 p = ' num2str(pk3FR)]) 
        