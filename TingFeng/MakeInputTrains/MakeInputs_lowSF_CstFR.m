clear 
close all

%% Step 1 - Define the statistics and generate spiketrains

%******************* USER DEFINED PARAMETERS *****************************************
fs      = 10000;           % Sampling rate (Hz)
T       = 2;               % Length of trains (sec)
FR = 10;                   % target Firing rate (Hz)
dbin = 0.010;              % binning window to compute similarity (sec)
Nbins = [2:2:20];          % number of occupied bins for each spiketrains. N = length(Nbins) is the number of spiketrains 
%******************* USER DEFINED PARAMETERS *****************************************

Ns = FR*T;  % total number of spikes in a spiketrain
bins = [0:dbin:T-dbin];
MaxDelay = dbin./2;
timeD = -MaxDelay:1./fs:MaxDelay-1./fs; % create a pool of numbers in the right range of delay values, between -dbin/2 and dbin/2, with a resolution defined by the sampling rate
% Nbins = [1:Ns] ;  % number of bins to fill with spikes
Burst = Ns./Nbins; % range of burstiness. The larger, the burstier.

%% create spiketrains
N = length(Nbins); %number of spiketrains
for n = 1:N % for each number of occupied bins value
 
        Burstiness = Ns/Nbins(n); 
        NspkInteger = fix(Burstiness); %rounds to nearest integer towards 0
        Nmiss = rem(Ns,Nbins(n)); %remainder of the euclidean division: Ns./Nbins = fix(Ns./Nbins) + remainder./Nbins

        OccBin{n} = repmat(NspkInteger,Nbins(n),1); %vector with number of spikes in each occupied bin
        
        if Nmiss > 0
            for i = 1:Nmiss %redistribute the missing spike among the bins already filled
                Iredist = randi(Nbins(n),1,1);
                OccBin{n}(Iredist) = OccBin{n}(Iredist)+1;
                clear Iredist
            end
        end

        if sum(OccBin{n})== Ns
            'good'
            num2str(n)
        else
            return
        end
        
         %associate a filled bin to an index in bins, take the center and
         %distribute spiketimes given their delays
         Ispk = sort(randperm(length(bins),Nbins(n))); %randomly selected unique indices of bins to be filled with spikes, sorted in ascending order
         CenterT = round(bins(Ispk)+dbin./2, fs); %center time of each spike, given with fs resolution
            clear j
            for j= 1:length(CenterT)
            D = randperm(length(timeD), OccBin{n}(j)); % randomly select as many unique indices of timeD as spikes need to be added
            Delay = timeD(D)'; % sample indices from the pool of possible delays in a [-dbin/2,dbin/2] window
            SpkT{n}{j} = repmat(CenterT(j),length(OccBin{n}(j)),1) + Delay;  % create a vector of the new spiketimes
    %         clear Delay D
            end

        SpkTimes{n} = cat(1, SpkT{n}{:});
        SpkCount{n} = spikecounts2matrix(bins, SpkTimes{n});
        FRspk(n) = sum(SpkCount{n})./T;
        tNbins(n) = (Nbins(n));
        mNbins(n) = (length(SpkCount{n}(SpkCount{n}>0)));
        tBurst(n) = (Burstiness); % target burstiness
        mBurst(n) = (sum(SpkCount{n})./mNbins(n)); % measured burstiness  (to compare with the target burstiness)
            if tBurst(n) - mBurst(n) > 1 ||tBurst(n) - mBurst(n) < -1
                return
            end
    %     clear lambda Nbins OccBin Ispk CenterT
    end

    tBurstall = tBurst;
    mBurstall = mBurst;
    

    % measure similarity
    for i = 1:length(SpkCount)
        for j = 1:length(SpkCount)
            [NDPmat(i,j), SFmat(i,j)] = NDP_SF(SpkCount{i}, SpkCount{j}); % compute NDP and SF
             Rdummy = corrcoef(SpkCount{i}, SpkCount{j}); % compute Pearson's coef of correlation
             Rmat(i,j) = Rdummy(1,2); clear Rdummy;
             
             mBursti(i,j) = mBurst(i);
             mBurstj(i,j) = mBurst(j);
             mNbinsi(i,j) = mNbins(i);
             mNbinsj(i,j) = mNbins(j);

             % track FR difference
             if sum(SpkCount{j}) <= sum(SpkCount{i}) 
                   difFR(i,j) = (sum(SpkCount{i})- sum(SpkCount{j}))./T; 
             else  
                   difFR(i,j) = (sum(SpkCount{j})- sum(SpkCount{i}))./T;
             end

             % track burstiness difference
             if mBurst(j) <= mBurst(i) 
                   difBurst(i,j) = mBurst(i) - mBurst(j); 
             else  
                   difBurst(i,j) = mBurst(j) - mBurst(i);
             end
             
             if mNbins(j) <= mNbins(i) 
                   difNbins(i,j) = mNbins(i) - mNbins(j); 
             else  
                   difNbins(i,j) = mNbins(j) - mNbins(i);
             end
        end    
    end

    [ Xcos, ndp_mean, ndp_sem ] = SymMat2List(NDPmat);
    [ Xsf, sf_mean, sf_sem ] = SymMat2List(SFmat);
    [ Xr, r_mean, r_sem ] = SymMat2List(Rmat);
    [ XdifFR, difFR_mean, difFR_sem ] = SymMat2List(difFR);
    [ XdifBurst, difBurst_mean, difBurst_sem ] = SymMat2List(difBurst);
    [ XdifNbins, difNbins_mean, difNbins_sem ] = SymMat2List(difNbins);
    [XmNbinsi, im, isem] = SymMat2List(mNbinsi);
    [XmNbinsj, jm, jsem] = SymMat2List(mNbinsj);

% make spiketrains from spiketimes with fs resolution + time vector (with
% fs resol) and associated labels
dt  = 1./fs;            % time step
time= 0:dt:T-dt;        % time axis

SpkRasters = spiketimes2matrix(time, SpkTimes); %transforms cells of spiketimes into a matrix of rasters, each column being a spiketrain. 

collabels{1} = 'Time, s';
for n = 1:N
    collabels{n+1}  = num2str(n);
end   

%% Display
figure('units', 'inch', 'pos', [1 1 8 10])

subplot(311) %rasters
displayflag                     = 0;
colord                          = lines(N); 
[spikex, spikey]                = makedisplayrasters(SpkTimes, displayflag, colord);
displaydisplayrasters(spikex, spikey, colord)
titletxt = {['minFR = ' num2str(min(FRspk)) 'maxFR' num2str(max(FRspk))]; ['R = ' num2str(r_mean) ', maxR = ' num2str(max(Xr)) ', minR = ' num2str(min(Xr))]};
title(titletxt, 'fontname', 'times', 'fontsize', 14)

subplot(312) % matrix of SF
imagesc(SFmat); hold on
axis square
colormap(hot); 
caxis([0 1]);
colorbar;
title('SF')

subplot(313) %
scatter(XdifNbins, Xsf)
grid off
axis square
xlabel('difference of # spiking bins')
ylabel('SF')
title('How burstiness differences impact SF')

cd('/Users/antoine/Google Drive/LabReasearch/Protocols/PatchProtocols&Config/SFvaryingBurstiness');
print('10_2sInputTrains_10Hz_SFrange_Burstinessraster_fig', '-depsc', '-painters');
%% Export Data
defunits = get(0, 'defaultfigureunits');
set(0, 'defaultfigureunits', 'pix')
   ButtonName = questdlg('Do you want to export these trains?', ...
                         'Export Data?', ...
                         'Yes', 'No', 'No');
set(0, 'defaultfigureunits', defunits)

switch ButtonName
    case 'Yes'
        
        basepath = '/Users/antoine/Google Drive/LabReasearch/jonesLab_git/jonesLab_code/Projects/PatSep_Antoine_2017_oldNamefile/MakeCorrelatedTrains/SF';
        cd(basepath)
        basename = ['2sInputTrains_' num2str(mean(FRspk)) 'Hz_SFrange_Burstiness'];
        [f, p] = uiputfile('*', 'Base name to export: (without extension)', basename);
        % Write binary rasters for Axograph import
        writeascii([time', SpkRasters], collabels, [p f 'binary.txt'] )
        % Write figures
        print([p f 'raster_fig.tif'], '-dtiff', '-r150')
        % Write mat-file
        save([p f 'raster_fig.mat'])
    case 'No'
        close all
end        

disp('DONE')








