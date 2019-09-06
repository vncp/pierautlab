%% Use MVD algorithm to generate Poisson spiketrains with a given correlation coeff and min SF
% Uses the toolbox described in Macke et al. 2009 (Neural Computation 21, 397-423)
clear all
close all

%% Step 1 - Define the statistics and generate spiketrains

%******************* USER DEFINED PARAMETERS *****************************************
fs      = 10000;           % Sampling rate (Hz)
T       = 2;               % Length of trains (sec)

N       = 10;               % Number of spiketrains to generate
corval  = 0.5;            % Cross-train Correlation (uniform for simplicity, but does not have to be)
error = 0.2;               % max difference between min and max correlation
minSF = 0.4;
% minFR = 5; %Hz
% maxFR = 30; %Hz

% set mean rates for the N spiketrains
% mnrate = randi([minFR, maxFR], N, 1); % randomly selected rates 
mnrate = [6:3:33]'; %[5 10 17 25 30]

dbin = 0.010;              % binsize to compute SF
%******************* USER DEFINED PARAMETERS *****************************************


% parameters
dt  = 1./fs;            % time step
time= 0:dt:T-dt;        % time axis
bins = [0:dbin:T];        % time vector for spikecounts2matrix (sec)

nogood = 1;
count = 0;
while nogood
count = count + 1;

% Set up the covariance matrix
mu  =   mnrate ./ fs;       	% Mean probability per bin
v       = mu;               % Variance for each train (this needs to be of order mu, or correlations disappear!)
M   =   fs .* T;              % Length of trains (# of time bins)
Sigxy   = sqrt(v)*sqrt(v'); % Matrix of variance products (sigma_x * sigma_y)
C       = corval .* Sigxy;  % Covariance matrix (for off-diagonals only)
C(eye(N)==1) = v;           % Put the variances on the diagonal

% Sample from a Poisson distribution
M   = fs .* T;          % Length of trains (# of time bins)
X = sampleCovPoisson(mu,C,M); 
X(X>1) = 1; % For some reason, some values are greater than 1, so truncate.

% get spiketimes and downsample spiketrain to have a raster of binsize = dbin. 
    collabels{1} = 'Time, s';
    for i = 1:N
        spiketimes{i} = time( X(i,:) > 0 ); % get stimtimes (
        SpkCount{i} = spikecounts2matrix(bins, spiketimes{i});
        collabels{i+1}  = num2str(i);
    end
    clear i
    
% measure Similarity    
    for i = 1:N % for each spiketrain    
        for j = 1:N % for each spiketrain
             [NDPmat(i,j), SFmat(i,j)] = NDP_SF(SpkCount{i}, SpkCount{j}); % compute NDP and SF
             Rdummy = corrcoef(SpkCount{i}, SpkCount{j}); % compute Pearson's coef of correlation
             Rmat(i,j) = Rdummy(1,2); clear Rdummy;
             
             % track FR of each spiketrain
             mFRi(i,j) = sum(SpkCount{i})./T;
             mFRj(i,j) = sum(SpkCount{j})./T;
             % track the difference of firing rate btw spiketrains (such that it is always positive)
            if sum(SpkCount{j}) <= sum(SpkCount{i}) 
               difFR(i,j) = (sum(SpkCount{i})- sum(SpkCount{j}))./T; 
            else  
               difFR(i,j) = (sum(SpkCount{j})- sum(SpkCount{i}))./T;
            end
             
        end  
    end
    
[ Xsf, sf_mean, sf_sem ] = SymMat2List(SFmat);
[ Xr, r_mean, r_sem ] = SymMat2List(Rmat);
[ XdifFR, difFR_mean, difFR_sem ] = SymMat2List(difFR);
[ XmFRi, mFRi_mean, mFRi_sem ] = SymMat2List(mFRi);
[ XmFRj, mFRj_mean, mFRj_sem ] = SymMat2List(mFRj);

    if min(Xsf) <= minSF && abs(min(Xr)- max(Xr)) <= error
        nogood = 0;
    end

end   

%% Display
figure('units', 'inch', 'pos', [1 1 8 10])

subplot(311) %rasters
displayflag                     = 0;
colord                          = lines(N); 
[spikex, spikey]                = makedisplayrasters(spiketimes, displayflag, colord);
displaydisplayrasters(spikex, spikey, colord)
titletxt = {['minFR = ' num2str(min(XmFRi)) 'maxFR' num2str(max(XmFRi))] ; ['R = ' num2str(r_mean) ', maxR = ' num2str(max(Xr)) ', minR = ' num2str(min(Xr))]};
title(titletxt, 'fontname', 'times', 'fontsize', 14)

subplot(312) % matrix of SF
imagesc(SFmat); hold on
axis square
colormap(hot); 
caxis([0 1]);
colorbar;
title('SF')

subplot(313) %
scatter(XdifFR, Xsf)
grid off
axis square
xlabel('difference of firing rate')
ylabel('SF')
title('How firing rate impacts SF')

cd('/Users/antoine/Google Drive/LabReasearch/Protocols/PatchProtocols&Config/SFvaryingFR')
print('10-2sInputTrains_R0.75_0.6to0.85_SFrange_PoissonFRraster_fig2.eps', '-depsc', '-painters')
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
        basename = ['2sInputTrains_R' num2str(r_mean) 'SFrange_PoissonFR'];
        [f, p] = uiputfile('*', 'Base name to export: (without extension)', basename);
        % Write binary rasters for Axograph import
        writeascii([time', X'], collabels, [p f 'binary.txt'] )
        % Write figures
        print([p f 'raster_fig.tif'], '-dtiff', '-r150')
        % Write mat-file
        save([p f 'raster_fig.mat'])
    case 'No'
        close all
end        

disp('DONE')