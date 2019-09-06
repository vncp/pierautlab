%% Use MVD algorithm to generate correlated Poisson spiketrains.
% Uses the toolbox described in Macke et al. 2009 (Neural Computation 21, 397-423)
clear all
close all

%% Step 1 - Define the statistics and generate spiketrains

%******************* USER DEFINED PARAMETERS *****************************************
fs      = 10000;            % Sampling rate (Hz)
T       = 2;               % Length of trains (sec)
N       = 5;               % Number of spiketrains to generate

mnrate  = 30.*ones(N,1); % Mean rate for each train (spike per sec)
                            % NOTE - sampleCovPoisson sems to produce
                            % slightly lower rates than the specified mu,
                            % can correct by making mu slightly higher.
corval  = 0.25;            % Cross-train Correlation (uniform for simplicity, but does not have to be)
%******************* USER DEFINED PARAMETERS *****************************************

% Set up the covariance matrix
mu  =   mnrate ./ fs;       	% Mean probability per bin
v       = mu;               % Variance for each train (this needs to be of order mu, or correlations disappear!)
M   =   fs .* T;              % Length of trains (# of time bins)
Sigxy   = sqrt(v)*sqrt(v'); % Matrix of variance products (sigma_x * sigma_y)
C       = corval .* Sigxy;  % Covariance matrix (for off-diagonals only)
C(eye(N)==1) = v;           % Put the variances on the diagonal



% Sample from a Poisson distribution
M   = fs .* T;          % Length of trains (# of time bins)
dt  = 1./fs;            % time step
time= 0:dt:T-dt;        % time axis



nogood = 1;
count = 0;
while nogood
    count = count + 1;
    disp(['Attempt #' num2str(count)])
    X = sampleCovPoisson(mu,C,M); 
    % For some reason, some values are greater than 1, so truncate.
    X(find(X)) = 1;

    mnrateHat = mean(sum(X,2)./T)';  % estimate mean rate
    mnrateSD  = std(sum(X,2)./T)';
    CorrHat = corr(X');     % estimate correlation
    CorrHat = mean( CorrHat( find(triu(CorrHat, 1) ) ) );
    CorrHat = mean( CorrHat(:) );
    CorrHatSD = std( CorrHat );
    disp([mnrateHat CorrHat])
    if [ abs(mean(mnrateHat)-mean(mnrate))./mean(mnrate) < 0.01 ] & [ abs(CorrHat-corval)./corval < 0.001 ]
        nogood = 0;
    end
end
    

% Get spiketimes & ISIs
dmn = 0.1.* 1./mean(mnrate);
ISI = {};
ISIbins = 0 : dmn : 100.*dmn;
ISIhist = [];
collabels{1} = 'Time, s';
for n = 1:N
    spiketimes{n} = time( find( X(n,:) ) );
    ISI{n}          = diff( spiketimes{n} );
    ISIhist(:, n)   = hist( ISI{n} , ISIbins);
    collabels{n+1}  = num2str(n);
end   

% Fit an exponential, as a test of Poisson-ness
x = repmat(ISIbins', 1, N);
[x, ind] = sort(x(:));
y = ISIhist(:);
y = y(ind);
guess = [100, 1./mean(mnrate), 1];
figure
guess = fminsearch('fit1exp', guess, [], x, y, 0);
[amp, tau, const] = deal(guess(1), guess(2), guess(3));
fit  =  amp .* exp(-x./tau) + const;
close(gcf)





% Display
figure('units', 'inch', 'pos', [1 1 8 10])
subplot(211)
displayflag                     = 0;
colord                          = lines(N); 
[spikex, spikey]                = makedisplayrasters(spiketimes, displayflag, colord);
displaydisplayrasters(spikex, spikey, colord)
titletxt = {['Counted Mean Rate = ' num2str(mnrateHat, '%1.2f') '/sec']; ['Mean Correlation = ' num2str(CorrHat, '%1.3f')]};
title(titletxt, 'fontname', 'times', 'fontsize', 14)

subplot(212)
plot(ISIbins, ISIhist, '.', 'markersize', 6); hold on
plot(x, fit, 'k-', 'linewidth', 2); 
text( mean(ISIbins), 0.5.*max(ISIhist(:)), ['ISI Fitted Mean Rate = ' num2str(1./tau, '%1.3f') '/sec'], 'fontname', 'times', 'fontsize', 14); hold off
box off
xlabel('ISI (sec)')
ylabel('Count')



% Export Data
defunits = get(0, 'defaultfigureunits');
set(0, 'defaultfigureunits', 'pix')
   ButtonName = questdlg('Do you want to export these trains?', ...
                         'Export Data?', ...
                         'Yes', 'No', 'No');
set(0, 'defaultfigureunits', defunits)



switch ButtonName
    case 'Yes'
        basepath = [ num2str(mnrateHat, '%1.0f') 'Hz_R' num2str(CorrHat, '%1.3f')];
        basename = [ num2str(mnrateHat, '%1.0f') 'Hz_R' num2str(CorrHat, '%1.3f') '_CorrelatedPoisson.'];
        mkdir(basepath)
        parentdir = cd;
        cd(basepath)
        [f, p] = uiputfile('*', 'Base name to export: (without extension)', basename);
        % Write binary rasters for Axograph import
        writeascii([time', X'], collabels, [p f 'binary.txt'] )
        % Write figures
        print([p f 'raster_fig.tif'], '-dtiff', '-r150')
        % Write mat-file
        save([p f 'raster_fig.mat'])
        cd(parentdir)
    case 'No'
        close all
end        

disp('DONE')





    

