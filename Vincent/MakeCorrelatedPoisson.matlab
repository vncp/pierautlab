%& Pieraut Lab - Recoded program for MATLAB
% Uses MVD algorithm to generate correlated Poisson
% Uses toolbox described in Macke et al.


clear all
close all

%% Define parameters to generate spiketrains

% User Defined Parameters
fs		= 10000; 			% Sampling rate (Hz)
T 		= 2;				% Length of trains (sec)
N 		= 5;				% Number of spiketrains
mnrate	= 30.*ones(N, 1);	% Mean rate for each train (spikes/s)
							% ones() creates array of all 1's
							% sampleCovPoisson() creates slightly lower rates, change mu up to compensate if needed
corval	= 0.25;				% Cross-train correlation (uniform)

% Calculated Parameters
mu		= mnrate./fs;		% Mean probability per bin
v		= mu;				% Variance for each train
M 		= fs .* T;			% Length of trains (frequency*timee)
Sigxy 	= sqrt(v)*sqrt(v');	% Matrix variance products (sigma.x*sigma.y)
C 		= corval.*Sigxy;	% Covariance matrix (off-diagonals)
C(eye(N)==1)	= v;		% Variances into diagonal of matrix

% Sample from a Poisson distribution
% Length of trains (# of time bins)
dt 		= 1./fs;			% time step size
time 	= 0:dt:T-dt;		% time axis (0 to T-dt with iterations of dt)


nogood	= 1;
count = 0;
while nogood 								% infinite loop (while true; only stops on break statement)
	count++;	
	disp(['Attempt #' num2str(count)]) 		% displays attempts to terminal
	X = sampleCovPoisson (mu, C, M);		% Samples from multivariate correlated binary and Poisson random variables
	X(find(X)) = 1;							% Truncates values to 1 because some are greater for some reason

	mnrateHat = mean(sum(X,2)./T)';			% Estimates mean rate
	mnrateSD = std(sum(X,2)./T)';			% Standard deiation
	CorrHat = corr(X');						% Estimates correlation
	CorrHat = mean(Corrhat(find(triu(CorrHat,1))));	% Gives means of the upper part of the triangular matrix
	CorrHat = mean (CorrHat(:));			% Gives mean of diagonals
	CorrHatSD = std(CorrHat);				% Standard deviation of CorrHat matrix
	disp([mnrateHat CorrHat])				% prints to terminal the estimated mean rate and the means of CorrHat

	if[abs(mean(mnrateHat)-mean(mnrate))./mean(mnrate)<0.01] & [abs(CorrHat-corval)./corval < 0.001]
		break
	end
end

% Spiketimes and ISIs
dmn = 0.1 .* 1 ./ mean(mnrate);
ISI = {};
ISIbins = 0:dmn:100.*dmn;
ISIhist = [];
collabels{1} = 'Time, s';

for n = 1:N
	spiketimes{n}	= time(find(X(n,:)));
	ISI{N}			= diff(spiketimes{n});
	ISIhist(:,n)	= hist(ISI{n} , ISIbins);
	collabels{n+1}	= num2str(n);
end

% Fits data into an exponential, to test against Poisson distribution
x = repmat(ISIbins',1,N);
[x, ind] = sort(x(:));;
y = ISIhist(:);
y = y(ind);
guess = [100, 1./mean(mnrate), 1];
figure
guess = fminsearch('fit1exp', guess, [], x, y, 0);		%
[amp, tau, const] = deal(guess(1),guess(2),guess(3));
fit = amp .* exp(-x./tau) + const;
close(gcf)

%Data Display
figure('units', 'inch', 'pos' [1 1 0 10])
subplot(211)
displayflag			= 0;
colord				= lines(N);

