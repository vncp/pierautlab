function [spiketimes, spikelocs, peaktimes, peaklocs] = detectspikes(time, data, signaltype, threshtype, thresh, peakflag, displayflag);

% [spiketimes, spikelocs, peaktimes, peaklocs] = detectspikes(time, data, signaltype, threshtype, thresh, peakflag, displayflag);
% 
% Detect the times of threshold crossings (POSITIVE ONLY, for negative crossings, data must be preprocesed appropriately to make crossing positive).
% TIME is an array of sequential time values, one for each point in DATA. 
% DATA is an array of data values, one for each point in TIME. Each column is a sweep.
% SIGNALTYPE is either 'raw' or 'gradient', and determines whether the raw
%   signal or the first derivative will be used for detection.
% THRESHTYPE is either 'direct' or 'SD' and determines whether an absolute amplitude 
%     or a multiple of the StDev of the signal will be used for detection.
% THRESH is the threshold of detection, expressed in units that depend on THRESHTYPE. 
% PEAKFLAG is either 0 or 1 and determines whether the times of peaks are also found (takes time). If 0, PEAKLOCS and PEAKTIMES will be empty
% DISPLAYFLAG is either 0 or 1 and specifies whether to display.
% SPIKETIMES is an array listing the time of each threshold crossing.
% SPIKELOCS is an index into the TIME or DATA arrays coresponding to the
%   spiketimes
% PEAKTIMES is the approximate time of peaks following threshold crossing
% PEAKLOCS is an array of indices corresponding to the peaks
%
% -MJ 2009
% 
% example:
% 
% time    = 1:1000;
% data    = filter( [1 2 2 1], 1, [1./rand(length(time), 1)].^0.5 );
% signaltype  = 'gradient';
% threshtype  = 'SD';
% thresh      = 2;
% peakflag    = 1;
% displayflag = 1;
% [spiketimes, spikelocs, peaktimes, peaklocs] = detectspikes(time, data, signaltype, threshtype, thresh, peakflag, displayflag)





len = length(time);

peaklocs = {};
peaktimes = {};


mintime = min(time);            % correct for if time doesn't start at zero
time    = time - mintime;



for sweep = 1:size(data,2)
    st = cputime;
    tic
    disp('  ')
    disp('Preprocessing data')
    if strcmp( signaltype, 'gradient' )
        disp('Computing gradient')
        signal = gradient( data(:, sweep) );


    
    else
        signal = data(:, sweep);
    end

    if strcmp( threshtype, 'SD' )
        truethresh = mean( signal ) + thresh.*std( signal );
    else
        truethresh = thresh;
    end
    toc



    %   Detect events
    tic
    disp('  ')
    disp('Finding events')
        over                            = find( signal >= truethresh );
        ideal                           = 0.*signal;
        ideal(over)                     = 1; 
        diffideal                       = [0; diff(ideal)]; 
        idealstarts                     = diffideal;
        idealstarts(idealstarts<1)      = 0;

    %     Get spiketimes
        spikelocs{sweep}    = find(diffideal==1);       % Indices of THRESHOLD CROSSING
        spiketimes{sweep}   = time(spikelocs{sweep});   % Time of THRESHOLD CROSSING
    toc   


    if peakflag == 1
    tic
    disp('  ')
    disp('Finding peaks')
        if length(spiketimes{sweep}) > 0
            for s = 1:length( spiketimes{sweep}) % Approximate the time of the peak 
                if s < length(spiketimes{sweep})
                    nxt = spikelocs{sweep}(s+1);
                else
                    nxt = length(time);
                end    

                [mx, indx]          = max( data( spikelocs{sweep}(s):nxt, sweep));
                peaklocs{sweep}(s)  = [indx+spikelocs{sweep}(s)-1]';
                peaktimes{sweep}(s) = [time(peaklocs{sweep}(s))]';
            end
        else
            % Enter a [] if this sweep has no spikes
            peaklocs{sweep} = [];
            peaktimes{sweep} = [];
        end    
    toc    
    end




    if displayflag
    tic
    disp('  ')
    disp('Displaying events')
    tic
        subplot(2, 1, 1)
            plot( time, data(:, sweep), 'k-', spiketimes{sweep}, data(spikelocs{sweep}, sweep), 'b.');  box off
            if peakflag && length(peaktimes{sweep})>0
                hold on
                plot( peaktimes{sweep}, data(peaklocs{sweep}, sweep), 'ro'); hold off
            end    
            set(gca, 'xlim', time([1, end]), 'ylim', 1e-3.*[-120 80])
            legend('data', 'crossing', 'peak', 'location', 'northeast')
        subplot(2, 1, 2)
            if strcmp( signaltype, 'gradient' ) & strcmp( threshtype, 'SD' )
                plot(time, signal, 'b-', time([1, end]), mean(signal).*[1 1], 'r-', time([1,end]), truethresh.*[1 1], 'r--'); box off
                set(gca, 'xlim', time([1, end]))
                legend('derivative', 'mean', ['mean+' num2str(thresh) '*SD'], 'location', 'northeast')
            else
                plot(time, signal, 'b-'); hold on; box off
                plot(time([1,end]), truethresh.*[1 1], 'r--'); hold off
                set(gca, 'xlim', time([1, end]))
                legend('data', 'threshold', 'location', 'northeast')
            end
    toc
    end    



    % undo time correction
    time = time + mintime;
    spiketimes{sweep}  = spiketimes{sweep}  + mintime;
    if peakflag
        peaktimes{sweep}   = peaktimes{sweep}   + mintime;
    end
    
    elt = cputime-st;
    disp('  ')
    disp( ['Total time for ' num2str(len) ' data points was ' num2str(elt) ' seconds - Approx ' num2str(elt./len) ' sec/pt.']);


end







