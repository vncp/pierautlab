function Y = spikecounts2matrix(time, spiketimes)



disp('Converting spiketime lists to rasters')

% Convert to cell if not already
notcell = 0;
if ~iscell(spiketimes)
    notcell = 1;
    spiketimes = {spiketimes};
end

nchans = length(spiketimes);
dt = mean(diff(time));

Y = zeros(length(time), nchans);
for chan = 1:nchans
    disp(['Channel ' num2str(chan) ' of ' num2str(nchans)])
    drawnow
    stimes = spiketimes{chan};
        if isempty(stimes)==1
            Y(:,chan) = zeros(size(time)); %in case of no firing (no AP detected)
        else
        Y(:, chan) = histc(stimes, time);
        end
end    

%display exemple
% Y = spikecounts2matrix([-100:100:1100], times1);
% Converting spiketime lists to binary rasters
% Channel 1 of 1
% figure; bar([-100:100:1100], Y, 'histc')
% figure; bar([-100:100:1100], Y)