% [Bcap, Binf] = burstinessbin(SpkCntVector)
%
% burstiness for binned spiketrains: B = nb of bins that could be occupied / nb of occupied bins
%
% argument: vector of spike counts in each bin. 
%output: 
% Binf is occupancy (number of spikes per bin)
% Bcap is compactness: 1 - the proportion of occupied bins (1 - occupied bins / total number of bins). ~1 is most compact (1 cannot be
% attained. 1-1/totbins is max), 0 is the true minimum and means all bins
% were occupied, i.e. least compact.

% Antoine Madar February 2018

function [Bcap, Binf] = burstinessbin(SpkCntVector)

    nspk = sum(SpkCntVector); % total number of spike in spiketrain
    maxBins = length(SpkCntVector);
    
    if nspk == 0 % when no spikes in spiketrain
     Binf = NaN;
     Bcap = NaN;

    else
%         % compute number of bins that could be occupied
% 
%         if nspk > maxBins
%             potentialBins = maxBins; % this caps the value of B at length(SpkCntVector).
%         else
%             potentialBins = nspk;
%         end

        % compute nb of occupied bins
        OccupBins = length(SpkCntVector(SpkCntVector > 0)); %number of bins that have at least a spike

        %compute Burstiness
        Binf = nspk./OccupBins;
        
        % compute compactness index (independent of nspk, i.e. indpdt of FR)
        Bcap = 1 - OccupBins./maxBins; 
        
    end

end


