% [NDP, SF] = NDP_SF(SpiketrainA, SpiketrainB)
%
% NDP and SF of a pair of binned rasters (spiketimes converted to counts using spikecounts2matrix with a given binsize)
%
% Antoine Madar April 2017 

function [NDP, SF] = NDP_SF(SpiketrainA, SpiketrainB) 

Dot = dot(SpiketrainA, SpiketrainB); %Dot product of column vectors A and B = ||A||*||B||*cos(theta), with theta = angle between vect A and B
Norm = norm(SpiketrainA).*norm(SpiketrainB); % ||A||x||B|| = magnitude of the dot product
NDP = Dot./Norm; % cos(theta) = normalized dot product

% compute ||A||./||B|| = scaling factor between the 2 vectors. It's a function of how different the firing rate in each bin is.
% To get a symetric matrix, always compute ||A||/||B|| with ||A|| < ||B||.
% This way, 0 < scaling factor < 1
% When A or B = 0 (i.e when there is no spike in one of the rasters), make the ratio = 0 
    if norm(SpiketrainA) > 0 && norm(SpiketrainB) > 0 %if both A and B have a norm > 0, i.e at least 1 spike in both of the rasters. 
        if norm(SpiketrainA) <= norm(SpiketrainB) 
           SF = norm(SpiketrainA)./norm(SpiketrainB); 
        else  
           SF = norm(SpiketrainB)./norm(SpiketrainA);
        end
    else
           SF = 0;
           
    end
    
    % when SF is 0, because either A, B or both don't have any spikes, make
    % it NaN to exclude it from further analysis
    if SF == 0
        SF = NaN;
    end

% compute Burstiness similarity
% N_spk
% BurstA = ./length(SpiketrainA)

end
