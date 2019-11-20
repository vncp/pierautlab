function [Within, Between] = SimMat_WandB(nStims, nrepeats, mat)
% [Within, Between] = SimMat_WandB(nStims, nrepeats, mat)
%
% takes a matrix of pairwise coefficients belonging to different groups and extract the Within and Between comparison groups elements, excluding the diagonal. 
% Within coeffs are coming from different repetitions of the same input groups (i.e. Resp1 vs Resp1, Resp2 vs Resp2, etc...). 
% Between coeffs are coming from the non-Within comparison groups (i.e. Resp1 vs Resp2, Resp1 vs Resp3, etc...)
% In a spiketrain similarity analysis, it provides, respectively, the Sw, i.e. the ability of
% the recorded system to reproduce the same output from repetitions of the same input, and the Sb, i.e. the actual output similarity.  
%
% Arguments:
% - nStims: number of groups (in the case of spiketrain similarity analysis, it is the number of input spiketrains in a protocol)
% - nrepeats: number of repetition of each input trace.
% - mat : similarity matrix of pairwise coefficient organized such that
%         all output sweeps coming from the same input trace are grouped together. 
%
% Outputs: 
% - Within is a Struct containing the following fields:
%           - List: column vector listing all the coef values comparing spiketrains from same input trace for all input traces.
%           - ListMean and ListSEM: average Sw and corresponding SEM (from values in List)
%           - ListGPmean: vector of length nStims listing the (nrepeats^2) mean values of Sw for each stim trace parent. 
%           - ListGP: a cell containing nStims column vectors listing the pairwise coeffs for each group
%
% - Between is a Struct containing the following fields:           
%           - List: column vector listing all the ( ((nStims^2 - nStims)./2)*nrepeats^2 ) coef values  comparing spiketrains coming from different input traces.
%           - ListMean and ListSEM: average Sb (i.e. Soutput) and corresponding SEM (from values in List)
%           - ListGPmean: vector of length ((nStims^2 - nStims)./2) listing the mean values of Sb for each non-within comparison group of outputs (i.e. Resp1 vs Resp2, Resp1 vs Resp3, etc...)   
%           - ListGP: a cell containing ((nStims^2 - nStims)./2) column vectors listing the pairwise coeffs for each group           - 
%
% Antoine Madar Septembre 2017

% Within
    for n = 1:nStims
            pts1 = 1 + (n-1).*nrepeats;
            pts2 = pts1+nrepeats-1;
            dummy1 = mat(pts1:pts2, pts1:pts2);
            nanmask = ones(size(dummy1));
            nanmask = triu(nanmask,1)./triu(nanmask,1);
            x = nanmask.*triu(dummy1, 1);
            WListGP{n} = x(~isnan(x)); % takes all the values of the matrix dummy1 to remove NaN from NaN mask and make it a column vector.
            WListGPmean(n)=mean(WListGP{n});
            clear dummy1 x nanmask
    end
    
    WList = cat(1,WListGP{:});
    WListMean    = mean( WList );
    WListSEM     = std( WList )./sqrt(length(WList));
    
    % package in a struct
    Within.List = WList;
    Within.ListMean = WListMean;
    Within.ListSEM = WListSEM;
    Within.ListGPmean = WListGPmean'; %' so that it's a column
    Within.ListGP = WListGP;
    
% Between
        clear n
    for n = 1:(nStims-1) % for nStims-1 group rows
        row1 = 1 + (n-1).*nrepeats; % start row 
        row2 = row1 + nrepeats - 1; % end row
        for m = n+1:nStims % for group columns that do not correspond to the input group of the selected row.
        col1 = 1 + (m-1).*nrepeats ; % start column
        col2 = col1 + nrepeats - 1; % end column
        BListGP{n,m} = mat(row1:row2, col1:col2); % takes a square (nrepeats x nrepeats) corresponding to an individual comparison group in the matrix mat (e.g. Resp1 vs Resp2, then Resp1 vs Resp2, ..., Resp1 vs Resp(nStims))
        BListGP{n,m} = BListGP{n,m}(~isnan(BListGP{n,m})); %removes NaN values (there shouldn't be any. If there were, it would create an error later in the similarity analysis) AND makes it a vector
        BmatGPmean(n,m) = mean(BListGP{n,m}); % warning: this creates artificial zeros in the lower part of the matrix and on the diagonal.
        end
    end
    BListGP = BListGP(~cellfun('isempty',BListGP)); %remove the empty cells generated when organizing the cell array into an n x m matrix
    
    % remove the zeros that were artificially created on the lower triangle + diag of BmatGPmean (but not the ones that might be in the upper triangle, as they are real values) + make it a column vector.
    nanmask = ones(size(BmatGPmean));
    nanmask = triu(nanmask,1)./triu(nanmask,1);
    x = nanmask.*triu(BmatGPmean, 1);
    BListGPmean = x(~isnan(x)); % takes all the values of the matrix dummy1 to remove NaN from NaN mask and make it a column vector.     
    
    BList = cat(1,BListGP{:});
    BListMean    = mean( BList );
    BListSEM     = std( BList )./sqrt(length(BList));

    % package in a struct
    Between.List = BList;
    Between.ListMean = BListMean;
    Between.ListSEM = BListSEM;
    Between.ListGPmean = BListGPmean;
    Between.ListGP = BListGP;
    
    

    
    
    