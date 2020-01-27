% take a symetric matrix X and return the list of unique values in the upper (or lower, doesn't matter since X is symetric) triangle, excluding the diagonal. 
% Note: It keeps NaN values
%
% Antoine Madar, 2018/11/6

function [list, average, sem] = SymMat2ListNaN(X) 
    
nanmask = ones(size(X));
nanmask = triu(nanmask,1)./triu(nanmask,1); 

% list = nanmask.*triu(X, 1);
list = X(~isnan(nanmask));
average = mean( list,'omitnan' );
sem  = std( list )./sqrt(length(list));
end