% take a symetric matrix X and return the list of unique values in the upper (or lower, doesn't matter since X is symetric) triangle, excluding the diagonal. 
%It also removes NaN values
%
% Antoine Madar, April 2017

function [list, average, sem] = SymMat2List(X) 
    
nanmask = ones(size(X));
nanmask = triu(nanmask,1)./triu(nanmask,1);
list = nanmask.*triu(X, 1);
average = mean( list,'omitnan' );
sem  = std( list )./sqrt(length(list));
list = list(~isnan(list));

end