function [EST, p_Rb, Fstat_Rb] = PolFitandFtest (Npol, X2, Between10msCat, X2gc, BetweenCat, displayflag) 
% Ftest based on a polynomial regression
%
% Arguments:
% - Npol: degree of the polynomial regression. 1 for a linear fit, 2 for a parabole, etc...
% - X2: column-vector containing the x-axis values of dataset for all points
% - Between10msCat: column-vector containing the y-axis values of dataset for all points
% - X2gc: column-vector containing the x-axis  values of CONTROL dataset for all points
% - BetweenCat: column-vector containing the y-axis values of CONTROL dataset for all points
% - displayflag: 0 or 1. 1 plots the F-test result from the function "Fdist_display" 
%
% Outputs:
% - EST is a structure with 3 fields "ctrl", "other" and "pool". It contains the y-values of the estimated parameters of the fitted polynomial for each dataset and the pooled dataset. 
%   To plot, use "plot(Absc, polyval(EST.ctrl, Absc),'k');" where Absc can be, for ex, Absc = 0:0.01:1;  
% - p_Rb: pvalue of the Ftest
% - Fstat_Rb: Fvalue of the Ftest
%
% Antoine Madar, October 2017

%pool of Btw Corr for GC and data set to compare, to check if they are from same distrib
PoolX = cat(1,X2,X2gc);
PoolY = cat(1, Between10msCat, BetweenCat); 

% 2nd degree polynomial fitting using polyfit

    % Rb : polyfit for data set to compare to GC
ESTbtw = polyfit( X2, Between10msCat, Npol);%ESTbtw contains the parameters (coeffs) of the polynomial. S contains degrees of freedom S.df and the norm of residual sum of squares S.normr 
SSbtw = sum((Between10msCat - polyval(ESTbtw,X2)).^2); % residual sum of squares 
dfbtw = length(X2) - length(ESTbtw); % degrees of freedom = Number of data points - number of parameters

    % Rb: polyfit for GC
ESTgcB = polyfit( X2gc, BetweenCat, Npol);%ESTbtw contains the 3 parameters (coeffs) of the polynomial. S contains degrees of freedom S.df and the norm of residual sum of squares S.normr 
SSbtwGC = sum((BetweenCat - polyval(ESTgcB,X2gc)).^2); % residual sum of squares 
dfbtwGC = length(X2gc) - length(ESTgcB);

    % Rb pool polyfit 
ESTpool = polyfit( PoolX, PoolY, Npol);
SSpool = sum((PoolY - polyval(ESTpool,PoolX)).^2); 
dfpool = length(PoolY) - length(ESTpool);

EST.ctrl = ESTgcB;
EST.other = ESTbtw;
EST.pool = ESTpool;

% F-test to assess significance of separate fittings versus pooled (combined data sets).
    
  % Rb 
SSsep = SSbtw + SSbtwGC;
dfsep = dfbtw + dfbtwGC;

Fstat_Rb = ((SSpool - SSsep)/(dfpool - dfsep))/(SSsep/dfsep);

df1 = dfpool-dfsep;
df2 = dfsep;
Pdes = 0.05;

%perform test (and plot result if displayflag == 1)
[F, pF, p_Rb] = Ftest_Antoine(df1, df2, Fstat_Rb, Pdes, displayflag); 
