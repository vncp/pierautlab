function [ MAT, MAT4] = pvalmat( comparison )
%pvalmat takes the output comparison matrix of the multcompare function
%(specifically the 2 first columns that state the groups compared and the
%6th column with p-values) and make a triangular matrix MAT with the groups
%compared as the X and Y axes, and the p-values as each corresponding
%elements in the lower triangle, zeros on diagonal and upper triangle. MAT4 is the same as MAT but with 0 when p>0.1, 1 when
%0.1>p>0.05, 2 when 0.01<p<0.05, 3 when 0.001<p<0.01, 4 when p<0.001.
%
%
% WARNING: only works on version of Matlab 2014b or posterior. (before, the
% multcompare doesn't have the 6th column containing p-values, and you have
% to compute them yourself)
%
%
%Example of how to plot a figure where MAT4 is displayed with 
% 0=white [1 1 1], 1=grey [0.5 0.5 0.5], 2=yellow [1 1 0], 3=orange [1 0.5 0], 4=red [1 0 0] (i.e all colors correspond to significant p-values) 
%
    %Anova to compare GC and other dataset at same input corr
    % FSGCOutputR = cat(1,Between10msCat, BetweenCat); %Y = LIST OF DATA POINTS
    % FSGCInputR = cat(1,X2,X2gc); %group1 (list of input corr for each data point)
    % 
    % for k=1:length(Between10msCat)
    %     BlockGroupFS{k,1}= 'FS';
    % end
    % for l=1:length(BetweenCat)
    %     BlockGroupGC{l,1}='GC';
    % end
    % BlockGroupFSGC = cat(1,BlockGroupFS,BlockGroupGC); %group2 (cell array)
    % G_FSGC{1} = FSGCInputR;
    % G_FSGC{2}= BlockGroupFSGC;
    % 
    % [P,ANOVATAB,STATS]=anovan( FSGCOutputR, G_FSGC, 'model', 'interaction'); % also plots a table with pvalues of 2 variables and interactions of variables
%     CompFSGC = multcompare(STATS, 'dimension', [1,2], 'display', 'off'); %dimensions allows to compare between groups   
%     [MAT_FSGC, MAT4_FSGC ] = pvalmat(CompFSGC);
%
% 
% figure
% subplot(1,2,1) % plot matrix of pvalues.
% 
%     linesx = [ [0 10]'    [0 10]'      [0 10]'      [0 10]' [0 10]'  [0 10]'  [0 10]' [0 10]'  [0 10]'  [0 10]' ];
%     linesy = [ [1 1]'  [2 2]'    [3 3]'    [4 4]' [5 5]' [6 6]' [7 7]' [8 8]' [9 9]' [10 10]' ] + 0.5;
% 
%     ytickpos    = 1:5;
%     yticklabs   = {'GC-0.56','GC-0.65', 'GC-0.74', 'GC-0.84', 'GC-0.88'};
%     xticklabs   = {'FS-0.56','FS-0.65', 'FS-0.74', 'FS-0.84', 'FS-0.88'};
% 
%     PVALcmap = [1 1 1; 0.7 0.7 0.7; 1 1 0; 1 0.5 0; 1 0 0]; % white, light grey, yellow, orange, red 
%     MAT4zoom = MAT4_FSGC(6:10,1:5); %take only the comparaison between GC (rows->yaxis) and FS (columns->xaxis)
%     MAT4zoom = [MAT4zoom, zeros(size(MAT4zoom,1),1), zeros(size(MAT4zoom,1),1)]; %add 2 columns. The last column will contain the colormap, in order to plot the right colors on the rest of the matrix when there is only 1 or 2 values in it. 
%     MAT4zoom(2,end) = 1;
%     MAT4zoom(3,end) = 2;
%     MAT4zoom(4,end) = 3;
%     MAT4zoom(5,end) = 4;
%     imagesc(MAT4zoom); hold on
%     xlim([0.5,5.5]); % exclude the meaningless columns (that were just here to have the right colormap
%     set(gca, 'ytick', ytickpos, 'yticklabel', yticklabs) 
%     set(gca, 'xtick', ytickpos, 'xticklabel', xticklabs)
%     
%     ylabel('GC-InputR', 'rot', 90)
%     xlabel('FS-InputR', 'rot', 0)
%     plot(2.*linesx, linesy, 'k-', 'linewidth', 1)
%     plot(linesy, 2.*linesx, 'k-', 'linewidth', 1)
%     axis square
%     colormap(gca, PVALcmap);
%     titlelabelCol = {'FS vs GC on same input corr: Anova, Post-hoc Tukey-Kramer p-values:'};
%     title(titlelabelCol)
%     axbox1 = get(gca, 'PlotBoxAspectRatio'); %to have 2 subplots of same size, I get the dimensions of the first subplot axes
% 
% subplot(1,2,2); % plot colormap with legend
% 
%     imagesc(MAT4zoom(:,end));
%     colormap(gca, PVALcmap);
%     ytickpos2    = 1:5;
%     yticklabs2   = {'ns <=> p>0.1', 'borderline <=> 0.05<p<0.1', '* <=> 0.01<p<0.05', '** <=> 0.001<p<0.01', '*** <=> p<0.001'};
%     set(gca, 'ytick', ytickpos2, 'yticklabel', yticklabs2);
%     axbox2 = get(gca, 'PlotBoxAspectRatio');                     
%     set(gca, 'PlotBoxAspectRatio', [axbox1(1), axbox2(2), axbox2(3)]) ; % set the second plot axes to same x dim as first subplot


%Antoine Madar, December 2015

NumOfGP = max(comparison(:,2));
MAT = zeros(NumOfGP, NumOfGP);
MAT4 = MAT;

for n = 1:length(comparison(:,1))
    pdummy = comparison(n,6);
    MAT(comparison(n,2),comparison(n,1)) = comparison(n,6);
        if pdummy>0.1
            MAT4(comparison(n,2),comparison(n,1)) = 0;
        elseif pdummy<0.1 && pdummy>0.05
            MAT4(comparison(n,2),comparison(n,1)) = 1;
        elseif pdummy<0.05 && pdummy>0.01
            MAT4(comparison(n,2),comparison(n,1)) = 2;
        elseif pdummy<0.01 && pdummy>0.001
            MAT4(comparison(n,2),comparison(n,1)) = 3;
        elseif pdummy<0.001
            MAT4(comparison(n,2),comparison(n,1)) = 4;
        end
end








