function [time, group, S] = parse_axograph(varargin);

% [time, group, S] = parse_axograph
% [time, group, S] = parse_axograph(filespec)
% [time, group, S] = parse_axograph(filespec, displayflag)
% 
% 
% Read an AxographX file into a struct S using
% read_axograph(). Then parse it into appropriate groups.
% 
% FILESPEC is a string containing the path and filename.
% 
% DISPLAYFLAG is either 0 or 1
% 
% If no arguments are provided, a file dialog will appear, and the data
% will be displayed in groups after parsing.
% 
% NOTE: Grouping is based on the original COLUMN NAMES in the Axograph
% file, which may be different than the grouping displayed in the Axograph window.
% 
% You then have a double array called 'time', and a cell array called 'group'. 
% Each element of 'group' is a double array containing all the sweeps in that group.
% Each individual trace can be accessed by calling 'group{group#}(rowindx, colindx)'.
% NOTE THE USE OF {} TO INDEX THE CELL ELEMENT, BUT USE OF () TO INDEX THE DOUBLE ARRAY.
% For example, these two statements are equivalent, and yield ALL columns from group #1.
%     x = group{1}
%     x = group{1}(:,:)
% In contrast, this statement yields only columns #3,5 and 7 from group #1:
%     x = group{1}(:,[3,5,7])
% This statement yields only ROWS 1-10 from group #1:
%     x = group{1}(1:10, :)
% etc...
% 
% MJones, 2010


% File selection dialog, if unspecified
if nargin < 1
    [filename, pathname] = uigetfile('*.*', 'Load an Axograph file:');
    filespec = [pathname, filename];
end  

if nargin >= 1
    filespec = varargin{1};
end   

if nargin == 2
   displayflag = varargin{2};
else
   displayflag = 1;
end



% Call the function that reads the data (read_axograph()) and load all the data into variable S (of type STRUCT). 
% Look at Command Window to see the contents of STRUCT S
S = read_axograph(filespec)


if sum(S.timeGroup) == 1                            % If there is only one unique group holding 'Time (s)'
    time = S.columnData{S.timeGroup( S.timeGroup )};	% Then put that group into the double array "time"
    datagroups = find(~S.timeGroup);                % Find the ones that *aren't* the time array
    
    group = {};                                 % Initialize the cell array to hold the separate groups
    for n = 1:length(datagroups)                % Loop through the data groups
        g = find( S.groupID==datagroups(n) );   % Find all the traces that match the current data groupID
        group{n} = cat( 2, S.columnData{ g } ); % Concatenate these traces together and store in the cell array "group"
    end
end



% Plot each group in a separate subplot, with the same labels as AxoGraph etc:
if displayflag
    figure('units', 'cent', 'pos', [2 10 20 20], 'color', 'w')
    for groupnum = 1:length(group)
        subplot( length(group), 1, groupnum, 'Fontname', 'times', 'fontsize', 14)
            plot( time, group{groupnum}(:,:) )
            title( ['GROUP #' num2str(groupnum)] , 'Fontname', 'times', 'fontsize', 18 )
            xlabel( S.groupNames( S.timeGroup ), 'Fontname', 'times', 'fontsize', 18 )
            ylabel( S.groupNames( datagroups(groupnum) ) , 'Fontname', 'times', 'fontsize', 18 )   
            box off; zoom on
    end
end






