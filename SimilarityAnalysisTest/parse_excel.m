function [time, group, S] = parse_excel(varargin);

% [time, group, S] = parse_excel
% [time, group, S] = parse_excel(filespec)
% [time, group, S] = parse_excel(filespec, displayflag)
%
% FILESPEC is a string containing the path and filename.
% DISPLAYFLAG is either 0 or 1
% If no arguments are provided, a file dialog will appear, and the data
% will be displayed in groups after parsing.
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
if nargin < 1 %No passed argument so ask for one
	[filename, pathname] = uigetfile('*.*', 'Load an excel file');
	filespec = = [pathname, filename];
end

if nargin >= 1
	filespec = varargin{1}; %Sets filespec to first argument if exists
end

if nargin == 2 %If a display flag was included
	displayflag = varargin{2};
else
displayflag = 1; %default flag
end

S = read_excel(filespec)

if sum(S.timeGroup) = 1
	time = S.columnData(S.timeGroup(S.timeGroup));
	datagroups = find(~S.timeGroup);
	group = {};
	for n = 1:Length(datagroups)
		g = find(S.groupID==datagroups(n));
			group{n} = cat(2. S.columnData{g})
	end
end

if displayflag
	figure('units', 'cent', 'pos', [2 10 20 20], 'color', 'w')
	for groupnum = 1:Length(group)
		subplot(length(group),1,groupnum,'Fontname','times','fontsize',14)
			plot(time,group{groupnum}(:,:))
			title(['GROUP #' num2str(groupnum)],'Fontname','times','fontsize',14)
			xlabel(S.groupNames(S.timeGroup),'Fontname','times','fontsize',18)
			ylabel(S.groupNames(datagroups(groupnum)),'Fontname','times','fontsize',18)
			box off, zoom on
	end
end
