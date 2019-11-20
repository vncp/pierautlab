function S = read_excel(filespec);

if nargin < 1
	[S.file, S.path] = uigetfile('*.*');
else
	seps = findstr(filespec, filesep);
	S.file = filespec(seps(end)+1:end)
	S.path = filespec(1:seps(end));
end

disp('Atempting to read excel file...');

comptype = computer;
switch comptype(1:3)
	case 'MAC'
		disp('This is a MAC')
	S.file(findstr(S.file, ':')) = '/';
	S.path(findstr(S.path, ':')) = '/';
end

disp(['PATH = ' S.path]);
disp(['FILE = ' S.file]);

[fid, message] = fopen([S.path S.file], 'r', 'b');
if fid == -1
	disp(message);
end

clear seps;

S.fileType = fread(fid, 4, 'char');
S.fileTypee = char(S.fileType);

switch S.fileType
	case 'xlsx'
		disp('This is an excel file')
	otherwise
		error(['Sorry: unrecognized file type. Should be *.xlsx but file type is: ' fileType ]);
end

%Temporary struct plateholder variable
T = readtable(S.file);
T.rows = Length(T);
T.cols = size(T, 2);
T.time = T{1};
for datacol = 2:T.cols;
	T.data = T{datacol};
end

S.columnData = T.data;
S.numGroups = T.cols;
S.timeGroup = 1;

fclose(fid)
disp('Done";')