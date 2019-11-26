function S = read_axograph(filespec);

% ------------------------------------------------------------------------
% Read and parse an AxoGraph binary file.
%
% S is a big struct with header info as well as data.
% 
% Versions:
% 12-12-2000 - Original program by Mathew V. Jones.
% 04-01-2006 - Speed significantly increased by removing unnecessary calls
%               to fseek() - Marcel Paz Goldschen-Ohm
% 02-07-2007 - Updated to read Axograph X file format - MVJones
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% FOR ORIGINAL FILE FORMAT DEFINITIONS, INFO and and SAMPLE I/O CODE, see
% AxoGraph_ReadWrite.h in the Developer Documentation folder of the AxoGraph installation.
% NOTE: There are some errors in those descriptions, so see also files in "Axograph X CORRECT I/O routines"
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% AxoGraph 4.6 Standard Format
% ============================
%
% Header:
% -------
%   Byte	Type		Contents
%   ----    ----        --------
%   0		OSType		AxoGraph file header identifier = 'AxGr'.
%   4		Integer		AxoGraph file format version number = 1.
%   6		Integer		Number of data columns to follow.
% 
% Each column:
% ------------
%   Byte	Type		Contents
%   ----    ----        --------
%   0		Longint		Number of data points in the column.
%   4		String[80]	Column title.
%   84		Real*4		1st data point.
%   88		Real*4		2nd data point.
%   ..		..			....
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% AxoGraph 4.6 Digitized Format
% =============================
%
% Header:
% -------
%   Byte	Type		Contents
%   ----    ----        --------
%   0		OSType		AxoGraph file header identifier = 'AxGr'.
%   4		Integer		AxoGraph file format version number = 2.
%   6		Integer		Number of columns to follow.
% 
% First column (X data):
% ----------------------
%   Byte	Type		Contents
%   ----    ----        --------
%   0		Longint		Number of data points in each of the subsequent columns.
%   4		String[80]	Column title.
%   84		Real*4		Sample interval.
% 
% Each subsequent column (Y data):
% --------------------------------
%   Does NOT start on byte 96 as you might suspect, but actually starts
%   on byte 100. Why, I have no idea.
%
%   Byte	Type		Contents
%   ----    ----        --------
%   0		Longint		Number of data points in the column.
%   4		String[80]	Column title.
%   84		Real*4		Scale factor.
%   88		Integer*2	1st Data point.
%   90		Integer*2	2nd Data point.
%   ..		...			....
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% AxoGraph X Data File Format
% ===================================
% 
% Header
% ------
% Byte	Type		Contents
% 0		char[4]		AxoGraph file header identifier = 'AxGx' - same as filename extension
% 4		long		AxoGraph X file format ID = a number between 3 (earliest version) and 6 (latest version)
% 8		long		Number of columns to follow
% 
% 
% Each column
% ----------------------
% Byte	Type		Contents
% 0		long		Number of points in the column ( columnPoints )
% 4		long		Column type 
% 8		long		Length of column title in bytes (Unicode - 2 bytes per character)
% 12		char*		Column title (Unicode 2 byte per char) - S.I. units should be in brackets e.g. 'Current (pA)'
% ??		??			Byte offset depends on length of column title string. 
% ..		...			Numeric type and layout depend on the column type
% ..		...			etc.
% 
% Six column types are supported...
% 	4: short 
% 	5: long
% 	6: float
% 	7: double
% 	9: 'series'
% 	10: 'scaled short'
% 
% In the first four column types, data is stored as a simple array of the corresponding type.
% The 'scaled short' column type stores data as a 'double' scaling factor and offset, and a 'short' array.
% The 'series' column type stores data as a 'double' first value and a 'double' increment.
% ------------------------------------------------------------------------

%#########################################################################
% Open in big-endian (Mac-compatible) format, read only.
%#########################################################################

disp('Loading Axograph binary file...');

if nargin < 1
    [S.file, S.path] = uigetfile('*.*');
else
	seps = findstr(filespec, '/');
    disp(['DEBUGGER: ' filesep])
	S.file = filespec(seps(end)+1:end);
	S.path = filespec(1:seps(end));
end

% NOTE: The new Matlab does not understand ":" as a correct file or path separator (idiots). 
% Need to convert all ":" to "/". 
% Note: This creates problems on PCs (idiots).
% Therefore need to handle them differently.
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
if fid == -1; disp(message); end

clear seps;

% Make sure it's an Axograph file.

S.fileType = fread(fid, 4, 'char');
S.fileType = char(S.fileType)';

switch S.fileType
	case 'AxGr'
        disp('This is an AxoGraph 4.x file'); 
        S.AxGrVersion = '4.x'; 
        idformat = 'int16';
	case 'axgx'
        disp('This is an AxoGraph X file (file type = axgx)'); 
        S.AxGrVersion = 'X'; 
        idformat = 'int32';
	case 'axgx'
        disp('This is an AxoGraph X file (file type = axgd)'); 
        S.AxGrVersion = 'X'; 
        idformat = 'int32';
    otherwise 
        error(['Sorry: This is an unrecognized file type. Should be "AxGr" or "axgx" (or "axgd"?). Instead, byte 4 says it is: ' fileType ]);

end 
    
 



% Get file format.

S.fileFormat = fread(fid, 1, idformat);
disp(['File Format ID# = ' num2str(S.fileFormat)])

if S.fileFormat == 1
    S.FileDesc = 'Axograph 4.x Standard File'; 
	disp(['FILE FORMAT = ' S.FileDesc]);
    numcolsbyteoff = 6;
elseif	S.fileFormat == 2
    S.FileDesc = 'Axograph 4.x Digitized File'; 
	disp(['FILE FORMAT = ' S.FileDesc]);
    numcolsbyteoff = 6;
elseif	S.fileFormat >= 3 & S.fileFormat <= 6
    S.FileDesc = 'Axograph X File'; 
	disp(['FILE FORMAT = ' S.FileDesc]);
    numcolsbyteoff = 8;
else
% 	sound(sin(1:1024))
	error(['Sorry: ' num2str(S.fileFormat) ' is an unrecognized file format. Should be 1 (Standard) or 2 (Digitized).']);
end






%#########################################################################
% Read data.
%#########################################################################

fseek(fid, numcolsbyteoff, -1);
S.numColumns = fread(fid, 1, idformat);

disp(['NUMBER OF COLUMNS = ' num2str(S.numColumns)]);

progressIndicator = waitbar(0, ['Reading ' num2str(S.numColumns) ' sweeps...']);

% AxoGraph 4.x Standard format.
if S.fileFormat == 1
    for aColumn = 1 : S.numColumns
        S.columnPts(aColumn)            = fread(fid, 1, 'int32');
        S.columnType(aColumn)           = 6;
		S.columnTitle{aColumn}          = char( fread(fid, 80, 'char')');
        S.columnTitleLength(aColumn)    = length(S.columnTitle{aColumn});          
        S.columnTypeName{aColumn}       = 'float';
        S.columnTypeFormat{aColumn}     = 'float32';
        S.columnScale(aColumn)          = 1;          
        S.columnData{aColumn}           = fread(fid, S.columnPts(aColumn), 'real*4');
        waitbar(aColumn / S.numColumns)
    end

% AxoGraph 4.x Digitized format.
elseif S.fileFormat == 2
    aColumn = 1;
		S.columnPts(aColumn)            = fread(fid, 1, 'int32');           % display(['Num Pts = ' num2str(S.columnPts(aColumn))]);
        S.columnType(aColumn)           = 6;
        S.columnTitle{aColumn}          = char( fread(fid, 80, 'uchar=>uchar')');  % display(['Title = ' S.columnTitle{aColumn}]);
        S.columnTitleLength(aColumn)    = length(S.columnTitle{aColumn});          
        S.columnTypeName{aColumn}       = 'float';
        S.columnTypeFormat{aColumn}     = 'float32';
        S.columnScale(aColumn)          = fread(fid, 1, 'real*4');          % display(['Sample Interval = ' num2str(S.columnScale(aColumn))]);
        S.columnData{aColumn}           = S.columnScale(aColumn) .* [1:S.columnPts(aColumn)]';
		waitbar(aColumn / S.numColumns)

    fseek(fid, 100, 'bof');
    for aColumn = 2 : S.numColumns
		S.columnPts(aColumn)            = fread(fid, 1, 'int32');           % display(['Num Pts = ' num2str(S.columnPts(aColumn))]);
        S.columnType(aColumn)           = 6;
        S.columnTitle{aColumn}          = char(fread(fid, 80, 'uchar=>uchar')');  % display(['Title = ' S.columnTitle{aColumn}]);
        S.columnTitleLength(aColumn)    = length(S.columnTitle{aColumn});          
        S.columnScale(aColumn)          = fread(fid, 1, 'real*4');          % display(['Scale Factor = ' num2str(S.columnScale(aColumn))]);
        S.columnTypeName{aColumn}       = 'float';
        S.columnTypeFormat{aColumn}     = 'float32';
        S.columnData{aColumn}           = S.columnScale(aColumn) .* fread(fid, S.columnPts(aColumn), 'integer*2');
		waitbar(aColumn / S.numColumns)
    end
   
% AxoGraph X format
elseif	S.fileFormat >= 3 & S.fileFormat <= 6
    for aColumn = 1 : S.numColumns
		S.columnPts(aColumn)            = fread(fid, 1, 'int32');           % display(['Num Pts = ' num2str(S.columnPts(aCoumn))]);
		S.columnType(aColumn)           = fread(fid, 1, 'int32')';   
		S.columnTitleLength(aColumn)    = fread(fid, 1, 'int32');          % In bytes, 2 bytes per character
        str                             = char( fread(fid, S.columnTitleLength(aColumn), 'char') );          % In bytes, 2 bytes per character
        S.columnTitle{aColumn}          = str(2:2:end)';                                                    % Get rid of dopey blanks
        switch S.columnType(aColumn)
%         Read data according to its format and length
%             NOTE: I THINK this will work for all column types, but have
%             only tested type 9 and 10
            case 4
                S.columnTypeName{aColumn}	= 'short';
                S.columnTypeFormat{aColumn}	= 'int16';
                S.columnScale(aColumn)      = 1;          
                S.columnData{aColumn}    	= fread(fid, S.columnPts(aColumn), S.columnTypeFormat{aColumn});
            case 5
                S.columnTypeName{aColumn}	= 'long';
                S.columnTypeFormat{aColumn}	= 'int32';
                 S.columnScale(aColumn)   	= 1;          
                 S.columnData{aColumn}    	= fread(fid, S.columnPts(aColumn), S.columnTypeFormat{aColumn});
            case 6
                S.columnTypeName{aColumn}	= 'float';
                S.columnTypeFormat{aColumn}	= 'float';
                S.columnScale(aColumn)      = 1;          
                S.columnData{aColumn}    	= fread(fid, S.columnPts(aColumn), S.columnTypeFormat{aColumn});
            case 7
                S.columnTypeName{aColumn}   = 'double';
                S.columnTypeFormat{aColumn} = 'real*8';
                S.columnScale(aColumn)      = 1;          
                S.columnData{aColumn}    	= fread(fid, S.columnPts(aColumn), S.columnTypeFormat{aColumn});
            case 9
                S.columnTypeName{aColumn}   = 'series';
                S.columnTypeFormat{aColumn} = 'real*8';
                seriesinfo                  = fread(fid, 2, S.columnTypeFormat{aColumn});
                S.columnScale(aColumn)      = 1;          
                S.columnData{aColumn}       = [ seriesinfo(1) + seriesinfo(2).*[1:S.columnPts(aColumn)] ]';
            case 10
                S.columnTypeName{aColumn}   = 'scaled short';
                S.columnTypeFormat{aColumn} = {'real*8'; 'real*8'; 'int16'};
                scaling                     = fread(fid, 1, S.columnTypeFormat{aColumn}{1});
                offset                      = fread(fid, 1, S.columnTypeFormat{aColumn}{2});
                data                        = fread(fid, S.columnPts(aColumn), S.columnTypeFormat{aColumn}{3});
                S.columnScale(aColumn)      = 1;          
                S.columnData{aColumn}       = scaling.*data + offset;
        end

		waitbar(aColumn / S.numColumns)

    end
    
    
    
    
end


%#########################################################################
% Detect Groups by unique column titles.
%#########################################################################
[S.groupNames, I, S.groupID] = unique( S.columnTitle );
S.numGroups = length( S.groupNames );
S.timeGroup = strcmp( S.groupNames, 'Time (s)');
if sum(S.timeGroup) < 1 % Check whether the time group has a different name
    S.timeGroup = strcmp( S.groupNames, 'Time, sec');
end
if sum(S.timeGroup) < 1 % Check whether the time group has a different name
    S.timeGroup = strcmp( S.groupNames, 'Time, s');
end

%#########################################################################
% Clean up.
%#########################################################################

close(progressIndicator);
fclose(fid);
disp('... Done.');

%#########################################################################











