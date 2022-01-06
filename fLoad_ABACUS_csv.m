%fLoad_ABACUS_csv.m
%An .m file that imports ABACUS-formatted data into MATLAB. 
%It attributes each column in the header (separated by commas) to its affiliated data column. 
%Units & ancillary information, which is in parentheses after the header name in the first row, is ignored. 
%The contents of the first row are printed in the command line window for bookkeeping.
%Use this file to read ABACUS-formatted .csv data files into MATLAB.


%%List name of file. The program assumes that it is .csv
%name='AB_AWS_Suggestion'

%%Read in numeric data (i.e. everything except the column headers in the 1st row) and label it matrix "x";
cmd1=sprintf('x = csvread(''%s.csv'',1,0);',name);eval(cmd1);

%%Open file to read in 1st line (column headers).
cmd2=sprintf('fid = fopen(''%s.csv'',''r'');',name);eval(cmd2);

%%Read in first line of file and name it Line1. 
while fid==3
    Line1 = fgetl(fid);
    if ~ischar(Line1),   break,   end

    %%Display 1st line, including units
    disp(Line1)

    %%Trick MATLAB into exiting loop after first line
    fid=fid+1;
end
fid=fid-1;

%%Close file.
fclose(fid);

%%Associate columns in matrix "x" with column headers in "Line1"
if Line1(length(Line1))==','
    Line1=Line1(1:length(Line1)-1);
end

Head=[];col=1;
for j=1:length(Line1)
%%Temporarily store data in Line1 until a comma is reached.
    temp=Line1(j);
    if temp~=','
%%Save temporary data as a Header until a comma is reached.
    Head=[Head temp];
    else
        
%%If units (in parentheses) are present, remove this information for the purpose of naming the column.
n=1:length(Head);
z=n(Head=='(');
if length(z)>0
Head=Head(1:z-1);
end

%%Associate each column in x (except the last) with the appropriate Header.
    cmd = sprintf('%s = x(:,%d);', Head,col); eval(cmd);
    Head=[];
    
    %%move on to next column
    col=col+1;
    end
end

%%Name the last column, which is not followed by a comma.
n=1:length(Head);
z=n(Head=='(');
if length(z)>0
Head=Head(1:z-1);
end
cmd = sprintf('%s = x(:,%d);', Head, col);eval(cmd);

clear Head Line1 j temp fid col z cmd cmd1 cmd2 n ans x

%% Save output as a .mat file if required.
%cmd=sprintf('save %s',name);eval(cmd);
%clear cmd name
