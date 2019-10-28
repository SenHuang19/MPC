function ExtractSimData_sentest(test_season)

%----- Generates data for the model validation -----------
%- Usage: test_season='fall_'; %"season_" is default for one loc
                               % season=fall,summer,winter
                               %"season_loc_" is other loc like loc="sf"
%-        ExtractSimData(test_season)

%-----Bottom floor-------------------------------------------

% Adding folders path
addpath('C:\Users\vchinde\Documents\Box Sync\SenJournalPaper\coeff\coeff\Bot_floor\raw_data');

fid = fopen(strcat(test_season,'zone4-1.csv'));
databot1= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Botfloor1 = ExtractZone_info(databot1);

fid = fopen(strcat(test_season,'zone4-2.csv'));
databot2= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Botfloor2 = ExtractZone_info(databot2);

fid = fopen(strcat(test_season,'zone4-3.csv'));
databot3= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Botfloor3 = ExtractZone_info(databot3);

fid = fopen(strcat(test_season,'zone4-4.csv'));
databot4= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Botfloor4 = ExtractZone_info(databot4);

fid = fopen(strcat(test_season,'zone4-5.csv'));
databot5= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Botfloor5 = ExtractZone_info(databot5);


%-----------Mid floor --------------------------------------------------
addpath('C:\Users\vchinde\Documents\Box Sync\SenJournalPaper\coeff\coeff\Mid_floor\raw_data');

fid = fopen(strcat(test_season,'zone4-1.csv'));
datamid1= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Midfloor1 = ExtractZone_info(datamid1);

fid = fopen(strcat(test_season,'zone4-2.csv'));
datamid2= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Midfloor2 = ExtractZone_info(datamid2);

fid = fopen(strcat(test_season,'zone4-3.csv'));
datamid3= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Midfloor3 = ExtractZone_info(datamid3);

fid = fopen(strcat(test_season,'zone4-4.csv'));
datamid4= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Midfloor4 = ExtractZone_info(datamid4);

fid = fopen(strcat(test_season,'zone4-5.csv'));
datamid5= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Midfloor5 = ExtractZone_info(datamid5);

%----------- Top floor -------------------------------------------------
addpath('C:\Users\vchinde\Documents\Box Sync\SenJournalPaper\coeff\coeff\Top_floor\raw_data');

fid = fopen(strcat(test_season,'zone4-1.csv'));
datatop1= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Topfloor1 = ExtractZone_info(datatop1);

fid = fopen(strcat(test_season,'zone4-2.csv'));
datatop2= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Topfloor2 = ExtractZone_info(datatop2);

fid = fopen(strcat(test_season,'zone4-3.csv'));
datatop3= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Topfloor3 = ExtractZone_info(datatop3);

fid = fopen(strcat(test_season,'zone4-4.csv'));
datatop4= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Topfloor4 = ExtractZone_info(datatop4);

fid = fopen(strcat(test_season,'zone4-5.csv'));
datatop5= textscan(fid,'%C%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
Topfloor5 = ExtractZone_info(datatop5);

Tamb=Topfloor5.Tamb;


QintTotal=[Botfloor1.Qint Botfloor2.Qint Botfloor3.Qint Botfloor4.Qint Botfloor5.Qint ...
    Midfloor1.Qint Midfloor2.Qint Midfloor3.Qint Midfloor4.Qint Midfloor5.Qint...
   Topfloor1.Qint Topfloor2.Qint Topfloor3.Qint Topfloor4.Qint Topfloor5.Qint];

filename=strcat(test_season,'Ext_data');
save(filename,'Tamb','QintTotal')

CTtemp=[Botfloor1.t1, Botfloor2.t1,Botfloor3.t1,Botfloor4.t1,Botfloor5.t1,...
    Midfloor1.t1,Midfloor2.t1,Midfloor3.t1,Midfloor4.t1,Midfloor5.t1,...
    Topfloor1.t1,Topfloor2.t1,Topfloor3.t1,Topfloor4.t1,Topfloor5.t1];


m_i= [Botfloor1.mf Botfloor2.mf Botfloor3.mf Botfloor4.mf Botfloor5.mf  ...
    Midfloor1.mf Midfloor2.mf Midfloor3.mf Midfloor4.mf Midfloor5.mf  ...
   Topfloor1.mf Topfloor2.mf Topfloor3.mf Topfloor4.mf Topfloor5.mf];


Prh_i= [Botfloor1.rh Botfloor2.rh Botfloor3.rh Botfloor4.rh Botfloor5.rh ...
    Midfloor1.rh Midfloor2.rh Midfloor3.rh Midfloor4.rh Midfloor5.rh...
   Topfloor1.rh Topfloor2.rh Topfloor3.rh Topfloor4.rh Topfloor5.rh];

filename1=strcat(test_season,'Model_valid');
save(filename1,'CTtemp','m_i','Prh_i','QintTotal','Tamb')
