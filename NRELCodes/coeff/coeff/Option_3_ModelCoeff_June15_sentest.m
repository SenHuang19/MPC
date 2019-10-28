clc
clearvars
close all


%% -------------Option 3------------------------


% % X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states 
% % U = [ m_i';Prh_i'];
% % D = [Tamb(i_sch);Qint(i_sch,:)';1];

tresolution='\1min';  % 5,60
timeresbot=strcat('.\Bot_floor',tresolution);   
                                    
BB=load(strcat(timeresbot,'\bot_floor.mat'));

timeresmid=strcat('.\Mid_floor',tresolution);

M=load(strcat(timeresmid,'\mid_floor.mat'));

timeresmid=strcat('.\Top_floor',tresolution);
T=load(strcat(timeresmid,'\top_floor.mat'));

% Bottom floor
 
Bf1.Summer.Office=BB.zone4_1_summer_result;
Bf2.Summer.Office=BB.zone4_2_summer_result;
Bf3.Summer.Office=BB.zone4_3_summer_result;
Bf4.Summer.Office=BB.zone4_4_summer_result;
Bf5.Summer.Office=BB.zone4_5_summer_result;

% % Middle floor

Mf1.Summer.Office=M.zone4_1_summer_result;
Mf2.Summer.Office=M.zone4_2_summer_result;
Mf3.Summer.Office=M.zone4_3_summer_result;
Mf4.Summer.Office=M.zone4_4_summer_result;
Mf5.Summer.Office=M.zone4_5_summer_result;

% % Top floor

Tf1.Summer.Office=T.zone4_1_summer_result;
Tf2.Summer.Office=T.zone4_2_summer_result;
Tf3.Summer.Office=T.zone4_3_summer_result;
Tf4.Summer.Office=T.zone4_4_summer_result;
Tf5.Summer.Office=T.zone4_5_summer_result;


%% ------------------------ A matrix
A.Summer.Office = diag([double(Bf1.Summer.Office.a2);double(Bf2.Summer.Office.a2);double(Bf3.Summer.Office.a2);double(Bf4.Summer.Office.a2);double(Bf5.Summer.Office.a2);...
                        double(Mf1.Summer.Office.a2);double(Mf2.Summer.Office.a2);double(Mf3.Summer.Office.a2);double(Mf4.Summer.Office.a2);double(Mf5.Summer.Office.a2);...
                        double(Tf1.Summer.Office.a2);double(Tf2.Summer.Office.a2);double(Tf3.Summer.Office.a2);double(Tf4.Summer.Office.a2);double(Tf5.Summer.Office.a2)]);

%% ------------------------B matrix 

%B_m_summerB=diag([double(Bf1.Summer.Office.a5),double(Bf2.Summer.Office.a5),double(Bf3.Summer.Office.a5),double(Bf4.Summer.Office.a5),double(Bf5.Summer.Office.a5),...
          %double(Mf1.Summer.Office.a5),double(Mf2.Summer.Office.a5),double(Mf3.Summer.Office.a5),double(Mf4.Summer.Office.a5),double(Mf5.Summer.Office.a5),...
          %double(Tf1.Summer.Office.a5),double(Tf2.Summer.Office.a5),double(Tf3.Summer.Office.a5),double(Tf4.Summer.Office.a5),double(Tf5.Summer.Office.a5)]);

B_m_summerB=diag([double(Bf1.Summer.Office.a3),double(Bf2.Summer.Office.a3),double(Bf3.Summer.Office.a3),double(Bf4.Summer.Office.a3),double(Bf5.Summer.Office.a3),...
          double(Mf1.Summer.Office.a3),double(Mf2.Summer.Office.a3),double(Mf3.Summer.Office.a3),double(Mf4.Summer.Office.a3),double(Mf5.Summer.Office.a3),...
          double(Tf1.Summer.Office.a3),double(Tf2.Summer.Office.a3),double(Tf3.Summer.Office.a3),double(Tf4.Summer.Office.a3),double(Tf5.Summer.Office.a3)]);
      
B_reheat_summerB= diag([double(Bf1.Summer.Office.a4),double(Bf2.Summer.Office.a4),double(Bf3.Summer.Office.a4),double(Bf4.Summer.Office.a4),double(Bf5.Summer.Office.a4),...
          double(Mf1.Summer.Office.a4),double(Mf2.Summer.Office.a4),double(Mf3.Summer.Office.a4),double(Mf4.Summer.Office.a4),double(Mf5.Summer.Office.a4),...
          double(Tf1.Summer.Office.a4),double(Tf2.Summer.Office.a4),double(Tf3.Summer.Office.a4),double(Tf4.Summer.Office.a4),double(Tf5.Summer.Office.a4)]);

B.Summer.Office=[B_m_summerB B_reheat_summerB];


%% ------------ E matrix

E_amb_sb=[double(Bf1.Summer.Office.a1);double(Bf2.Summer.Office.a1);double(Bf3.Summer.Office.a1);double(Bf4.Summer.Office.a1);double(Bf5.Summer.Office.a1);...
          double(Mf1.Summer.Office.a1);double(Mf2.Summer.Office.a1);double(Mf3.Summer.Office.a1);double(Mf4.Summer.Office.a1);double(Mf5.Summer.Office.a1);...
          double(Tf1.Summer.Office.a1);double(Tf2.Summer.Office.a1);double(Tf3.Summer.Office.a1);double(Tf4.Summer.Office.a1);double(Tf5.Summer.Office.a1)];

E_int_sb=diag([double(Bf1.Summer.Office.a5),double(Bf2.Summer.Office.a5),double(Bf3.Summer.Office.a5),double(Bf4.Summer.Office.a5),double(Bf5.Summer.Office.a5),...
          double(Mf1.Summer.Office.a5),double(Mf2.Summer.Office.a5),double(Mf3.Summer.Office.a5),double(Mf4.Summer.Office.a5),double(Mf5.Summer.Office.a5),...
          double(Tf1.Summer.Office.a5),double(Tf2.Summer.Office.a5),double(Tf3.Summer.Office.a5),double(Tf4.Summer.Office.a5),double(Tf5.Summer.Office.a5)]);

E_const_sb=[double(Bf1.Summer.Office.a0);double(Bf2.Summer.Office.a0);double(Bf3.Summer.Office.a0);double(Bf4.Summer.Office.a0);double(Bf5.Summer.Office.a0);...
          double(Mf1.Summer.Office.a0);double(Mf2.Summer.Office.a0);double(Mf3.Summer.Office.a0);double(Mf4.Summer.Office.a0);double(Mf5.Summer.Office.a0);...
          double(Tf1.Summer.Office.a0);double(Tf2.Summer.Office.a0);double(Tf3.Summer.Office.a0);double(Tf4.Summer.Office.a0);double(Tf5.Summer.Office.a0)];   
      
E.Summer.Office=[E_amb_sb E_int_sb E_const_sb]; 

 

%Xnew=A*X+B*U+E*D;

% A1.Summer.Baltimore=A.Summer.Baltimore.^5;
% 
% 
% B1.Summer.Baltimore=A.Summer.Baltimore.^4*B.Summer.Baltimore+A.Summer.Baltimore.^3*B.Summer.Baltimore+A.Summer.Baltimore.^2*B.Summer.Baltimore+A.Summer.Baltimore*B.Summer.Baltimore+B.Summer.Baltimore;
% 
% E1.Summer.Baltimore=A.Summer.Baltimore.^4*E.Summer.Baltimore+A.Summer.Baltimore.^3*E.Summer.Baltimore+A.Summer.Baltimore.^2*E.Summer.Baltimore+A.Summer.Baltimore*E.Summer.Baltimore+E.Summer.Baltimore;

%% Chiller Coefficients

chillerresol=strcat('.\chiller',tresolution);

load(strcat(chillerresol,'\chiller.mat'));

%load('.\chiller\1min\chiller.mat')

Chiller.Summer.Office = summer_chiller_result;

%% Boiler Coefficients


boilerresol=strcat('.\boiler',tresolution);

load(strcat(boilerresol,'\boiler.mat'));

%load('.\boiler\1min\boiler.mat')

Boiler.Summer.Office = summerboiler_result;
%% Fan Coefficients

load('.\fan\fan.mat')

Fan.Summer.Office.Fan1 = summer_fan1_result;
Fan.Summer.Office.Fan2 = summer_fan2_result;
Fan.Summer.Office.Fan3 = summer_fan3_result;
Fan.Summer.Office.Fan4 = summer_fan4_result;


save('Simparamdata1min','Chiller','Fan','Boiler','A','B','E')
