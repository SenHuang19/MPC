clc
clearvars
% close all



                                    


Option_3_ModelCoeff_June15_sentest

load('summer_Model_valid.mat')
%%-----------------Zone temperatures--------------------------------%%

% X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states
%U = [ m_i';Prh_i'];
%D = [Tamb(i_sch);Qint(i_sch,:)'];

m=1440*4;
n=1440*5;

QintTotal=QintTotal(m:n,:);
m_i=m_i(m:n,:);
%m_i(:,6:10)=m_i(:,6:10)./10;
Prh_i=Prh_i(m:n,:);
%Prh_i(:,6:10)=Prh_i(:,6:10)./10;
Tamb=Tamb(m:n);
Tpredict(:,1)=CTtemp(m,:)';
%% zone temperature simulation

for i_sch=1:length(QintTotal)
    
    U = [m_i(i_sch,:)';Prh_i(i_sch,:)'];
    D = [Tamb(i_sch);QintTotal(i_sch,:)';1];
    
    
        Tpredict(:,i_sch+1)=A.Summer.Office*Tpredict(:,i_sch)+B.Summer.Office*U+E.Summer.Office*D;
        %Tpredict(:,i_sch+1)=A.Summer.Office*CTtemp(i_sch,:)'+B.Summer.Office*U+E.Summer.Office*D; %one-step ahead prediction
    
    
end

%--------------Basement and Bottom floor zones ------------------------
figure(1)
plot(Tpredict(1,:),'r')
hold on
plot(CTtemp(m:n,1),'b')
xlabel('time')
ylabel('Botfloor1')
legend('Predict','actual')
figure(2)
plot(Tpredict(2,:),'r')
hold on
plot(CTtemp(m:n,2),'b')
xlabel('time')
ylabel('BotFloor2 ')
%legend('actual','setpoint')
figure(3)
plot(Tpredict(3,:),'r')
hold on
plot(CTtemp(m:n,3),'b')
xlabel('time')
ylabel('BotFloor3 ')
figure(4)
plot(Tpredict(4,:),'r')
hold on
plot(CTtemp(m:n,4),'b')
xlabel('time')
ylabel('BotFloor4 ')
%legend('actual','setpoint')
figure(5)
plot(Tpredict(5,:),'r')
hold on
plot(CTtemp(m:n,5),'b')
xlabel('time')
ylabel('BotFloor5 ')

%legend('actual','setpoint')

%-------------Middle floor zones -------------------------
figure(6)
plot(Tpredict(6,:),'r')
hold on
plot(CTtemp(m:n,6),'b')
xlabel('time')
ylabel('MidFloor1 temperature')
legend('Predict','actual')
figure(7)
plot(Tpredict(7,:),'r')
hold on
plot(CTtemp(m:n,7),'b')
xlabel('time')
ylabel('MidFloor2 ')
%legend('actual','setpoint')
figure(8)
plot(Tpredict(8,:),'r')
hold on
plot(CTtemp(m:n,8),'b')
xlabel('time')
ylabel('MidFloor3 ')
figure(9)
plot(Tpredict(9,:),'r')
hold on
plot(CTtemp(m:n,9),'b')
xlabel('time')
ylabel('MidFloor4 ')
%legend('actual','setpoint')
figure(10)
plot(Tpredict(10,:),'r')
hold on
plot(CTtemp(m:n,10),'b')
xlabel('time')
ylabel('MidFloor5 ')

%---------------Top floor zones ----------------------------
figure(11)
plot(Tpredict(11,:),'r')
hold on
plot(CTtemp(m:n,11),'b')
xlabel('time')
ylabel('TopFloor1 temperature')
legend('Predict','actual')
figure(12)
plot(Tpredict(12,:),'r')
hold on
plot(CTtemp(m:n,12),'b')
xlabel('time')
ylabel('TopFloor2 ')
%legend('actual','setpoint')
figure(13)
plot(Tpredict(13,:),'r')
hold on
plot(CTtemp(m:n,13),'b')
xlabel('time')
ylabel('TopFloor3')
%legend('actual','setpoint')
figure(14)
plot(Tpredict(14,:),'r')
hold on
plot(CTtemp(m:n,14),'b')
xlabel('time')
ylabel('TopFloor4 ')
%legend('actual','setpoint')
figure(15)
plot(Tpredict(15,:),'r')
hold on
plot(CTtemp(m:n,15),'b')
xlabel('time')
ylabel('TopFloor5 ')
%legend('actual','setpoint')

