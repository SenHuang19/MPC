% This script simulates a supervisory controller for an HVAC system. The
% controller computes the zone temperature set-point based on the current
% time and the outdoor dry-bulb temperature. The building is simulated by
% EnergyPlus. This simulation is the same as that implemented in
% simple.mdl, but uses plain Matlab code instead of Simulink.
%
% This script illustrates the usage of class mlepProcess in the MLE+
% toolbox for feedback control which involves whole building simulation.
% It has been tested with Matlab R2009b and EnergyPlus 6.0.0.
%
% This example is taken from an example distributed with BCVTB version
% 0.8.0 (https://gaia.lbl.gov/bcvtb).
%
% This script is free software.
%
% (C) 2010-2015 by Truong Nghiem and Willy Bernal (nghiem@seas.upenn.edu, Willy.BernalHeredia@nrel.gov)
%
% CHANGES:
%   2012-04-23  Fix an error with E+ 7.0.0: Matlab must read data from E+
%               before sending any data to E+.
%   2011-07-13  Update to new version of MLE+ which uses mlepInit for
%               system settings.
%   2011-04-28  Update to new version of MLE+ which uses Java for running
%               E+.
%   2015-07-30  Update to new version of MLE+ which uses the Matlab system 
%               call for running E+. Updated BCVTB libraries to 1.5. 
%               

tic

%% Create an mlepProcess instance and configure it

ep = mlepProcess;
% ep.arguments = {'Building', 'USA_CA_San.Francisco.Intl.AP.724940_TMY3'};
%ep.arguments = {'B_3rdCoSim_v3_FanOn', 'USA_VA_Sterling-Washington.Dulles.Intl.AP.724030_TMY3'};
ep.arguments = {'B_Office_20181228', 'USA_MD_Baltimore-Washington.Intl.AP.724060_TMY3'};
ep.acceptTimeout = 6000;

VERNUMBER = 2;  % version number of communication protocol (2 for E+ 7.2.0)

%% Start EnergyPlus cosimulation
[status, msg] = ep.start;  

if status ~= 0
    error('Could not start EnergyPlus: %s.', msg);
end

[status, msg] = ep.acceptSocket;

if status ~= 0
    error('Could not connect to EnergyPlus: %s.', msg);
end

%% The main simulation loop

t_end = 2;
t_EDC = 1;
deltaT = 60*60;  % time step
step_per_day = 24*3600/deltaT;
kStep = 1;  % current simulation step
MAXSTEPS = step_per_day*t_end+1;  % first arguement is the time step (hour =1, 15min =4) max simulation time = 4 days

% eta_fce = 0.5;
% Params_Inputs = [];

% logdata stores set-points, outdoor temperature, and zone temperature at
% each time step.
% logdata = zeros(MAXSTEPS, 4);
logdata_SP = [];
logdata_output = [];
logdata_opt = [];
InitialStates_vec = [];
InitialStates = [];
%load('PRBS_data.mat')
%RandSPGen =(24-23.5).*rand(MAXSTEPS,1) + 23.5;
%SP_data=dataPRBS(300:300+8761);
% SP_data=[repmat([20 20 20.3 20.3],1,ceil((MAXSTEPS/(4*10)))),...   
%     repmat([20.4 20.4 20.7 20.7],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([20.8 20.8 21.1 21.1],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([21.2 21.5 21.5 21.5],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([21.6 21.9 21.9 21.9],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([22 22 22.3 22.3],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([22.4 22.4 22.7 22.7],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([22.8 22.8 23.1 23.1],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([23.2 23.5 23.5 23.5],1,ceil((MAXSTEPS/(4*10)))),...
%     repmat([23.6 23.9 23.9 23.9],1,ceil((MAXSTEPS/(4*10))))]; 
% randomization command BB=AA(randperm(length(AA)));
%SP_data=[repmat([20 20 20 20],1,24),repmat([20 20 20 20],1,8),repmat([20 20 20.5 20.5],1,10),repmat([20 20 20 20],1,6),repmat([20 20 21 21],1,8760/4)];


%SP_data=repmat([20 20 21 21],1,8760/4);

X0_init = 22*ones(15,1);  % '15' stands for the dimension of the state space model

InitialStates = [];
InitialStates.X = X0_init;

tic
while kStep <= MAXSTEPS    
    % Read a data packet from E+
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Parse it to obtain building outputs
    [flag, eptime, outputs] = mlepDecodePacket(packet);
    if flag ~= 0, break; end
        
    % BEGIN Compute next set-points
    if kStep < step_per_day*t_EDC
        T_float = [];
        for i = 1:16
            T_float = [T_float,22,22];
        end
        SP = [0,0.1*ones(1,4),T_float];
        OptSchedule = [];
        OptStates = [];
    else
        if ~(kStep == step_per_day*t_EDC)
            if strcmp(SolverStatus.Status,'Solved') || strcmp(SolverStatus.Status,'Inaccurate/Solved') % if not solved, use last step's initial values
                InitialStates.X = outputs(1:15)';  % '15' stands for the dimension of the state space model
            end
        end
        
        [OptSchedule,OptStates,SolverStatus] = func_optimization_v1(kStep,InitialStates);
        
        SP = zeros(1,34);
        SP(1) = OptSchedule.Prh(1);
        SP(2:5) = 0.1;
        SP(6:7) = 22;
        for j_zone = 1:5
            SP(8+(j_zone-1)*2:8+(j_zone*2)-1) = OptStates.x(1,j_zone);
            SP(18+(j_zone-1)*2:18+(j_zone*2)-1) = OptStates.x(1,5+j_zone);
            SP(28+(j_zone-1)*2:28+(j_zone*2)-1) = OptStates.x(1,10+j_zone);
        end

        
%         T_float = [];
%         for i = 1:16
%             T_float = [T_float,20,24];
%         end
%         SP = [0,0.1*ones(1,4),T_float];
        
%         T_float = [];
%         for i = 1:16
%             T_float = [T_float,22,22];
%         end
%         SP = [0,0.1*ones(1,4),T_float];

    end
    
    logdata_OptSchedule{kStep,1} = OptSchedule;
    logdata_OptStates{kStep,1} = OptStates;
    
    
    % END Compute next set-points
    
    % Write to inputs of E+
    ep.write(mlepEncodeRealData(VERNUMBER, 0, (kStep-1)*deltaT, SP));    

    % Save to logdata

%     logdata(kStep, :) = [SP outputs];
    logdata_SP = [logdata_SP;SP];
    
%     logdata(kStep, :) = [SP outputs];
    logdata_output = [logdata_output;outputs];
    
    kStep = kStep + 1;
end

% Stop EnergyPlus
ep.stop;

disp(['Stopped with flag ' num2str(flag)]);

% figure(1)
% yyaxis left
% plot(logdata_output(:,29),'-*')
% xlabel('Time ','fontweight','bold')
% ylabel('Heating Coil Heating Rate','fontweight','bold')
% yyaxis right
% plot(logdata_output(:,76),'-o')
% %plot(logdata_output(:,107),'-o')
% ylabel('Zone Air flow rate','fontweight','bold')
% 
% figure(2)
% plot(logdata_output(:,1),'-*')
% xlabel('Time','fontweight','bold')
% ylabel('Zone Air Temperature (Basement)','fontweight','bold')
% hold on
% stairs(logdata_SP(:,6))
% ylabel('Zone Air Temperature setpoint (Basement)','fontweight','bold')

%Testcaseresult
% Remove unused entries in logdata
% kStep = kStep - 1;
% if kStep < MAXSTEPS
%     logdata((kStep+1):end,:) = [];
% end

%% Post process

% for i_zone = 1:15
%     figure
%     title(strcat('Zone ',num2str(i_zone)))
%     subplot(3,1,1)
%     plot(logdata_OptSchedule{24,1}.m_z(:,i_zone))
%     legend('m_z')
%     subplot(3,1,2)
%     plot(logdata_OptSchedule{24,1}.Prh(:,i_zone))
%     legend('Prh')
%     subplot(3,1,3)
%     plot(logdata_OptStates{24,1}.x(:,i_zone))
%     legend('Tz')
% end


% % ==========FLAGS==============
% % Flag	Description
% % +1	Simulation reached end time.
% % 0	    Normal operation.
% % -1	Simulation terminated due to an unspecified error.
% % -10	Simulation terminated due to an error during the initialization.
% % -20	Simulation terminated due to an error during the time integration.
% 
toc