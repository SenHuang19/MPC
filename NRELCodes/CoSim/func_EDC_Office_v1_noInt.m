% clearvars
% clc
% close all

function [OptSchedule,OptStates,SolverStatus] = func_EDC_v1_noInt(Flags,Prices,Inputs,InitialStates,Params)

% run the EDC for the next 24-hour period, need the following as inputs

% Flags
% Prices
% Inputs
% Initial conditions
% Params

%% Flags
Flag_Plot = Flags.Flag_Plot;
Flag_SingleRun = Flags.Flag_SingleRun;
Flag_FC = Flags.Flag_FC;
Flag_Battery = Flags.Flag_Battery;
Flag_Thermal = Flags.Flag_Thermal;
Flag_Timestep = Flags.Flag_Timestep;
Flag_AS = Flags.Flag_AS;


%% Prices
c_e = Prices.c_e;
c_ng = Prices.c_ng;
c_as = Prices.c_as;
c_up = Prices.c_up;
c_down = Prices.c_down;


%% Initial conditions for building and equipment status
X0 = InitialStates.X;
Pfc_f0 = InitialStates.Pfc_f;
Pfc_h0 = InitialStates.Pfc_h;
Eb0 = InitialStates.Eb;
Eth0 = InitialStates.Eth;
w0 = InitialStates.w;
Pf_total = zeros(1,1);
Pcc_total = zeros(1,1);
Ppre_total = zeros(1,1);


%% Inputs
T_amb = Inputs.T_amb;
Qint = Inputs.Qint;
Pe_nonHVAC = Inputs.Pe_nonHVAC;
P_HotWater = Inputs.P_HotWater;


%% Params
eta_d = Params.eta_d;
eta_c = Params.eta_c;
eta_th_d = Params.eta_th_d;
eta_th_c = Params.eta_th_c;
eta_th_loss = Params.eta_th_loss;
eta_fc_e = Params.eta_fc_e;
eta_fc_h = Params.eta_fc_h;
kas = Params.kas;
Ras_max = Params.Ras_max*Flag_AS;

T_high = Params.T_high;
T_low = Params.T_low;
Eb_max = Params.Eb_max;
Eth_max = Params.Eth_max;
Pfc_f_max = Params.Pfc_f_max;
Pfc_f_min = Params.Pfc_f_min;
Pb_max = Params.Pb_max;
Pth_max = Params.Pth_max;

m_min = Params.m_min;
m_max = Params.m_max;
Prh_max = Params.Prh_max;
FC_ramp_up = Params.FC_ramp_up;
FC_ramp_down = Params.FC_ramp_down;
min_up = Params.min_up;
min_down = Params.min_down;


%% ROM parameters
A_mat = Params.A_mat;
B_mat = Params.B_mat;
E_mat = Params.E_mat;
Fan_mat = Params.Fan_mat;
Chiller_mat = Params.Chiller_mat;
Boiler_mat = Params.Boiler_mat;

Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(3,1); % supply air temperature
Zoneparam.T_approx = 22*ones(3,1); % approximation for room temperature for complexity
N_AHU = 3;
alpha_oa = 0.1;


%% Simulation constants
if Flag_Timestep == 0 
    dt = 1*60; % Time step for discrete thermal dynamics model, seconds
elseif Flag_Timestep == 1
    dt = 5*60;
end
dt_h = 60*60; % Time step for scheduling period
N_sch = 24; % Number of scheduling period
N_state = 3600/dt*N_sch; % Number of thermal dynamics model steps
N_perhour = dt_h/dt; % Number of thermal dynamics model steps per hour


%% Other params calculation

Tmix_vec = alpha_oa*T_amb + (1-alpha_oa)*Zoneparam.T_approx(1,1); % Assuming OA ratio and Tapprox are the same for all zones


%% EDC optimization
cvx_begin

% Decision variables:
% zone temperatue set-point, battery charge, battery discharge, supply
% air flow rate, FC fuel consumption set-point, FC heat output
% set-point, total reheat, AS capacity
variables u(N_sch,7) m_z(N_sch,15) Prh(N_sch,15) Ras_z(N_sch,15)
% u =[Pb_in,Pb_out,Pfc_f_ref,Pfc_h_ref,Ras,Pth_in,Pth_out]
% m_z: Zone airflow rate, '15' indicates 15 zones, similar for reheat (Prh)


% State variables:
% zone temperature, battery SOC
variables x(N_sch,15,N_perhour) Eb(N_sch,1,N_perhour) Eth(N_sch,1,N_perhour) % Three dimensional vector (second column represents 'size')
%('15' includes zones (Basement,Bottomfloor(6),Middlefloor(6),Topfloor(6))

% Interger variables:
% FC on/off, FC start up, FC shut down
% variable w(N_sch+max(min_up,min_down)) binary
% variable w_up(N_sch+max(min_up,min_down)) binary
% variable w_down(N_sch+max(min_up,min_down)) binary


expression Pcc_total(N_sch)
expression Ppre_total(N_sch)
expression Pf_total(N_sch)
expression Pcc_i(3)
expression J_i(N_sch)

% calculate cost
for i = 1:N_sch
    
    % Total cost incurred from various sources (electricity, natural gas, fuel cell, AS)
    % Electricity cost
    
    %% Chiller cost
  
    m_z_power = m_z;
    m_z_power(:,6:10) = m_z_power(:,6:10)*10;
    Prh_power = Prh;
    Prh_power(:,6:10) = Prh_power(:,6:10)*10;
    
    Tmix_i = Tmix_vec(i);
    for n_f = 1:N_AHU
        if Tmix_i >= Zoneparam.Ts(n_f)
            Pcc_i(n_f,1) = sum(m_z_power(i,5*n_f-4:5*n_f))*Zoneparam.Cp*(Tmix_i-Zoneparam.Ts(n_f,1));
        else
            Pcc_i(n_f,1) = 0;
        end
    end
    
    Pcc_total(i,1) = sum(Pcc_i);   % addition of the CAV system chiller power consumption
    
    Pch_i = Chiller_mat.d0 + Chiller_mat.d1*T_amb(i) + Chiller_mat.d2*(Pcc_total(i,1));  % Regression equation for chiller

    %%  Fan
    for n_f = 1:N_AHU
        eval(strcat('cfan_0 = Fan_mat.Fan',num2str(n_f),'.c0;'));
        eval(strcat('cfan_1 = Fan_mat.Fan',num2str(n_f),'.c1;'));
        eval(strcat('cfan_2 = Fan_mat.Fan',num2str(n_f),'.c2;'));
        Pf_i(n_f,1) = cfan_0 + cfan_1*sum(m_z_power(i,5*n_f-4:5*n_f)) + cfan_2*power(sum(m_z_power(i,5*n_f-4:5*n_f)),2); %m_s(i,n_f).^2;  % m_i is the supply air flow rate 
    end
    
    Pf_total(i,1) = sum(Pf_i);
    
    %
    Pfc_e_i = u(i,3)*eta_fc_e; % FC elec output
    
    P_e_i = Pch_i + Pf_total(i,1) + u(i,1) - Pfc_e_i - u(i,2) + Pe_nonHVAC(i,1); % grid power
    
    J_e_i = c_e(i)*P_e_i;
    
    Ppre_total(i,1) = sum(Prh_power(i,:),2) + P_HotWater(i,1);
    
    Pcc_boiler(i,1) = Boiler_mat.e0 + Boiler_mat.e1*T_amb(i) + Boiler_mat.e2*(Ppre_total(i,1)-u(i,4)+u(i,6)-u(i,7));
    
    P_ng_i = u(i,3) + Pcc_boiler(i,1);
    J_gas_i = c_ng(i)*P_ng_i;
    
    % AS payment
    J_as_i = -c_as(i)*u(i,5);
    
    % FC on/off
%     w_up_i = w_up(i);
%     w_down_i = w_down(i);
%     J_fc_com_i = c_up*w_up_i + c_down*w_down_i;
    
    % Total cost
    J(i,1) = J_e_i + J_gas_i + J_as_i;
    
end


minimize ( sum(J) )


subject to

for i_sch = 1:N_sch  % (nS=24)
    
    m_i = m_z(i_sch,:);
    Prh_i = Prh(i_sch,:);
    Pb_in_i = u(i_sch,1);
    Pb_out_i = u(i_sch,2);
    Pfc_f_ref_i = u(i_sch,3);
    Pfc_h_ref_i = u(i_sch,4);
    Ras_i = u(i_sch,5);
    Ras_z_i = Ras_z(i_sch,:)';
    Pth_in_i = u(i_sch,6);
    Pth_out_i = u(i_sch,7);
%     w_i = w(i_sch);
    
    if i_sch == 1
        Pfc_f_ref_iminus = Pfc_f0;
        Pfc_h_ref_iminus = Pfc_h0;
%         w_iminus = w0;
    else
        Pfc_f_ref_iminus = u((i_sch-1),3);
        Pfc_h_ref_iminus = u((i_sch-1),4);
%         w_iminus = w(i_sch-1);
    end
    
    Ras_i == sum(Ras_z_i);
    
    %Limits
    0 <= Pfc_h_ref_i <= Pfc_f_ref_i*eta_fc_h;
    0 <= sum(Prh(i_sch,:)) + Pth_in_i - Pth_out_i - Pfc_h_ref_i; 
    0 <= Ras_z_i <= Ras_max;
    m_min <= m_i';
%     m_i' <= m_max;
    Pfc_f_min*1 <= Pfc_f_ref_i <= Pfc_f_max*1;
    -FC_ramp_down <= Pfc_f_ref_i - Pfc_f_ref_iminus <= FC_ramp_up;
    0 <= Pb_in_i <= Pb_max;
    0 <= Pb_out_i <= Pb_max;
    0 <= Pth_in_i <= Pth_max*1;
    0 <= Pth_out_i <= Pth_max*1;
    Prh_i >= 0;
%     Prh_i' <= Prh_max; 
    
    % FC commitment
%     for i_up = 1:min_up % FC minimum up time
%         -w_iminus + w_i - w(i_sch+i_up) <= 0;
%     end
%     for i_down = 1:min_down % FC minimum down time
%         w_iminus - w_i + w(i_sch+i_down) <= 1;
%     end
%     
%     -w_iminus + w_i - w_up(i_sch) <= 0; % FC start up
%     w_iminus - w_i - w_down(i_sch) <= 0; % FC shut down
  
    
    % States in faster time step
    
    for i_state = 1:N_perhour
        
        T_i = x(i_sch,:,i_state)';
        
        Eb_i = Eb(i_sch,1,i_state);
        Eth_i = Eth(i_sch,1,i_state);
        
        if i_sch == 1 && i_state == 1
            X = X0;
            Eb_iminus = Eb0;
            Eth_iminus = Eth0;
        elseif i_state == 1
            X = x(i_sch-1,:,N_perhour)';
            Eb_iminus = Eb(i_sch-1,1,N_perhour);
            Eth_iminus = Eth(i_sch-1,1,N_perhour);
        else
            X = x(i_sch,:,i_state-1)';
            Eb_iminus = Eb(i_sch,1,i_state-1);
            Eth_iminus = Eth(i_sch,1,i_state-1);
        end
        
        U = [m_i';Prh_i'];  % size(U)=51=[15;15]
        D = [T_amb(i_sch);Qint(i_sch,:)';1]; %size(D)=16=[1;15]
        
        if Flag_Timestep == 0
            Xnew = A_mat*X+B_mat*U+E_mat*D; % X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states
        elseif Flag_Timestep == 1
            Xnew = A1*X+B1*U+E1*D; % X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states
        end

        Eb_new = Eb_iminus + dt*eta_c*Pb_in_i - dt*1/eta_d*Pb_out_i;
        Eth_new = Eth_iminus + dt*eta_th_c*Pth_in_i - dt*1/eta_th_d*Pth_out_i;
        
        Xnew == T_i;
        Eb_new == Eb_i;
        Eth_new == Eth_i;
     
        T_low + kas*Ras_z_i <= T_i <= T_high - kas*Ras_z_i;
        0 <= Eb_i <= Eb_max;
        -Eth_max <= Eth_i <= Eth_max;
        
    end
    
%    [Xnew(1:6);Xnew(8:12);Xnew(14:18)]' == 22; % temperature reach set-point
%     Xnew == 22;
    
end

Eb_new == Eb0; % final battery SOC
Eth_new == Eth0;

cvx_end


OptSchedule = [];
OptSchedule.u = u;
OptSchedule.m_z = m_z;
OptSchedule.Prh = Prh;

OptStates = [];
OptStates.x = x;
OptStates.Eb = Eb;
OptStates.Eth = Eth;
% OptStates.w = w;

SolverStatus = [];
SolverStatus.Status = cvx_status;
SolverStatus.OptVal = cvx_optval;

%% Post process

% % Calculate costs and other variables
% 
% for i = 1:N_sch
% 
%     % Chiller cost
%   
%     m_z_power = m_z;
%     m_z_power(:,6:10) = m_z_power(:,6:10)*10;
%     
%     Tmix_i = Tmix_vec(i);
%     for n_f = 1:N_AHU
%         if Tmix_i >= Zoneparam.Ts(n_f)
%             Pcc_i(n_f,1) = sum(m_z_power(i,5*n_f-4:5*n_f))*Zoneparam.Cp*(Tmix_i-Zoneparam.Ts(n_f,1));
%         else
%             Pcc_i(n_f,1) = 0;
%         end
%     end
%     
%     Pcc_total(i,1) = sum(Pcc_i);   % addition of the CAV system chiller power consumption
%     
%     Pch_i = Chiller_mat.d0 + Chiller_mat.d1*T_amb(i) + Chiller_mat.d2*(Pcc_total(i,1));  % Regression equation for chiller
% 
%     Pch_vec(i,1) = Pch_i;
% 
%     %  Fan
%     for n_f = 1:N_AHU
%         eval(strcat('cfan_0 = Fan_mat.Fan',num2str(n_f),'.c0;'));
%         eval(strcat('cfan_1 = Fan_mat.Fan',num2str(n_f),'.c1;'));
%         eval(strcat('cfan_2 = Fan_mat.Fan',num2str(n_f),'.c2;'));
%         Pf_i(n_f,1) = cfan_0 + cfan_1*sum(m_z_power(i,5*n_f-4:5*n_f)) + cfan_2*power(sum(m_z_power(i,5*n_f-4:5*n_f)),2); %m_s(i,n_f).^2;  % m_i is the supply air flow rate 
%     end
%     
%     Pf_total(i,1) = sum(Pf_i);
%     
%     %
%     Pfc_e_i = u(i,3)*eta_fc_e; % FC elec output
%     
%     P_e_i = Pch_i + Pf_total(i,1) + u(i,1) - Pfc_e_i - u(i,2); % grid power
%     P_e_vec(i,1) = P_e_i;
%     
%     J_e_i = c_e(i)*P_e_i;
%     J_e_vec(i,1) = J_e_i;
%     
%     Ppre_total(i,1) = sum(Prh(i,:),2);
%     Ppre_boiler(i,1) = Ppre_total(i,1) - u(i,4) + u(i,6) - u(i,7);
%     
%     Pcc_boiler(i,1) = Boiler_mat.e0 + Boiler_mat.e1*T_amb(i) + Boiler_mat.e2*(Ppre_total(i,1) - u(i,4) + u(i,6) - u(i,7));
%     
%     P_ng_i = u(i,3) + Pcc_boiler(i,1);
%     J_gas_i = c_ng(i)*P_ng_i;
%     J_gas_vec(i,1) = J_gas_i;
%     
%     % AS payment
%     J_as_i = -c_as(i)*u(i,5);
%     J_as_vec(i,1) = J_as_i;
%     
%     % FC on/off
% %     w_up_i = w_up(i);
% %     w_down_i = w_down(i);
% %     J_fc_com_i = c_up*w_up_i + c_down*w_down_i;
%     
%     % Total cost
%     J(i,1) = J_e_i + J_gas_i + J_as_i;
%     
% end
% 
% 
% % Plotting
% 
% if Flag_Plot
%     
%     tvec_sch = (0:N_sch)';
%     tvec_state = (0:N_sch*N_perhour)'/N_perhour;
%     
%    %---------------------Temperatures------------------------------ 
%    
%     figure
%     for i =1:5
%         subplot(3,2,i)
%         yyaxis left
%         plot(tvec_sch(2:end),x(:,i,end))
% %         yyaxis right
% %         plot(tvec_sch(2:end),T_amb)
% %         ylabel('Ambient')
%         legend('T')
%     end
%     
%     figure
%     for i =1:5
%         subplot(3,2,i)
%         yyaxis left
%         plot(tvec_sch(2:end),x(:,5+i,end))
% %         yyaxis right
% %         plot(tvec_sch(2:end),T_amb)
% %         ylabel('Ambient')
%         legend('T')
%     end
%     
%     figure
%     for i =1:5
%         subplot(3,2,i)
%         yyaxis left
%         plot(tvec_sch(2:end),x(:,10+i,end))
% %         yyaxis right
% %         plot(tvec_sch(2:end),T_amb)
% %         ylabel('Ambient')
%         legend('T')
%     end
%     
%   %---------------------- flow rates and reheats------------------------  
%       figure
%       for i =1:5
%  
%           subplot(3,2,i)
%           yyaxis left
%           plot(tvec_sch(1:end-1),m_z(:,i))
%           hold on
%           plot(tvec_sch(1:end-1),m_min(i)*ones(24,1),'--r')
%           hold on
%           plot(tvec_sch(1:end-1),m_max(i)*ones(24,1),'--b')
%           ylabel('flow rate')
%           yyaxis right
%           plot(tvec_sch(1:end-1),Prh(:,i))
%           ylabel('Reheat')
%           legend('flow rate','min flow','max flow','reheat')
%         
%       end
%       
%       figure
%       for i =1:5
%  
%           subplot(3,2,i)
%           yyaxis left
%           plot(tvec_sch(1:end-1),m_z(:,5+i))
%           hold on
%           plot(tvec_sch(1:end-1),m_min(5+i)*ones(24,1),'--r')
%           hold on
%           plot(tvec_sch(1:end-1),m_max(5+i)*ones(24,1),'--b')
%           ylabel('flow rate')
%           yyaxis right
%           plot(tvec_sch(1:end-1),Prh(:,5+i))
%           ylabel('Reheat')
%         legend('flow rate','min flow','max flow','reheat')
%       end
%     
%        figure
%       for i =1:5
%  
%           subplot(3,2,i)
%           yyaxis left
%           plot(tvec_sch(1:end-1),m_z(:,10+i))
%           hold on
%           plot(tvec_sch(1:end-1),m_min(10+i)*ones(24,1),'--r')
%           hold on
%           plot(tvec_sch(1:end-1),m_max(10+i)*ones(24,1),'--b')
%           ylabel('flow rate')
%           yyaxis right
%           plot(tvec_sch(1:end-1),Prh(:,10+i))
%           ylabel('Reheat')
%           legend('flow rate','min flow','max flow','reheat')
%         
%       end
%       
%       
%   
%     % FC
%     Pb_in = u(:,1);
%     Pb_out = u(:,2);
%     Pfc_f_ref = u(:,3);
%     Pfc_h_ref = u(:,4);
%     Ras = u(:,5);
%     Pth_in = u(:,6);
%     Pth_out = u(:,7);
%     
%     Eb_vec = [];
%     for i = 1:N_sch
%         Eb_vec = [Eb_vec;reshape(Eb(i,1,:),[N_perhour,1])];
%     end
%     
%     Eth_vec = [];
%     for i = 1:N_sch
%         Eth_vec = [Eth_vec;reshape(Eth(i,1,:),[N_perhour,1])];
%     end
%     
%     figure
%     subplot(2,1,1)
%     plot(tvec_sch(1:end-1),Pb_in)
%     hold on 
%     plot(tvec_sch(1:end-1),Pb_out)
%     legend('Pb in','Pb out')
%     grid on
%     subplot(2,1,2)
%     plot(tvec_state(2:end),Eb_vec)
%     legend('Eb')
%     grid on
%     
%     figure
%     plot(tvec_sch(1:end-1),Pfc_f_ref)
%     hold on
%     plot(tvec_sch(1:end-1),Pfc_h_ref)
%     plot(tvec_sch(1:end-1),sum(Prh,2))
%     legend('Pfc f','Pfc h','P h')
%     grid on
%     
%     figure
%     plot(tvec_sch(1:end-1),Pch_vec)
%     hold on
%     plot(tvec_sch(1:end-1),Pf_total)
%     plot(tvec_sch(1:end-1),Pcc_boiler)
%     plot(tvec_sch(1:end-1),P_e_vec)
%     ylabel('Power')
%     legend('Chiller','Fan','Boiler','Grid')
%     
%     figure
%     plot(tvec_sch(1:end-1),Pcc_total)
%     hold on
%     plot(tvec_sch(1:end-1),Ppre_total)
%     plot(tvec_sch(1:end-1),sum(Qint,2))
%     ylabel('Power')
%     legend('Cooling coil','Reheart','Q int')
%     
%     figure
%     plot(tvec_sch(1:end-1),J_e_vec)
%     hold on
%     plot(tvec_sch(1:end-1),J_gas_vec)
%     plot(tvec_sch(1:end-1),J_as_vec)
%     legend('Elec','NG','AS')
%     
%     figure
%     plot(tvec_sch(2:end),T_amb)
%     legend('Ambient')
%     
%     figure
%     subplot(2,1,1)
%     plot(tvec_sch(1:end-1),Pth_in)
%     hold on 
%     plot(tvec_sch(1:end-1),Pth_out)
%     legend('Pth in','Pth out')
%     grid on
%     subplot(2,1,2)
%     plot(tvec_state(2:end),Eth_vec)
%     legend('Eth')
%     grid on
%     
%     figure
%     subplot(2,1,1)
%     plot(Pch_vec)
%     hold on
%     plot(Pcc_total)
%     legend('ch','cc')
%     subplot(2,1,2)
%     plot(Ppre_boiler)
%     hold on
%     plot(Pcc_boiler)
%     legend('rh','boiler')
%     
% end
% 
% 
% 
