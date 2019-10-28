function [OptSchedule,OptStates,SolverStatus] = func_optimization_v1(tshift,InitialStates)

% tshift = 1;
% X0_init = 22*ones(15,1);
% InitialStates = [];
% InitialStates.X = X0_init;

season = 'Summer';
location = 'Baltimore';
building = 'Office';
elec_price = 'High';
natgas_price = 'Low';
mismatch = 'Zero';

season = season;
location = location;
building = building; % 'Office', 'Hotel'
c_e_type = elec_price;
c_ng_type = natgas_price;
Flag_Mismatch = mismatch;


%% Load inputs and parameters

load('..\coeff\coeff\Simparamdata60min.mat');
load('..\coeff\coeff\summer_Model_valid.mat');
load OfficePricingInfo.mat;
if strcmp(location,'Baltimore')
    load Inputs_OfficeROM_sizing_Baltimore.mat;
elseif strcmp(location,'SanFrancisco')
    load Inputs_OfficeROM_sizing_SanFrancisco.mat;
end

m_min = Params_office_sizing(:,2);
m_min(6:10) = m_min(6:10)/10;
m_max = Params_office_sizing(:,1);
m_max(6:10) = m_max(6:10)/10;
Prh_max = Params_office_sizing(:,3);
Prh_max(6:10) = Prh_max(6:10)/10;

A_mat = A.Summer.Office;
B_mat = B.Summer.Office;
E_mat = E.Summer.Office;
Fan_mat = Fan.Summer.Office;
Chiller_mat = Chiller.Summer.Office;
Boiler_mat = Boiler.Summer.Office;


%% Simulation constants

dt_control = 60*60; % Controller time step, sec
dt_model = 60*60; % Building model time step, sec 
pred_horizon = 24; % Prediction horizon, h

N_control = pred_horizon*3600/dt_control;
N_model = pred_horizon*3600/dt_model;
N_model_per_controlstep = dt_control/dt_model;


%% Prices
%     c_e_type = 'Low';
%     c_ng_type = 'High';

eval(strcat('c_e_vec = ElecPrice.',c_e_type,';')); % $/kwh
eval(strcat('c_ng_vec_mmbtu = GasPrice.',c_ng_type,';')); % $/mmbtu
c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kwh

c_e = c_e_vec(tshift:tshift+23);
c_ng = c_ng_vec(tshift:tshift+23);


%% Parameters (tunable)

T_high = 24*ones(15,1); % temperature bounds
T_low = 20*ones(15,1);

%% Initial conditions
X0 = InitialStates.X;

%% Inputs

Tamb_downsample = downsample(Tamb,(dt_control/60));
T_amb = Tamb_downsample(tshift:(tshift+23));
QintTotal_downsample = downsample(QintTotal,(dt_control/60));
Qint = QintTotal_downsample(tshift:(tshift+23),:);

%% Run the EDC

Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(3,1); % supply air temperature
Zoneparam.T_approx = 22*ones(3,1); % approximation for room temperature for complexity
N_AHU = 3;
alpha_oa = 0.1;


%% Other params calculation

Tmix_vec = alpha_oa*T_amb + (1-alpha_oa)*Zoneparam.T_approx(1,1); % Assuming OA ratio and Tapprox are the same for all zones

%% Optimization 

cvx_begin

% Decision variables:
variables m_z(N_control,15) Prh(N_control,15) 

% State variables:
variables x(N_control,15,N_model_per_controlstep) 

expression Pcc_total(N_control)
expression Ppre_total(N_control)
expression Pf_total(N_control)
expression Pcc_i(3)
expression J_i(N_control)

% calculate cost
for i = 1:N_control
    
    % Chiller cost
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
    Pcc_total(i,1) = sum(Pcc_i); % Cooling at the cooling coil
    
    Pch_i = Chiller_mat.d0 + Chiller_mat.d1*T_amb(i) + Chiller_mat.d2*(Pcc_total(i,1));  % Regression equation for chiller
    Pch_vec(i,1) = Pch_i;
    
    %  Fan
    for n_f = 1:N_AHU
        eval(strcat('cfan_0 = Fan_mat.Fan',num2str(n_f),'.c0;'));
        eval(strcat('cfan_1 = Fan_mat.Fan',num2str(n_f),'.c1;'));
        eval(strcat('cfan_2 = Fan_mat.Fan',num2str(n_f),'.c2;'));
        Pf_i(n_f,1) = cfan_0 + cfan_1*sum(m_z_power(i,5*n_f-4:5*n_f)) + cfan_2*power(sum(m_z_power(i,5*n_f-4:5*n_f)),2); %m_s(i,n_f).^2;  % m_i is the supply air flow rate 
    end
    Pf_total(i,1) = sum(Pf_i);
   
    % Grid power
    P_e_i = Pch_i + Pf_total(i,1); 
    P_e_vec(i,1) = P_e_i;

    J_e_i = c_e(i)*P_e_i;
    
    % Heating
    Ppre_total(i,1) = sum(Prh_power(i,:),2);
    
    Pcc_boiler(i,1) = Boiler_mat.e0 + Boiler_mat.e1*T_amb(i) + Boiler_mat.e2*(Ppre_total(i,1));
    
    P_ng_i = Pcc_boiler(i,1);
    J_gas_i = c_ng(i)*P_ng_i;
    
    % Total cost
    J(i,1) = J_e_i + J_gas_i;
    
end


minimize ( sum(J) )


subject to

for i_control = 1:N_control  % (nS=24)
    
    m_i = m_z(i_control,:);
    Prh_i = Prh(i_control,:);
   
    %Limits
    m_min <= m_i';
%     m_i' <= m_max;
    Prh_i >= 0;
%     Prh_i' <= Prh_max; 
    
    % States in faster time step
    
    for i_model = 1:N_model_per_controlstep
        
        T_i = x(i_control,:,i_model)';
        
        if i_control == 1 && i_model == 1
            X = X0;
        elseif i_model == 1
            X = x(i_control-1,:,N_model_per_controlstep)';
        else
            X = x(i_control,:,i_model-1)';
        end
        
        U = [m_i';Prh_i'];  
        D = [T_amb(i_control);Qint(i_control,:)';1]; 
        
        Xnew = A_mat*X+B_mat*U+E_mat*D;

        Xnew == T_i; % Temperature dynamics const.
      
        T_low <= T_i <= T_high; % Temperature range const.
        
    end

end

cvx_end


OptSchedule = [];
OptSchedule.m_z = m_z;
OptSchedule.Prh = Prh;

OptStates = [];
OptStates.x = x;

SolverStatus = [];
SolverStatus.Status = cvx_status;
SolverStatus.OptVal = cvx_optval;


