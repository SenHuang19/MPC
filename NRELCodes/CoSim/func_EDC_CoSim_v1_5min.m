function [OptSchedule,OptStates,SolverStatus] = func_EDC_CoSim_v1(tshift,InitialStates)

scenario_id = 2;
season = 'Summer';
location = 'Baltimore';
building = 'Office';
elec_price = 'High';
natgas_price = 'Low';
mismatch = 'Zero';

scenarioID = scenario_id;
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

dt = 5*60;
dt_h = 60*60; % Time step for scheduling period
N_sch = 3600/dt*24; % Number of scheduling period

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

Tamb_downsample = downsample(Tamb,12);
T_amb = Tamb_downsample(tshift:(tshift+23));
QintTotal_downsample = downsample(QintTotal,12);
Qint = QintTotal_downsample(tshift:(tshift+23),:);

%% Run the EDC

Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(3,1); % supply air temperature
Zoneparam.T_approx = 22*ones(3,1); % approximation for room temperature for complexity
N_AHU = 3;
alpha_oa = 0.1;


%% Other params calculation

Tmix_vec = alpha_oa*T_amb + (1-alpha_oa)*Zoneparam.T_approx(1,1); % Assuming OA ratio and Tapprox are the same for all zones


%% EDC optimization
cvx_begin

    % Decision variables:
    variables m_z(N_sch,15) Prh(N_sch,15)

    % State variables:
    % zone temperature
    variables x(N_sch,15)

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

        P_e_i = Pf_total(i,1); % grid power

        J_e_i = c_e(i)*P_e_i;

        Ppre_total(i,1) = sum(Prh_power(i,:),2);

        Pcc_boiler(i,1) = Boiler_mat.e0 + Boiler_mat.e1*T_amb(i) + Boiler_mat.e2*(Ppre_total(i,1));

        P_ng_i = Pcc_boiler(i,1);
        J_gas_i = c_ng(i)*P_ng_i;

        % Total cost
        J(i,1) = J_e_i + J_gas_i;

    end


    minimize ( sum(J) )


    subject to

    for i_sch = 1:N_sch  % (nS=24)

        m_i = m_z(i_sch,:);
        Prh_i = Prh(i_sch,:);

        %Limits
        m_min <= m_i';
%         m_i' <= m_max;
        Prh_i >= 0;
%         Prh_i' <= Prh_max; 

        T_i = x(i_sch,:)';
        
        if i_sch == 1
            X = X0;
        else
            X = x(i_sch-1,:)';
        end
        
        U = [m_i';Prh_i'];  % size(U)=51=[15;15]
        D = [T_amb(i_sch);Qint(i_sch,:)';1]; %size(D)=16=[1;15]
        
%         if Flag_Timestep == 0
%             Xnew = A_mat*X+B_mat*U+E_mat*D; % X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states
%         elseif Flag_Timestep == 1
%             Xnew = A1*X+B1*U+E1*D; % X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states
%         end
        
        Xnew = A_mat*X+B_mat*U+E_mat*D; % X=[Bot_floor_zones,Mid_floor_zones,Top_floor_zones] total = 15 states
        
        Xnew == T_i;
        
        T_low <= T_i <= T_high ;

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





