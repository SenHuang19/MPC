clearvars
close all

load SimuData_MPC_1h.mat

for i_zone = 1:2
    figure
    title(strcat('MPC Zone ',num2str(i_zone)))
    subplot(3,1,1)
    plot(logdata_OptSchedule{24,1}.m_z(:,i_zone))
    legend('m_z')
    subplot(3,1,2)
    plot(logdata_OptSchedule{24,1}.Prh(:,i_zone))
    legend('Prh')
    subplot(3,1,3)
    plot(logdata_OptStates{24,1}.x(:,i_zone))
    legend('Tz')
end

tshift = 24;

location = 'Baltimore';
elec_price = 'High';
natgas_price = 'Low';
c_e_type = elec_price;
c_ng_type = natgas_price;

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

dt = 60*60;
dt_h = 60*60; % Time step for scheduling period
N_sch = 24; % Number of scheduling period
N_state = 3600/dt*N_sch; % Number of thermal dynamics model steps
N_perhour = dt_h/dt; % Number of thermal dynamics model steps per hour

eval(strcat('c_e_vec = ElecPrice.',c_e_type,';')); % $/kwh
eval(strcat('c_ng_vec_mmbtu = GasPrice.',c_ng_type,';')); % $/mmbtu
c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kwh

c_e = c_e_vec(tshift:tshift+23);
c_ng = c_ng_vec(tshift:tshift+23);

m_z_power = logdata_output(25:end,16:30)/0.1;
Prh_power = logdata_output(25:end,31:45);

alpha_oa = 0.1;

Tamb_downsample = downsample(Tamb,60);
T_amb = Tamb_downsample(tshift:(tshift+23));
QintTotal_downsample = downsample(QintTotal,60);
Qint = QintTotal_downsample(tshift:(tshift+23),:);



Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(3,1); % supply air temperature
Zoneparam.T_approx = 22*ones(3,1); % approximation for room temperature for complexity
N_AHU = 3;
alpha_oa = 0.1;


Tmix_vec = alpha_oa*T_amb + (1-alpha_oa)*Zoneparam.T_approx(1,1); % Assuming OA ratio and Tapprox are the same for all zones


for i = 1:24

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

J_overall_MPC = sum(J)

%%

clearvars

load SimuData_Default_1h.mat

for i_zone = 1:2
    figure
    title(strcat('Default Zone ',num2str(i_zone)))
    subplot(3,1,1)
    plot(logdata_output(25:end,i_zone+15)/0.1)
    legend('m_z')
    subplot(3,1,2)
    plot(logdata_output(25:end,i_zone+30))
    legend('Prh')
    subplot(3,1,3)
    plot(logdata_output(25:end,i_zone))
    legend('Tz')
end

tshift = 24;

location = 'Baltimore';
elec_price = 'High';
natgas_price = 'Low';
c_e_type = elec_price;
c_ng_type = natgas_price;

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

dt = 60*60;
dt_h = 60*60; % Time step for scheduling period
N_sch = 24; % Number of scheduling period
N_state = 3600/dt*N_sch; % Number of thermal dynamics model steps
N_perhour = dt_h/dt; % Number of thermal dynamics model steps per hour

eval(strcat('c_e_vec = ElecPrice.',c_e_type,';')); % $/kwh
eval(strcat('c_ng_vec_mmbtu = GasPrice.',c_ng_type,';')); % $/mmbtu
c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kwh

c_e = c_e_vec(tshift:tshift+23);
c_ng = c_ng_vec(tshift:tshift+23);

m_z_power = logdata_output(25:end,16:30)/0.1;
Prh_power = logdata_output(25:end,31:45);

alpha_oa = 0.1;

Tamb_downsample = downsample(Tamb,60);
T_amb = Tamb_downsample(tshift:(tshift+23));
QintTotal_downsample = downsample(QintTotal,60);
Qint = QintTotal_downsample(tshift:(tshift+23),:);



Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(3,1); % supply air temperature
Zoneparam.T_approx = 22*ones(3,1); % approximation for room temperature for complexity
N_AHU = 3;
alpha_oa = 0.1;


Tmix_vec = alpha_oa*T_amb + (1-alpha_oa)*Zoneparam.T_approx(1,1); % Assuming OA ratio and Tapprox are the same for all zones


for i = 1:24

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

J_overall_Default = sum(J)

%%

clearvars

load SimuData_Default_1h_Single.mat

% for i_zone = 6:10
%     figure
%     title(strcat('Default Zone ',num2str(i_zone)))
%     subplot(3,1,1)
%     plot(logdata_output(25:end,i_zone+15)/0.1)
%     legend('m_z')
%     subplot(3,1,2)
%     plot(logdata_output(25:end,i_zone+30))
%     legend('Prh')
%     subplot(3,1,3)
%     plot(logdata_output(25:end,i_zone))
%     legend('Tz')
% end

tshift = 24;

location = 'Baltimore';
elec_price = 'High';
natgas_price = 'Low';
c_e_type = elec_price;
c_ng_type = natgas_price;

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

dt = 60*60;
dt_h = 60*60; % Time step for scheduling period
N_sch = 24; % Number of scheduling period
N_state = 3600/dt*N_sch; % Number of thermal dynamics model steps
N_perhour = dt_h/dt; % Number of thermal dynamics model steps per hour

eval(strcat('c_e_vec = ElecPrice.',c_e_type,';')); % $/kwh
eval(strcat('c_ng_vec_mmbtu = GasPrice.',c_ng_type,';')); % $/mmbtu
c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kwh

c_e = c_e_vec(tshift:tshift+23);
c_ng = c_ng_vec(tshift:tshift+23);

m_z_power = logdata_output(25:end,16:30)/0.1;
Prh_power = logdata_output(25:end,31:45);

alpha_oa = 0.1;

Tamb_downsample = downsample(Tamb,60);
T_amb = Tamb_downsample(tshift:(tshift+23));
QintTotal_downsample = downsample(QintTotal,60);
Qint = QintTotal_downsample(tshift:(tshift+23),:);



Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(3,1); % supply air temperature
Zoneparam.T_approx = 22*ones(3,1); % approximation for room temperature for complexity
N_AHU = 3;
alpha_oa = 0.1;


Tmix_vec = alpha_oa*T_amb + (1-alpha_oa)*Zoneparam.T_approx(1,1); % Assuming OA ratio and Tapprox are the same for all zones


for i = 1:24

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

J_overall_Default_Single = sum(J)