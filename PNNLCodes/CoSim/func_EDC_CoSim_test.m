function [Optimal_Temp_Ctrl] = func_EDC_CoSim_test(ST, c_e, c_ng, T_out, Q_int, T_ini)

% ST = 5; % sampling time: 1, 5, 60 min
sample_time = strcat(num2str(ST), 'min');
num_samples = 24*60/ST; % number of samples

num_zone = 16; % number of zones

%% load the coefficients a_0 ... a_5
top_floor = load(strcat('Top_floor/', sample_time, '/top_floor.mat'));
mid_floor = load(strcat('Mid_floor/', sample_time, '/mid_floor.mat'));
bot_floor = load(strcat('Bot_floor/', sample_time, '/bot_floor.mat'));

a_0 = zeros(num_zone, 1);
a_1 = zeros(num_zone, 1);
a_2 = zeros(num_zone, 1);
a_3 = zeros(num_zone, 1);
a_4 = zeros(num_zone, 1);
a_5 = zeros(num_zone, 1);

for i = 1 : num_zone
    if  1 <= i && i <= 5 % top_floor
        a_0(i) = eval(strcat('top_floor.zone4_', num2str(i), '_summer_result.a0'));
        a_1(i) = eval(strcat('top_floor.zone4_', num2str(i), '_summer_result.a1'));
        a_2(i) = eval(strcat('top_floor.zone4_', num2str(i), '_summer_result.a2'));
        a_3(i) = eval(strcat('top_floor.zone4_', num2str(i), '_summer_result.a3'));
        a_4(i) = eval(strcat('top_floor.zone4_', num2str(i), '_summer_result.a4'));
        a_5(i) = eval(strcat('top_floor.zone4_', num2str(i), '_summer_result.a5'));
    end
    if  6 <= i && i <= 10 % mid_floor
        a_0(i) = eval(strcat('mid_floor.zone4_', num2str(i-5), '_summer_result.a0'));
        a_1(i) = eval(strcat('mid_floor.zone4_', num2str(i-5), '_summer_result.a1'));
        a_2(i) = eval(strcat('mid_floor.zone4_', num2str(i-5), '_summer_result.a2'));
        a_3(i) = eval(strcat('mid_floor.zone4_', num2str(i-5), '_summer_result.a3'));
        a_4(i) = eval(strcat('mid_floor.zone4_', num2str(i-5), '_summer_result.a4'));
        a_5(i) = eval(strcat('mid_floor.zone4_', num2str(i-5), '_summer_result.a5'));
    end
    if  11 <= i && i <= 15 % bot_floor
        a_0(i) = eval(strcat('bot_floor.zone4_', num2str(i-10), '_summer_result.a0'));
        a_1(i) = eval(strcat('bot_floor.zone4_', num2str(i-10), '_summer_result.a1'));
        a_2(i) = eval(strcat('bot_floor.zone4_', num2str(i-10), '_summer_result.a2'));
        a_3(i) = eval(strcat('bot_floor.zone4_', num2str(i-10), '_summer_result.a3'));
        a_4(i) = eval(strcat('bot_floor.zone4_', num2str(i-10), '_summer_result.a4'));
        a_5(i) = eval(strcat('bot_floor.zone4_', num2str(i-10), '_summer_result.a5'));
    end
    if  i == 16 % basement
        a_0(i) = eval(strcat('bot_floor.zone', num2str(i-10), '_summer_result.a0'));
        a_1(i) = eval(strcat('bot_floor.zone', num2str(i-10), '_summer_result.a1'));
        a_2(i) = eval(strcat('bot_floor.zone', num2str(i-10), '_summer_result.a2'));
        a_3(i) = eval(strcat('bot_floor.zone', num2str(i-10), '_summer_result.a3'));
        a_4(i) = eval(strcat('bot_floor.zone', num2str(i-10), '_summer_result.a4'));
        a_5(i) = eval(strcat('bot_floor.zone', num2str(i-10), '_summer_result.a5'));
    end
end


%% Read the initial zone temperature, the predicted outdoor temperature and internal heat gain

% T_ini = zeros(num_zone, 1); % initial zone temperature
% 
% T_out_data = csvread('../data/Top_floor/summer_zone4-1.csv', 1, 1, [1 1 60*24 1]);
% Q_int_data = zeros(60*24, num_zone);
% for i = 1 : num_zone
%     if  1 <= i && i <= 5 % top_floor
%         T_ini(i) = csvread(strcat('../data/Top_floor/summer_zone4-', num2str(i), '.csv'), 1, 3, [1 3 1 3]);
%         Q_int_data(:, i) = csvread(strcat('../data/Top_floor/summer_zone4-', num2str(i), '.csv'), 1, 2, [1 2 60*24 2]);
%     end
%     if  6 <= i && i <= 10 % mid_floor
%         T_ini(i) = csvread(strcat('../data/Mid_floor/summer_zone4-', num2str(i-5), '.csv'), 1, 3, [1 3 1 3]);
%         Q_int_data(:, i) = csvread(strcat('../data/Mid_floor/summer_zone4-', num2str(i-5), '.csv'), 1, 2, [1 2 60*24 2]);
%     end
%     if  11 <= i && i <= 15 % bot_floor
%         T_ini(i) = csvread(strcat('../data/Bot_floor/summer_zone4-', num2str(i-10), '.csv'), 1, 3, [1 3 1 3]);
%         Q_int_data(:, i) = csvread(strcat('../data/Bot_floor/summer_zone4-', num2str(i-10), '.csv'), 1, 2, [1 2 60*24 2]);
%     end
%     if  i == 16 % basement
%         T_ini(i) = csvread(strcat('../data/Bot_floor/summer_zone', num2str(i-10), '.csv'), 1, 3, [1 3 1 3]);
%         Q_int_data(:, i) = csvread(strcat('../data/Bot_floor/summer_zone', num2str(i-10), '.csv'), 1, 2, [1 2 60*24 2]);
%     end
% end
% 
% % T_out_avrgsample = zeros(num_samples, 1);
% % Q_int_avrgsample = zeros(num_samples, num_zone);
% % for j = 1 : ST : 60*24 % calculate the average outdoor temperature and internal heat gain during each sampling period
% %     k = (j+ST-1)/ST;
% %     T_out_avrgsample(k) = sum(T_out_data(j:j+ST-1))/ST;
% %     for i = 1 : num_zone
% %         Q_int_avrgsample(k, i) = sum(Q_int_data(j:j+ST-1, i))/ST;
% %     end
% % end
% 
% T_out_downsample = downsample(T_out_data, ST);
% Q_int_downsample = downsample(Q_int_data, ST);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optimization engine

% mpc parameters
% tshift = 1;

N_schd = num_samples;
N_zone = num_zone;
N_AHU = 4;

% load price (hourly price)
% load OfficePricingInfo.mat;
% c_e_vec = eval(strcat('ElecPrice.Low;')); % $/kWh
% c_ng_vec_mmbtu = eval(strcat('GasPrice.High;')); % $/mmbtu
% c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kWh
% 
% c_e_vec = kron(c_e_vec, ones(60/ST, 1));
% c_ng_vec = kron(c_ng_vec, ones(60/ST, 1));
% c_e = c_e_vec(tshift:tshift+N_schd-1);
% c_ng = c_ng_vec(tshift:tshift+N_schd-1);

% load boiler, chiller, and fan
load(strcat('boiler/', sample_time, '/boiler.mat'));
load(strcat('chiller/', sample_time, '/chiller.mat'));
load(strcat('fan/fan.mat'));

load limits.mat;
% m_min = Params_office_sizing(:,2);
% m_min(6:10) = m_min(6:10)/10;
% m_max = Params_office_sizing(:,1);
% m_max(6:10) = m_max(6:10)/10;
% Prh_max = Params_office_sizing(:,3);
% Prh_max(6:10) = Prh_max(6:10)/10;

% temperature bounds
T_low = 20*ones(N_zone, 1);
T_hgh = 24*ones(N_zone, 1);

% ourdoor temperature and internal heat gain
% T_out = (T_out_downsample(tshift:tshift+N_schd-1))';
% Q_int = (Q_int_downsample(tshift:tshift+N_schd-1,:))';

Zoneparam.Cp = 1000; % specific heat of air
Zoneparam.Ts = f2c(55)*ones(N_AHU,1); % supply air temperature
Zoneparam.T_approx = 22*ones(N_AHU,1); % approximation of room temperature

alpha_oa = 0.1;
Tmix_vec = alpha_oa*T_out + (1-alpha_oa)*Zoneparam.T_approx(1); % same OA ratio and T_approx for all zones

cvx_begin

    % Decision variables and expressions
    
    variable m_z(N_zone, N_schd)
    variable Prh(N_zone, N_schd)
    variable T_z(N_zone, N_schd)

    expression Pcc_i(N_AHU)
    expression Pcc_total(N_schd)
    
    expression Pf_i(N_AHU)
    expression Pf_total(N_schd)
    
    expression Ppre_total(N_schd)
    expression Pcc_boiler(N_schd)
    
    expression J(N_schd)
    
    m_z_power = m_z;
    m_z_power(6:10,:) = m_z_power(6:10,:) * 10; % mid-floor
    
    Prh_power = Prh;
    Prh_power(6:10,:) = Prh_power(6:10,:) * 10; % mid-floor

    % cost function formulation
    
    for i = 1 : N_schd

        %%  chiller
        
        Tmix_i = Tmix_vec(i);
        for n_f = 1 : N_AHU
            if  1 <= n_f && n_f <= 3
                if  Tmix_i >= Zoneparam.Ts(n_f)
                    Pcc_i(n_f) = sum(m_z_power(5*n_f-4:5*n_f,i)) * Zoneparam.Cp * (Tmix_i-Zoneparam.Ts(n_f));
                else
                    Pcc_i(n_f) = 0;
                end
            end
            if  n_f == 4
                if  Tmix_i >= Zoneparam.Ts(n_f)
                    Pcc_i(n_f) = m_z_power(16,i) * Zoneparam.Cp * (Tmix_i-Zoneparam.Ts(n_f));
                else
                    Pcc_i(n_f) = 0;
                end
            end
        end
        Pcc_total(i) = sum(Pcc_i);
        Pch_i = summer_chiller_result.d0 + summer_chiller_result.d1*T_out(i) + summer_chiller_result.d2*Pcc_total(i);

        %%  fan
        
        for n_f = 1 : N_AHU
            cfan_0 = eval(strcat('summer_fan', num2str(n_f), '_result.c0;'));
            cfan_1 = eval(strcat('summer_fan', num2str(n_f), '_result.c1;'));
            cfan_2 = eval(strcat('summer_fan', num2str(n_f), '_result.c2;'));
            if  1 <= n_f && n_f <= 3
                Pf_i(n_f) = cfan_0 + cfan_1*sum(m_z_power(5*n_f-4:5*n_f,i)) + cfan_2*power(sum(m_z_power(5*n_f-4:5*n_f,i)),2);
            end
            if  n_f == 4
                Pf_i(n_f) = cfan_0 + cfan_1*m_z_power(16,i); % + cfan_2*power(m_z_power(16,i),2);
            end
        end
        Pf_total(i) = sum(Pf_i);

        P_e_i = Pch_i + Pf_total(i); % grid power
        J_e_i = c_e(i) * P_e_i * ST/60;

        %%  boiler
        
        Ppre_total(i) = sum(Prh_power(:,i));
        Pcc_boiler(i) = summerboiler_result.e0 + summerboiler_result.e1*T_out(i) + summerboiler_result.e2*Ppre_total(i);

        J_gas_i = c_ng(i) * Pcc_boiler(i) * ST/60;

        % Total cost
        J(i) = J_e_i + J_gas_i;

    end

    minimize (sum(J))

    subject to

        for i_sch = 1 : N_schd
            
            m_min <= m_z(:, i_sch); % <= m_max;
            zeros(N_zone, 1) <= Prh(:, i_sch); % <= Prh_max;
            
            if  i_sch == 1
                T_z(:, i_sch) == a_0 + a_1*T_out(i_sch) + a_2.*T_ini' + a_3.*m_z(:,i_sch) + a_4.*Prh(:,i_sch) + a_5.*Q_int(:,i_sch);
            else
                T_z(:, i_sch) == a_0 + a_1*T_out(i_sch) + a_2.*T_z(:,i_sch-1) + a_3.*m_z(:,i_sch) + a_4.*Prh(:,i_sch) + a_5.*Q_int(:,i_sch);
            end

            T_low <= T_z(:, i_sch) <= T_hgh;

        end

cvx_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OptSchedule = [];
% OptSchedule.m_z = m_z;
% OptSchedule.Prh = Prh;
% 
% OptStates = [];
% OptStates.T_z = T_z;
% 
% SolverStatus = [];
% SolverStatus.Status = cvx_status;
% SolverStatus.OptVal = cvx_optval;

Optimal_Temp_Ctrl = T_z(:,1);