% run a scenario for 7 days
% at hourly time step, EDC produces Tset, call ROM to simulate the "real"
% building for one hour

function run_scenario_parallel(scenario_id, season, location, building, elec_price, natgas_price, mismatch)
if ~ischar(scenario_id)
    scenario_id = num2str(scenario_id);
end

if ismac
    Unix = 1;
elseif isunix
    Unix = 1;
elseif ispc
    Unix = 0;
else
    disp('Platform not supported')
end

scenarioID = scenario_id;
season = season;
location = location;
building = building; % 'Office', 'Hotel'
c_e_type = elec_price;
c_ng_type = natgas_price;
Flag_Mismatch = mismatch;

tic

N_steps = 168; % 168 for 7-day simulation

%% Flags
Flag_Plot = 0;
Flag_SingleRun = 1;
Flag_FC = 1;
Flag_Battery = 1;
Flag_Thermal = 1;
Flag_Timestep = 0; % 0 for 1 min, 1 for 5 min
Flag_Integer = 0;
Flag_AS = 0;
%     Flag_Mismatch = 'Zero'; % 'Zero', 'Low', 'Mid', 'High'
Flag_Perfect = 0;

%% Load inputs and parameters
%     season = 'Winter';
%     location = 'Baltimore';
%     scenarioID = '1'; % will create a folder and save the results as .mat

if strcmp(building,'Office')
    
    load Inputs_OfficeROM_Params.mat;
    load OfficeActualInfo.mat;
    load OfficeForecastInfo.mat;
    load OfficePricingInfo.mat;
    if strcmp(location,'Baltimore')
        load Inputs_OfficeROM_sizing_Baltimore.mat;
    elseif strcmp(location,'SanFrancisco')
        load Inputs_OfficeROM_sizing_SanFrancisco.mat;
    end
    load HotWaterEnergyInfo.mat
    
    m_min = Params_office_sizing(:,2);
    m_min(6:10) = m_min(6:10)/10;
    m_max = Params_office_sizing(:,1);
    m_max(6:10) = m_max(6:10)/10;
    Prh_max = Params_office_sizing(:,3);
    Prh_max(6:10) = Prh_max(6:10)/10;
    
    eval(strcat('A_mat = A.',season,'.',location,';'));
    eval(strcat('B_mat = B.',season,'.',location,';'));
    eval(strcat('E_mat = E.',season,'.',location,';'));
    eval(strcat('Fan_mat = Fan.',season,'.',location,';'));
    eval(strcat('Chiller_mat = Chiller.',season,'.',location,';'));
    eval(strcat('Boiler_mat = Boiler.',season,'.',location,';'));
    
    
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
    
    %% Noise
    
    if strcmp(Flag_Mismatch,'Zero')
        noise_sd_temperature = 0;
    elseif strcmp(Flag_Mismatch,'Low')
        noise_sd_temperature = 0.1;
    elseif strcmp(Flag_Mismatch,'Mid')
        noise_sd_temperature = 0.5;
    elseif strcmp(Flag_Mismatch,'High')
        noise_sd_temperature = 1;
    end
    
    %% Prices
    %     c_e_type = 'Low';
    %     c_ng_type = 'High';
    
    eval(strcat('c_e_vec = ElecPrice.',c_e_type,';')); % $/kwh
    eval(strcat('c_ng_vec_mmbtu = GasPrice.',c_ng_type,';')); % $/mmbtu
    c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kwh
    
    c_as_perday = 0.03*ones(N_sch,1)*1e-3; % ancillary service price
    c_as_perday(N_sch*2/3:N_sch*5/6) = 0.05*1e-3;
    c_as_vec = repmat(c_as_perday,10,1);
    c_up = 3; % FC start-up cost
    c_down = 3; % FC shut-down cost
    
    %% Parameters (tunable)
    eta_d = 0.9; % battery discharging efficiency
    eta_c = 0.9; % battery charging efficiency
    eta_th_d = 0.9; % thermal storage discharging efficiency
    eta_th_c = 0.9; % thermal storage charging efficiency
    eta_th_loss = 0.05; % thermal storage static loss
    eta_fc_e = 0.5; % electrical efficiency of the FC
    eta_fc_h = 0.35; % heat efficiecny of the FC
    kas = 1e-4; % coefficent to calculate temperature buffer from ancillary service
    Ras_max = 1e4; % max AS capacity
    
    T_high = 24*ones(15,1); % temperature bounds
    T_low = 20*ones(15,1);
    Eb_max = 1e4*3600*Flag_Battery; % battery energy capacity
    Eth_max = 1e4*3600*0.5*Flag_Thermal; % thermal storage energy capacity
    Pfc_f_max = 1e4*Flag_FC; % FC fuel consumption capacity
    Pfc_f_min = 0.2*Pfc_f_max*Flag_FC; % FC fuel consumption minimum operating point
    Pb_max = 1e4*Flag_Battery; % battery power capacity
    Pth_max = 1e4*0.5*Flag_Thermal; % thermal storage power capacity
    
    % season = 'Winter';
    % location = 'Baltimore';
    % if strcmp(location,'Baltimore')
    %     load ROM_June14_office_sizing_Baltimore.mat
    % elseif strcmp(location,'SanFrancisco')
    %     load ROM_June14_office_sizing_SF.mat
    % end
    % m_min = Params_office_sizing(:,2);
    % m_min(6:10) = m_min(6:10)/10;
    % m_max = Params_office_sizing(:,1)*2;
    % m_max(6:10) = m_max(6:10)/10;
    % Prh_max = Params_office_sizing(:,3);
    %
    % load ROM_June14_office.mat
    % eval(strcat('A_mat = A.',season,'.',location,';'));
    % eval(strcat('B_mat = B.',season,'.',location,';'));
    % eval(strcat('E_mat = E.',season,'.',location,';'));
    % eval(strcat('Fan_mat = Fan.',season,'.',location,';'));
    % eval(strcat('Chiller_mat = Chiller.',season,'.',location,';'));
    % eval(strcat('Boiler_mat = Boiler.',season,'.',location,';'));
    
    FC_ramp_up = 0.25*Pfc_f_max*Flag_FC; % FC ramp rate limits
    FC_ramp_down = 0.25*Pfc_f_max*Flag_FC;
    
    min_up = 2; % FC minimum up time
    min_down = 2; % FC minimum down time
    
    
    
    %% Params for the "real" ROM
    % load Data_temp.mat;
    
    %% Initial conditions
    X0_init = 22*ones(15,1);  % '15' stands for the dimension of the state space model
    Pfc_f0_init = 0;
    Pfc_h0_init = 0;
    Eb0_init = Eb_max*0.5*Flag_Battery;
    Eth0_init = Eth_max*0.5*Flag_Thermal;
    w0_init = 0;
    
    for i = 1:N_steps
        
        %% Initial conditions
        if i == 1
            X0 = X0_init;  % '15' stands for the dimension of the state space model
            Pfc_f0 = Pfc_f0_init;
            Pfc_h0 = Pfc_h0_init;
            Eb0 = Eb0_init;
            Eth0 = Eth0_init;
            w0 = w0_init;
        else
            if strcmp(SolverStatus.Status,'Solved') || strcmp(SolverStatus.Status,'Inaccurate/Solved') % if not solved, use last step's initial values
                X0 = ROM_Results.x(:,end);  % '15' stands for the dimension of the state space model
                Pfc_f0 = OptSchedule.u(1,3);
                Pfc_h0 = OptSchedule.u(1,4);
                Eb0 = OptStates.Eb(1,1,end);
                Eth0 = OptStates.Eth(1,1,end);
                if Flag_Integer == 1
                    w0 = OptStates.w(1);
                else
                    w0 = w0_init;
                end
            end
        end
        
        
        %% Inputs
        
        %     eval(strcat('load ROM_June14_office_vali_inputs_',season,'.mat;'));
        eval(strcat('T_amb_vec = DryBulbWeatherActual.',season,'.',location,';'));
        eval(strcat('Qint_cell = InternalLoadActual.',season,'.',location,';'));
        Qint_vec = [Qint_cell.Peri_bot1,Qint_cell.Peri_bot2,Qint_cell.Peri_bot3,Qint_cell.Peri_bot4,Qint_cell.Core_bot,...
            Qint_cell.Peri_mid1,Qint_cell.Peri_mid2,Qint_cell.Peri_mid3,Qint_cell.Peri_mid4,Qint_cell.Core_mid,...
            Qint_cell.Peri_top1,Qint_cell.Peri_top2,Qint_cell.Peri_top3,Qint_cell.Peri_top4,Qint_cell.Core_top];
        eval(strcat('Pe_nonHVAC_vec = NonHVACElectricActual.',season,'.',location,'*1e3;'));
        eval(strcat('P_HotWater_vec = HotWaterEnergyActual.',building,'.',season,'.',location,';'));
        
        eval(strcat('T_amb_cell = DryBulbWeatherForecast.',season,'.',location,';'));
        eval(strcat('Qint_cell_Forecast = InternalLoadForecast.',season,'.',location,';'));
        
        tshift = i;
        T_amb = T_amb_vec(tshift:(tshift+23));
        Qint = Qint_vec(tshift:(tshift+23),:);
        Pe_nonHVAC = Pe_nonHVAC_vec(tshift:(tshift+23));
        P_HotWater = P_HotWater_vec(tshift:(tshift+23));
        T_amb_Forecase = T_amb_cell(i,:)';
        Qint_Forecast = [Qint_cell_Forecast.Peri_bot1(i,:)',Qint_cell_Forecast.Peri_bot2(i,:)',Qint_cell_Forecast.Peri_bot3(i,:)',Qint_cell_Forecast.Peri_bot4(i,:)',Qint_cell_Forecast.Core_bot(i,:)'...
            Qint_cell_Forecast.Peri_mid1(i,:)',Qint_cell_Forecast.Peri_mid2(i,:)',Qint_cell_Forecast.Peri_mid3(i,:)',Qint_cell_Forecast.Peri_mid4(i,:)',Qint_cell_Forecast.Core_mid(i,:)'...
            Qint_cell_Forecast.Peri_top1(i,:)',Qint_cell_Forecast.Peri_top2(i,:)',Qint_cell_Forecast.Peri_top3(i,:)',Qint_cell_Forecast.Peri_top4(i,:)',Qint_cell_Forecast.Core_top(i,:)'];
        Pe_nonHVAC_Forecast = Pe_nonHVAC;
        P_HotWater_Forecast = P_HotWater;
        
        % If EDC perfect knowledge
        if Flag_Perfect
            T_amb_Forecase = T_amb;
            Qint_Forecast = Qint;
        end
        
        
        %% Noise
        
        noise_temperature = randn(15,60).*noise_sd_temperature/sqrt(60);
        
        %% Collect arguments for EDC
        
        InitialStates = [];
        InitialStates.X = X0;
        InitialStates.Pfc_f = Pfc_f0;
        InitialStates.Pfc_h = Pfc_h0;
        InitialStates.Eb = Eb0;
        InitialStates.Eth = Eth0;
        InitialStates.w = w0;
        
        Inputs_EDC = [];
        Inputs_EDC.T_amb = T_amb_Forecase;
        Inputs_EDC.Qint = Qint_Forecast;
        Inputs_EDC.Pe_nonHVAC = Pe_nonHVAC_Forecast;
        Inputs_EDC.P_HotWater = P_HotWater_Forecast;
        
        Params = [];
        Params.eta_d = eta_d;
        Params.eta_c = eta_c;
        Params.eta_th_d = eta_th_d;
        Params.eta_th_c = eta_th_c;
        Params.eta_th_loss = eta_th_loss;
        Params.eta_fc_e = eta_fc_e;
        Params.eta_fc_h = eta_fc_h;
        Params.kas = kas;
        Params.Ras_max = Ras_max;
        
        Params.T_high = T_high;
        Params.T_low = T_low;
        Params.Eb_max = Eb_max;
        Params.Eth_max = Eth_max;
        Params.Pfc_f_max = Pfc_f_max;
        Params.Pfc_f_min = Pfc_f_min;
        Params.Pb_max = Pb_max;
        Params.Pth_max = Pth_max;
        
        Params.m_min = m_min;
        Params.m_max = m_max;
        Params.Prh_max = Prh_max;
        Params.FC_ramp_up = FC_ramp_up;
        Params.FC_ramp_down = FC_ramp_down;
        Params.min_up = min_up;
        Params.min_down = min_down;
        
        Params.season = season;
        Params.location = location;
        
        Params.A_mat = A_mat;
        Params.B_mat = B_mat;
        Params.E_mat = E_mat;
        Params.Fan_mat = Fan_mat;
        Params.Chiller_mat = Chiller_mat;
        Params.Boiler_mat = Boiler_mat;
        
        Params.noise_temperature = noise_temperature;
        
        Prices = [];
        Prices.c_e = c_e_vec(tshift:tshift+23);
        Prices.c_ng = c_ng_vec(tshift:tshift+23);
        Prices.c_as = c_as_vec(tshift:tshift+23);
        Prices.c_up = c_up;
        Prices.c_down = c_down;
        
        Flags = [];
        Flags.Flag_Plot = Flag_Plot;
        Flags.Flag_SingleRun = Flag_SingleRun;
        Flags.Flag_FC = Flag_FC;
        Flags.Flag_Battery = Flag_Battery;
        Flags.Flag_Thermal = Flag_Thermal;
        Flags.Flag_Timestep = Flag_Timestep;
        Flags.Flag_AS = Flag_AS;
        
        
        %% Run the EDC
        
        if Flag_Integer == 1
            [OptSchedule,OptStates,SolverStatus] = func_EDC_Office_v1(Flags,Prices,Inputs_EDC,InitialStates,Params);
        else
            [OptSchedule,OptStates,SolverStatus] = func_EDC_Office_v1_noInt(Flags,Prices,Inputs_EDC,InitialStates,Params);
        end
        
        %% Run ROM
        
        Inputs_ROM = [];
        if strcmp(SolverStatus.Status,'Solved') || strcmp(SolverStatus.Status,'Inaccurate/Solved')
            Inputs_ROM.T_set = OptStates.x(1,:,end)';
        else
            Inputs_ROM.T_set = ones(15,1)*22;
        end
        if i == 1
            Inputs_ROM.X0 = X0;
        else
            Inputs_ROM.X0 = ROM_Results.x(:,end);
        end
        Inputs_ROM.T_amb = T_amb(1);
        Inputs_ROM.Qint = Qint(1,:);
        ROM_Results = func_runROM_Office_1hour(Params,Inputs_ROM);
        if ~(strcmp(ROM_Results.SolverStatus,'Solved') || strcmp(ROM_Results.SolverStatus,'Inaccurate/Solved'))
            ROM_Results.x(:,end) = Inputs_ROM.T_set;
        end
        
        %% Run ROM default case
        
        Inputs_ROM = [];
        Inputs_ROM.T_set = ones(15,1)*22;
        if i == 1
            Inputs_ROM.X0 = X0;
        else
            Inputs_ROM.X0 = ROM_Results_default.x(:,end);
        end
        Inputs_ROM.T_amb = T_amb(1);
        Inputs_ROM.Qint = Qint(1,:);
        ROM_Results_default = func_runROM_Office_1hour(Params,Inputs_ROM);
        if ~(strcmp(ROM_Results_default.SolverStatus,'Solved') || strcmp(ROM_Results_default.SolverStatus,'Inaccurate/Solved'))
            ROM_Results_default.x(:,end) = Inputs_ROM.T_set;
        end
        
        %% Record data
        ROM_Results_output{i,1} = ROM_Results;
        ROM_Results_default_output{i,1} = ROM_Results_default;
        OptSchedule_output{i,1} = OptSchedule;
        OptStates_output{i,1} = OptStates;
        SolverStatus_output{i,1} = SolverStatus;
        noise_temperature_output{1,1} = noise_temperature;
        
            
        %% Post process
        
        % Create folder and save results
        if Unix == 1
            eval(strcat('mkdir ./Results/Scenario_',scenarioID))
            eval(strcat('save(''./Results/Scenario_',scenarioID,'/ScenarioResults.mat'')'))
        else
            eval(strcat('mkdir .\Results\Scenario_',scenarioID))
            eval(strcat('save(''.\Results\Scenario_',scenarioID,'\ScenarioResults.mat'')'))
        end
        
    end
   
elseif strcmp(building,'Hotel')
    
    load Inputs_HotelROM_Params.mat;
    % load Inputs_Scenario.mat;
    load HotelActualInfo.mat;
    load HotelForecastInfo.mat;
    load HotelPricingInfo.mat;
    if strcmp(location,'Baltimore')
        load Inputs_HotelROM_sizing_Baltimore.mat
    elseif strcmp(location,'SanFrancisco')
        load Inputs_HotelROM_sizing_SanFrancisco.mat
    end
    load HotWaterEnergyInfo.mat
    
    m_min = Params_hotel_sizing(:,2);
    m_min(9) = m_min(9)/168;
    m_min(10) = m_min(10)/4;
    m_min(11) = m_min(11)/11;
    m_max = Params_hotel_sizing(:,1)*10;
    m_max(9) = m_max(9)/168;
    m_max(10) = m_max(10)/4;
    m_max(11) = m_max(11)/11;
    Prh_max = Params_hotel_sizing(:,3)*10;
    Prh_max(9) = Prh_max(9)/168;
    Prh_max(10) = Prh_max(10)/4;
    Prh_max(11) = Prh_max(11)/11;
    
    
    eval(strcat('A_mat = A.',season,'.',location,';'));
    eval(strcat('B_mat = B.',season,'.',location,';'));
    eval(strcat('E_mat = E.',season,'.',location,';'));
    eval(strcat('Fan_mat = Fan.',season,'.',location,';'));
    eval(strcat('Chiller_mat = Chiller.',season,'.',location,';'));
    eval(strcat('Boiler_mat = Boiler.',season,'.',location,';'));
    
    
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
    
    %% Noise
    
    if strcmp(Flag_Mismatch,'Zero')
        noise_sd_temperature = 0;
    elseif strcmp(Flag_Mismatch,'Low')
        noise_sd_temperature = 0.1;
    elseif strcmp(Flag_Mismatch,'Mid')
        noise_sd_temperature = 0.5;
    elseif strcmp(Flag_Mismatch,'High')
        noise_sd_temperature = 1;
    end
    
    %% Prices
%     c_e_type = 'Low';
%     c_ng_type = 'High';
    
    eval(strcat('c_e_vec = ElecPrice.',c_e_type,';')); % $/kwh
    eval(strcat('c_ng_vec_mmbtu = GasPrice.',c_ng_type,';')); % $/mmbtu
    c_ng_vec = c_ng_vec_mmbtu/293.07; % $/kwh
    
    c_as_perday = 0.03*ones(N_sch,1)*1e-3; % ancillary service price
    c_as_perday(N_sch*2/3:N_sch*5/6) = 0.05*1e-3;
    c_as_vec = repmat(c_as_perday,10,1);
    c_up = 3; % FC start-up cost
    c_down = 3; % FC shut-down cost
    
    %% Parameters (tunable)
    eta_d = 0.9; % battery discharging efficiency
    eta_c = 0.9; % battery charging efficiency
    eta_th_d = 0.9; % thermal storage discharging efficiency
    eta_th_c = 0.9; % thermal storage charging efficiency
    eta_th_loss = 0.05; % thermal storage static loss
    eta_fc_e = 0.5; % electrical efficiency of the FC
    eta_fc_h = 0.35; % heat efficiecny of the FC
    kas = 1e-4; % coefficent to calculate temperature buffer from ancillary service
    Ras_max = 1e4; % max AS capacity
    
    T_high = 24*ones(15,1); % temperature bounds
    T_low = 20*ones(15,1);
    Eb_max = 1e4*3600*Flag_Battery; % battery energy capacity
    Eth_max = 1e4*3600*0.5*Flag_Thermal; % thermal storage energy capacity
    Pfc_f_max = 1e4*Flag_FC; % FC fuel consumption capacity
    Pfc_f_min = 0.2*Pfc_f_max*Flag_FC; % FC fuel consumption minimum operating point
    Pb_max = 1e4*Flag_Battery; % battery power capacity
    Pth_max = 1e4*0.5*Flag_Thermal; % thermal storage power capacity
    
    % season = 'Winter';
    % location = 'Baltimore';
    % if strcmp(location,'Baltimore')
    %     load ROM_June14_office_sizing_Baltimore.mat
    % elseif strcmp(location,'SanFrancisco')
    %     load ROM_June14_office_sizing_SF.mat
    % end
    % m_min = Params_office_sizing(:,2);
    % m_min(6:10) = m_min(6:10)/10;
    % m_max = Params_office_sizing(:,1)*2;
    % m_max(6:10) = m_max(6:10)/10;
    % Prh_max = Params_office_sizing(:,3);
    %
    % load ROM_June14_office.mat
    % eval(strcat('A_mat = A.',season,'.',location,';'));
    % eval(strcat('B_mat = B.',season,'.',location,';'));
    % eval(strcat('E_mat = E.',season,'.',location,';'));
    % eval(strcat('Fan_mat = Fan.',season,'.',location,';'));
    % eval(strcat('Chiller_mat = Chiller.',season,'.',location,';'));
    % eval(strcat('Boiler_mat = Boiler.',season,'.',location,';'));
    
    FC_ramp_up = 0.25*Pfc_f_max*Flag_FC; % FC ramp rate limits
    FC_ramp_down = 0.25*Pfc_f_max*Flag_FC;
    
    min_up = 2; % FC minimum up time
    min_down = 2; % FC minimum down time
    
    
    
    %% Params for the "real" ROM
    % load Data_temp.mat;
    
    %% Initial conditions
    X0_init = 22*ones(15,1);  % '15' stands for the dimension of the state space model
    Pfc_f0_init = 0;
    Pfc_h0_init = 0;
    Eb0_init = Eb_max*0.5*Flag_Battery;
    Eth0_init = Eth_max*0.5*Flag_Thermal;
    w0_init = 0;
    
    for i = 1:N_steps
        
        %% Initial conditions
        if i == 1
            X0 = X0_init;  % '15' stands for the dimension of the state space model
            Pfc_f0 = Pfc_f0_init;
            Pfc_h0 = Pfc_h0_init;
            Eb0 = Eb0_init;
            Eth0 = Eth0_init;
            w0 = w0_init;
        else
            if strcmp(SolverStatus.Status,'Solved') || strcmp(SolverStatus.Status,'Inaccurate/Solved') % if not solved, use last step's initial values
                X0 = ROM_Results.x(:,end);  % '15' stands for the dimension of the state space model
                Pfc_f0 = OptSchedule.u(1,3);
                Pfc_h0 = OptSchedule.u(1,4);
                Eb0 = OptStates.Eb(1,1,end);
                Eth0 = OptStates.Eth(1,1,end);
                if Flag_Integer == 1
                    w0 = OptStates.w(1);
                else
                    w0 = w0_init;
                end
            end
        end
        
        
        %% Inputs
        
        %     eval(strcat('load ROM_June14_office_vali_inputs_',season,'.mat;'));
        eval(strcat('T_amb_vec = DryBulbWeatherActual.',season,'.',location,';'));
        eval(strcat('Qint_cell = InternalLoadActual.',season,'.',location,';'));
        Qint_vec = [Qint_cell.Basement,Qint_cell.Storage,Qint_cell.Retail2,Qint_cell.Retail1,Qint_cell.Mechroom,...
            Qint_cell.Lobby,Qint_cell.Cafe,Qint_cell.Laundry,mean(Qint_cell.Guestroom2,2),Qint_cell.Corridor2,...
            mean(Qint_cell.Guestroom3,2),Qint_cell.Banquet,Qint_cell.Corridor3,Qint_cell.Dining,Qint_cell.Kitchen];
        eval(strcat('Pe_nonHVAC_vec = NonHVACElectricActual.',season,'.',location,'*1e3;'));
        eval(strcat('P_HotWater_vec = HotWaterEnergyActual.',building,'.',season,'.',location,';'));
        
        eval(strcat('T_amb_cell = DryBulbWeatherForecast.',season,'.',location,';'));
        eval(strcat('Qint_cell_Forecast = InternalLoadForecast.',season,'.',location,';'));
        
        tshift = i;
        T_amb = T_amb_vec(tshift:(tshift+23));
        Qint = Qint_vec(tshift:(tshift+23),:);
        Pe_nonHVAC = Pe_nonHVAC_vec(tshift:(tshift+23));
        P_HotWater = P_HotWater_vec(tshift:(tshift+23));
        T_amb_Forecase = T_amb_cell(i,:)';
        Qint_cell_Forecast_Guestroom2 = (Qint_cell_Forecast.Guestroom2.Room1(i,:)' + Qint_cell_Forecast.Guestroom2.Room2(i,:)' + Qint_cell_Forecast.Guestroom2.Room3(i,:)' + ...
            Qint_cell_Forecast.Guestroom2.Room4(i,:)' + Qint_cell_Forecast.Guestroom2.Room5(i,:)' + Qint_cell_Forecast.Guestroom2.Room6(i,:)')/6;
        Qint_cell_Forecast_Guestroom3 = (Qint_cell_Forecast.Guestroom3.Room1(i,:)' + Qint_cell_Forecast.Guestroom3.Room2(i,:)' + Qint_cell_Forecast.Guestroom3.Room3(i,:)')/3;
        Qint_Forecast = [Qint_cell_Forecast.Basement(i,:)',Qint_cell_Forecast.Storage(i,:)',Qint_cell_Forecast.Retail2(i,:)',Qint_cell_Forecast.Retail1(i,:)',Qint_cell_Forecast.Mechroom(i,:)'...
            Qint_cell_Forecast.Lobby(i,:)',Qint_cell_Forecast.Cafe(i,:)',Qint_cell_Forecast.Laundry(i,:)',Qint_cell_Forecast_Guestroom2,Qint_cell_Forecast.Corridor2(i,:)'...
            Qint_cell_Forecast_Guestroom3,Qint_cell_Forecast.Banquet(i,:)',Qint_cell_Forecast.Corridor3(i,:)',Qint_cell_Forecast.Dining(i,:)',Qint_cell_Forecast.Kitchen(i,:)'];
        Pe_nonHVAC_Forecast = Pe_nonHVAC;
        P_HotWater_Forecast = P_HotWater;
        
        
        % If EDC perfect knowledge
        if Flag_Perfect
            T_amb_Forecase = T_amb;
            Qint_Forecast = Qint;
        end
        
        
        %% Noise
        
        noise_temperature = randn(15,60).*noise_sd_temperature/sqrt(60);
        
        %% Collect arguments for EDC
        
        InitialStates = [];
        InitialStates.X = X0;
        InitialStates.Pfc_f = Pfc_f0;
        InitialStates.Pfc_h = Pfc_h0;
        InitialStates.Eb = Eb0;
        InitialStates.Eth = Eth0;
        InitialStates.w = w0;
        
        Inputs_EDC = [];
        Inputs_EDC.T_amb = T_amb_Forecase;
        Inputs_EDC.Qint = Qint_Forecast;
        Inputs_EDC.Pe_nonHVAC = Pe_nonHVAC_Forecast;
        Inputs_EDC.P_HotWater = P_HotWater_Forecast;
        
        Params = [];
        Params.eta_d = eta_d;
        Params.eta_c = eta_c;
        Params.eta_th_d = eta_th_d;
        Params.eta_th_c = eta_th_c;
        Params.eta_th_loss = eta_th_loss;
        Params.eta_fc_e = eta_fc_e;
        Params.eta_fc_h = eta_fc_h;
        Params.kas = kas;
        Params.Ras_max = Ras_max;
        
        Params.T_high = T_high;
        Params.T_low = T_low;
        Params.Eb_max = Eb_max;
        Params.Eth_max = Eth_max;
        Params.Pfc_f_max = Pfc_f_max;
        Params.Pfc_f_min = Pfc_f_min;
        Params.Pb_max = Pb_max;
        Params.Pth_max = Pth_max;
        
        Params.m_min = m_min;
        Params.m_max = m_max;
        Params.Prh_max = Prh_max;
        Params.FC_ramp_up = FC_ramp_up;
        Params.FC_ramp_down = FC_ramp_down;
        Params.min_up = min_up;
        Params.min_down = min_down;
        
        Params.season = season;
        Params.location = location;
        
        Params.A_mat = A_mat;
        Params.B_mat = B_mat;
        Params.E_mat = E_mat;
        Params.Fan_mat = Fan_mat;
        Params.Chiller_mat = Chiller_mat;
        Params.Boiler_mat = Boiler_mat;
        
        Params.noise_temperature = noise_temperature;
        
        Prices = [];
        Prices.c_e = c_e_vec(tshift:tshift+23);
        Prices.c_ng = c_ng_vec(tshift:tshift+23);
        Prices.c_as = c_as_vec(tshift:tshift+23);
        Prices.c_up = c_up;
        Prices.c_down = c_down;
        
        Flags = [];
        Flags.Flag_Plot = Flag_Plot;
        Flags.Flag_SingleRun = Flag_SingleRun;
        Flags.Flag_FC = Flag_FC;
        Flags.Flag_Battery = Flag_Battery;
        Flags.Flag_Thermal = Flag_Thermal;
        Flags.Flag_Timestep = Flag_Timestep;
        Flags.Flag_AS = Flag_AS;
        
        
        %% Run the EDC
        
        if Flag_Integer == 1
            [OptSchedule,OptStates,SolverStatus] = func_EDC_Hotel_v1(Flags,Prices,Inputs_EDC,InitialStates,Params);
        else
            [OptSchedule,OptStates,SolverStatus] = func_EDC_Hotel_v1_noInt(Flags,Prices,Inputs_EDC,InitialStates,Params);
        end
        
        %% Run ROM
        
        Inputs_ROM = [];
        if strcmp(SolverStatus.Status,'Solved') || strcmp(SolverStatus.Status,'Inaccurate/Solved')
            Inputs_ROM.T_set = OptStates.x(1,:,end)';
        else
            Inputs_ROM.T_set = ones(15,1)*22;
        end
        if i == 1
            Inputs_ROM.X0 = X0;
        else
            Inputs_ROM.X0 = ROM_Results.x(:,end);
        end
        Inputs_ROM.T_amb = T_amb(1);
        Inputs_ROM.Qint = Qint(1,:);
        ROM_Results = func_runROM_Hotel_1hour(Params,Inputs_ROM);
        if ~(strcmp(ROM_Results.SolverStatus,'Solved') || strcmp(ROM_Results.SolverStatus,'Inaccurate/Solved'))
            ROM_Results.x(:,end) = Inputs_ROM.T_set;
        end
        
        %% Run ROM default case
        
        Inputs_ROM = [];
        Inputs_ROM.T_set = ones(15,1)*22;
        if i == 1
            Inputs_ROM.X0 = X0;
        else
            Inputs_ROM.X0 = ROM_Results_default.x(:,end);
        end
        Inputs_ROM.T_amb = T_amb(1);
        Inputs_ROM.Qint = Qint(1,:);
        ROM_Results_default = func_runROM_Hotel_1hour(Params,Inputs_ROM);
        if ~(strcmp(ROM_Results_default.SolverStatus,'Solved') || strcmp(ROM_Results_default.SolverStatus,'Inaccurate/Solved'))
            ROM_Results_default.x(:,end) = Inputs_ROM.T_set;
        end
        
        %% Record data
        ROM_Results_output{i,1} = ROM_Results;
        ROM_Results_default_output{i,1} = ROM_Results_default;
        OptSchedule_output{i,1} = OptSchedule;
        OptStates_output{i,1} = OptStates;
        SolverStatus_output{i,1} = SolverStatus;
        noise_temperature_output{1,1} = noise_temperature;
        
        %% Post process
        
        % Create folder and save results
        if Unix == 1
            eval(strcat('mkdir ./Results/Scenario_',scenarioID))
            eval(strcat('save(''./Results/Scenario_',scenarioID,'/ScenarioResults.mat'')'))
        else
            eval(strcat('mkdir .\Results\Scenario_',scenarioID))
            eval(strcat('save(''.\Results\Scenario_',scenarioID,'\ScenarioResults.mat'')'))
        end
        
    end
    
end

toc



