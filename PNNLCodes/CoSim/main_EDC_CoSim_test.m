clear; clc;

ST = 5; % sampling time: 1, 5, 60 min
sample_time = strcat(num2str(ST), 'min');
num_samples = 24*60/ST; % number of samples

num_zone = 16; % number of zones

top_floor = load(strcat('../Top_floor/', sample_time, '/top_floor.mat'));
mid_floor = load(strcat('../Mid_floor/', sample_time, '/mid_floor.mat'));
bot_floor = load(strcat('../Bot_floor/', sample_time, '/bot_floor.mat'));

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

CoSim_Input = ws2struct();
[OptSchedule,OptStates,SolverStatus] = func_EDC_CoSim_test(CoSim_Input);