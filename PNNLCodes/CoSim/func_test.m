function [T_set, sol_status] = func_test(ST, T_out, T_ini, Q_int, eng_price, gas_price)

ST_size = size(ST)

T_out_size = size(T_out)

T_ini_size = size(T_ini)

Q_int_size = size(Q_int)

eng_price_size = size(eng_price)

gas_price_size = size(gas_price)

T_set = 23*ones(16,1);

sol_status = 'solved';