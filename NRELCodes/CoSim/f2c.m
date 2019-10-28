function c = f2c(f)
% Convert temperature in F to Kelvin
% 1 input: f - temperature in F, can be scalar or vector
% 1 output: c - temperature in C, same size as input

c = (f-32).*5./9;