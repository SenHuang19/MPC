function f = c2f(c)
% Convert temperature in Kelvin to F
% 1 input: k - temperature in C, can be scalar or vector
% 1 output: f - temperature in F, same size as input
%#eml

f = c.*9./5 + 32;