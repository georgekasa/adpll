
%https://people.engr.tamu.edu/spalermo/ecen620/2007_design_procedure_adpll_kratyuk_tcas2.pdf
% A Design Procedure for All-Digital Phase-Locked
%Loops Based on a Charge-Pump Phase-Locked-Loop Analogy
clc
clear
Fref = 48e6;
Fdco = 1.20e9;
Tref = 1/Fref;
Ndiv = Fdco/Fref;
tdc_res = 55e-12; % TDC time resolution
Kdco = 700e3;
phaseMargin = 40.0; % In degrees should be at least over 60 in verilog
% Unity gain bandwidth
Wugbw = 2*pi*0.2e6;
Ndiv = Fdco/Fref;
% Zero frequency calculation
Wz = Wugbw/tan(2*pi*phaseMargin/360);
% Equivalent Icp calculation
Icp = Tref/tdc_res;
% Equivalent resistance R & C
R = (2.0*pi*Ndiv/(Icp*Kdco)) * ((Wz^2)/sqrt(Wz^2 + (Wugbw)^2));
C = tan(2.0*pi*phaseMargin/360.0) / (R*Wugbw);

alpha = R - Tref / (2*C); %eq 11
beta = Tref / C; % eq 12

disp(alpha)
disp(beta)

[rounded_alpha, closest_value_alpha] = closest_power_of_two_inverse(alpha);
fprintf("The closest power Alpha of 2^-x for %f is 2^-%d = %f\n", alpha, rounded_alpha, closest_value_alpha);



[rounded_beta, closest_value_beta] = closest_power_of_two_inverse(beta);
fprintf("The closest power beta of 2^-x for %f is 2^-%d = %f\n", beta, rounded_beta, closest_value_beta);


function [rounded_x, closest_value] = closest_power_of_two_inverse(number)
    if number <= 0
        error("The number must be positive.");
    end

    % Calculate x = -log2(number)
    x = -log2(number);

    % Round x to the nearest integer
    rounded_x = round(x);

    % Calculate the closest 2^-x
    closest_value = 2^(-rounded_x);
end
