#REFERENCES

#1) BOOK of Staszewski's
#2) for vco spurs Supply-Insensitive Frequency Synthesis for an LDO-Free Powering Scheme in SoCs
#3) dtc/tdc spurs =DTC and TDC IC Design for Ultra-Low-Power ADPLL Master of Science Thesis Peng Chen

clc
clear
pkg load control

Freference = 26e6;
Fdco = 1.8e9;
Tref = 1.0/Freference;
Tdco = 1.0/Fdco;
FCW = Fdco/Freference;
FCW_frac = abs(floor(FCW)-FCW);

Kpush_vco = 1e9;
Aripple = 20e-3;#vp2p
fripple = 10e9;
INL_pp = 0.1;


ref_pn = -130;
%TDC & DCO spec.
dtc_res = 40e-12;
tdc_res = 20e-12;
kdco = 20e3;
%Filter parameters
rho_filter = 2^(-15);
alpha_filter = 2^(-7);
lamda1 = 2^(-3);
lamda2 = lamda1;
lamda3 = lamda1;
lamda4 = 2^(-4);


f = logspace(1, 10, 100); % Frequency range from 10 Hz to 10 GHz
start_index_jitter = find(f >= 1e3, 1);
end_index_jitter = find(f <= 1e9, 1, 'last');


s = 1j * 2 * pi * f; % Complex frequency variable
%Quantization noise TDC page 121 book, eq 4.58
L_tdc_quant = (((2*pi)^2)/12) * ((tdc_res/Tdco)^2) *Tref;


%tranfer Function for high order ADPLL (IIR + PI)

%Hol page 137 book, eq 4.95
Hol_numerator = (rho_filter*Freference^2).*(s./((rho_filter*Freference)/alpha_filter) +1);
Hol_denominator = s.*s;%1*s2
Hol_tf = Hol_numerator./Hol_denominator;
H_closed_tdc_NoIIR = L_tdc_quant.*Hol_tf./(1+Hol_tf);
%Hol = (alpha * f_ref) ./ s + (ro * f_ref^2) ./ s.^2;

% H_iir1 Filter
H_iir1_filter_num = 1+s./Freference ;
H_iir1_filter_den = (s./(lamda1 .* Freference)) +1;
H_iir1_filter_tf = H_iir1_filter_num./ H_iir1_filter_den;

% H_iir2 Filter
H_iir2_filter_num = 1+s./Freference ;
H_iir2_filter_den = (s./(lamda2 * Freference)) +1;
H_iir2_filter_tf = H_iir2_filter_num./ H_iir2_filter_den;

% H_iir3 Filter
H_iir3_filter_num = 1+s./Freference ;
H_iir3_filter_den = (s./(lamda3 * Freference)) +1;
H_iir3_filter_tf = H_iir3_filter_num./ H_iir3_filter_den;

% H_iir4 Filter
H_iir4_filter_num = 1+s./Freference ;
H_iir4_filter_den = (s./(lamda4 * Freference)) +1;
H_iir4_filter_tf = H_iir4_filter_num./ H_iir4_filter_den;

% Multiply all the transfer functions with IIR filters
Hol_total_iir = Hol_tf .*(H_iir1_filter_tf .* H_iir2_filter_tf .* H_iir3_filter_tf .* H_iir4_filter_tf);
H_closed_tdc = Hol_total_iir./(1+Hol_total_iir);
% CLosed loop transfer function of the TDC multplied with Quant noise eq 4.97
H_closed_quant = L_tdc_quant.*(abs(H_closed_tdc).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%DCO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DCO quantization noise (without Sigma-Delta dithering)
%extremely important to compare this with the "natural" PN of the DCO itself!!!!
%eq 4.72
L_dco_quant = (1/12).*((kdco./f).^2).*(1/Freference).*sinc(kdco./Freference)^2;


%Transfer Function of DCO eq 4.98 NO MASH
Hclosed_dco = 1.0./(1.0 + Hol_total_iir);
Hclosed_dco_quant = (abs(Hclosed_dco).^2).*L_dco_quant;


%Transfer Function of reference eq 4.96
Hclosed_ref = FCW.*Hol_total_iir./(1.0 + Hol_total_iir);
reference_pn_lin = 10.^(ref_pn./10);
Hclosed_ref_noise = reference_pn_lin.*abs(Hclosed_ref).^2;

adpll_total_pn = H_closed_quant + Hclosed_dco_quant + Hclosed_ref_noise;

%Noises shaped by the H(s)
%L_tdc = 10*log10(L_tdc_quant);
%PN_tdc_noIIR = 20*log10(abs(H_closed_tdc_NoIIR));
PN_tdc_IIR = 10.*log10(abs(H_closed_quant));
%L_dco = 10*log10(L_dco_quant);
PN_dco_IIR = 10.*log10(abs(Hclosed_dco_quant));
PN_ref_IIR = 10.*log10(abs(Hclosed_ref_noise));
% Total ADPLL noise PSD

adpll_pn = 10.*log10(adpll_total_pn);

figure;

hold on
%semilogx(f, PN_tdc_noIIR, 'b', 'LineWidth', 1.5);
semilogx(f, PN_tdc_IIR, 'b', 'LineWidth', 1.5);
%semilogx(f, L_dco, 'y', 'LineWidth', 1.5);
semilogx(f, PN_dco_IIR, 'g', 'LineWidth', 1.5);
semilogx(f, PN_ref_IIR, 'r', 'LineWidth', 1.5);
semilogx(f, adpll_pn, 'black', 'LineWidth', 1.5);
xlabel('Frequency Offset (Hz)');
ylabel('Phase Noise (dBc/Hz)');
grid on;

legend({'TDC PN (IIR)', 'DCO PN (IIR)', 'Reference PN', 'ADPLL PN'}, 'Location', 'northeast');




#jitter razavi pll page 50-51, eq 2.32
% Find the indices corresponding to frequency range from 1 kHz to 1 GHz

integrated_noise_adpll_radians = sqrt(2*trapz(f(start_index_jitter:end_index_jitter),
                  adpll_total_pn(start_index_jitter:end_index_jitter)));#jrms this is in radians

integrated_noise_adpll_sec = Tdco*integrated_noise_adpll_radians/(2*pi);
disp(["Integrated jitter: ", num2str(integrated_noise_adpll_sec)]);


#spur of vco pushing
spur_vco = 20*log10(Kpush_vco*Aripple/(4*fripple));
%spurs from TDC from quantization
spur_loc_TDC = (FCW_frac*Tdco/tdc_res)*Freference;

%spurs DTC/TDC INL
%As shown in Figure 2-6, IN Lpp means when DTC digital input control code is 16, INL is 0.1
%LSB; when the control code is 48, INL is -0.1 LSB. After DTC, the delay FREF is 2.5 ps away, dtc_res = 25PS
L_dtc = (0.25*pi^2)*((INL_pp*dtc_res)/Tdco)^2;
L_dtc_power = 10*log10(L_dtc); %eq 2.25, to check eq 2.45 if is almost the same

L_tdc = (0.25*pi^2)*((INL_pp*tdc_res)/Tdco)^2;
L_tdc_power = 10*log10(L_tdc); %eq 2.25

disp(["VCO pushing Power: ", num2str(spur_vco), " dBc"]);
disp(["DTC Power: ", num2str(L_dtc_power), " dBc"]);
disp(["TDC Power: ", num2str(L_tdc_power), " dBc"]);


print("finish")














