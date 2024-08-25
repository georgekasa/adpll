clc
clear
pkg load control

Freference = 26e6;
Fdco = 1.8e9;
Tref = 1.0/Freference;
Tdco = 1.0/Fdco;
FCW = Fdco/Freference;

%TDC & DCO spec.
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

% CLosed loop transfer function of the TDC multplied with Quant noise eq 4.97
H_closed_tdc = L_tdc_quant.*Hol_total_iir./(1+Hol_total_iir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%DCO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DCO quantization noise (without Sigma-Delta dithering)
%extremely important to compare this with the "natural" PN of the DCO itself!!!!
%eq 4.72
L_dco_quant = (1/12).*((kdco./f).^2).*(1/Freference).*sinc(kdco./Freference)^2;


%Transfer Function of DCO eq 4.98 NO MASH
Hclosed_dco = 1.0./(1.0 + Hol_total_iir)
Hclosed_dco_quant = Hclosed_dco.*L_dco_quant;


%Transfer Function of reference eq 4.96
Hclosed_ref = FCW.*Hol_total_iir./(1.0 + Hol_total_iir)


%Noises
L_tdc = 10*log10(L_tdc_quant);
PN_tdc_noIIR = 20*log10(abs(H_closed_tdc_NoIIR));
PN_tdc_IIR = 20*log10(abs(H_closed_tdc));
L_dco = 10*log10(L_dco_quant);
PN_dco_IIR = 20*log10(abs(Hclosed_dco_quant));
PN_ref_IIR = 20*log10(abs(Hclosed_ref));
% Total ADPLL noise PSD
adpll_pn = PN_tdc_IIR + PN_dco_IIR + PN_ref_IIR;


figure;
semilogx(f, L_tdc, 'g', 'LineWidth', 1.5);%tdc quantiz
hold on
%semilogx(f, PN_tdc_noIIR, 'b', 'LineWidth', 1.5);
semilogx(f, PN_tdc_IIR, 'b', 'LineWidth', 1.5);
%semilogx(f, L_dco, 'y', 'LineWidth', 1.5);
semilogx(f, PN_dco_IIR, 'y', 'LineWidth', 1.5);
semilogx(f, PN_ref_IIR, 'r', 'LineWidth', 1.5);
semilogx(f, adpll_pn, 'black', 'LineWidth', 1.5);
xlabel('Frequency Offset (Hz)');
ylabel('Phase Noise (dBc/Hz)');
grid on;



%tdc_pn_cl = tdc_pn_lin.*abs(H_cl_tdc).^2;


print("finish")














