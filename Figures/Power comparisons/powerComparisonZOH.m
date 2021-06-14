clear;

verbose = false;
% for i=1:20
%%  ----------  Circuit parameters  --------------
% number_of_bits = 1+i;
number_of_bits = 6;
input_BW=256;
V_thermal = 25.27*10^-3;
k = 1.38064852e-23;
F_per_M = 0.001025;
T = 293;
sine_amplitude = .001;
sine_bias = 0;
Vref = 2;
noise_quantization = sqrt(Vref^2/12*2^(-2*number_of_bits));
C_mismatch_parameter = 3.4878e-09;

ADC_input = database_load('EEGdata_ch1.mat', 1)'/10;

% Clock & sine timing settings
clk_freq = 2.1*(number_of_bits+1)*input_BW;
clk_period = 1/clk_freq;

%% Circuit non-idealities

non_idealities_on = 1;
% Sample time jitter
sampling_jitter_on = non_idealities_on;
sample_time_jitter = 16e-12;

% Cap sizing
minimum_technology_cap = 1.995e-15;
minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
minimum_matching_cap = getCapSizes1(number_of_bits,C_mismatch_parameter, minimum_technology_cap, F_per_M);
minimum_cap = max([minimum_technology_cap,minimum_sampling_cap, minimum_matching_cap]);

% Switching noise
switch_noise_on = non_idealities_on;
Cs = minimum_cap*10;
Cf = minimum_cap;
sampling_gain=Cs/Cf;
switching_noise = k*T/Cf;


% Comparator offset/hysteresis & noise
comparator_offset_on = non_idealities_on;
comparator_high_offset = 0.2*10e-4;
comparator_low_offset = -.1*10e-4;
comparator_noise_on = non_idealities_on;
comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;

%% LNA
LNA_noise_on = non_idealities_on;
LNA_gain = 1000;
LNA_bandwidth = 6*input_BW; %% in Hz
LNA_NEF = 1.08;
gmoverid = 20; %%I=150nA

% Determine minimum current if bandwidth-constrained
gm_min_1 = LNA_bandwidth*2*pi*minimum_cap;
Id_min_1 = gm_min_1*gmoverid;
LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if slewrate-constrained
SR_required = Vref*clk_freq; %Minimum current
Id_min_2 = SR_required*minimum_cap;
gm_min_2 = Id_min_2*gmoverid;
LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if noise-constrained
LNA_input_referred_noise_3 = noise_quantization/(LNA_gain*sampling_gain);
Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

% Determine which of the three cases limits LNA performance
I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
if (I_LNA==Id_min_3 && Id_min_3>130e-09)
    warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
end

LNA_input_referred_noise = min([LNA_input_referred_noise_1,LNA_input_referred_noise_2,LNA_input_referred_noise_3]); % Take maximum current of 3 required, noise will be minimum due to NEF formula
LNA_SR = I_LNA/(minimum_cap);

if(sine_amplitude>Vref/(2*sampling_gain*LNA_gain))
    warning("Careful! Input signal might be saturating ADC conversion range.\n");
end

%% Power calculations
C_0_DAC = minimum_cap;
P_dac = (2^number_of_bits*clk_freq*C_0_DAC)/(number_of_bits+1) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));

V_eff = 0.1;
V_fs = Vref;
C_load_sar = minimum_cap;
P_comp = 2*number_of_bits*log(2)*clk_freq*number_of_bits/(number_of_bits+1)*V_eff*V_fs*C_load_sar;

alpha_logic = 0.4;
C_logic = .77*1e-15;
P_logic = alpha_logic * (2*number_of_bits+1)*8*C_logic*Vref^2*(clk_freq-clk_freq/(number_of_bits+1));

P_tot_SAR = (P_dac + P_comp + P_logic);

V_thermal = 25.27*10^-3;
I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA = Vref*I_LNA;

I_SH = 2*sampling_gain*C_logic*Vref*(clk_freq/(number_of_bits+1));
P_SH = Vref * I_SH;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = clk_freq/(number_of_bits+1) * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR+P_LNA + P_transmission_full_signal + P_SH;
P_LNA_percentage = 100*P_LNA/ (P_tot);
P_SAR_percentage =  100*P_tot_SAR/ (P_tot);
P_comp_percentage = 100*P_comp/P_tot;
P_dac_percentage = 100*P_dac/P_tot;
P_SH_percentage = 100*P_SH/P_tot;
P_logic_percentage = 100*P_logic/P_tot;
P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;




%% ZOH calculation
%% Power calculations
nb_ZOH = 182;
nb_samples = 384;
C_0_DAC_ZOH = minimum_cap;
P_dac_ZOH = (2^number_of_bits*clk_freq*C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));
P_dac_ZOH = P_dac_ZOH * (nb_samples - nb_ZOH)/nb_samples; %% DAC inactive during ZOH cycles  Incorporate leakage though!!

C_load_sar_ZOH = minimum_cap;
E_comp_ZOH = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp_ZOH = (clk_freq-clk_freq/(number_of_bits/3))*(((nb_samples-nb_ZOH)/nb_samples * E_comp_ZOH) + (nb_ZOH/nb_samples)*2/(number_of_bits+3)*E_comp_ZOH); %include leakage power

P_logic_ZOH = alpha_logic * ( (2*number_of_bits+5)*8*C_logic + 3*C_logic) * Vref^2 * (clk_freq-clk_freq/(number_of_bits+3)); %include leakage power  2n+1=> 2n+5 (2 extra FF for additional cycles, 2 extra for delaying outputs), + 3C_logic for ZOH gates

P_tot_SAR_ZOH = (P_dac+P_comp+P_logic);


I_LNA_ZOH = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA_ZOH = Vref*I_LNA;

I_SH_ZOH = 2*sampling_gain*minimum_cap*Vref*clk_freq/(number_of_bits+3);
P_SH_ZOH = Vref * I_SH;

% Transmission power


P_transmission_full_signal_ZOH = (1-nb_ZOH/nb_samples)*clk_freq/(number_of_bits+3) * number_of_bits * transmission_power_per_bit;


P_tot_ZOH = P_tot_SAR_ZOH+P_LNA_ZOH + P_transmission_full_signal_ZOH + P_SH_ZOH;
P_LNA_percentage_ZOH = 100*P_LNA_ZOH/ (P_tot);
P_SAR_percentage_ZOH =  100*P_tot_SAR_ZOH/ (P_tot);
P_comp_percentage_ZOH = 100*P_comp_ZOH/P_tot;
P_dac_percentage_ZOH = 100*P_dac_ZOH/P_tot;
P_SH_percentage_ZOH = 100*P_SH_ZOH/P_tot;
P_logic_percentage_ZOH = 100*P_logic_ZOH/P_tot;
P_transmission_full_percentage_ZOH = 100*P_transmission_full_signal_ZOH/P_tot;


figure
hold on
bar(1e9*[P_SH, P_SH_ZOH; P_logic, P_logic_ZOH; P_dac, P_dac_ZOH;  P_comp, P_comp_ZOH;]);
set(gca, 'XTick', 1:4, 'XTickLabels', { "S&H", "Logic", "DAC",  "Comparator"});
ylabel("Power consumption (nW)")
plot_paper
legend(["Regular", "ZOH"], "Location", "northwest");
axes('Position',[.6 .6 .3 .3])
box on
bar(1e6*[P_LNA, P_LNA_ZOH; P_tot_SAR, P_tot_SAR_ZOH; P_transmission_full_signal, P_transmission_full_signal_ZOH; ]);
set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "Others", "Transmission"});
ylabel("Power consumption (\muW)")
title("With TX and LNA")
plot_paper

