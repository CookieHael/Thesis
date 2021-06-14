clear;

verbose = false;
% for i=1:20
%%  ----------  Circuit parameters  --------------
% number_of_bits = 1+i;
number_of_bits = 6;
sine_amplitude = .001;
sine_bias = 0;
Vref = 2;
noise_quantization = Vref/(sqrt(12)*2^(number_of_bits));
V_thermal = 25.27*10^-3;
C_mismatch_parameter = 3.4878e-09;
F_per_M = 0.001025;
   
% ADC_input = load('/Users/mac/Documents/MATLAB/Thesis/ECG data/100m (0).mat').x;
% ADC_input = ADC_input(2:400);
ADC_input = database_load('EEGdata_ch1.mat', 1)'/10;



%% Clock & sine timing settings
input_BW = 256; %Set to 256Hz
nb_cycles_per_sample = number_of_bits + 3;
clk_freq = 2.1*(nb_cycles_per_sample)*input_BW;
clk_period = 1/clk_freq;

%% Circuit non-idealities

non_idealities_on = 1;
% Sample time jitter
sampling_jitter_on = non_idealities_on;
sample_time_jitter = 16e-12;

% Switching noise
switch_noise_on = non_idealities_on;
k = 1.38064852e-23;
T = 293;
minimum_technology_cap = 1.995e-15;
minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
minimum_matching_cap = getCapSizes(number_of_bits,C_mismatch_parameter, minimum_technology_cap, F_per_M);
minimum_cap = max([minimum_technology_cap,minimum_sampling_cap, minimum_matching_cap]);
Cs = minimum_cap*10;
Cf = minimum_cap;
sampling_gain=Cs/Cf;
switching_noise = k*T/Cf;


%%% Comparator offset/hysteresis & noise
comparator_offset_on = non_idealities_on;
comparator_high_offset = 0.2*10e-4;
comparator_low_offset = -.1*10e-4;
comparator_noise_on = non_idealities_on;
comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;


%% LNA
LNA_noise_on = 1;
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
Id_min_2 = SR_required*2*minimum_cap;
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
LNA_input_referred_noise = min([LNA_input_referred_noise_1, LNA_input_referred_noise_2, LNA_input_referred_noise_3]);
LNA_SR = I_LNA/(2*minimum_cap);


if(sine_amplitude>Vref/(2*sampling_gain*LNA_gain))
    fprintf("Careful! Input signal might be saturating ADC conversion range.\n");
end

%% ZOH settings
% ZOH compression
ZOH_delta = 2*Vref/(2^number_of_bits);

%%% ----------------  Testing  ----------------

%% Simulate

nb_samples = length(ADC_input);
sim_out = sim('ZOHADC', (nb_samples+1)*clk_period*nb_cycles_per_sample);
nb_ZOH = sim_out.yout.signals(2).values(end-1) %2nd last value because sim needs to run slightly longer, but that sample may get ZOH'd


%% Result quality calculations
% sndr = sinad_ADC(sim_out.yout.signals(1).values, nb_cycles_per_sample, clk_freq);
% rms_out = calculate_rms(1000*reshape_output(squeeze(sim_out.yout.signals(3).values), nb_cycles_per_sample), reshape_output(sim_out.yout.signals(1).values, nb_cycles_per_sample))
% sndr_minus_ideal = sndr - 6.02 * number_of_bits - 1.76 + 20 * log(Vref/((2 * sine_amplitude) * LNA_gain * sampling_gain));

%% Power calculations
C_0_DAC = minimum_cap;
P_dac = (2^number_of_bits*clk_freq*C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));
P_dac = P_dac * (nb_samples - nb_ZOH)/nb_samples; %% DAC inactive during ZOH cycles  Incorporate leakage though!!

V_eff = 0.1;
V_fs = Vref;
C_load_sar = minimum_cap;
E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp = (clk_freq-clk_freq/(number_of_bits/3))*(((nb_samples-nb_ZOH)/nb_samples * E_comp) + (nb_ZOH/nb_samples)*2/(number_of_bits+3)*E_comp); %include leakage power

alpha_logic = 0.4;
C_logic = .77*1e-15;
P_logic = alpha_logic * ( (2*number_of_bits+5)*8*C_logic + 3*C_logic) * Vref^2 * (clk_freq-clk_freq/(number_of_bits+3)); %include leakage power  2n+1=> 2n+5 (2 extra FF for additional cycles, 2 extra for delaying outputs), + 3C_logic for ZOH gates

P_tot_SAR = (P_dac+P_comp+P_logic);

V_thermal = 25.27*10^-3;
I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA = Vref*I_LNA;

I_SH = 2*sampling_gain*minimum_cap*Vref*clk_freq/(number_of_bits+3);
P_SH = Vref * I_SH;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = (1-nb_ZOH/nb_samples)*clk_freq/(number_of_bits+3) * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR+P_LNA + P_transmission_full_signal + P_SH;
P_LNA_percentage = 100*P_LNA/ (P_tot);
P_SAR_percentage =  100*P_tot_SAR/ (P_tot);
P_comp_percentage = 100*P_comp/P_tot;
P_dac_percentage = 100*P_dac/P_tot;
P_SH_percentage = 100*P_SH/P_tot;
P_logic_percentage = 100*P_logic/P_tot;
P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;


figure
hold on
plot(10000*ADC_input(1:end))
plot(sim_out.yout.signals(1).values(3:end))
hold off
xlabel("Sample index")
ylabel("Output voltage (V)")
legend("Input signal", "Output signal")
plot_paper

% figure
% hold on
% bar([P_LNA_percentage; P_SH_percentage; P_comp_percentage; P_dac_percentage;P_logic_percentage]);
% set(gca, 'XTick', 1:5, 'XTickLabels', {"LNA", "S&H", "Comparator", "DAC", "Logic"});
% ylabel("Relative power draw percentages without tranmission")
% title("Total power: " + P_tot*1e6 + "\muW")
% axes('Position',[.7 .7 .2 .2])
% plot_paper
% box on
% bar([P_LNA_percentage; P_SAR_percentage; P_transmission_full_percentage; ]);
% set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "SAR", "Transmission"});
% title("With TX")
% plot_paper


% figure
% hold on
% bar([P_SH_percentage; P_dac_percentage;P_logic_percentage; P_comp_percentage]);
% set(gca, 'XTick', 1:4, 'XTickLabels', { "S&H", "DAC", "Logic", "Comparator"});
% ylabel("Relative power draw percentages without Tx & LNA")
% plot_paper
% 
% title("Total power: " + P_tot*1e6 + "\muW")
% axes('Position',[.7 .7 .2 .2])
% box on
% bar([P_LNA_percentage; P_SAR_percentage; P_transmission_full_percentage; ]);
% set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "SAR", "Transmission"});
% title("With TX and LNA")
% plot_paper


%% Comparison plot
nb_ZOH=0;
clk_freq = 2.1*input_BW*(number_of_bits+1);
C_0_DAC = minimum_cap;
P_dac_c = (2^number_of_bits*clk_freq*C_0_DAC)/(number_of_bits+1) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));

V_eff = 0.1;
V_fs = Vref;
C_load_sar = minimum_cap;
P_comp_c = 2*number_of_bits*log(2)*clk_freq*number_of_bits/(number_of_bits+1)*V_eff*V_fs*C_load_sar;

alpha_logic = 0.4;
C_logic = .77*1e-15;
P_logic_c = alpha_logic * (2*number_of_bits+1)*8*C_logic*Vref^2*(clk_freq-clk_freq/(number_of_bits+1));

P_tot_SAR_c = (P_dac + P_comp + P_logic);

V_thermal = 25.27*10^-3;
I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA_c = Vref*I_LNA;

I_SH = 2*sampling_gain*C_logic*Vref*(clk_freq/(number_of_bits+1));
P_SH_c = Vref * I_SH;

P_transmission_full_signal_c = (1-nb_ZOH/nb_samples)*clk_freq/(number_of_bits+3) * number_of_bits * transmission_power_per_bit;


P_tot_c = P_tot_SAR_c + P_LNA + P_transmission_full_signal_c + P_SH;
P_SAR_percentage_c =  P_SAR_percentage*P_tot_SAR/P_tot_SAR_c;
P_comp_percentage_c = P_comp_percentage*P_comp/P_comp_c;
P_dac_percentage_c = P_dac_percentage*P_dac/P_dac_c;
P_logic_percentage_c = P_logic_percentage*P_logic/P_logic_c;
P_tot_c_percentage = 100*P_tot/P_tot_c;
P_transmission_full_percentage_c = P_transmission_full_percentage*P_transmission_full_signal/P_transmission_full_signal_c;


figure
hold on
bar(1e6*[ P_SH, P_SH; P_dac_c, P_dac; P_logic_c, P_logic; P_comp_c,P_comp ; ]);
set(gca, 'XTick', 1:4, 'XTickLabels', {"S&H", "DAC", "Logic", "Comparator"});
ylabel("Power consumption (\muW)")
legend("No ZOH", "ZOH", 'Location','northwest')
axes('Position',[.7 .7 .2 .2])
box on
bar(1e6*[ P_tot_c, P_tot;P_transmission_full_signal_c, P_transmission_full_signal; P_LNA, P_LNA]);
set(gca, 'XTick', 1:3, 'XTickLabels', {"Total", "Transmission", "LNA", "Others"});
title("With TX")

