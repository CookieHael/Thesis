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
I_LNA_reg = I_LNA;
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



%% CS power
k = 1.38064852e-23;
T = 293;
V_thermal = 25.27*10^-3;
minimum_technology_cap = 1.995 * 1e-15;
F_per_M = 0.001025;
area_per_C = 1/ F_per_M;  %1/ (F/m^2)
C_mismatch_parameter = 3.4878e-09;
C_logic = 1*1e-15;



noise_multiplier = 1;

%% Process parameters
Vref = 2;
gmoverid = 20; %%I=150nA
V_eff = 1/gmoverid;
fmax_ADC = 10e6;   %% Determine from simulation!


% ADC_input_OS =  database_load('upsampled_EEG', epoch)'/10;
%data = load('100m (0).mat').x/10;


%% System settings
number_of_bits = 6;
input_BW = 256; %This should be the Nyquist frequency, and input should be oversampled wrt this
signal_peak_amplitude = .00002;
signal_bias = 0;
nb_MAC = 384;
cap_ratio = 15;
nb_channels = 150;
tertiary_matrix = false;
nb_activations = 2;
cap_gain_error = 1e-2;

non_idealities_on = 1;
noise_on = 1;
leakage_on = 1;

if ~noise_on
    warning("Noise currently off!")
end

if noise_multiplier>1 
    warning("Noise multiplier not one!")
end

if ~leakage_on
    warning("Leakage currently  off!")
end

noise_quantization = sqrt((Vref)^2/12*2^(-2*number_of_bits));


%% Clock & sine timing settings

clk_freq = 2.1*input_BW;
clk_period = 1/clk_freq;

%% Sample time jitter
sampling_jitter_on = noise_on;
sample_time_jitter = 16e-12;

%% Capacitor sizing
minimum_sampling_cap = 2 * 12*k*T*2^(2*number_of_bits)/Vref^2;  % Double to increase noise performance
[minimum_match_cap_DAC, minimum_match_cap_MAC] = getCapSizes(number_of_bits,C_mismatch_parameter, minimum_technology_cap, cap_ratio, cap_gain_error, F_per_M);
MAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_MAC]);
DAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_DAC]);

%% MAC non-idealities & noise
C_sense_size_multiplier = 2.5;
C1 = C_sense_size_multiplier*MAC_cap;
C2 = C1*cap_ratio; 
MM_incomplete_transfer_coeff = C1/(C1+C2);
if (non_idealities_on)
    MM_attenuation_coeff = C2/(C1+C2);  %Badly chosen name: this is voltage attenuation by charge sharing architecture
else
    MM_attenuation_coeff = 1;  %Actually more like 1-leakage
end
MM_sampling_noise = noise_multiplier*noise_on*((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2));  % Formula with square root is sigma, model uses variance

%% Leakage
leak_current = leakage_on *1*10^(-15);

%% Sample & hold (with gain)
switch_noise_on = noise_on;
Cs_unit_cap = 4/(MM_incomplete_transfer_coeff);
Cf_unit_cap = 2;
sampling_gain = Cs_unit_cap/Cf_unit_cap;
switching_noise = noise_multiplier*k*T/(Cf_unit_cap*minimum_sampling_cap);


%% Comparator offset/hysteresis & noise
comparator_offset_on = 1;
comparator_high_offset = 655e-6;
comparator_low_offset = comparator_high_offset;
comparator_noise_on = 1;
comparator_noise_rms = noise_multiplier* 1*k*T/C_logic * (2/3) ;


%% LNA
LNA_noise_on = 1;
LNA_gain = 1000;
LNA_bandwidth = 3*input_BW; %% in Hz
LNA_NEF = 1.08;

% Determine minimum current if bandwidth-constrained
gm_min_1 = LNA_gain*LNA_bandwidth*2*pi*nb_activations*MAC_cap;
Id_min_1 = gm_min_1*gmoverid;
LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if slewrate-constrained
SR_required = Vref*clk_freq; %Minimum current
Id_min_2 = SR_required*nb_activations*MAC_cap;
LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if noise-constrained
LNA_input_referred_noise_3 = 0.5*noise_quantization/(LNA_gain*sampling_gain*MM_incomplete_transfer_coeff);
Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

% Determine which of the three cases limits LNA performance
I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
if (I_LNA==Id_min_1 && Id_min_3>130e-09)
    warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
end

LNA_in_rms = min([LNA_input_referred_noise_1,LNA_input_referred_noise_2,LNA_input_referred_noise_3]); % Take maximum current of 3 required, noise will be minimum due to NEF formula
LNA_SR = I_LNA/(nb_activations*MAC_cap);
pd = makedist('Normal', 0, 1*LNA_in_rms);
rng(420699);
LNA_noise_vector = [pd.random(nb_MAC,1);0]';


%% Total gain check
total_gain = LNA_gain*sampling_gain*MM_incomplete_transfer_coeff;
if(signal_peak_amplitude*total_gain>Vref/2)
    warning("Careful! Input signal might be saturating ADC conversion range.\n");
end


%% ZOH settings

ZOH_delta = 0*Vref/(2^number_of_bits);
if (ZOH_delta>0)
    fprintf("Reminder: ZOH is on, signal quality may be degraded.\n");
end


%% Determine minimum required number of ADCs

OCR = nb_activations*(number_of_bits+3);
nb_ADC = ceil(OCR*clk_freq/fmax_ADC);


%% Power calculations without ZOH
C_0_DAC_CS = DAC_cap;
E_dac = (2^(number_of_bits) * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(signal_bias)^2 - (1/2)^number_of_bits*Vref*(signal_bias));
P_dac_CS = E_dac * nb_channels*clk_freq/(nb_MAC) * (nb_channels);

alpha_logic = 0.4;
P_logic_SAR_CS = alpha_logic * (2*number_of_bits*8*C_logic) * Vref^2 * nb_channels/(nb_MAC)*clk_freq*number_of_bits/(number_of_bits+1); %include leakage power
P_logic_phi_CS = (log2(nb_MAC) + 1) * nb_MAC * 8*C_logic * Vref^2 * clk_freq;  %alpha logic assumed 1 for shift register
P_logic_CS = P_logic_SAR_CS+P_logic_phi_CS;

P_switch_CS = 4*C_logic* (2*nb_activations)* Vref^2*clk_freq;


C_load_sar_CS = 4*C_logic;  %Assume latch used for input to SAR
E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp_CS = nb_channels/(nb_MAC) * clk_freq * number_of_bits/(number_of_bits+1) * E_comp;

P_tot_SAR_CS = (P_dac+P_comp+P_logic);

P_LNA_CS = Vref*I_LNA;

I_SH_CS = 2*sampling_gain*minimum_sampling_cap*Vref*clk_freq*nb_channels/(nb_MAC);
P_SH_CS = Vref * I_SH;

% Transmission power

P_transmission_full_signal_CS = clk_freq/nb_MAC * nb_channels * number_of_bits * transmission_power_per_bit;


P_tot_CS = P_tot_SAR_CS + P_LNA_CS + P_SH_CS + P_transmission_full_signal_CS+P_switch_CS;
% P_LNA_percentage = 100*P_LNA/P_tot;
% P_SH_percentage = 100*P_SH/P_tot;
% P_comp_percentage = 100*P_comp/P_tot;
% P_dac_percentage = 100*P_dac/P_tot;
% P_SAR_percentage = 100*P_tot_SAR/P_tot;
% P_switch_percentage = 100*P_switch/P_tot;
% P_logic_percentage = 100*P_logic/P_tot;
% P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;


figure
hold on
bar(1e9*[P_SH, P_SH_CS; P_logic, P_logic_CS; P_dac, P_dac_CS;  P_comp, P_comp_CS;]);
set(gca, 'XTick', 1:4, 'XTickLabels', { "S&H", "Logic", "DAC",  "Comparator"});
ylabel("Power consumption (nW)")
plot_paper
legend(["Regular", "CS"], "Location", "northwest");
axes('Position',[.6 .6 .3 .3])
box on
bar(1e6*[P_LNA, P_LNA_CS; P_tot_SAR, P_tot_SAR_CS; P_transmission_full_signal, P_transmission_full_signal_CS; ]);
set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "Others", "Transmission"});
ylabel("Power consumption (\muW)")
title("With TX and LNA")
plot_paper

