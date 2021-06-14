function [P_tot,total_C_unit] = getSystemCharacteristics(number_of_bits,nb_channels, nb_MAC, nb_activations, cap_ratio, leakage_size_boost, sampling_gain)
%GETSYSTEMCHARACTERISTICS Get system power and area consumption for certain
%set point

%%  Constants 
k = 1.38064852e-23;
T = 293;
V_thermal = 25.27*10^-3;
minimum_technology_cap = 1.995 * 1e-15;
F_per_M = 0.001025;
area_per_C = 1/ F_per_M;  %1/ (F/m^2)
C_mismatch_parameter = 3.4878e-09;
C_logic = 1*1e-15;
noise_multiplier  = 1;

%% Process parameters
Vref = 2;
gmoverid = 20; %%I=150nA
V_eff = 1/gmoverid;
fmax_ADC = 10e6;   %% Determine from simulation!



%% System settings
input_BW = 256; %This should be the Nyquist frequency, and input should be oversampled wrt this
signal_peak_amplitude = .00002;
signal_bias = 0;
tertiary_matrix = false;
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
C_sense_size_multiplier = leakage_size_boost;
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
Cf_unit_cap = 2;
Cs_unit_cap = Cf_unit_cap*sampling_gain;
switching_noise = noise_multiplier*k*T/(Cf_unit_cap*minimum_sampling_cap);


%% Comparator offset/hysteresis & noise
comparator_offset_on = noise_on;
comparator_high_offset = 8e-3;
comparator_low_offset = comparator_high_offset;
comparator_noise_on = noise_on;
comparator_noise_rms = noise_multiplier* 1*k*T* 2/3 / C_logic;


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

%% Determine minimum required number of ADCs

OCR = nb_activations*(number_of_bits+3);
nb_ADC = ceil(OCR*clk_freq/fmax_ADC);


%% Power calculations without ZOH
C_0_DAC = DAC_cap;
E_dac = (2^(number_of_bits) * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(signal_bias)^2 - (1/2)^number_of_bits*Vref*(signal_bias));
P_dac = E_dac * nb_channels*clk_freq/(nb_MAC) * (nb_channels);

alpha_logic = 0.4;
P_logic_SAR = alpha_logic * (2*number_of_bits*8*C_logic) * Vref^2 * nb_channels/(nb_MAC)*clk_freq*number_of_bits/(number_of_bits+1); %include leakage power
P_logic_phi = (log2(nb_MAC) + 1) * nb_MAC * 8*C_logic * Vref^2 * clk_freq;  %alpha logic assumed 1 for shift register
P_logic = P_logic_SAR+P_logic_phi;

P_switch = 4*C_logic* (2*nb_activations)* Vref^2*clk_freq;

V_fs = Vref;
C_load_sar = 4*C_logic;  %Assume latch used for input to SAR
E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp = nb_channels/(nb_MAC) * clk_freq * number_of_bits/(number_of_bits+1) * E_comp;

P_tot_SAR = (P_dac+P_comp+P_logic);

P_LNA = Vref*I_LNA;

I_SH = 2*sampling_gain*minimum_sampling_cap*Vref*clk_freq*nb_channels/(nb_MAC);
P_SH = Vref * I_SH;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = clk_freq/nb_MAC * nb_channels * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR + P_LNA + P_SH + P_transmission_full_signal+P_switch;

%% Calculate total capacitance/area required - trade-off multi-ADC VS OC-ADC
%Total C & overclocking ratio

min_nb_ADC = ceil(nb_activations*(number_of_bits+3)/fmax_ADC);
C_MAC = C_sense_size_multiplier*nb_activations + C_sense_size_multiplier*cap_ratio*nb_channels;
C_SAR = 2^number_of_bits * (min_nb_ADC);
total_C = C_MAC*MAC_cap + C_SAR*DAC_cap + (Cs_unit_cap + Cf_unit_cap)*MAC_cap*min_nb_ADC;
OCR = nb_activations*(number_of_bits+3)/(min_nb_ADC);
        
total_C_unit = round(total_C/minimum_technology_cap);
area_C = total_C*area_per_C;

end

