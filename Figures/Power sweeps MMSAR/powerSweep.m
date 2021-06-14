clear;
verbose = false;
tic

%%  Constants 
k = 1.38064852e-23;
T = 293;
V_thermal = 25.27*10^-3;
minimum_technology_cap = 1.995 * 1e-15;
area_per_C = 0.001025;  %F/m^2
C_mismatch_parameter = 3.4878e-09;
C_logic = .77*1e-15;

%% Process parameters
Vref = 2;
gmoverid = 20; %%I=150nA
V_eff = 1/gmoverid;
fmax_ADC = 10e6;   %% Determine from simulation!


P_dac_array = NaN(7, 100, 100);
P_tot_SAR_array = NaN(7, 100, 100);
P_switch_array = NaN(7, 100, 100);
P_logic_phi_array = NaN(7, 100, 100);
P_logic_array = NaN(7, 100, 100);
P_comp_array = NaN(7, 100, 100);
P_LNA_array = NaN(7, 100, 100);
P_SH_array = NaN(7, 100, 100);
P_transmission_full_signal_array = NaN(7, 100, 100);
P_tot_array = NaN(7, 100, 100);
total_C_array = NaN(7, 100, 100);


for nb_MAC = round(linspace(50, 512, 100))
    for nb_channels = round(linspace(50, nb_MAC, 10*round(512/50)))
        parfor number_of_bits = 4:10

%% System settings
input_BW = 256; %This should be the Nyquist frequency, and input should be oversampled wrt this
signal_peak_amplitude = .00002;
signal_bias = 0;
cap_ratio = 15;
nb_activations = 2;
cap_gain_error = 1e-2;
noise_on=1;
noise_quantization = sqrt((Vref)^2/12*2^(-2*number_of_bits));

%% Clock & sine timing settings

clk_freq = 2.1*input_BW;
clk_period = 1/clk_freq;

%% Capacitor sizing
minimum_sampling_cap = 2 * 12*k*T*2^(2*number_of_bits)/Vref^2;  %half of value for noise reasons
[minimum_match_cap_DAC, minimum_match_cap_MAC] = getCapSizes(number_of_bits,C_mismatch_parameter, minimum_technology_cap, cap_ratio, cap_gain_error, area_per_C);
MAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_MAC]);
DAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_DAC]);

C_sense_size_multiplier=1;
C1 = C_sense_size_multiplier*MAC_cap;
C2 = C1*cap_ratio; 
MM_incomplete_transfer_coeff = C1/(C1+C2);

%% Sample & hold (with gain)
switch_noise_on = noise_on;
Cs_unit_cap = 20/(MM_incomplete_transfer_coeff)/2.5;
Cf_unit_cap = 20;
sampling_gain = Cs_unit_cap/Cf_unit_cap;
switching_noise = k*T/(Cf_unit_cap*minimum_sampling_cap);

%% LNA
LNA_noise_on = 1;
LNA_gain = 500;
LNA_bandwidth = 3*input_BW; %% in Hz
LNA_NEF = 1.08;

% Determine minimum current if bandwidth-constrained
gm_min_1 = LNA_bandwidth*2*pi*nb_activations*MAC_cap;
Id_min_1 = gm_min_1*gmoverid;
LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if slewrate-constrained
SR_required = Vref*clk_freq; %Minimum current
Id_min_2 = SR_required*nb_activations*MAC_cap;
gm_min_2 = Id_min_2*gmoverid;
LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if noise-constrained
LNA_input_referred_noise_3 = 0.5*noise_quantization/(LNA_gain*sampling_gain*MM_incomplete_transfer_coeff);
Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

% Determine which of the three cases limits LNA performance
I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
% if (I_LNA==Id_min_3 && Id_min_3>130e-09)
%     warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
% end

%% Power calculation
nb_ZOH=0;

C_0_DAC = DAC_cap;
E_dac = (2^(number_of_bits) * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(signal_bias)^2 - (1/2)^number_of_bits*Vref*(signal_bias));
P_dac = E_dac * nb_channels * nb_activations*(number_of_bits+3) * clk_freq/(nb_MAC) * (nb_channels-nb_ZOH);
P_dac_array(number_of_bits-3, nb_MAC, nb_channels) = P_dac;


alpha_logic = 0.4;
P_logic_SAR = alpha_logic * ( (2*number_of_bits+5)*8*C_logic) * Vref^2 * nb_channels/(nb_MAC)*clk_freq*number_of_bits/(number_of_bits+3) * nb_activations*(number_of_bits+3); %include leakage power  2n+1=> 2n+5 (2 extra FF for additional cycles, 2 extra for delaying outputs), + 3C_logic for ZOH gates
P_logic_ZOH = alpha_logic * C_logic * 17 * Vref^2 * nb_channels/(nb_MAC)*clk_freq * 2/(number_of_bits+3) * nb_activations*(number_of_bits+3); % 3 * AND + 1 * OR + 1 * NOT assuming AND&OR are 4*C_logic
P_logic_phi = (log2(nb_MAC) + 1) * nb_MAC * 8*C_logic * Vref^2 * clk_freq;  %alpha logic assumed 1 for shift register
P_logic = P_logic_SAR + P_logic_ZOH + P_logic_phi;
P_logic_phi_array(number_of_bits-3, nb_MAC, nb_channels) = P_logic_phi;
P_logic_array(number_of_bits-3, nb_MAC, nb_channels) = P_logic;

P_switch = 4*C_logic* (nb_MAC+nb_activations)* Vref^2*clk_freq;
P_switch_array(number_of_bits-3, nb_MAC, nb_channels) = P_switch;

V_fs = Vref;
C_load_sar = 4*C_logic;  %Assume latch used for input to SAR
E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp = nb_channels/(nb_MAC) * clk_freq * (number_of_bits+3) * ( (nb_channels-nb_ZOH)/nb_MAC * E_comp + 3/(number_of_bits+3) * (1-(nb_channels-nb_ZOH)/nb_MAC) * E_comp);
P_comp_array(number_of_bits-3, nb_MAC, nb_channels) = P_comp;


P_tot_SAR = (P_dac+P_comp+P_logic);
P_tot_SAR_array(number_of_bits-3, nb_MAC, nb_channels) = P_tot_SAR;

P_LNA = Vref*I_LNA;
P_LNA_array(number_of_bits-3, nb_MAC, nb_channels) = P_LNA;

I_SH = 2*sampling_gain*minimum_sampling_cap*Vref*clk_freq*nb_channels/(nb_MAC);
P_SH = Vref * I_SH;
P_SH_array(number_of_bits-3, nb_MAC, nb_channels) = P_SH;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = clk_freq/nb_MAC * (nb_channels-nb_ZOH) * number_of_bits * transmission_power_per_bit;
P_transmission_full_signal_array(number_of_bits-3, nb_MAC, nb_channels) = P_transmission_full_signal;


P_tot = P_tot_SAR + P_LNA + P_SH + P_transmission_full_signal+P_switch;
P_tot_array(number_of_bits-3, nb_MAC, nb_channels) = P_tot;

%% Calculate total capacitance/area required - trade-off multi-ADC VS OC-ADC
%Total C & overclocking ratio

min_nb_ADC = ceil(nb_activations*(number_of_bits+3)/fmax_ADC);
total_C_unit = zeros(4,1);
area_C = zeros(4,1);
    
C_MAC = C_sense_size_multiplier*nb_activations + C_sense_size_multiplier*cap_ratio*nb_channels; 
C_SH = min_nb_ADC * (1 + sampling_gain);
C_SAR = 2^number_of_bits * min_nb_ADC;

total_C = C_MAC*MAC_cap + C_SH*minimum_sampling_cap + C_SAR*DAC_cap + (Cs_unit_cap + Cf_unit_cap)*MAC_cap;
total_C_array(number_of_bits-3, nb_MAC, nb_channels) = total_C;
OCR = nb_activations*(number_of_bits+3)/(min_nb_ADC);

        end
    end
end
toc
save("Output", 'P_dac_array', 'P_tot_SAR_array', 'P_switch_array', 'P_switch_array', 'P_logic_array', 'P_comp_array', 'P_LNA_array', 'P_SH_array', 'P_transmission_full_signal_array', 'P_tot_array', 'total_C_array');

