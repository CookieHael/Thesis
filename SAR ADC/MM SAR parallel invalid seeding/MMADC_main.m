clear;
verbose = false;
tic
% for i=1:20
%%  ----------  Constants & circuit settings  --------------
% number_of_bits = 1+i;
number_of_bits = 6;
input_BW = 800; %This should be the Nyquist frequency, and input should be oversampled wrt this
input_sampling_rate = 8000;
sine_amplitude = .001;
sine_bias = 0;
Vref = 2;
k = 1.38064852e-23;
T = 293;

%% Inputs
%ADC_input = load('/Users/mac/Documents/MATLAB/Thesis/SAR ADC/MM SAR/Neural data/upsampled_EEG.mat').output;
ADC_input = load('ADC_input.mat').in;

%% Sensing matrix settings

nb_MAC = 384;
cap_ratio = 35;
nb_channels = 150;
nb_activations = 2;
% [sensing_matrix_large, sensing_matrix_corrected] = generateCorrectedSRBM(nb_channels, nb_MAC, nb_activations, cap_ratio);
sensing_matrix_large = load("sensing150.mat").sensing_matrix_large;  %%Capacitor ratio was 35
sensing_matrix_corrected = load("sensing150.mat").sensing_matrix_corrected;

% sensing_matrix_large = load("sensing192.mat").Phi;
% sensing_matrix_corrected = load("sensing192.mat").corrected_Phi;
if( nb_MAC ~= size(sensing_matrix_large, 2) || nb_channels ~= size(sensing_matrix_large,1))
    error("Wrong matrix dimensions");
end
%% Clock & sine timing settings

%nyquist_frequency = sine_frequency*number_of_bits*2;
clk_freq = 2.1*(nb_MAC+1)/nb_MAC*input_BW;
clk_period = 1/clk_freq;

%% Circuit non-idealities

noise_quantization = Vref^2/12*2^(-2*number_of_bits);
non_idealities_on = 1;
noise_on = 1;

% Sample time jitter
sampling_jitter_on = noise_on;
sample_time_jitter = 16e-12;


% Capacitors
minimum_technology_cap = 1.995 * 1e-15;
minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
minimum_cap = max(minimum_technology_cap,minimum_sampling_cap);

% MAC non-idealities & noise
C_sense_size = 1;
C1 = C_sense_size*minimum_cap;
C2 = C1*cap_ratio; 
MM_incomplete_transfer_coeff = C1/(C1+C2);
if (non_idealities_on)
    MM_leakage_per_clk_cycle = C2/(C1+C2);  %Badly chosen name: this is voltage attenuation by charge sharing architecture
else
    MM_leakage_per_clk_cycle = 1;  %Actually more like 1-leakage
end

MM_sampling_noise = noise_on*((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2));  % Formula with square root is sigma, model uses variance


% Sample & hold switching noise
switch_noise_on = noise_on;
Cs = minimum_cap*10;
Cf = minimum_cap;
sampling_gain = Cs/Cf;
switching_noise = k*T/Cf;


%%% Comparator offset/hysteresis & noise
comparator_offset_on = noise_on;
comparator_high_offset = 0.2*10e-4;
comparator_low_offset = -.1*10e-4;
comparator_noise_on = noise_on;
comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;


% LNA
LNA_noise_on = noise_on;
LNA_bandwidth = 2.5*input_BW; %% in Hz
% LNA_input_referred_noise = 3.5*10^-6;
% LNA_gain = 100;
% LNA_NEF = 2;
LNA_input_referred_noise = sqrt(pi/2*LNA_bandwidth*(194*10^-9)^2);
LNA_gain = 10^(36/20)*12*1.5849;
LNA_NEF = 1.08;
if(sine_amplitude>Vref/(2*sampling_gain*LNA_gain*MM_incomplete_transfer_coeff))
    warning("Careful! Input signal might be saturating ADC conversion range.\n");
end

% Opamp




%% ZOH

ZOH_delta = 0*Vref/(2^number_of_bits);


%% Simulate
tic
cs_out = zeros(nb_channels,1);
i=nb_channels;
j=1;
sim_counter = 1;
sim_number=0;
nb_samples = 2;

while(i>0)
    fprintf("Simulation count: " + sim_counter + ".\n");
    sim_counter = sim_counter + 1;
    uncommentBlocks(20);
    if (i<20)
        commentBlocks(i+1,20);
    end
    
    sensing_matrix = sensing_matrix_large(j:j+min(19,i-1),:);
    
    sim_out = sim('MMADC', nb_samples*clk_period*(nb_MAC+1));
    cs_out(j) = sim_out.yout.signals(1).values(2);
    for l = 1:min(i-1,19)
        cs_out(j+l) = sim_out.yout.signals(l+3).values(2);
    end
    j = j + min(i,20);
    i = i - min(20,i);
    nb_ZOH = sim_out.yout.signals(2).values(length(sim_out.yout.signals(3).values));
end
toc

%% Compressive sensing reconstruction      
if (non_idealities_on)
    A = sensing_matrix_corrected*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
else
    A = sensing_matrix_large*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
end
recovery = BSBL_BO(A, cs_out, 1:15:nb_MAC, 0, 'prune_gamma',-1, 'max_iters',10);
   
recovered_signal = (idct(recovery.x))';
recovered_signal = recovered_signal(2:nb_MAC);
ADC_input_scaled = 1000*ADC_input(1:383);
recovered_signal_scaled = recovered_signal*2/(max(recovered_signal)-min(recovered_signal));
r = lowpass(recovered_signal, 0.7);
recovered_signal_scaled_lowpass = r*2/(max(r)-min(r));
rms_out = calculateRMS(ADC_input_scaled, recovered_signal);
rms_out_scaled = calculateRMS(ADC_input_scaled, recovered_signal_scaled);
rms_out_scaled_lowpass = calculateRMS(ADC_input_scaled, recovered_signal_scaled_lowpass);


%% Power calculations
C_0_DAC = minimum_cap;
E_dac = (2^(number_of_bits) * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));
P_dac = E_dac * nb_channels*clk_freq/(nb_MAC+1) * (nb_channels-nb_ZOH);

alpha_logic = 0.4;
C_logic = .77*1e-15;
P_logic_SAR = alpha_logic * ( (2*number_of_bits+5)*8*C_logic) * Vref^2 * nb_channels/(nb_MAC+1)*clk_freq*number_of_bits/(number_of_bits+3); %include leakage power  2n+1=> 2n+5 (2 extra FF for additional cycles, 2 extra for delaying outputs), + 3C_logic for ZOH gates
P_logic_ZOH = alpha_logic * C_logic * 17 * Vref^2 * nb_channels/(nb_MAC+1)*clk_freq * 2/(number_of_bits+3); % 3 * AND + 1 * OR + 1 * NOT assuming AND&OR are 4*C_logic
P_logic = P_logic_SAR + P_logic_ZOH;

V_eff = 0.1;
V_fs = Vref;
C_load_sar = 4*C_logic;  %Assume latch used for input to SAR
E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp = nb_channels/(nb_MAC+1) * clk_freq * ( (nb_channels-nb_ZOH)/nb_MAC * E_comp + 3/(number_of_bits+3) * (1-(nb_channels-nb_ZOH)/nb_MAC) * E_comp);

P_tot_SAR = (P_dac+P_comp+P_logic);

V_thermal = 25.27*10^-3;
I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA = Vref*I_LNA;


% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = clk_freq * (nb_channels/nb_MAC) * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR + P_LNA + P_transmission_full_signal;
P_LNA_percentage = round(100*P_LNA/P_tot,2);
P_comp_percentage = round(100*P_comp/P_tot,2);
P_dac_percentage = round(100*P_dac/P_tot,2);
P_SAR_percentage = round(100*P_tot_SAR/P_tot,2);
P_transmission_full_percentage = round(100*P_transmission_full_signal/P_tot,2);


%% Calculate total capacitance/area required - trade-off multi-ADC VS OC-ADC
%Total C & overclocking ratio
total_C_unit = zeros(5,1);
OCR = zeros(5,1);

for nb_ADC = 1:5

    C_MAC = C_sense_size*nb_activations + C_sense_size*cap_ratio*nb_channels;

    C_SH = nb_ADC * (1 + sampling_gain);

    C_SAR = 2^number_of_bits * nb_ADC;

    total_C_unit(nb_ADC) = C_MAC + C_SH + C_SAR;
    %area_C = total_C*area_per_C;
    OCR(nb_ADC) = nb_activations*(number_of_bits+3)/nb_ADC;
end

total_C = total_C_unit*minimum_cap;


%% Figure plotting
% Reconstruction
figure
subplot(3,1,1)
hold on
plot(recovered_signal)
plot(1000*ADC_input)
hold off
legend("Output values", "Input values");
title("Non-scaled, RMS = " + rms_out)
subplot(3,1,2)
hold on
plot(recovered_signal_scaled)
plot(1000*ADC_input)
hold off
legend("Output values", "Input values");
title("Scaled, RMS = " + rms_out_scaled)
subplot(3,1,3)
hold on
plot(recovered_signal_scaled_lowpass)
plot(1000*ADC_input)
hold off
legend("Output values", "Input values");
title("Scaled+lowpassed, RMS = " + rms_out_scaled_lowpass)

% Power & area
figure
hold on
subplot(311)
bar([P_LNA_percentage; P_SAR_percentage; P_transmission_full_percentage]);
title("Power per block in percentage of total. Total power = " + P_tot);
set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "SAR", "Transmission"});
subplot(312)
plot(total_C_unit)
title("Area per configuration");
ylabel("Total amount of unit caps");
xlabel("Number of ADCs");
set(gca, 'XTick', 1:5);
subplot(313)
plot(OCR)
set(gca, 'XTick', 1:5);
xlabel("Number of ADCs");
ylabel("Overclocking ratio");
title("Overclocking ratio for given number of ADCs")


%% Helper functions
function sndr_out = sinad_ADC(input_vector, nb_cycles_per_sample, clk_freq)
    y = zeros( round( size(input_vector,1) / (2*nb_cycles_per_sample) ), 1);
    i=1;
    for j = 1:size(y,1)
        y(j) = input_vector(i);
        i = i + 2*nb_cycles_per_sample;
    end
   % sinad(y, clk_freq/(nb_cycles_per_sample));
    sndr_out = sinad(y, clk_freq/(nb_cycles_per_sample));
end

function out = reshape_output(input, nb_cycles_per_sample)
    nb_skipped = 10*nb_cycles_per_sample; %% Skip the first few samples. These may be way off due to startup behaviour of LNA.
    out = zeros(round((size(input,1))/(2*nb_cycles_per_sample))-nb_skipped,1);
    i=nb_skipped;
    for j = 1:size(out,1)
        out(j) = input(i);
        i = i + 2*nb_cycles_per_sample;
    end
end

