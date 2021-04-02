clear;
verbose = false;

% for i=1:20
%%%  ----------  Circuit parameters  --------------
% number_of_bits = 1+i;
number_of_bits = 6;
sine_freq = 20000;
sine_amplitude = .001;
sine_bias = 0;
Vref = 2;
k = 1.38064852e-23;
T = 293;

ADC_input = load('/Users/mac/Documents/MATLAB/Thesis/ECG data/100m (0).mat').x;


%% MAC settings

nb_MAC = 10;

   
%% Clock & sine timing settings
sine_frequency_rad = 2 * pi * sine_freq;
%nyquist_frequency = sine_frequency*number_of_bits*2;
nb_cycles_per_sample = number_of_bits + 4 + nb_MAC;
clk_freq = 2.1*(nb_cycles_per_sample)*sine_frequency_rad;
clk_period = 1/clk_freq;


%% Circuit non-idealities

noise_quantization = Vref^2/12*2^(-2*number_of_bits);
non_idealities_on = 0;

% Sample time jitter
sampling_jitter_on = non_idealities_on;
sample_time_jitter = 16e-12;


% Capacitors
minimum_technology_cap = 1.995 * 1e-15;
minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
minimum_cap = max(minimum_technology_cap,minimum_sampling_cap);

% MAC
C1 = minimum_cap;
C2 = (nb_MAC)*minimum_cap;    %% V_accumulated = sum of C1/(C1+C2)*V(i) = C1/C2 * V(i) * nb_MAC => C2=nb_MAC*C1 => V_accumulated = V(i) assuming V(i) didn't change too  much
MM_sampling_noise = ((k*T * C1 / (C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2))^1/2;
MM_incomplete_transfer_coeff = C1/(C1+C2);
MM_leakage_per_clk_cycle = 1;

% Sample&hold switching noise
switch_noise_on = non_idealities_on;
Cs = minimum_cap * 10;
Cf = minimum_cap;
sampling_gain = Cs/Cf;
switching_noise = k*T/Cf;


%%% Comparator offset/hysteresis & noise
comparator_offset_on = non_idealities_on;
comparator_high_offset = 0.2*10e-4;
comparator_low_offset = -.1*10e-4;
comparator_noise_on = non_idealities_on;
comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;


% LNA
LNA_noise_on = non_idealities_on;
LNA_bandwidth = 40000; %% in Hz
LNA_input_referred_noise = 3.5*10^-6;
LNA_gain = 100;
LNA_NEF = 2;

if(sine_amplitude>Vref/(2*sampling_gain*LNA_gain))
    fprintf("Careful! Input signal might be saturating ADC conversion range.\n");
end


%% Compression
% ZOH compression

ZOH_delta = Vref/(2^number_of_bits);

%%% ----------------  Testing  ----------------

%% Simulate

nb_samples = 3000;
sim_out = sim('MAC_ADC', nb_samples*clk_period*(nb_cycles_per_sample));
nb_ZOH = sim_out.yout.signals(2).values(length(sim_out.yout.signals(3).values))


%% SNDR calculation
%sndr = sinad_ADC(sim_out.yout.signals(1).values, nb_cycles_per_sample, clk_freq);
rms_out = calculate_rms(1000*reshape_output(squeeze(sim_out.yout.signals(3).values), nb_cycles_per_sample), reshape_output(sim_out.yout.signals(1).values, nb_cycles_per_sample))
%sndr_minus_ideal = sndr - 6.02 * number_of_bits - 1.76 + 20 * log(Vref/((2 * sine_amplitude) * LNA_gain * sampling_gain));

%% Power calculations
C_0_DAC = minimum_cap;
P_dac = (2^(number_of_bits) * clk_freq * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));
P_dac = P_dac * (nb_samples - nb_ZOH)/nb_samples; %% DAC inactive during ZOH cycles  Incorporate leakage though!!

V_eff = 0.1;
V_fs = Vref;
C_load_sar = minimum_cap;
P_comp = 2*number_of_bits*log(2)*(clk_freq/(number_of_bits+3))*V_eff*V_fs*C_load_sar;
P_comp = ((nb_samples-nb_ZOH)/nb_samples * P_comp) + (1-((nb_samples-nb_ZOH)/nb_samples))*3/(number_of_bits+2)*P_dac; %include leakage power

alpha_logic = 0.4;
C_logic = .77*1e-15;
P_logic = alpha_logic * ( (2*number_of_bits+5)*8*C_logic + 3*C_logic) * Vref^2 * (clk_freq-clk_freq/(number_of_bits+3)); %include leakage power  2n+1=> 2n+5 (2 extra FF for additional cycles, 2 extra for delaying outputs), + 3C_logic for ZOH gates

P_tot_SAR = (P_dac+P_comp+P_logic)*1e9;

V_thermal = 25.27*10^-3;
I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA = Vref*I_LNA*10^9;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = 10^9*clk_freq/(number_of_bits+3) * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR + P_LNA + P_transmission_full_signal;
P_LNA_percentage = round(100*P_LNA/ (P_tot),2);
P_comp_percentage = round(100*P_comp/P_tot,2);
P_dac_percentage = round(100*P_dac/P_tot,2);
P_transmission_full_percentage = round(100*P_transmission_full_signal/P_tot,2);


%% Print results numbers
% fprintf("ADC resolution = " + number_of_bits + " bits.\n");
% fprintf("Calculated SNDR is " + sndr + " dB. This is " + sndr_minus_ideal + "dB from ideal. \n");
% if (verbose)
%     fprintf("Quantization noise = " + 10^6 * noise_quantization + "uVrms.\n");
%     fprintf("Minimal capacitor value = " + round(minimum_sampling_cap,4) + "fF.\n");
%     fprintf("Chosen sampling capacitor =  "+ Cs*10^9 +"nF. If larger than minimum, this lowers ADC power efficiency.\n")
% end
% fprintf("Total power draw =  " + P_tot + " nW.\n");
% if (verbose)
%     fprintf("LNA power draw = " + P_LNA + " nW, or " + P_LNA_percentage + "%% of the total power draw.\n");
%     fprintf("Comparator power draw = " + P_comp + " nW, or " + P_comp_percentage + "%% of the total power draw.\n");
%     fprintf("DAC power draw = " + P_dac + " nW, or " + P_dac_percentage + "%% of the total power draw.\n");
% end
% fprintf("Transmission power draw = " + P_transmission_full_signal + " nW, or " + P_transmission_full_percentage + "%% of the total power draw.\n");
% fprintf("Truncated signal's SNDR = " + sndr_compression1 + " dB, with compression ratio of " + compression1_ratio + "%%.\nTotal compression 1 power draw = " + P_tot_compression1 + "nW, or " + round(100*P_tot_compression1/P_tot, 1) + "%% of uncompressed total power.\n");
% end

%% Sample time jitter sweep

% sampling_jitter_on = 1;
% sample_time_jitters = [0 linspace(1e-5, 1e-8, 50)];
% 
% SNDR = zeros(size(sample_time_jitters));
% 
% for j = 1:size(sample_time_jitters,2)
%     sample_time_jitter = sample_time_jitters(j);
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(sample_time_jitters, SNDR)
% exportgraphics(f, 'sampling_time_jitter_sweep.pdf');
% sampling_jitter_on = 0;


%% Switching noise sweep

% switch_noise_on = 1;
% k = 1.38064852e-23;
% T = 293;
% nb_points = 50;
% Cs_array = linspace(1e-15, 1e-14, nb_points);
% Cf_array = Cs_array/sampling_gain;
% nb_points = size(Cs_array, 2);
% SNDR = zeros(1,nb_points);
% 
% for j = 1:nb_points
%     Cs = Cs_array(j);
%     Cf = Cf_array(j);
%     sampling_gain=Cs/Cf;
%     sim_out = sim('MMADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(Cs_array, SNDR);
% exportgraphics(f, 'switching_noise_sweep.pdf');
% switch_noise_on = 0;


%% Comparator offset

% comparator_offset_on = 1;
% comparator_high_offset_array = linspace(1e-5, 5e-1, 50);
% comparator_low_offset_array = comparator_high_offset_array/2;
% nb_points = size(comparator_high_offset_array,2);
% SNDR = zeros(1,nb_points);
% 
% for j = 1:nb_points
%     comparator_high_offset = comparator_high_offset_array(j);
%     comparator_low_offset = comparator_low_offset_array(j);
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(comparator_high_offset_array, SNDR);
% exportgraphics(f, 'comparator_offset_sweep.pdf');
% comparator_offset_on = 0;

% %% LNA noise
% 
% nb_points = 50;
% LNA_noise_on = 1;
% LNA_bandwidth = 10000;
% LNA_array = linspace(0, 5*1e-3, nb_points);
% SNDR = zeros(1,nb_points);
% 
% for j = 1:nb_points
%     LNA_input_reffered_noise = LNA_array(j);
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(comparator_high_offset_array, SNDR);
% exportgraphics(f, 'LNA_noise_sweep.pdf');
% LNA_noise_on = 0;
%     

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

function rms_out = calculate_rms(output_values, input_values)
hold on
plot(output_values)
plot(input_values)
hold off

if (length(input_values) ~= length(output_values))
    rms_out = -inf;
    return
end
input_values = input_values(2:length(input_values));
output_values = output_values(1:length(output_values)-1);
rms_out = sqrt(  sum( (input_values-output_values).^2 )  / length(output_values) );
end
