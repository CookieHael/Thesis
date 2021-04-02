clear;

verbose = false;
SNDRi = zeros(12,1);
comppoweri = zeros(12,1);
sarpoweri = zeros(12,1);
LNApoweri = zeros(12,1);
totpoweri = zeros(12,1);
dacpoweri = zeros(12,1);
logicpoweri = zeros(12,1);
transmissionpoweri = zeros(12,1);

for i=1:12
    %%%  ----------  Circuit parameters  --------------
    number_of_bits = 1+i;
    sine_freq = 200;
    sine_amplitude = .001;
    sine_bias = 0;
    Vref = 2;
    noise_quantization = Vref^2/12*2^(-2*number_of_bits);


    % Clock & sine timing settings
    sine_frequency_rad = 2*pi*sine_freq;
    %nyquist_frequency = sine_frequency*number_of_bits*2;
    clk_freq = 2.1*(number_of_bits+1)*sine_frequency_rad;
    clk_period = 1/clk_freq;


    % Sample time jitter
    sampling_jitter_on = 1;
    sample_time_jitter = 16e-12;

    % Switching noise
    switch_noise_on = 1;
    k = 1.38064852e-23;
    T = 293;
    minimum_technology_cap = 1e-15;
    minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
    minimum_cap = max(minimum_technology_cap,minimum_sampling_cap);
    Cs = minimum_cap*10;
    Cf = minimum_cap;
    sampling_gain=Cs/Cf;
    switching_noise = k*T/Cf;


    %%% Comparator offset/hysteresis & noise
    comparator_offset_on = 0;
    comparator_high_offset = 0.2*10e-4;
    comparator_low_offset = -.1*10e-4;
    comparator_noise_on = 1;
    comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;


    % LNA
    LNA_noise_on=1;
    LNA_bandwidth = 400;
    LNA_input_referred_noise = 3.5*10^-6;
    LNA_gain = 100;
    LNA_NEF = 2;

    if(sine_amplitude>Vref/(2*sampling_gain*LNA_gain))
        fprintf("Careful! Input signal might be saturating ADC conversion range.\n");
    end

    % ZOH compression
    ZOH_delta = .05;
    
    % Truncation compression
    number_of_bits_transmission = 2;

    %%% ----------------  Testing  ----------------

    %% Simulate
    sim_out = sim('ADC', 20000*clk_period);


    %% SNDR calculation
    sndr = sinad_ADC(sim_out.yout.signals(1).values, number_of_bits, clk_freq);
    sndr_minus_ideal = sndr - 6.02*number_of_bits - 1.76 + 20*log(Vref/((2*sine_amplitude)*LNA_gain*sampling_gain))
    SNDRi(i,1) = sndr;

    %% Power calculations
    C_0_DAC = minimum_cap;
    P_dac = 2/3 * (2^(number_of_bits))*clk_freq*C_0_DAC*Vref^2/(number_of_bits+1);
    %P_dac = ((2^(number_of_bits))*clk_freq*C_0_DAC)/(number_of_bits+1) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias+sine_amplitude)^2 - (1/2)^number_of_bits*Vref*(sine_bias+sine_amplitude));
    
    V_eff = 0.1;
    V_fs = Vref;
    C_load_sar = minimum_cap;
    P_comp = 2*number_of_bits*log(2)*(clk_freq-clk_freq/(2*number_of_bits+2))*V_eff*V_fs*C_load_sar;

    comppoweri(i,1) = P_comp;

    alpha_logic = 0.4;
    C_logic = .77*1e-15;
    P_logic = alpha_logic * (2*number_of_bits+1)*8*C_logic*Vref^2*(clk_freq-clk_freq/(number_of_bits+1));
    logicpoweri(i,1) = P_logic;

    P_tot_SAR = (P_dac+P_comp+P_logic);
    sarpoweri(i,1) = P_tot_SAR;

    V_thermal = 25.27*10^-3;
    I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
    P_LNA = Vref*I_LNA;
    LNApoweri(i,1) = P_LNA;

    % Transmission power

    transmission_power_per_bit = 10^-9;
    P_transmission_full_signal = clk_freq/(number_of_bits+2) * number_of_bits * transmission_power_per_bit;
    transmissionpoweri(i,1) = P_transmission_full_signal;

    P_tot = P_tot_SAR + P_LNA + P_transmission_full_signal;
    totpoweri(i,1) = P_tot;
    P_LNA_percentage = round(100*P_LNA/ (P_tot),2);
    P_comp_percentage = round(100*P_comp/P_tot,2);
    P_dac_percentage = round(100*P_dac/P_tot,2);
    P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;

    % Compressed signals

%     sndr_compression1 = sinad_ADC(sim_out.yout.signals(2).values, number_of_bits, clk_freq);
%     P_transmission_compression1 = 10^9*clk_freq/(number_of_bits+1) * number_of_bits_transmission * transmission_power_per_bit;
%     compression1_ratio = round(100*(number_of_bits-number_of_bits_transmission)/number_of_bits,2);   %Voorstelling? Relatief vs absoluut
%     P_tot_compression1 = P_LNA+P_tot_SAR + P_transmission_compression1;


    %% Print results numbers
    fprintf("ADC resolution = " + number_of_bits + " bits.\n");
    fprintf("Calculated SNDR is " + sndr + " dB. This is " + sndr_minus_ideal + "dB from ideal. \n");
    if (verbose)
        fprintf("Quantization noise = " + 10^6 * noise_quantization + "uVrms.\n");
        fprintf("Minimal capacitor value = " + round(minimum_sampling_cap,4) + "fF.\n");
        fprintf("Chosen sampling capacitor =  "+ Cs*10^9 +"nF. If larger than minimum, this lowers ADC power efficiency.\n")
    end
    fprintf("Total power draw =  " + P_tot + " nW.\n");
    if (verbose)
        fprintf("LNA power draw = " + P_LNA + " nW, or " + P_LNA_percentage + "%% of the total power draw.\n");
        fprintf("Comparator power draw = " + P_comp + " nW, or " + P_comp_percentage + "%% of the total power draw.\n");
        fprintf("DAC power draw = " + P_dac + " nW, or " + P_dac_percentage + "%% of the total power draw.\n");
    end
     fprintf("Transmission power draw = " + P_transmission_full_signal + " nW, or " + P_transmission_full_percentage + "%% of the total power draw.\n");
%     fprintf("Truncated signal's SNDR = " + sndr_compression1 + " dB, with compression ratio of " + compression1_ratio + "%%.\nTotal compression 1 power draw = " + P_tot_compression1 + "nW, or " + round(100*P_tot_compression1/P_tot, 1) + "%% of uncompressed total power.\n");

    fprintf("Run: " +i);
end

hold on
set(gca,'fontsize', 25);
set(gca, 'YScale', 'log')
scatter(SNDRi, totpoweri, 150,'filled','DisplayName', 'total power');
scatter(SNDRi, comppoweri, 150, 'filled','DisplayName', 'comp power');
scatter(SNDRi, sarpoweri, 150, 'filled','DisplayName', 'sar power');
scatter(SNDRi, LNApoweri, 150, 'filled','DisplayName', 'lna power');
scatter(SNDRi, dacpoweri, 150, 'filled','DisplayName', 'dac power');
scatter(SNDRi, logicpoweri, 150, 'filled','DisplayName', 'logic power');
scatter(SNDRi, transmissionpoweri, 150, 'filled','DisplayName', 'transmission power');
legend
plot(SNDRi, totpoweri);
plot(SNDRi, comppoweri);
plot(SNDRi, sarpoweri);
plot(SNDRi, LNApoweri);
plot(SNDRi, dacpoweri);
plot(SNDRi, logicpoweri);
plot(SNDRi, transmissionpoweri);
title("Power vs SNDR (non-ideal)");
xlabel("SNDR (dB)");
ylabel("Power (W)");
hold off

function sndr_out = sinad_ADC(input_vector, number_of_bits, clk_freq)
y = zeros(round((size(input_vector,1)-2)/(2*(number_of_bits+1))),1);
i=1;
for j = 1:size(y,1)
    y(j) = input_vector(i);
    i= i + 2 + 2*number_of_bits;
end
%sinad(y, clk_freq/(number_of_bits+1))
sndr_out = sinad(y, clk_freq/(number_of_bits+1));
end
