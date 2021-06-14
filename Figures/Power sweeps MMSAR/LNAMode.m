
clear;
I1 = zeros(12,1);
I2 = zeros(12,1);
I3 = zeros(12,1);
I_LNA = zeros(12,1);

for i=2:12
    %%  Constants 
    k = 1.38064852e-23;
    T = 293;
    V_thermal = 25.27*10^-3;
    minimum_technology_cap = 1.995 * 1e-15;
    F_per_M = 0.001025;
    area_per_C = 1/ F_per_M;  %1/ (F/m^2)
    C_mismatch_parameter = 3.4878e-09;
    C_logic = .77*1e-15;


    %% Process parameters
    Vref = 1.1;
    gmoverid = 20; %%I=150nA
    V_eff = 1/gmoverid;
    fmax_ADC = 10e6;   %% Determine from simulation!

    %% System settings
    number_of_bits = i;
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


    noise_quantization = sqrt((Vref)^2/12*2^(-2*number_of_bits));


    %% Clock & sine timing settings

    clk_freq = 2.1*input_BW;
    clk_period = 1/clk_freq;

    %% Capacitor sizing
    minimum_sampling_cap = 2 * 12*k*T*2^(2*number_of_bits)/Vref^2;  % Double to increase noise performance
    [minimum_match_cap_DAC, minimum_match_cap_MAC] = getCapSizes(number_of_bits,C_mismatch_parameter, minimum_technology_cap, cap_ratio, cap_gain_error, F_per_M);
    MAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_MAC]);
    DAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_DAC]);

    %% MAC non-idealities & noise
    C_sense_size_multiplier = 5;
    C1 = C_sense_size_multiplier*MAC_cap;
    C2 = C1*cap_ratio; 
    MM_incomplete_transfer_coeff = C1/(C1+C2);
    
    %% Sample & hold (with gain)
    switch_noise_on = noise_on;
    Cs_unit_cap = 1.1/2*40/(MM_incomplete_transfer_coeff);
    Cf_unit_cap = 20;
    sampling_gain = Cs_unit_cap/Cf_unit_cap;

    %% LNA
    LNA_noise_on = 1;
    LNA_gain = 1000;
    LNA_bandwidth = 3*input_BW; %% in Hz
    LNA_NEF = 1.08;

    % Determine minimum current if bandwidth-constrained
    gm_min_1 = LNA_gain*LNA_bandwidth*2*pi*nb_activations*MAC_cap;
    I1(i) = gm_min_1*gmoverid;
    LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*I1(i)/(pi*V_thermal*4*k*T*LNA_bandwidth)));

    % Determine minimum current if slewrate-constrained
    SR_required = Vref*clk_freq; %Minimum current
    I2(i) = SR_required*nb_activations*MAC_cap;
    LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*I2(i)/(pi*V_thermal*4*k*T*LNA_bandwidth)));

    % Determine minimum current if noise-constrained
    LNA_input_referred_noise_3 = 0.5*noise_quantization/(LNA_gain*sampling_gain*MM_incomplete_transfer_coeff);
    I3(i) = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

    % Determine which of the three cases limits LNA performance
    I_LNA(i) = max([I1(i), I2(i), I3(i)]);

end

figure
hold on
plot(I1)
plot(I2)
plot(I3)
plot(I_LNA, "--*")
hold off
set(gca, "YScale", "log")
ylabel("LNA current requirement (A)")
xlabel("ADC resolution (bit)")
legend("BW-limited", "Slewrate-limited", "Noise-limited", "Location", "northwest")
plot_paper