clear;

verbose = false;
SNDRi = zeros(12,1);
SNDRideal = zeros(12,1);

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
    LNA_noise_on=0;
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
    sndr_minus_ideal = sndr - 6.02*number_of_bits - 1.76 + 20*log(Vref/((2*sine_amplitude)*LNA_gain*sampling_gain));
    SNDRi(i,1) = sndr;
    SNDRideal(i,1) = 6.02*number_of_bits + 1.76 - 20*log(Vref/((2*sine_amplitude)*LNA_gain*sampling_gain));

    fprintf("Run: " +i);
end

C1=ones(12,1);
C1 = C1*[0, 0.4470, 0.7410];
C2 = ones(12,1);
C2 = C2*[0.8500, 0.3250, 0.0980];
hold on
set(gca,'fontsize', 18);
scatter(linspace(2,13, 12),SNDRi, 150, C1 , 'filled','DisplayName', 'Noisy ADC: sources as given before, LNA noiseless');
scatter(linspace(2,13, 12),SNDRideal,150, C2 , 'filled','DisplayName', 'Ideal SNDR');
legend
plot(linspace(2,13, 12),SNDRi, 'Color', [0, 0.4470, 0.7410]);
plot(linspace(2,13, 12),SNDRideal,'Color', [0.8500, 0.3250, 0.0980]);
title("Ideal SNDR vs noisy ADC with noiseless LNA");
ylabel("SNDR (dB)");
xlabel("Internal number of bits");
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