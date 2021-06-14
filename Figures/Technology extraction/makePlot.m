clear;

%% gm/ID
gmoverid=load("gmoverid.mat").gmoverid;
currents=load("currents.mat").currents;

figure
hold on
plot(1e6*currents, gmoverid)
label = plot(0.13, 20, 'o');
xlabel("Drain current (\muA)")
ylabel("g_m/I_d")
plot_paper

%% leakage

leakcurrent = 1e12*load("leakcurrent.mat").leakcurrent;
times = load("times.mat").times;

figure
hold on
plot(times, abs(leakcurrent))
ylabel("Leakage current (pA)")
xlabel("time")
set(gca, 'YScale', 'log')
plot_paper