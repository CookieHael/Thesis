% Load saved figures
c=hgload('Coefficients75chAllOn.fig');
k=hgload('CoeffErrors75chAllOn.fig');
% Prepare subplots
figure
h(1)=subplot(2,1,1, "LineStyle", "None");
ylabel="Signal (V)";
xlabel = "Coefficient index";
legend=(["Ideal coefficients", "Output coefficients", "Corrected coefficients"]);
plot_paper
h(2)=subplot(2,1,2, "LineStyle", "None");
legend=("Error on corrected");
ylabel="Signal (V)";
xlabel = "Coefficient index";
plot_paper
% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1));
copyobj(allchild(get(k,'CurrentAxes')),h(2));
% Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')