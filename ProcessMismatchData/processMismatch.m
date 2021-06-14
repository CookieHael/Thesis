clear;
sigma = [3.056, 14.81, 23.38, 41.25, 2.565, 1.785, 225.3, 578.5, 1575, 570.8] * 1e-15;  % sigma's are sigma(4.605*R*C) => sigma(c) = sigma's/(4.605*R)
% sigma = [4.624, 3.692, 13.02, 0.5768, 5.273, 242.4, 388.8, 574.4, 1193]*1e-15;
capSizes = [1, 2, 4, 6, 8, 12, 16, 20, 24, 28] * 1.995 * 1e-15;
capAreas = capSizes/.001025;

deltaC = sigma/(abs(log(.01))*1e8);
deltaCOverC = deltaC./capSizes;
linearSpread = deltaCOverC.*(capAreas.^(3/2));


%Get coefficients of a line fit through the data.
coefficients = polyfit(1./sqrt(capSizes), deltaCOverC, 1);
%Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(1./sqrt(capSizes)), max(1./sqrt(capSizes)), 1000);
%Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
%Plot everything.
plot(1./sqrt(capSizes), deltaCOverC, 'b.', 'MarkerSize', 15); % Plot training data.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
grid on;


fprintf("Fitted process parameter: "+ coefficients(1) + ".\n");
% plot(1./sqrt(capSizes), deltaSigmaOverC)
