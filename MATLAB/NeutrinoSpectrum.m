% Define the file name
filename = 'ngc1068_spectrum_95.txt';

purple = [0.7, 0.5, 0.9];
% Read the data from the file, skipping the header row
data = readmatrix(filename, 'Delimiter', '\t', 'NumHeaderLines', 1);

% Extract the columns into separate variables
energy = data(:,1); % Energy [GeV]
flux_best_fit = data(:,2); % E^2.flux (best-fit) [TeVcm^(-2)s^(-1)]
flux_95_low = data(:,3); % E^2.flux (95%-low) [TeVcm^(-2)s^(-1)]
flux_95_up = data(:,4); % E^2.flux (95%-up) [TeVcm^(-2)s^(-1)]

% Create figure
figure;
hold on;

% Fill the area between confidence limits in a light purple shade
fill([energy; flipud(energy)], [flux_95_up; flipud(flux_95_low)], purple, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Set log scale for both axes
set(gca, 'XScale', 'log', 'YScale', 'log');

% Plot the best-fit line in solid purple
loglog(energy, flux_best_fit, 'Color', purple, 'LineWidth', 2);

% Plot the confidence limits as dashed purple lines
loglog(energy, flux_95_up, '--', 'Color', purple, 'LineWidth', 1.5);
loglog(energy, flux_95_low, '--', 'Color', purple, 'LineWidth', 1.5);

% Labels and formatting
xlabel('E [GeV]');
ylabel('E^2 \Phi [TeV cm^{-2} s^{-1}]');


legend({'95% Confidence Interval', 'Best-Fit', '95% Confidence levels'}, 'Location', 'SouthWest');
