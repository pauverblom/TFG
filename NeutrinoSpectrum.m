% Define the file name
filename = 'ngc1068_spectrum_95.txt';

% Read the data from the file, skipping the header row
data = readmatrix(filename, 'Delimiter', '\t', 'NumHeaderLines', 1);

% Extract the columns into separate variables
energy = data(:,1); % Energy [GeV]
flux_best_fit = data(:,2); % E^2.flux (best-fit) [TeVcm^(-2)s^(-1)]
flux_95_low = data(:,3); % E^2.flux (95%-low) [TeVcm^(-2)s^(-1)]
flux_95_up = data(:,4); % E^2.flux (95%-up) [TeVcm^(-2)s^(-1)]

% Display a portion of the data
disp('First few rows of the data:');
disp(table(energy(1:5), flux_best_fit(1:5), flux_95_low(1:5), flux_95_up(1:5), ...
    'VariableNames', {'Energy_GeV', 'Flux_BestFit', 'Flux_95Low', 'Flux_95Up'}));

plot(energy, flux_best_fit)
legend('Energy GeV', 'E^2 flux')