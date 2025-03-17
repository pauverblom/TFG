clear;
clc;

% Load and process the Eichmann et al. data
eichmann_et_al_data = readmatrix("eichmannetal.txt"); % E^2 * Phi (eV s^-1 cm^-2) vs. E (eV)
differential_flux_data = [eichmann_et_al_data(:,1), eichmann_et_al_data(:,2) ./ (eichmann_et_al_data(:,1).^2)];

% Create background photons energies using logspace (ensure logarithmic scale by using log10)
background_photons_energies = logspace(log10(min(eichmann_et_al_data(:,1))), log10(max(eichmann_et_al_data(:,1))), 1000);
nested_integral_data = [background_photons_energies', nestedintegralfunc(background_photons_energies, differential_flux_data)']; % (d/R)^2 scaling factor 

% Test photon energies and constants
mass_electron = 0.511 * 10^6; % Electron mass in eV
test_photons_energies = logspace(log10(mass_electron^2/max(eichmann_et_al_data(:,1))) + 0.1, 14, 1000); % Energies of our original photons in eV
rates = zeros(size(test_photons_energies));
epsilonmax = max(differential_flux_data(:,1)); % Max background photon energy in eV

Radius = 7 * 10^(-5); % Approximate radius of the corona in pc.
thomson_cross_section = 6.65 * 10^(-24); % Thomson cross section in cm^2

% Define functions for beta and cross-section times s
beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);
crossection_times_s = @(s) (3 / 4) .* mass_electron.^2 .* ...
    ((3 - beta(s).^4) .* log((1 + beta(s)) ./ (1 - beta(s))) - 2 .* beta(s) .* (2 - beta(s).^2));

% Loop over test photon energies to compute rates
for i = 1:length(test_photons_energies)
    E = test_photons_energies(i);
    % Use interp1 with extrapolation default (returning 0) for out-of-bound values
    integrand = @(s) crossection_times_s(s) .* ...
        interp1(nested_integral_data(:,1), nested_integral_data(:,2), s./(4 .* E), 'linear', 0);
    
    % Integration limits: from 4*m_e^2 to 4*E*epsilonmax
    rates(i) = 1/(8 * E^2) * thomson_cross_section * ...
        integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax) * 3.086 * 10^18; 
end

% Plot Optical Depth
figure;
loglog(test_photons_energies, rates * Radius)
grid on;
title('Optical depth as a function of the energy of the photon')
xlabel('E [eV]');
ylabel('\tau_{\gamma\gamma}');

% Plot Interaction Length
figure;
loglog(test_photons_energies, 1./rates)
grid on;
hold on;
yline(Radius, 'Label', 'Radius of the Corona', 'LineStyle', '--')
title('Interaction length as a function of the energy of the photon')
xlabel('E [eV]');
ylabel('\lambda [pc]');

% Save interaction length data
lengths = 1./rates;
interactionLength = [test_photons_energies; lengths];
writematrix(interactionLength, 'interactionLengthData.txt');

% Nested integral function
function value = nestedintegralfunc(epsilon_min, data)
 
    scaling_factor = 2.1e22;
    
    x = data(:,1); % Energies of the background photons
    y = data(:,2); % Differential flux for said energies

    % Force x and y to be column vectors
    x = x(:);
    y = y(:);
    
    value = zeros(size(epsilon_min)); % Initialize output array

    for i = 1:length(epsilon_min)
        % Find the first index where x is >= the current epsilon_min value
        idx = find(x >= epsilon_min(i), 1);
        if isempty(idx)
            value(i) = 0; % If no values match, return 0
        else
            x_subset = x(idx:end);
            y_subset = y(idx:end);
            % Ensure both subsets are column vectors
            x_subset = x_subset(:);
            y_subset = y_subset(:);
            
            % Calculate the integrand: n / ε²
            n_over_epsilon_squared = (y_subset ./ (x_subset.^2)) * scaling_factor * 4 / 3e10;
            
            % Only integrate if we have at least 2 points
            if numel(x_subset) < 2
                value(i) = 0;
            else
                value(i) = trapz(x_subset, n_over_epsilon_squared);
            end
        end
    end
end
