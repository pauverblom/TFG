clear;
clc;

eichmann_et_al_data = readmatrix("eichmannetal.txt");

differential_flux_data = [eichmann_et_al_data(:,1),eichmann_et_al_data(:,2)./(eichmann_et_al_data(:,1).^2)]; 

background_photons_energies = logspace(-5, 10, 600);

nested_integral_data = [background_photons_energies', nestedintegralfunc(background_photons_energies, differential_flux_data)'];


energies = logspace(2, 14, 1000); % Energies of our original photons in eV
rates = zeros(size(energies)); % Initialize interaction rates for different energies
epsilonmax = max(differential_flux_data(:,1)); % Max background photon energy in eV
mass_electron = 0.511 * 10 ^ 6; % Electron mass in eV

Radius = 7*10^(-5); % Approximate radius of the corona in pc.
thomson_cross_section = 6.65*10^(-24); % Thomson cross section in cm^2

beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);

crossection_times_s = @(s) (3 / 4) .* mass_electron .^ 2 .* ((3 - (beta(s)).^4) .* log((1 + beta(s))./(1 - beta(s))) - 2 .* beta(s) .* (2 - (beta(s)).^2));


for i = 1:length(energies)
    E = energies(i);
    integrand = @(s) crossection_times_s(s) .* interp1(nested_integral_data(:, 1), nested_integral_data(:, 2), s./(4 .* E), 'linear');
    rates(i) = 1/(8 * E^2) * thomson_cross_section * integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax) * 3.086 * 10 ^ 18; % Rates, multiplied by 3.086e18 to convert to pc^-1
end


figure;
loglog(energies, rates * Radius)
grid on;
title('Optical depth as a function of the energy of the photon')
xlabel('E [eV]');
ylabel('\tau_{\gamma\gamma}');


figure;
loglog(energies, 1./(rates))
grid on;
hold on;
yline(Radius, 'Label','Radius of the Corona','LineStyle','--' )
title('Interaction length as a function of the energy of the photon')
xlabel('E [eV]');
ylabel('\lambda [pc]');

lengths = 1./(rates);

interactionLength = [energies(1,:); lengths(1,:)];

writematrix(interactionLength, 'interactionLengthData.txt');


function value = nestedintegralfunc(epsilon_min, data)

    scaling_factor = 2.1e22; % (d/R)^2 scaling factor 
    x = data(:,1); % Energies of the background photons
    y = data(:,2); % Differential Flux for said energies

    value = zeros(size(epsilon_min)); % Initialize output array

    for i = 1:length(epsilon_min)
        idx = find(x >= epsilon_min(i), 1);
        if isempty(idx)
            value(i) = 0; % If no values match, return 0
        else
            x_subset = x(idx:end);
            n_over_epsilon_squared = y(idx:end) ./ x_subset.^2 * scaling_factor * 4 / 3e10; % n
            value(i) = trapz(x_subset, n_over_epsilon_squared, 1);
        end
    end
end


