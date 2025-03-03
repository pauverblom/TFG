clear;
clc;

energies = logspace(3.42, 14, 10000); % Energies of our original photons
rates = zeros(size(energies)); % Initialize interaction rates for different energies
epsilonmax = 10 ^ 8; % Max background photon energy in eV
mass_electron = 0.511 * 10 ^ 6; % Electron mass in eV

Radius = 7*10^(-5); % Approximate radius of the corona in pc.
thomson_cross_section = 6.65*10^(-24); % Thomson cross section in cm^2

beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);

crossection_times_s = @(s) (3 / 4) * mass_electron ^ 2 .* ((3 - (beta(s)).^4) .* log((1 + beta(s))./(1 - beta(s))) - 2 * beta(s) .* (2 - (beta(s)).^2));


for i = 1:length(energies)
    E = energies(i);
    integrand = @(s) crossection_times_s(s) .* integralfunc(s/(4 * E));
    rates(i) = 1/(8 * E^2) * thomson_cross_section * integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax) * 3.086 * 10 ^ 18; % Rates, multiplied by 3.086e18 to convert to pc^-1
end


figure;
loglog(energies, rates * Radius)
title('Optical depth as a function of the energy of the photon')
xlabel('E [eV]');
ylabel('\tau_{\gamma\gamma}');


figure;
loglog(energies, 1./(rates * Radius))
hold on;
yline(Radius, 'Label','Radius of the Corona','LineStyle','--' )
title('Interaction length as a function of the energy of the photon')
xlabel('E [eV]');
ylabel('\lambda [pc]');

lengths = 1./(rates * Radius);

interactionLength = [energies(1,:); lengths(1,:)];

writematrix(interactionLength, 'interactionLengthData.csv');

function y = integralfunc(epsilon_min)
     % Initialize output as a zero array of the same size as epsilon
    y = zeros(size(epsilon_min));
    scaling_factor = 1.9*10^22;
    for i = 1:numel(epsilon_min)
        e = epsilon_min(i); % Extract each element separately
        if all(e < 4.3*10^-2)
            y(i) = 4/3*(10^(-0.3)*(4.3*10^(-2) - e) + 6.75*10^(-3)) * scaling_factor;
        elseif all(4.3*10^-2 < e & e < 2.56*10^5)
            y(i) = 4/3*(-(10^(-6.1))/3.25*((2.56*10^5)^(-3.25)-e^(-3.25))+2.353*10^(-25)) * scaling_factor;
        elseif all(2.56*10^5 < e)
            y(i) = 4/3*(-(10^25)/9*((10^8)^(-9)-e^(-9))) * scaling_factor;
        else
            y(i) = 0; % Default value outside range
        end
    end
end

