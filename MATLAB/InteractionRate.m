clear;
clc;

energies = logspace(3.42, 14, 1800); % Energies of our original photons
rates = zeros(size(energies)); % Initialize interaction rates for different energies
epsilonmax = 10^8; % Max background photon energy in eV
mass_electron = 0.511*10^6; % Electron mass in eV

beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);

crossection_times_s = @(s) 2 * pi * (1/137)^2 * ((3 - (beta(s)).^2) * log((1 + beta(s))/(1 - beta(s))) - 2 * beta(s) + (2 - (beta(s)).^2));



for i = 1:length(energies)
    E = energies(i);
    integrand = @(s) crossection_times_s(s) .* integralfunc(s/(4 * E));
    rates(i) = 1/(8 * E^2) * integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax);
end

loglog(energies, rates)


function y = integralfunc(epsilon_min)
     % Initialize output as a zero array of the same size as epsilon
    y = zeros(size(epsilon_min));

    for i = 1:numel(epsilon_min)
        e = epsilon_min(i); % Extract each element separately
        if all(e < 4.3*10^-2)
            y(i) = 1/3*(10^(-0.3)*(4.3*10^(-2) - e) + 6.75*10^(-3));
        elseif all(4.3*10^-2 < e & e < 2.56*10^5)
            y(i) = 1/3*(-(10^(-6.1))/3.25*((2.56*10^5)^(-3.25)-e^(-3.25))+2.353*10^(-25));
        elseif all(2.56*10^5 < e)
            y(i) = 1/3*(-(10^25)/9*((10^8)^(-9)-e^(-9)));
        else
            y(i) = 0; % Default value outside range
        end
    end
end