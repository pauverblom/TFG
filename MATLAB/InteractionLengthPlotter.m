clear all;
%% Global Physical Constants and Parameters (remain unchanged)
% All values are in eV, parsecs (pc), or derived units
mass_electron = 0.511e6;                        % Electron rest mass [eV]
schwarzschild_radius = 1.4e-6;                  % Schwarzschild radius [pc]
distance = 10.1e6;                              % Distance to NGC 1068 [pc]

corona_mult_values = [1, 10, 50, 100];


% Load and process the Eichmann et al. data
eichmann_et_al_data = readmatrix("eichmannetal.txt"); % E^2 * Phi (eV s^-1 cm^-2) vs. E (eV)


% Background photon data (independent of the parameters below)
background_photons_energies = logspace(log10(min(eichmann_et_al_data(:,1))), log10(max(eichmann_et_al_data(:,1))), 1000);
differential_flux_data = [eichmann_et_al_data(:,1), eichmann_et_al_data(:,2) ./ (eichmann_et_al_data(:,1).^2)];
thomson_cross_section = 6.65e-24; % Thomson cross section in cm^2

% Pre-calculate nested integral data (does not change with our parameters)
nested_integral_data = [background_photons_energies', nestedintegralfunc(background_photons_energies, differential_flux_data)'];

% Test photon energies and constants (also independent of our parameters)
problem_photons_energies = logspace(log10(mass_electron^2/max(eichmann_et_al_data(:,1))) + 0.1, 14, 1000); % in eV
epsilonmax = max(differential_flux_data(:,1)); % Max background photon energy in eV

% Define functions for beta and cross-section times s
beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);
crossection_times_s = @(s) (3/4) .* mass_electron.^2 .* ...
    ((3 - beta(s).^4) .* log((1 + beta(s))./(1 - beta(s))) - 2 .* beta(s) .* (2 - beta(s).^2));

depths(1, :) = problem_photons_energies;

for j = 1:length(corona_mult_values)
    corona_mult_value = corona_mult_values(j);
    corona_radius = schwarzschild_radius * corona_mult_value;
    scaling_factor = (distance / corona_radius)^2;                   % Flux scaling factor
    rates = zeros(size(problem_photons_energies));
        for i = 1:length(problem_photons_energies)
            E = problem_photons_energies(i);
            integrand = @(s) crossection_times_s(s) .* ...
                interp1(nested_integral_data(:,1), nested_integral_data(:,2), s/(4*E), 'linear', 0);
            % Integration limits: from 4*m_e^2 to 4*E*epsilonmax
            rates(i) = 1/(8 * E^2) * thomson_cross_section * ...
                integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax) * 3.086e18; % cm->pc conversion
        end
    
    % Pre-calculate neutrino spectrum and Fermi data quantities (independent of our parameters)
    %%
    % Plot Optical Depth
    depths(j + 1, :) = rates * corona_radius * scaling_factor;
    
    % Plot Interaction Length
    %figure;
    %loglog(problem_photons_energies, 1./rates / scaling_factor)
    %grid on;
    %hold on;
    %yline(corona_radius, 'Label', 'Radius of the Corona', 'LineStyle', '--')
    %title('Interaction length as a function of the energy of the photon')
    %xlabel('E [eV]');
    %ylabel('\lambda [pc]');

end

%%

depths = depths';

writematrix(depths, 'optical_depths_as_function_of_Rc.txt')

function value = nestedintegralfunc(epsilon_min, data)

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
        % Calculate the integrand: n / ε²
        n_over_epsilon_squared = (y_subset ./ (x_subset.^2)) * 4 / 3e10;
        % Only integrate if we have at least 2 points
        if numel(x_subset) < 2
            value(i) = 0;
        else
            value(i) = trapz(x_subset, n_over_epsilon_squared);
        end
    end
end
end
