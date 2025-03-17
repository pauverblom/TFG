
mass_electron = 0.511e6; % Electron mass in eV
corona_radius = 7e-5; % Radius of the corona in parsecs. If a particle is further than this radius, we consider it has escaped.
distance = 10.1e6; % Distance to NGC 1068 in pc.
particle_multiplicity = 20; % To increase the minimum amount of nums_particles for the minimum flux.
scaling_factor = (distance/corona_radius)^2;
alpha = 2;

% Load and process the Eichmann et al. data
eichmann_et_al_data = readmatrix("eichmannetal.txt"); % E^2 * Phi (eV s^-1 cm^-2) vs. E (eV)
neutrino_spectrum_data = readmatrix('ngc1068_spectrum_95.txt', 'Delimiter', '\t', 'NumHeaderLines', 1);
fermi_data = [6.0050E09, 5.4513E-13;6.0050E10,4.9875E-13;6.0050E08,1.2867E-12;2.0017E09,1.0198E-12;2.0017E08,2.2596E-12;1.0008E09,2.7675E-12;6.0050E09,7.8144E-13;6.0050E10,2.0034E-12;6.0050E08,1.8611E-12;2.0017E09,8.6998E-13;2.0017E08,2.6401E-12;6.0050E09,1.9979E-12;6.0050E10,2.9960E-12;6.0050E08,1.9921E-12;2.0017E09,1.7621E-12]; % Some Data from Fermi in barn



background_photons_energies = logspace(log10(min(eichmann_et_al_data(:,1))), log10(max(eichmann_et_al_data(:,1))), 1000);
differential_flux_data = [eichmann_et_al_data(:,1), eichmann_et_al_data(:,2) ./ (eichmann_et_al_data(:,1).^2)];


thomson_cross_section = 6.65 * 10^(-24); % Thomson cross section in cm^2



% Create background photons energies using logspace (ensure logarithmic scale by using log10)
nested_integral_data = [background_photons_energies', nestedintegralfunc(background_photons_energies, differential_flux_data)']; % (d/R)^2 scaling factor

% Test photon energies and constants
problem_photons_energies = logspace(log10(mass_electron^2/max(eichmann_et_al_data(:,1))) + 0.1, 14, 1000); % Energies of our problem photons in eV. Adding +0.1 at the beginning in order to avoid diverging integrals.
rates = zeros(size(problem_photons_energies));
epsilonmax = max(differential_flux_data(:,1)); % Max background photon energy in eV


% Define functions for beta and cross-section times s
beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);
crossection_times_s = @(s) (3 / 4) .* mass_electron.^2 .* ...
    ((3 - beta(s).^4) .* log((1 + beta(s)) ./ (1 - beta(s))) - 2 .* beta(s) .* (2 - beta(s).^2));

% Loop over test photon energies to compute rates
for i = 1:length(problem_photons_energies)
    E = problem_photons_energies(i);
    % Use interp1 with extrapolation default (returning 0) for out-of-bound values
    integrand = @(s) crossection_times_s(s) .* ...
        interp1(nested_integral_data(:,1), nested_integral_data(:,2), s./(4 .* E), 'linear', 0);

    % Integration limits: from 4*m_e^2 to 4*E*epsilonmax
    rates(i) = 1/(8 * E^2) * thomson_cross_section * ...
        integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax) * 3.086 * 10^18; % Multiplied by 3.086e18 to convert from cm to pc
end

% Save interaction length data
lengths = 1./rates;

% Nested integral function
interaction_length_data = [problem_photons_energies; lengths / scaling_factor];



neutrino_energies = neutrino_spectrum_data(:,1) * 10^9; % Energy of the neutrinos [eV]

E_squared_times_flux_neutrinos_best_fit = neutrino_spectrum_data(:,2) * 10^12; % E^2 * Flux (neutrinos) (best-fit) [eV cm^(-2)s^(-1)]

flux_neutrinos_best_fit = E_squared_times_flux_neutrinos_best_fit ./ (neutrino_energies.^2); % Divide by E^2 to get the differential flux [eV^(-1) cm^(-2)s^(-1)]

deltaE = (0.011570650269349 * neutrino_energies -0.000090594152042); % In order to have correct binning


fermi_energies = fermi_data(:,1); % Fermi energies in eV
E_squared_times_flux_fermi = fermi_data(:,2) * 6.242e11; % E^2 * Flux (fermi) [eV cm^(-2)s^(-1)] (1 Erg = 6.242e11 eV)


nums_particles = round(flux_neutrinos_best_fit .* deltaE / (min(flux_neutrinos_best_fit .*deltaE)) * particle_multiplicity); % Divide by (min(flux_neutrinos_best_fit .*deltaE)) and multiply by particle_multiplicity to have more test particles for the flux.

max_angle = pi;
final_energies = [];
multiplicity = [];  % initialize an array for weights

for k = 1:length(neutrino_energies)
    num_particles = nums_particles(k);
    for p = 1:num_particles
        energy = neutrino_energies(k)
        collisions = 0;
        position = zeros(3, 1);

        radius = corona_radius * rand()^(1/(alpha));
        position(:, 1) = sph2cart(2 * pi * rand(), acos(2 * rand() - 1), radius);
        % Random initial direction
        current_dir = randn(3, 1);
        current_dir = current_dir / norm(current_dir);  % Normalize to unit vector
        i = 1;

        while norm(position) < corona_radius

            collisions = collisions + 1;
            R = 1 / interp1(interaction_length_data(1, :), interaction_length_data(2, :), energy, 'linear'); % Interaction rate
            r = rand(); % Generate uniform random number in [0,1]
            interaction_length = -log(r) / R; % Sample random interaction length
            i = i + 1; % Increase the number of steps
            % Generate random local_step_direction within max_angle of [0; 0; 1]
            theta = acos(1 - (1 - cos(max_angle)) * rand());  % Polar angle uniformly distributed within spherical cap
            phi = 2 * pi * rand();  % Azimuthal angle uniformly distributed between 0 and 2*pi

            % Convert to Cartesian x,y,z coordinates (assuming initial direction along z-axis)

            local_step_direction = [sin(theta) * cos(phi);
                sin(theta) * sin(phi);
                cos(theta)];

            % Find rotation matrix that aligns [0,0,1] with current_direction

            rot_axis = cross([0; 0; 1], current_dir); % Axis of rotation
            angle = acos(dot([0; 0; 1], current_dir));  % Angle we need to rotate

            if norm(rot_axis) > 1e-6
                rot_axis = rot_axis / norm(rot_axis);

                % Rodrigues' rotation formula (inlined for speed)
                c = cos(angle);
                s = sin(angle);
                C = 1 - c;

                x = rot_axis(1);
                y = rot_axis(2);
                z = rot_axis(3);

                xy = x * y; % Precompute diagonal terms to increase speed (barely)
                xz = x * z;
                yz = y * z;

                rotMat = [x * x * C + c ,  xy * C - z * s, xz * C + y * s;
                    xy * C + z * s,  y * y * C + c , yz * C - x * s;
                    xz * C - y * s,  yz * C + x * s, z * z * C + c];
            else
                rotMat = eye(3); % If the deviation is too small, we can spare the computation
            end

            % Rotate local_step into our direction using R and multiply by interaction_length
            new_step = rotMat * local_step_direction * interaction_length;

            % Update position
            position(1, 1) = position(1, 1) + new_step(1, 1);
            position(2, 1) = position(2, 1) + new_step(2, 1);
            position(3, 1) = position(3, 1) + new_step(3, 1);
            % Update last direction
            if p == 1 && k == 200
                storedenergy(i - 1) = energy;
                storedinteractionlength(i - 1) = interaction_length;
                storedposition(:, i) = position(:, 1);
            end
            current_dir = new_step / norm(new_step);

            energy = energy * 0.5; % Multiply E by f = 0.5, assuming we are in the case that E * epsilon / m_e^2 ~ 1
        end
        final_energies(end+1) = energy; % Append the final energy
        multiplicity = [multiplicity; 2^collisions];  % record multiplicity for this particle
    end
end

E_min = min(final_energies);
E_max = max(final_energies);
E_bins = E_min; % Start bin edges

while E_bins(end) < E_max
    next_bin = E_bins(end) + (0.011570650269349 * E_bins(end) -0.000090594152042);
    if next_bin > E_max
        break;
    end
    E_bins = [E_bins, next_bin]; % Append next bin edge
end

[N, ~, bin] = histcounts(final_energies, E_bins);


weighted_counts = zeros(size(N));
for i = 1:length(final_energies)
    % find the bin index for this particle and add its weight
    idx = bin(i);
    if idx > 0
        weighted_counts(idx) = weighted_counts(idx) + multiplicity(i);
    end
end

E_centers = (E_bins(1:end-1) + E_bins(2:end)) / 2;
bin_widths = diff(E_bins);

dN_dE = weighted_counts ./ bin_widths;

simulated_flux = dN_dE * min(flux_neutrinos_best_fit .*deltaE) / 10;

flux_95_low = neutrino_spectrum_data(:,3)*10^12; % E^2.flux (95%-low) [TeVcm^(-2)s^(-1)]
flux_95_up = neutrino_spectrum_data(:,4)*10^12; % E^2.flux (95%-up) [TeVcm^(-2)s^(-1)]
purple = [0.7, 0.5, 0.9];




figure;
hold on;
fill([neutrino_energies; flipud(neutrino_energies)], [flux_95_up; flipud(flux_95_low)], purple, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

loglog(neutrino_energies, E_squared_times_flux_neutrinos_best_fit, 'Color', purple, 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
loglog(neutrino_energies, flux_95_up, '--', 'Color', purple, 'LineWidth', 1.5);
loglog(neutrino_energies, flux_95_low, '--', 'Color', purple, 'LineWidth', 1.5);
scatter(fermi_energies, E_squared_times_flux_fermi)

scatter(E_centers, simulated_flux .* E_centers .^ 2, 'filled','red')
xlabel('E [eV]');
ylabel('E^2 \Phi [eV cm^{-2} s^{-1}]');
title('Spectral Energy Distribution')


hold off;


x = 1:length(storedenergy); % Index values

% Create a color gradient based on the index
colors = spring(length(x)); % Use the 'jet' colormap or any other colormap
figure;

subplot(2, 2, [1, 3]); % Use a 2x2 grid, occupy the left side
hold on;
for i = 1:length(x)
    plot3(storedposition(1, i:i+1), storedposition(2, i:i+1), storedposition(3, i:i+1), 'Color', colors(i, :), 'LineWidth', 2)
end
[X, Y, Z] = sphere(50); % 50 controls the resolution

X = X * corona_radius;
Y = Y * corona_radius;
Z = Z * corona_radius;
surf(X,Y,Z, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor','white'); % Transparent sphere
axis equal;
view(30, 30); % Adjust these values to change the view
title('3D Position Plot');
xlabel('x [pc]')
ylabel('y [pc]')
zlabel('z [pc]')

hold off;


% Plot the line with changing colors
subplot(2, 2, 2); % Top-right subplot
hold on;
for i = 1:length(x)-1
    plot(x(i:i+1), storedenergy(i:i+1), 'Color', colors(i, :), 'LineWidth', 2);
end
title('Energy of the particle');
xlabel('Step');
ylabel('E [eV]');
yscale log;
hold off;

subplot(2, 2, 4); % Bottom-right subplot
hold on;
for i = 1:length(x)-1
    plot(x(i:i+1), storedinteractionlength(i:i+1), 'Color', colors(i, :), 'LineWidth', 2);
end
yline(corona_radius, 'Label','Radius of the Corona','LineStyle','--')
title('Interaction Length of the photon');
xlabel('Step');
ylabel('\lambda [pc]');
yscale log;
hold off;



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
        x_subset = x_subset(:);
        y_subset = y_subset(:);

        % Calculate the integrand: n / ε²
        %n_over_epsilon_squared = (y_subset ./ (x_subset.^2)) * scaling_factor * 4 / 3e10;
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
