clear;
%% Global Physical Constants and Parameters (remain unchanged)
% All values are in eV, parsecs (pc), or derived units
mass_electron = 0.511e6;                        % Electron rest mass [eV]
schwarzschild_radius = 1.4e-6;                  % Schwarzschild radius [pc]
distance = 10.1e6;                              % Distance to NGC 1068 [pc]
particle_multiplicity = 20;                     % Minimum particle multiplier
bin_slope = 0.011570650269349 + 0.25;
bin_offset = 0;




pd = makedist('Normal','mu',0.5,'sigma',0.01);  % Energy fraction normal distribution

% Load and process the Eichmann et al. data
eichmann_et_al_data = readmatrix("eichmannetal.txt"); % E^2 * Phi (eV s^-1 cm^-2) vs. E (eV)
neutrino_spectrum_data = readmatrix('ngc1068_spectrum_95.txt', 'Delimiter', '\t', 'NumHeaderLines', 1);
fermi_data = [5.4757E10, 4.7090E-13;
    5.4757E09,4.8478E-13;
    3.1622E11,1.0276E-12;
    1.7320E10,7.2498E-13;
    1.7320E09,9.5810E-13;
    5.4715E08,1.6786E-12;
    1.7320E08,1.3815E-12;
    7.0679E07,7.2834E-12];

% Background photon data (independent of the parameters below)
background_photons_energies = logspace(log10(min(eichmann_et_al_data(:,1))), log10(max(eichmann_et_al_data(:,1))), 1000);
differential_flux_data = [eichmann_et_al_data(:,1), eichmann_et_al_data(:,2) ./ (eichmann_et_al_data(:,1).^2)];
thomson_cross_section = 6.65e-24; % Thomson cross section in cm^2

% Pre-calculate nested integral data (does not change with our parameters)
nested_integral_data = [background_photons_energies', ...
    nestedintegralfunc(background_photons_energies, differential_flux_data)'];

% Test photon energies and constants (also independent of our parameters)
problem_photons_energies = logspace(log10(mass_electron^2/max(eichmann_et_al_data(:,1))) + 0.1, 14, 1000); % in eV
epsilonmax = max(differential_flux_data(:,1)); % Max background photon energy in eV

% Define functions for beta and cross-section times s
beta = @(s) sqrt(1 - 4 * mass_electron^2 ./ s);
crossection_times_s = @(s) (3/4) .* mass_electron.^2 .* ...
    ((3 - beta(s).^4) .* log((1 + beta(s))./(1 - beta(s))) - 2 .* beta(s) .* (2 - beta(s).^2));

% Pre-calculate neutrino spectrum and Fermi data quantities (independent of our parameters)
neutrino_energies = neutrino_spectrum_data(:,1) * 1e9; % [eV]
E_squared_times_flux_neutrinos_best_fit = neutrino_spectrum_data(:,2) * 1e12; % [eV cm^{-2} s^{-1}]
flux_neutrinos_best_fit = E_squared_times_flux_neutrinos_best_fit ./ (neutrino_energies.^2);
deltaE = (bin_slope * neutrino_energies + bin_offset);
fermi_energies = fermi_data(:,1); % [eV]
E_squared_times_flux_fermi = fermi_data(:,2) * 6.242e11; % [eV cm^{-2} s^{-1}]

%% Define parameter ranges to loop over
%corona_mult_values = [1,2,3,4,5,6,7,8,9,10,11,13,15,18,20,25,30,35,40,50,60,70,80,90,100];    % Example corona_multiplicity values
max_angle_values   = pi;     % Example max_angle values (in radians)
%alpha_values       = [0,0.01]; 
% Example alpha values

values = [0.1,57.55];
alpha_values = values(1);
corona_mult_values = values(2);


nums_particles = round(flux_neutrinos_best_fit .* deltaE / (min(flux_neutrinos_best_fit .* deltaE)) * particle_multiplicity);

% Compute total number of particles over all parameter combinations and energy bins
num_param_combinations = length(corona_mult_values) * length(max_angle_values) * length(alpha_values);
total_particles_per_combination = sum(nums_particles);
total_particles = num_param_combinations * total_particles_per_combination;

Count = 0;
Waitbar = waitbar(0, 'Initializing waitbar');
tic;


%% Loop over parameter combinations
for cm = corona_mult_values
    % Update corona parameters that depend on corona_multiplicity
    corona_multiplicity = cm;
    corona_radius = corona_multiplicity * schwarzschild_radius;      % [pc]
    scaling_factor = (distance / corona_radius)^2;                   % Flux scaling factor

    % Recompute rates and interaction length data (scaling_factor dependent)
    rates = zeros(size(problem_photons_energies));
    for i = 1:length(problem_photons_energies)
        E = problem_photons_energies(i);
        integrand = @(s) crossection_times_s(s) .* ...
            interp1(nested_integral_data(:,1), nested_integral_data(:,2), s/(4*E), 'linear', 0);
        % Integration limits: from 4*m_e^2 to 4*E*epsilonmax
        rates(i) = 1/(8 * E^2) * thomson_cross_section * ...
            integral(integrand, 4 * mass_electron^2, 4 * E * epsilonmax) * 3.086e18; % cm->pc conversion
    end
    lengths = 1 ./ rates;
    interaction_length_data = [problem_photons_energies; lengths / scaling_factor];
    
    for max_ang = max_angle_values
        max_angle = max_ang;
        
        for a = alpha_values
            alpha = a;
            
            %% Initialize arrays for simulation results (for each parameter combo)
            final_energies = [];
            multiplicity = [];  % Record weight (2^#collisions)
            escape_steps = [];
            initial_distances = [];
            
            % Calculate number of particles per energy bin
            %% Loop over neutrino energies and simulate particle propagation
            for k = 1:length(neutrino_energies)
                num_particles = nums_particles(k);
                for p = 1:num_particles

                   
                    energy = neutrino_energies(k);
                    collisions = 0;
                    % Set an initial position outside the corona and then randomize
                    % position = [corona_radius; corona_radius; corona_radius];
                    % Ensure the initial position is inside the corona:
                        % For radial distribution use: r = corona_radius * rand^(1/alpha)
                        % (If alpha==0, you might want to set r = 0 to avoid division by zero)
                    if alpha == 0
                        radius = 0;
                    else
                        radius = corona_radius * rand()^(1/alpha);
                    end
                    [x, y, z] = sph2cart(2*pi*rand(), acos(2*rand()-1), radius);
                    position = [x; y; z];
                    storedposition(:, 1) = position;
                    initial_distances(end+1) = norm(position);
                    
                    % Random initial direction (unit vector)
                    current_dir = randn(3,1);
                    current_dir = current_dir / norm(current_dir);
                    iStep = 1;
                    
                    % Propagate the particle until it escapes the corona
                    while norm(position) < corona_radius

                        collisions = collisions + 1;
                        % Get interaction rate from precomputed interaction_length_data
                        R = 1 / interp1(interaction_length_data(1,:), interaction_length_data(2,:), energy, 'linear', 0);
                        r = rand();
                        interaction_length = -log(r) / R; % Sampled step length
                        iStep = iStep + 1;
                        
                        % Generate a random step direction within max_angle of [0;0;1]
                        theta = acos(1 - (1 - cos(max_angle))*rand());
                        phi = 2 * pi * rand();
                        local_step_direction = [sin(theta)*cos(phi);
                                                sin(theta)*sin(phi);
                                                cos(theta)];
                        
                        % Rotate local step to align with current direction
                        rot_axis = cross([0; 0; 1], current_dir);
                        angle = acos(dot([0; 0; 1], current_dir));
                        if norm(rot_axis) > 1e-6
                            rot_axis = rot_axis / norm(rot_axis);
                            c = cos(angle); s = sin(angle); C = 1 - c;
                            x_a = rot_axis(1); y_a = rot_axis(2); z_a = rot_axis(3);
                            rotMat = [x_a*x_a*C + c,    x_a*y_a*C - z_a*s, x_a*z_a*C + y_a*s;
                                      x_a*y_a*C + z_a*s, y_a*y_a*C + c,     y_a*z_a*C - x_a*s;
                                      x_a*z_a*C - y_a*s, y_a*z_a*C + x_a*s, z_a*z_a*C + c];
                        else
                            rotMat = eye(3);
                        end
                        
                        new_step = rotMat * local_step_direction * interaction_length;
                        position = position + new_step;
                        if p == 1 && k == 200
                            storedenergy(collisions) = energy;
                            storedinteractionlength(collisions) = interaction_length;
                            storedposition(:, collisions + 1) = position(:, 1);
                        end
                        if norm(new_step) > 0
                            current_dir = new_step / norm(new_step);
                        end
                        % Update energy using the random energy fraction from pd
                        energy = energy * random(pd);
                    end
                    
                    escape_steps(end+1) = iStep;
                    final_energies(end+1) = energy;
                    multiplicity = [multiplicity; 2^collisions];
                end
                % After processing all particles for energy index k, update the waitbar.
                % Weight the progress by the number of particles processed.
                Count = Count + num_particles;
                elapsedTime = toc;
                Perc = Count / total_particles;
                estimatedRemaining = elapsedTime / Perc - elapsedTime;
                Hrs = floor(estimatedRemaining / 3600);
                Min = floor(mod(estimatedRemaining, 3600) / 60);
                waitbar(Perc, Waitbar, sprintf('%.1f%%, %02d h %02d min remaining', Perc*100, Hrs, Min));

            end
            
            %% Bin the final energies to compute simulated flux
            E_min = min(final_energies);
            E_max = max(final_energies);
            E_bins = E_min; % Start bin edge array
            while E_bins(end) < E_max
                next_bin = E_bins(end) + (bin_slope * E_bins(end) + bin_offset);
                if next_bin > E_max
                    break;
                end
                E_bins = [E_bins, next_bin];
            end
            [N, ~, bin] = histcounts(final_energies, E_bins);
            weighted_counts = zeros(size(N));
            for i = 1:length(final_energies)
                idx = bin(i);
                if idx > 0
                    weighted_counts(idx) = weighted_counts(idx) + multiplicity(i);
                end
            end
            E_centers = (E_bins(1:end-1) + E_bins(2:end)) / 2;
            bin_widths = diff(E_bins);
            dN_dE = weighted_counts ./ bin_widths;
            simulated_flux = dN_dE * min(flux_neutrinos_best_fit .* deltaE) / particle_multiplicity;
            
            %% Plotting results
            %figure('Visible','off');
            %hold on;
            % Plot the neutrino spectrum with uncertainty band
            % Note: flux_95_low/up are taken from columns 3 and 4 of neutrino_spectrum_data
            %flux_95_low = neutrino_spectrum_data(:,3) * 1e12;
            %flux_95_up = neutrino_spectrum_data(:,4) * 1e12;
            %purple = [0.7, 0.5, 0.9];
            %fill([neutrino_energies; flipud(neutrino_energies)], ...
            %    [flux_95_up; flipud(flux_95_low)], purple, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            %loglog(neutrino_energies, E_squared_times_flux_neutrinos_best_fit, 'Color', purple, 'LineWidth', 2);
            %loglog(neutrino_energies, flux_95_up, '--', 'Color', purple, 'LineWidth', 1.5);
            %loglog(neutrino_energies, flux_95_low, '--', 'Color', purple, 'LineWidth', 1.5);
            %scatter(fermi_energies, E_squared_times_flux_fermi);
            %scatter(E_centers, simulated_flux .* E_centers.^2, 'filled','red');
            %set(gca, 'XScale', 'log', 'YScale', 'log');
            %xlabel('E [eV]');
            %ylabel('E^2 \Phi [eV cm^{-2} s^{-1}]');
            %title(sprintf('$R_c = %d \\, R_s, \\quad \\theta_{\\max} = %.2f, \\quad \\alpha = %.1f$', ...
            %    corona_multiplicity, max_angle, alpha), 'Interpreter', 'latex');
            %hold off;

            % Save the figure with a filename that encodes the parameter values
            %fig_filename = sprintf('corona_%d_angle_%.2f_alpha_%.1f.fig', corona_multiplicity, max_angle, alpha);
            %fig_filepath = fullfile('figures', fig_filename);
            %fprintf('Image Produced: %s \n', fig_filename)
            %savefig(fig_filepath);
            
            %% SAVE THE DATA FOR SURFACE FITTING

            filepath = sprintf('data/corona_%.2f_angle_%.2f_alpha_%.2f.csv', corona_multiplicity, max_angle, alpha);            
            writematrix([E_centers; simulated_flux .* E_centers.^2]',filepath)
            fprintf('Data Saved: %s \n', filepath)



        end
    end
end

close(Waitbar)

% %%
% x = 1:length(storedenergy); % Index values
% 
% % Create a color gradient based on the index
% colors = spring(length(x)); % Use the 'jet' colormap or any other colormap
% figure;
% 
% subplot(2, 2, [1, 3]); % Use a 2x2 grid, occupy the left side
% hold on;
% for i = 1:length(x)
%     plot3(storedposition(1, i:i+1), storedposition(2, i:i+1), storedposition(3, i:i+1), 'Color', colors(i, :), 'LineWidth', 2)
% end
% [X, Y, Z] = sphere(50); % 50 controls the resolution
% 
% X = X * corona_radius;
% Y = Y * corona_radius;
% Z = Z * corona_radius;
% surf(X,Y,Z, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor','white'); % Transparent sphere
% axis equal;
% view(30, 30); % Adjust these values to change the view
% title('3D Position Plot');
% xlabel('x [pc]')
% ylabel('y [pc]')
% zlabel('z [pc]')
% 
% hold off;
% 
% 
% % Plot the line with changing colors
% subplot(2, 2, 2); % Top-right subplot
% hold on;
% for i = 1:length(x)-1
%     plot(x(i:i+1), storedenergy(i:i+1), 'Color', colors(i, :), 'LineWidth', 2);
% end
% title('Energy of the particle');
% xlabel('Step');
% ylabel('E [eV]');
% yscale log;
% hold off;
% 
% subplot(2, 2, 4); % Bottom-right subplot
% hold on;
% for i = 1:length(x)-1
%     plot(x(i:i+1), storedinteractionlength(i:i+1), 'Color', colors(i, :), 'LineWidth', 2);
% end
% yline(corona_radius, 'Label','Radius of the Corona','LineStyle','--')
% title('Interaction Length of the photon');
% xlabel('Step');
% ylabel('\lambda [pc]');
% yscale log;
% hold off;

%figure;

%histogram(escape_steps)

%figure;

%histogram(initial_distances)


%%


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
