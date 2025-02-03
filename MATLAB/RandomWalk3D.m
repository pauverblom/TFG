clc; clear; close all;

num_steps = 400;        % Number of steps
step_size = 1;          % Step size
num_particles = 2000;   % Number of particles
max_angle_values = nonLinspace(0, pi, 20, "exp10"); % We want to test more points closer to the origin, to resolve more detail. Hence the use of the nonLinspace function

avg_distances = zeros(size(max_angle_values));
expected_rmsd = sqrt(num_steps) * step_size;

% Iterate over different max_angle values
for t = 1:length(max_angle_values)
    max_angle = max_angle_values(t);
    final_positions = zeros(3, num_particles);
    
    % Simulate multiple particles for the same angle
    for p = 1:num_particles
        position = zeros(3, num_steps);

        %radius = rand()^(1/10);
        %position(:, 1) = sph2cart(2 * pi * rand(), acos(2 * rand() - 1), radius);

        % Random initial direction
        current_dir = randn(3, 1);
        current_dir = current_dir / norm(current_dir);  % Normalize to unit vector

        for i = 2:num_steps
            % Generate random local_step_direction within max_angle of [0; 0; 1]
            theta = acos(1 - (1 - cos(max_angle)) * rand());  % Polar angle uniformly distributed within spherical cap
            phi = 2 * pi * rand();  % Azimuthal angle uniformly distributed between 0 and 2*pi

            % Convert to Cartesian x,y,z coordinates (assuming initial direction along z-axis)

            local_step_direction = [sin(theta) * cos(phi);
                                    sin(theta) * sin(phi);
                                              cos(theta)];

            % Find rotation matrix that aligns [0,0,1] with current_direction

            axis = cross([0; 0; 1], current_dir); % Axis of rotation
            angle = acos(dot([0; 0; 1], current_dir));  % Angle we need to rotate
            
            if norm(axis) > 1e-6
                axis = axis / norm(axis);

                % Rodrigues' rotation formula (inlined for speed)
                c = cos(angle);
                s = sin(angle);
                C = 1 - c;

                x = axis(1); 
                y = axis(2); 
                z = axis(3);

                xy = x * y; % Precompute diagonal terms to increase speed (barely)
                xz = x * z;
                yz = y * z;

                R = [x * x * C + c ,  xy * C - z * s, xz * C + y * s;
                     xy * C + z * s,  y * y * C + c , yz * C - x * s;
                     xz * C - y * s,  yz * C + x * s, z * z * C + c];
            else
                R = eye(3); % If the deviation is too small, we can spare the computation
            end

            % Rotate local_step into our direction using R and multiply by step_size
            new_step = R * local_step_direction * step_size;

            % Update position
            position(1, i) = position(1, i - 1) + new_step(1, 1);
            position(2, i) = position(2, i - 1) + new_step(2, 1);
            position(3, i) = position(3, i - 1) + new_step(3, 1);

            % Update last direction
            current_dir = new_step / norm(new_step);
        end

        % Store final positions for a singular particles for a given max_angle
        final_positions(:, p) = position(:, num_steps);
    end
    % Compute average final distance for a given max_angle
    avg_distances(t) = mean(sqrt(sum(final_positions.^2, 1)));
end
%%
% Plot results
figure;
plot(max_angle_values, avg_distances, 'r-o', 'LineWidth', 1.5);
hold on;

% Check expected behavior at max_angle = pi
if any(max_angle_values == pi)
    fprintf('Expected RMSD: %.2f\n', expected_rmsd);
    fprintf('Measured RMSD at max_angle = pi: %.2f\n', avg_distances(end)); % Average distance at max_angle = pi
    
    % Plot separate histograms of final positions for max_angle = pi
    figure;
    histogram(final_positions(1, :), 'FaceColor', 'r', 'EdgeColor', 'k', 'NumBins', 20);
    title('Histogram of Final X Positions (\vartheta_{max} = π)');
    xlabel('Positions');
    ylabel('Frequency');
    grid on;
    hold on;
    fprintf('Mean X position: %.2f\n', mean(final_positions(1, :)));
    fprintf('X position variance: %.2f\n', var(final_positions(1, :)));
    
    histogram(final_positions(2, :), 'FaceColor', 'g', 'EdgeColor', 'k', 'NumBins', 20);
    fprintf('Mean Y position: %.2f\n', mean(final_positions(2, :)));
    fprintf('Y position variance: %.2f\n', var(final_positions(2, :)));
    
    histogram(final_positions(3, :), 'FaceColor', 'b', 'EdgeColor', 'k', 'NumBins', 20);
    fprintf('Mean Z position: %.2f\n', mean(final_positions(3, :)));
    fprintf('Z position variance: %.2f\n', var(final_positions(3, :)));
    legend('X positions', 'Y positions', 'Z positions');
end

figure;
plot3(position(1, :), position(2, :), position(3, :));
