clc; clear; close all;

num_step_values = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 4000, 6000, 8000, 10000];        % We are going to check the RMSD for a different number of steps
step_size = 1;          % Step size
num_particles = 5000;   % Number of particles
max_angle = pi; % We check only for max_angle = pi, since we should get the behaviour expected from a random walk

avg_distances = zeros(size(num_step_values));
expected_rmsd = zeros(size(num_step_values));

% Iterate over different num_step values
for t = 1:length(num_step_values)
    num_steps = num_step_values(t);
    expected_rmsd(t) = sqrt(num_steps) * step_size;
    final_positions = zeros(3, num_particles);

    % Simulate multiple particles
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

        % Store final positions
        final_positions(:, p) = position(:, num_steps);
    end
    % Compute average final distance
    avg_distances(t) = mean(sqrt(sum(final_positions.^2, 1)));

    if t == length(num_step_values)
        % Calculate mean displacement
        mean_displacement = mean(final_positions, 2);
        fprintf('Final mean displacement: [%.4f, %.4f, %.4f]\n', mean_displacement);
        
        % Calculate covariance matrix
        cov_matrix = cov(final_positions');
        disp('Final covariance matrix:');
        disp(cov_matrix);
    end
end
%%
% Plot RMSD results
figure;
plot(num_step_values, expected_rmsd, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(num_step_values, avg_distances, 'r--x', 'LineWidth', 1.5, 'MarkerSize', 8);
legend('Theoretical RMSD ($\sqrt{N}$)', 'Simulated RMSD', 'Location', 'northwest', 'Interpreter', 'latex');
xlabel('Number of Steps (N)');
ylabel('RMS Displacement');
title('3D Random Walk: RMSD vs Number of Steps');
grid on;

% Plot position histograms for largest N
figure;
subplot(1,3,1);
histogram(final_positions(1,:), 'FaceColor', 'r', 'EdgeColor', 'k');
title('X Positions');
xlabel('Position');
ylabel('Frequency');

subplot(1,3,2);
histogram(final_positions(2,:), 'FaceColor', 'g', 'EdgeColor', 'k');
title('Y Positions');
xlabel('Position');

subplot(1,3,3);
histogram(final_positions(3,:), 'FaceColor', 'b', 'EdgeColor', 'k');
title('Z Positions');
xlabel('Position');

sgtitle('Final Position Distributions (N = 10,000)');