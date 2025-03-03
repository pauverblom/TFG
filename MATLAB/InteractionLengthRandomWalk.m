clc; clear; close all;

 % Load the interaction length matrix

interaction_length_data = readmatrix('interactionLengthData.csv'); % Automatically detects delimiter

corona_radius = 7e-5;     % Radius of the corona in parsecs. If a particle is further than this radius, we consider it has escaped. 
%step_size = 1;          % Step size
num_particles = 200;    % Number of particles
energies = logspace(12,14,num_particles);
max_angle_values = linspace(pi, pi, 1); % Values of the max turning angle
%max_radius = 7e-5;
%alpha_values = [0.01, 0.1, 0.5, 1, 3, 5, 10];
alpha_values = 1;
n = 10000;              % Number of points per distribution

figure;
hold on; % Hold on so that multiple curves appear on the same figure

for q = 1:length(alpha_values)
    alpha = alpha_values(q);
    escape_steps = zeros(length(max_angle_values), num_particles);
    mean_escape_steps = zeros(length(max_angle_values), 1);

    % Iterate over different max_angle values
    for t = 1:length(max_angle_values)
        max_angle = max_angle_values(t);
        
        % Simulate multiple particles for the same angle
        for p = 1:num_particles
            position = zeros(3, 1);
            energy = energies(p);
            %radius = max_radius * rand()^(1/(alpha));   
            %position(:, 1) = sph2cart(2 * pi * rand(), acos(2 * rand() - 1), radius);
    
            % Random initial direction
            current_dir = randn(3, 1);
            current_dir = current_dir / norm(current_dir);  % Normalize to unit vector
    
            i = 1;
            
            while norm(position) < corona_radius
                interaction_length = interp1(interaction_length_data(1, :), interaction_length_data(2, :), energy, 'linear');
                i = i + 1; % Increase the number of steps
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
    
                % Rotate local_step into our direction using R and multiply by interaction_length
                new_step = R * local_step_direction * interaction_length;
    
                % Update position
                position(1, 1) = position(1, 1) + new_step(1, 1);
                position(2, 1) = position(2, 1) + new_step(2, 1);
                position(3, 1) = position(3, 1) + new_step(3, 1);
    
                % Update last direction
                current_dir = new_step / norm(new_step);
                energy = energy*0.95;
            end
    
            % Store final positions for a singular particles for a given max_angle
            escape_steps(t, p) = i - 1; % I subtract 1 since array indices start at one, but we have 1 step less.
            finalenergy(t, p) = energy;
        end
        mean_escape_steps(t, 1) = mean(escape_steps(t, :));
        % Compute average final distance for a given max_angle
    end
    %plot(max_angle_values, mean_escape_steps, 'DisplayName', ['\alpha = ' num2str(alpha)]);
end

plot(finalenergy)

% Label the axes and add the legend
%xlabel('Max Angle');
%ylabel('Mean Escape Steps');
%title('Mean Escape Steps vs. Max Angle for Different \alpha Values');
%legend('Location', 'best');
%hold off;