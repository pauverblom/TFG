
% Number of random points to generate
N = 1000;
max_angle = 0.7;

% Generate a random direction for the spherical cap
theta_d = acos(2 * rand() - 1);  % Random inclination
phi_d = 2 * pi * rand();  % Random azimuth

% Initialize arrays to store the points
phi_vals = zeros(1, N);
theta_vals = zeros(1, N);
x_vals = zeros(1, N);
y_vals = zeros(1, N);
z_vals = zeros(1, N);

% Generate points with rotated angles
for i = 1:N
    % Generate local angles around the +Z axis
    phi = 2 * pi * rand();  % Uniform distribution between 0 and 2*pi
    theta = acos(1 - (1 - cos(max_angle)) * rand());  % Uniform distribution within spherical cap

    % Rotate the angles to align with (theta_d, phi_d)
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);

    % Convert rotated angles
    theta_rot = acos(z * cos(theta_d) + x * sin(theta_d));  % New inclination
    phi_rot = atan2(y, x * cos(theta_d) - z * sin(theta_d)) + phi_d;  % New azimuth

    % Store values
    phi_vals(i) = mod(phi_rot, 2*pi);  % Ensure azimuth is in [0, 2π]
    theta_vals(i) = theta_rot;

    % Convert back to Cartesian coordinates
    x_vals(i) = sin(theta_rot) * cos(phi_rot);
    y_vals(i) = sin(theta_rot) * sin(phi_rot);
    z_vals(i) = cos(theta_rot);
end

% Plot points in 3D space
figure;
scatter3(x_vals, y_vals, z_vals, 10, 'filled');
title('Random Points in 3D Space (Oriented Spherical Cap)');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
grid on;

% Plot points in 2D space (phi vs. theta)
figure;
scatter(phi_vals, theta_vals, 10, 'filled');
title('Random Points in 2D Space (Phi vs. Theta)');
xlabel('Phi');
ylabel('Theta');
grid on;
