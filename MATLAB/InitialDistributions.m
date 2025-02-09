% Parameters
Rmax = 30;              % Maximum radius of the sphere
n = 10000;              % Number of points per distribution
alphas = [0.01, 0.1, 0.5, 1, 3, 5, 10];     % Different alpha values to explore
figure;
for i = 1:length(alphas)
    alpha = alphas(i);
    
    % Generate random numbers for the radial distribution
    u = rand(n, 1);
    % Inverse transform sampling: r = Rmax * u^(1/alpha)
    r = Rmax * u.^(1/alpha);
    
    % Generate uniformly distributed directions on the sphere.
    % One common method:
    theta = 2 * pi * rand(n, 1);       % azimuthal angle in [0, 2*pi]
    % For the polar angle, using the transformation that gives uniform points:
    phi = acos(2 * rand(n, 1) - 1);      % polar angle in [0, pi]
    
    % Convert spherical coordinates to Cartesian coordinates
    x = r .* sin(phi) .* cos(theta);
    y = r .* sin(phi) .* sin(theta);
    z = r .* cos(phi);
    
    % Plot the points in a subplot
    subplot(2, 4, i);
    scatter3(x, y, z, 1, r, 'filled', ''); % color points by z for visual depth
   
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(['\alpha = ' num2str(alpha)]);
end
sgtitle('Effect of \alpha on the 3D Distribution of Points');
