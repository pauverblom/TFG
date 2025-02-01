maxpoints = 10000;
position = zeros(3, maxpoints);

for i=1:maxpoints
    theta = acos(2 * rand() - 1);  % Random inclination
    phi = 2 * pi * rand();  % Random azimuth
    r = rand()^(1/3); % 1/3 to evenly distribute them inside a sphere
    position(1, i) = r*sin(theta)*cos(phi);
    position(2, i) = r*sin(theta)*sin(phi);
    position(3, i) = r*cos(theta);
end

plot3(position(1,:), position(2,:), position(3,:), '.')