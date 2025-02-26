
a = logspace(-6,8,100); % Lower bounds for ambient photon energy
b = 10^8;    % Upper bound for ambient photon energy

% Perform numerical integrations
for j = 1:100
    results(j) = integral(@phi_over_epsilon_squared, a(j), b);
end



% Display result
loglog(a,1/(3*10^10)*results)
hold on;
loglog(a, integralfunc(a))

% Function definition (same as before)
function y = phi_over_epsilon_squared(epsilon)
    % Initialize output as a zero array of the same size as epsilon
    y = zeros(size(epsilon));

    for i = 1:numel(epsilon)
        e = epsilon(i); % Extract each element separately
        if all(10^(-6) < e & e < 4.3*10^-2)
            y(i) = 10^9.7 * e^(4 - 2 - 2);
        elseif all(4.3*10^-2 < e & e < 2.56*10^5)
            y(i) = 10^3.9 * e^(-0.25 - 2 - 2);
        elseif all(2.56*10^5 < e & e < 10^8)
            y(i) = 10^35 * e^(-6 - 2 - 2);
        else
            y(i) = 0; % Default value outside range
        end
    end
end

function y = integralfunc(epsilon_min)
     % Initialize output as a zero array of the same size as epsilon
    y = zeros(size(epsilon_min));

    for i = 1:numel(epsilon_min)
        e = epsilon_min(i); % Extract each element separately
        if all(10^(-6) < e & e < 4.3*10^-2)
            y(i) = 1/3*(10^(-0.3)*(4.3*10^(-2)-e)+6.75*10^(-3));
        elseif all(4.3*10^-2 < e & e < 2.56*10^5)
            y(i) = 1/3*(-10^(-6.1)/3.25*((2.56*10^5)^(-3.25)-e^(-3.25))+2.353*10^(-25));
        elseif all(2.56*10^5 < e & e < 10^8)
            y(i) = 1/3*(-10^25/9*((10^8)^(-9)-e^(-9)));
        else
            y(i) = 0; % Default value outside range
        end
    end
end