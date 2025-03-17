neutrino_spectrum_data = readmatrix('ngc1068_spectrum_95.txt', 'Delimiter', '\t', 'NumHeaderLines', 1);

energies = neutrino_spectrum_data(:,1) * 10^9; % Energy [eV]

true_bin_sizes = diff(energies);
E_mid = (energies(1:end-1) + energies(2:end)) / 2;
p = polyfit(E_mid, true_bin_sizes, 1)

plot(E_mid, true_bin_sizes, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'True Bin Size'); hold on;
plot(E_mid, polyval(p, E_mid), 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Fit: ΔE = %.3fE + %.3f', p(1), p(2)));
xlabel('Energy (GeV)'); ylabel('Bin Width (eV)');
title('Linear Fit of Bin Width vs Energy');
legend('Location', 'best'); grid on;

% Display equation
fprintf('Fitted equation: ΔE = %.6f E + %.6f\n', p(1), p(2));