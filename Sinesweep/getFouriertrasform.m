function [] = getFouriertrasform(pVolt, ptime, pNamefile)

fprintf('Compute Fourier trasform for %s:\n', 'slow')
% compute the Fourier trasform of the signal
Y = fft(pVolt);
Fs = (1/ptime(2)).';      % Sampling frequency
T = ptime(2);             % Sampling period
L = length(ptime);        % Length of signal
fprintf('\tSampling frequency: %i \n\tLength of signal: %i\n', Fs, L);

% Compute the two-sided spectrum P2
P2 = abs(Y .* T);          % compute magnitude
P1 = P2(1:L/2+1);

% Define the frequency domain f and plot the single-sided amplitude
% spectrum P1
P1(2:end-1) = P1(2:end-1);
f = Fs * (0:(L/2))/L;

% plot the frequencies
figure();
plot(f,P1);
xlim([0 20])
grid on;
namefile = strcat('Fouriertrasform',pNamefile);
saveas(gcf, namefile, 'epsc')
end