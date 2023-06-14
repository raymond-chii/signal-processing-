%%%%Raymond
%% (a) Set the filter design parameters
n = 6; % filter order
f_lo = 10e3; % low frequency cutoff [Hz]
f_hi = 12e3; % high frequency cutoff [Hz]
r_p = 2; % passband ripple [dB]
f_s = 40e3; % sampling rate [Hz]
f_Ny = f_s/2; % Nyquist bandwidth [Hz]

% Design the filter and obtain b;a vectors
f_crit = [f_lo/f_Ny, f_hi/f_Ny]; % critical frequencies
[z, p, k] = cheby1(n, r_p, f_crit); % Chebyshev Type I filter design
[b, a] = zp2tf(z, p, k); % zero-pole to transfer function conversion

% Plot the pole-zero plot
figure;
zplane(z, p);
title('Chebyshev Type I Filter');

%% (b) Generate the frequency response vector
f = linspace(0, f_Ny, 10e4); 
f_dig = 2*pi*f/f_s;
H = freqz(b, a, f_dig); 

%% (c) Plot the magnitude and phase response
figure;
subplot(211);
H_db = 20*log10(abs(H)); % magnitude in dB
plot(f/1e3, H_db);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response');
xlim([0 f_Ny/1e3]);
ylim([-40 2]);
grid on;

subplot(212);
H_ph = rad2deg(unwrap(angle(H))); % phase in degrees
plot(f/1e3, H_ph);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Phase Response');
xlim([0 f_Ny/1e3]);
grid on;

%% (d) Find the frequencies where the gain falls below -30dB
a = find(H_db <= -30);
lvl1 = f(a(1));
lvl2 = f(a(end)+1);