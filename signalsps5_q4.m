%%question 4
C = 10e-9;
R = 1e3; 
N = 1e4;
f = linspace(0, 1e6, N); 
s = 2*pi*1i*f; 
Kmax = 398.1;
K0 = 2;
K1 = K0;
K2 = 0.5*K0 + 0.5*Kmax;
K3 = 0.2*K0 + 0.8*Kmax;

% transfer function
H1 = (K1/(R*C)*s)./(s.^2 + ((4-K1)/(R*C))*s + 2/(R^2*C^2));
H2 = (K2/(R*C)*s)./(s.^2 + ((4-K2)/(R*C))*s + 2/(R^2*C^2));
H3 = (K3/(R*C)*s)./(s.^2 + ((4-K3)/(R*C))*s + 2/(R^2*C^2));

% Plot magnitude and phase responses for each case
figure;
subplot(211);
plot(f/1e3, 20*log10(abs(H1)));
ylim([min(20*log10(abs(H1(f>=1e6)))), max(20*log10(abs(H1))) + 5]);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response for K1');
subplot(212);
plot(f/1e3, unwrap(angle(H1))*180/pi);
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
title('Phase Response for K1');

figure;
subplot(211);
plot(f/1e3, 20*log10(abs(H2)));
ylim([min(20*log10(abs(H2(f>=1e6)))), max(20*log10(abs(H2))) + 5]);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response for K2');
subplot(212);
plot(f/1e3, unwrap(angle(H2))*180/pi);
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
title('Phase Response for K2');

figure;
subplot(211);
plot(f/1e3, 20*log10(abs(H3)));
ylim([min(20*log10(abs(H3(f>=1e6)))), max(20*log10(abs(H3))) + 5]);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response for K3');
subplot(212);
plot(f/1e3, unwrap(angle(H3))*180/pi);
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
title('Phase Response for K3');

%%Case 1: K1 = K0
f = logspace(3,6,N);
s = 2*pi*1i*f;
K1 = K0;
H1 = (K1/(R*C)*s) ./ (s.^2 + (4-K1)/(R*C)*s + 2/(R^2*C^2));
figure;
subplot(2,1,1);
semilogx(f/1e3, 20*log10(abs(H1))); 
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Bode Plot for Case 1: K1 = K0');
grid on;
subplot(2,1,2);
semilogx(f/1e3, unwrap(angle(H1))*180/pi);
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
grid on;

%%Case 2: K2 = 0.5K0 + 0.5Kmax
K2 = 0.5*K0 + 0.5*Kmax;
H2 = (K2/(R*C)*s) ./ (s.^2 + (4-K2)/(R*C)*s + 2/(R^2*C^2));
figure;
subplot(2,1,1);
semilogx(f/1e3, 20*log10(abs(H2))); 
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Bode Plot for Case 2: K2 = 0.5K0 + 0.5Kmax');
grid on;
subplot(2,1,2);
semilogx(f/1e3, unwrap(angle(H2))*180/pi);
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
grid on;

%%Case 3: K3 = 0.2K0 + 0.8Kmax
K3 = 0.2*K0 + 0.8*Kmax;
H3 = (K3/(R*C)*s) ./ (s.^2 + (4-K3)/(R*C)*s + 2/(R^2*C^2));
figure;
subplot(2,1,1);
semilogx(f/1e3, 20*log10(abs(H3))); 
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Bode Plot for Case 3: K3 = 0.2K0 + 0.8Kmax');
grid on;
subplot(2,1,2);
semilogx(f/1e3, unwrap(angle(H3))*180/pi);
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
grid on;

