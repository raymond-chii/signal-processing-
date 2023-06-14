clc
clear
close all
%%%% Lei(Raymond) Chi signals ps08

%% q1

% Set the parameters
a = 0.544;
v = 5;
N = 1e6;

X = randn(N, 1);
t = trnd(v, N, 1);
U = rand(N, 1); 
cauchy = a*tan(pi*U); 

frac = sum(abs(X) < 1) / N;
frac_cauchy = sum(abs(cauchy) < 1) / N;
frac_t = sum(abs(t) < 1) / N;

figure;
histogram(X);
hold on;
xline(1, '--k');
xline(-1, '--k');
title(['Normal: ' num2str(frac*100) '% below 1']);

figure;
histogram(t);
hold on;
xline(1, '--k');
xline(-1, '--k');
title(['t(5): ' num2str(frac_cauchy*100) '% below 1']);

figure;
histogram(cauchy);
hold on;
xline(1, '--k');
xline(-1, '--k');
title(['Cauchy: ' num2str(frac_t*100) '% below 1']);

X_segments = reshape(X, [100000, 10]);
t_segments = reshape(t, [100000, 10]);
cauchy_segments = reshape(cauchy, [100000, 10]);

X_means = mean(X_segments)
t_means = mean(t_segments)
cauchy_means = mean(cauchy_segments)

%% q2 a)6

num = [1 4 0.2];
den = [1 -1.6 0.81];
[z, p, k] = tf2zpk(num, den);

all(abs(p) < 1)
all(abs(z) < 1)

figure;
zplane(z, p);
title('Pole-Zero Plot of H(z)');


%% q2 b)

% Constants
N = 10e4;
sigma = sqrt(2);

% 1
v = sigma*randn(N, 1);
x = filter(num, den, v);
r_x = zeros(1,7); 

% 2

for m = 0:6
    for n = 1:(N-m)
        r_x(m+1) = r_x(m+1) + x(n)*x(n+m);
    end
    r_x(m+1) = r_x(m+1)/(N-m);
end


% 3
figure;
stem(-6:6, [fliplr(r_x(2:end)), r_x])
xlabel('Lag (m)');
ylabel('Correlation');

% 4
R = toeplitz(r_x)

% 5
eigvals = eig(R)
% R is positive definite.

% 6
R2 = (1/(N-6))*R*R';
err = norm(R - R2)


%% c)

% 1
[s_est, w] = pwelch(x, hamming(512), 256, 512);
plot(w, s_est);
xlabel('Normalized Digital Radian Frequency');
ylabel('Power Spectral Density');
title('Welch Periodogram Method');

% 2
[~, idx] = max(s_est);
w0 = w(idx) % The frequency where the PSD has a peak is

% 3
b = [1, 4, 0.2];
a = [1, -1.6, 0.81];
[r, p, k] = residuez(b, a)
theta_p = angle(p(1))
disp(['The pole angle is approximately ', num2str(theta_p), ' radians.']);
disp(['The frequency where the PSD has a peak is approximately ', num2str(w0), ' radians.']);

%% d)

[a, varv] = aryule(x, 4);
v0 = randn(size(x));
x0 = filter(1, a, v0);
figure;
subplot(2,1,1);
stem(1:100, x0(1:100));
title('x0');
subplot(2,1,2);
stem(1:100, x(1:100));
title('x');






