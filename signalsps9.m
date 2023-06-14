clc
clear
close all
%%%% Lei(Raymond) Chi signals ps09


%% question 3

% Part A
 
E = [3 1; 2 1];
A0 = [-0.3 0.4; -0.4 -0.3];
eigval1 = eig(A0)
A = E * A0 * inv(E);
eigval2 = eig(A)

% Part B

s=sym('s');

g = inv(s*eye(2)-A); 
E = ilaplace(g);
matfunc = matlabFunction(E); 

% Part C
xd0 = [2; 1];
fs = 10;
dt = 1/fs;

xd = zeros(2, 100);
xd(:, 1) = xd0;
for n = 2:100
    
    xd(:, n) = matfunc(n * dt) * xd0

end 

% Part D
xde = zeros(2, 100);
xdm = zeros(2, 100);
xde(:, 1) = xd0; 
xdm(:, 1) = xd0; 
Ade = eye(size(A)) + dt * A;
Adm = eye(size(A)) + dt * A + ((dt)^2/2) * A^2;
for n = 1:99
    xde(:, n+1) = Ade * xde(:, n)
    xdm(:, n+1) = Adm * xdm(:, n)
end

% Part E

figure;
hold on;
plot(xd(1, :), xd(2, :), 'g');
plot(xde(1, :), xde(2, :), 'r');
plot(xdm(1, :), xdm(2,:), 'b');
xlabel('x1');
ylabel('x2');
hold off;


% Part F

e_xde = max(abs(xde - xd));
e_xdm = max(abs(xd - xdm));

max(e_xdm)
max(e_xde)

% Part G

% the error value suggested that xde(Midpoint) is better than xdm 
% The approximate error does not converge to zero 
% as time approaches zero, which gives the impression that the Euler technique 
% generates a higher error rate at smaller numbers. 








