clc
clear 
close all
%%%% Lei(Raymond) Chi signal processing ps3 Q3

%% question b
x = [1:7];
h = [3,2,1,4];
y = downsamp(x , 2);
figure
stem(y)
p = upsamp(x , 2);
figure
stem(p)

%% question c
g = upsamp(h , 2)
A = downsamp(conv(g , x) , 2)
B = conv(h , downsamp(x , 2))
%% question d
max_abs_diff = max(abs(A(1:end) - B(1:end)))
%since the maximun abs difference is zero means array a and b are the same
%% question e
C = upsamp(conv(h , x) , 2)
D = conv(g, upsamp(x , 2))
% i dont know why but im getting one extra 0 at the last column

%% question a
function y = downsamp(x, m)
    y = x(1:m:end);
end

function y = upsamp(x, m)
    N = length(x);
    y = zeros(1, N*m);
    y(1:m:N*m) = x;
end 