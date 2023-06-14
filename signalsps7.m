clc
clear
close all
%%%% Lei(Raymond) Chi signals ps07
% The Experiment 
%% Part I: Signal Generation 
M = 100;
K = 20;
N = 200;
PdB = [0, -2, -4];
PndB = 10;
var_ = 10.^(PdB./10);% var power 
var_noise = 10.^(PndB./10);
L = length(PdB);
S = L_vectors(M, K, L); % matrix S

A = data_matrix(L, M, N, S, var_, var_noise); % matrix A
%% Part II: Analysis
R = (1/N)*A*(A');

[eigvec,eigval0]= eig(R);
[eigval,idx]= sort(diag(eigval0),'descend'); 
eigvec= eigvec(:,idx);

[U, sval, V ] = svd(A);
sval = diag(sval);

Ul = U(:, 1:L);
figure;
stem(sval, 'filled');
xlabel('s_value Index');
ylabel('s_value');
title('s_values of A');

figure;
stem(eigval, 'filled');
xlabel('eig_val Index');
ylabel('eig_val');
title('Sorted eig_val of R');

sig_ratio = sval(3) / sval(4) % ratio of the third largest to the fourth 
% largest singular value
wavelength_ratio = eigval(3) / eigval(4) % ratio of the third largest to 
% the fourth largest eigenvalues

Ps = Ul*(Ul.');
Pn = eye(100)-Ps;
Rin = R^(-1);

Smusic = zeros(1,3);
for i = 1:L
    Smusic(i) = 1/((S(:,i).')*Pn*S(:,i));
end

Smvdr = zeros(1,3);
for j = 1:L
    Smvdr(j) = (S(:,j).')*Rin*S(:,j);
end

test_L = 20;
test_S = [S L_vectors(M, K, test_L)];

music_test = zeros(1, test_L+L);
mvdr_test = zeros(1, test_L+L);

for i = 1:test_L+3
    music_test(i) = 1/((test_S(:,i).')*Pn*test_S(:,i)); 
end
for i = 1:test_L+3
mvdr_test(i) = (test_S(:,i).')*Rin*test_S(:,i);
end

max_music_test = max(music_test)
median_music_test = median(music_test)
mean_music_test = mean(music_test)

max_mvdr_test = max(mvdr_test)
median_mvdr_test = median(mvdr_test)
mean_mvdr_test = mean(mvdr_test)

% base on the result, MVDR seems to work better at identifying the correct
% source vectors, since the number seems more compact rather than scattered

%% Part III: Try Again 
M = 100;
N = 50;
S_ = L_vectors(M, K, L);
A = data_matrix(L, M, N, S_, var_, var_noise);
R = (1/N)*A*(A.');

[eigvec,eigval0]= eig(R);
[eigval,idx]= sort(diag(eigval0),'descend'); 
eigvec= eigvec(:,idx);

[U, sval, V] = svd(A);
sval = diag(sval);
Ul = U(:, 1:L);

figure;
stem(sval, 'filled');
xlabel('s_value Index');
ylabel('s_value');
title('s_values of A_part3');
sig_ratio = sval(3) / sval(4) 

Ps = Ul*(Ul.');
Pn = eye(100)-Ps;

Smusic_part3 = zeros(1,3);
for i = 1:L
    Smusic_part3(i) = 1/((S_(:,i).')*Pn*S_(:,i));
end

test_S_part3 = [S L_vectors(M, K,test_L)];

music_test_part3 = zeros(1, test_L+L);

for j = 1:test_L+L
    music_test_part3(j) = 1/((test_S_part3(:,j).')*Pn*test_S_part3(:,j));
end

Smusic
Smusic_part3

music_test
music_test_part3


max_music_test_part3 = max(music_test_part3)
median_music_test_part3 = median(music_test_part3)
mean_music_test_part3 = mean(music_test_part3)

%The first MUSIC algorithm has better performance than the part3 one based
% on the maximum, median, and mean values of the test results.

%% Part IV: One Last Thing 
display_the_matrix = S*(S.')

% The matrix S*(S^T) is a square matrix with size L x L, shows the 
% relationship between different components of the source vectors in terms 
% of their covariance. It conveys how the source vectors are correlated 
% with each other, which is important in this problem because it helps us 
% understand the structure of the signals and how they are related to each other.


%% fucntions 
function S = L_vectors(M, K, L)
    S = zeros(M, L);
    for i = 1:L
        perm = randperm(M, K);
        S(perm, i) = 1/sqrt(K);
    end
    S = S./ vecnorm(S); % scalar 
end

function A = data_matrix(L, M, N, S, var_, var_noise)
    B = randn(L, N); %B = L x N
    for j = 1:L
        B(j, :) = sqrt(var_(j)).*B(j, :);
    end
    V = randn(M, N); %V = M x N
    for j = 1:N
        V(:,j) = sqrt(var_noise)*V(:, j);
    end
    A = S*B+(1/sqrt(M))*V;
end