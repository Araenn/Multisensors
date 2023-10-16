clc; clear; close all


f1 = 0.0020;

n = 1000;
t = (0:n-1)';
sigma = 0.05;

% sources
s1 = sin(2 * pi * f1 * t);
s2 = 5 * rand(n, 1);

% centrees reduites
s1 = (s1 - mean(s1));
s1 = s1 / std(s1);
s2 = (s2 - mean(s2));
s2 = s2 / std(s2);

% matrice de melange
m = [1 0.9; 0.95 1; 1 0.7; 0.2 1; 1 0.6];

% observations
bruit = sigma * randn(5, n);
s = [s1 s2]';
y = m * s + bruit;

% decorrelation
Ryb = cov(y');
[Qyb, Dyb] = eig(Ryb);
Dys = Dyb(4:5, 4:5);
Qys = Qyb(1:5, 4:5);

B = Dys^(-1/2) * Qys';

x = B * y;
Rx = cov(x');

% separation

tau = 0.9;
N = size(s, 1);
for i = 1:N
    for j = 1:N
        for k = 1:N
            for l = 1:N
                cum(i, j, k, l) = mean(x(i, :).*x(j,:).*x(k,:).*x(l,:)) - mean(x(i,:).*x(j,:))*mean(x(k,:).*x(l,:)) - mean(x(i,:).*x(k,:))*mean(x(j,:).*x(l,:)) - mean(x(i,:).*x(l,:))*mean(x(j,:).*x(k,:));
            end
        end
    end
end
Cx = cum(:, :, 1, 1);

choix = menu("Choisir entre corrélation et cumulant : ", "Corrélation", "Cumulant");
if choix == 1
    Tx = Rx - tau;
    T1 = Rx;
else
    Tx = Cx - tau;
    T1 = Cx;
end
% si vp de Tx distinctes 2 à 2

Tx = (Tx + Tx') / 2;
[U, D] = eig(Tx);
%z = U' * x;

% sinon

T2 = Tx;
[V, L] = eig(T2, T1);
z = V' * x;
DP = V' * B * m

figure(1)
subplot(2,1,1)
plot(t, s1)
grid()
title("Sinus")
subplot(2,1,2)
plot(t, s2)
grid()
title("Rand")


figure(2)
subplot(2,1,1)
plot(t, y(1,:))
grid()
title("y1")
subplot(2,1,2)
plot(t, y(2,:))
grid()
title("y2")

figure(3)
subplot(2,1,1)
plot(t, x(1,:))
grid()
title("x1")
subplot(2,1,2)
plot(t, x(2,:))
grid()
title("x2")

figure(4)
subplot(2,1,1)
plot(t, z(2,:))
grid()
title("z1")
subplot(2,1,2)
plot(t, z(1,:))
grid()
title("z2")