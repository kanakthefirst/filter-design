fs1 = 95;
fp1 = 100;
fp2 = 175;
fs2 = 180;

fs = 600;

ws1 = tan(pi*fs1/fs);
wp1 = tan(pi*fp1/fs);
wp2 = tan(pi*fp2/fs);
ws2 = tan(pi*fs2/fs);

w0_square = wp1*wp2;
b = wp2-wp1;

ws1lpf = (ws1^2 - w0_square)/(b*ws1);
ws2lpf = (ws2^2 - w0_square)/(b*ws2);

wp = 1;
ws = min(abs(ws1lpf), abs(ws2lpf));

delta = 0.15;
ep = sqrt(-1 + 1/(1-delta)^2);
es = sqrt(-1 + 1/delta^2);

k = wp/ws;% 0.8748
k1 = ep/es;% 0.0940

[K, K_] = ellipk(k);% [2.1850, 1.6775]
[K1, K1_] = ellipk(k1);% [1.5743, 3.7566]

N = (K1_/K1)/(K_/K);% 3.1080
N = ceil(N);% 4

k = ellipdeg(N, k1);% 0.9595
wsnew = 1/k;% 1.0422

L = floor(N/2);
r = mod(N,2);
i = (1:L)';

u = (2*i-1)/N; zeta_i = cde(u,k);
za = wp * 1i./(k*zeta_i) % filter zeros

v0 = -1i*asne(1i/ep, k1)/N;
pa = wp * 1i*cde(u-1i*v0, k) % filter poles
pa0 = wp * 1i*sne(1i*v0, k);
