fp1 = 85;
fs1 = 90;
fs2 = 130;
fp2 = 135;

fs = 425;

wp1 = tan(pi*fp1/fs);
ws1 = tan(pi*fs1/fs);
ws2 = tan(pi*fs2/fs);
wp2 = tan(pi*fp2/fs);

w0_square = wp1*wp2;
b = wp2-wp1;

ws1lpf = (b*ws1)/(w0_square-ws1^2);
ws2lpf = (b*ws2)/(w0_square-ws2^2);

wp = 1;
ws = min(abs(ws1lpf), abs(ws2lpf));

delta = 0.15;
ep = sqrt(-1 + 1/(1-delta)^2);
es = sqrt(-1 + 1/delta^2);

k = wp/ws;% 0.7895
k1 = ep/es;% 0.0940

[K, K_] = ellipk(k);% [1.9770, 1.7609]
[K1, K1_] = ellipk(k1);% [1.5743, 3.7566]

N = (K1_/K1)/(K_/K);% 2.6791
N = ceil(N);% 3

k = ellipdeg(N, k1);% 0.8571
wsnew = 1/k;% 1.1667

L = floor(N/2);
r = mod(N,2);
i = (1:L)';

u = (2*i-1)/N; zeta_i = cde(u,k);
za = wp * 1i./(k*zeta_i) % filter zeros

v0 = -1i*asne(1i/ep, k1)/N;
pa = wp * 1i*cde(u-1i*v0, k) % filter poles
pa0 = wp * 1i*sne(1i*v0, k)
