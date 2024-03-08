i = %i
pi = %pi
s = %s
z = %z
N = 3;
fs = 425;
ze = 1.2604*i;
po = -0.1153 + 0.9936*i;
p0 = -0.6232;

num = 1-s*2*real(1/ze)+s^2*abs(1/ze)^2;
den = (1-s*2*real(1/po)+s^2*abs(1/po)^2)*(1-s*real(1/p0));
Hlpf = num/den;
disp(Hlpf);

omega = 0:0.0001:2;
eval = horner(Hlpf, i*omega);

scf(0);clf(0);
xtitle('Magnitude Response of Elliptic Lowpass Filter', '$\Omega_L$','$|H_{lpf}(j\Omega_L)|$');
xgrid(color('grey'));
plotframe([0, 0, 1, 1]);
plot(omega, abs(eval));
plot([0, 2], [0.85, 0.85], 'm--');
plot([0, 2], [0.15, 0.15], 'm--');

scf(1);clf(1);
xtitle('Phase Response of Elliptic Lowpass Filter', '$\Omega_L$','Phase (in radians)');
xgrid(color('grey'));
plot(omega, atan(imag(eval), real(eval)));

Hbsf = horner(Hlpf, (0.8237552*s)/(s^2+1.1263572));
disp(Hbsf);
H = horner(Hbsf, (z-1)/(z+1));
disp(H);

omega = 0:0.0001:1;
eval = horner(H, exp(i*pi*omega));

scf(2);clf(2);
xtitle('Magnitude Response', 'Un-normalized Frequency (in kHz)','Magnitude');
xgrid(color('grey'));
plotframe([0, 0, 1, 1]);
plot(omega*fs/2, abs(eval));
plot([0, fs/2], [0.85, 0.85], 'm--');
plot([0, fs/2], [0.15, 0.15], 'm--');

scf(3);clf(3);
xtitle('Phase Response', 'Un-normalized Frequency (in kHz)','Phase (in radians)');
xgrid(color('grey'));
plot(omega*fs/2, atan(imag(eval), real(eval)));
