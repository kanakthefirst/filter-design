i = %i
pi = %pi
s = %s
z = %z
N = 4;
fs = 600;
ze = [1.0639*i, 1.7690*i];
po = [-0.0310 + 0.9995*i, -0.3536 + 0.7071*i];

Hlpf = 0.85;
for a=1:2
    Hlpf = Hlpf*(1-s*2*real(1/ze(a))+s^2*abs(1/ze(a))^2)/(1-s*2*real(1/po(a))+s^2*abs(1/po(a)^2));
end
disp(Hlpf);

omega = 0:0.0001:2;
eval = horner(Hlpf, i*omega);

scf(0);clf(0);
xtitle('Magnitude Response of Elliptic Lowpass Filter', '$\Omega_L$','$|H_{lpf}(j\Omega_L)|$');
xgrid(color('grey'));
plot(omega, abs(eval));
plot([0, 2], [0.85, 0.85], 'm--');
plot([0, 2], [0.15, 0.15], 'm--');

scf(1);clf(1);
xtitle('Phase Response of Elliptic Lowpass Filter', '$\Omega_L$','Phase (in radians)');
xgrid(color('grey'));
plot(omega, atan(imag(eval), real(eval)));

Hbpf = horner(Hlpf, (s^2+0.7524175)/(0.7258751*s));
disp(Hbpf);
H = horner(Hbpf, (z-1)/(z+1));
disp(H);

omega = 0:0.0001:1;
eval = horner(H, exp(i*pi*omega));

scf(2);clf(2);
xtitle('Magnitude Response', 'Un-normalized Frequency (in kHz)','Magnitude');
xgrid(color('grey'));
plot(omega*fs/2, abs(eval));
plot([0, fs/2], [0.85, 0.85], 'm--');
plot([0, fs/2], [0.15, 0.15], 'm--');

scf(3);clf(3);
xtitle('Phase Response', 'Un-normalized Frequency (in kHz)','Phase (in radians)');
xgrid(color('grey'));
plot(omega*fs/2, atan(imag(eval), real(eval)));
