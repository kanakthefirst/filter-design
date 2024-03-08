function y = eval(x)
    y = horner(H_fir, exp(i*x))
endfunction

pi = %pi;
i = %i;
z = %z;

f = [95, 100, 175, 180];
fs = 600;
w = f*2*pi/fs; //[0.9948377, 1.0471976, 1.8325957, 1.8849556]
wc1 = (w(2) + w(1))/2; //1.0210176
wc2 = (w(4) + w(3))/2; //1.8587757

delta = 0.15;
A = -20*log10(delta); //16.478175

if A > 50 then;
    alpha = 0.1102*(A-8.7);
elseif A >= 21 then
    alpha = 0.5842*(A-21)^(0.4) + 0.07886*(A-21);
else
    alpha = 0;
end
deltawt = w(2)-w(1); //0.0523599
len = 1 + (A - 8)/(2.285*deltawt); //71.862675
N = ceil((len - 1)/2); //36

b = alpha/N; //0
n = -N:N;
h_ideal = sinc(wc2*n)*wc2/pi - sinc(wc1*n)*wc1/pi;
h_fir = h_ideal .* window('kr', 2*N+1, b);
H_fir = poly(h_fir, 'z', 'coeff');
while abs(eval([-w(1), -w(4)])) > delta || abs(eval([-w(2), -w(3)])) < 1 - delta
    N = N + 1;
    n = -N:N;
    h_ideal = sinc(wc2*n)*wc2/pi - sinc(wc1*n)*wc1/pi;
    h_fir = h_ideal .* window('kr', 2*N+1, b);
    H_fir = poly(h_fir, 'z', 'coeff');
end
disp(N) //46
disp(h_fir)
omega = 0:0.0001:1;
data = eval(pi*omega);
[phase, mag] = phasemag(data);

scf(0);clf();xgrid(color('grey'));
xtitle('Magnitude of FIR Bandpass Filter', 'Un-normalized frequency (in kHz)', 'Magnitude');
plot(omega*fs/2, 10^(mag/20));
plot(f, abs(eval(-w)), 'rx');
plot([0, 300], delta, 'm--');
plot([0, 300], 1-delta, 'm--');

scf(1);clf();xgrid(color('grey'));
xtitle('Phase of FIR Bandpass Filter', 'Un-normalized frequency (in kHz)', 'Phase (in radians)');
plot(omega*fs/2, -phase/180*pi);
plot([f(2), f(2)], [0, -45], 'm--', [f(3), f(3)], [0, -45], 'm--')

scf(2);clf();xgrid(color('grey'));
xtitle('Time Domain Sequence', 'Time index n', 'h_fir[n]');
plot2d3(h_fir);
