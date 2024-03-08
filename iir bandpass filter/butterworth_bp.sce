//functions
function y=q(m)
    y = floor((m-1)/10)
endfunction
function y=r(m)
    y = m - 10*q(m)
endfunction
function y=BL(m)
    y = 10 + 5*q(m) + 13*r(m)
endfunction
function y=BH(m)
    y = BL(m) + 75
endfunction
function y = normalized(x)
    y = 2*%pi*x/fs
endfunction
function y = bilinear(x)
    y = tan(x/2)
endfunction
function y=freq_transform(w)
    y = (w^2 - w0square) / (b*w)
endfunction
function y=kth_root(k)
    y = wc*exp(%i*%pi*(0.5 + 0.5/n + k/n))
endfunction
function m=mag_2(r)
    m = 1 / (1 + (freq_transform(bilinear(r))/wc)^(2*n))
endfunction

//un-normalized discrete time specs for bandpass
//monotonic passband => butterworth approximation
m = 55;
discrete_time_wp1 = BL(m);//100
discrete_time_ws1 = discrete_time_wp1 - 5;//95
discrete_time_wp2 = BH(m);//175
discrete_time_ws2 = discrete_time_wp2 + 5;//180

//normalized digital specs for bandpass
fs = 600;
digital_ws1 = normalized(discrete_time_ws1);//0.9948377
digital_wp1 = normalized(discrete_time_wp1);//1.0471976
digital_wp2 = normalized(discrete_time_wp2);//1.8325957
digital_ws2 = normalized(discrete_time_ws2);//1.8849556

//analog specs for bandpass
analog_ws1 = bilinear(digital_ws1);//0.5429557
analog_wp1 = bilinear(digital_wp1);//0.5773503
analog_wp2 = bilinear(digital_wp2);//1.3032254
analog_ws2 = bilinear(digital_ws2);//1.3763819

//parameters for frequency transformation
w0square = analog_wp2*analog_wp1;//0.7524175 = (0.8674200)^2
b = analog_wp2 - analog_wp1;//0.7258751

//analog specs for lowpass
ws1L = freq_transform(analog_ws1);//-1.1611157
wp1L = freq_transform(analog_wp1);//-1.0000000
wp2L = freq_transform(analog_wp2);//1.0000000
ws2L = freq_transform(analog_ws2);//1.1430597

wp = 1;
ws = min(abs(ws1L), abs(ws2L));//1.1430597

//parameters for butterworth approximation
delta = 0.15;
d1 = 1 / (1 - delta)^2 - 1;//0.3840830
d2 = 1 / delta^2 - 1;//43.444444
n = ceil(log(d2/d1)/(2*log(ws/wp)));//= ceil(17.681654) = 18

t = 0.5;
wclo = wp / d1^(1/(2*n));//1.0269369
wchi = ws / d2^(1/(2*n));//1.0293682
wc = t*wclo + (1-t)*wchi;//1.0281525

//denominator of the analog lowpass transfer function
s = %s;
den = poly(1, "s", "coeff");

//a figure
mag_res = scf(2);
clf(2);
title('LHP Poles')
xlabel('$\mathbb{R}$');
ylabel('$j\mathbb{R}$');
xgrid(color('grey'));

//lhp poles
for k = 0:n-1
    den = den * (s - kth_root(k));
    plot(real(kth_root(k)), imag(kth_root(k)), 'rx');
end

angle = 0:0.001:2*%pi;
plot(wc*cos(angle), wc*sin(angle));

//analog lowpass transfer function
h_analog_lowpass = wc^n / real(den);
disp(h_analog_lowpass);
//analog bandpass transfer function
h_analog_bandpass = horner(h_analog_lowpass, (s^2 + w0square) / (b*s));
disp(h_analog_bandpass);
//discrete time bandpass transfer function
z = %z;
numfin = horner(h_analog_bandpass.num, (z-1)/(z+1)).num;
denfin = horner(h_analog_bandpass.den, (z-1)/(z+1)).num;
c = (1+z)^n;
h_discrete_time_bandpass = c*numfin/denfin;
disp(h_discrete_time_bandpass);
//a figure
mag_res = scf(0);
clf(0);
title('Magnitude Response')
xlabel('$\omega/\pi$');
ylabel('Magnitude');
xgrid(color('grey'));

//magnitude
for r = 0:0.001:1;
    mag(1000*r+1) = mag_2(%pi*r)
    eval(1000*r+1) = horner(h_analog_lowpass, %i*freq_transform(bilinear(%pi*r)))
end

plot(0:0.001:1, sqrt(mag));

f = [digital_ws1, digital_wp1, digital_wp2, digital_ws2]
for fr = f
    plot(fr/%pi, sqrt(mag_2(fr)), 'rx')
end

//a figure
phi_res = scf(1);
clf(1);
title('Phase Response')
xlabel('$\omega/\pi$');
ylabel('Phase (in radians)');
xgrid(color('grey'));

//phase
phi = atan(imag(eval), real(eval));

for i = 2:1001
    j = 0;
    if phi(i)-phi(i-1)>0 then
        j = j + 1;
        phi(i:1001) = phi(i:1001) - 2*%pi*j;
    end
end

plot(0:0.001:1, phi);
