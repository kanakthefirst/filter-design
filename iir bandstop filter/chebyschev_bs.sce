//functions
function y=q(m)
    y = floor((m-1)/10)
endfunction
function y=r(m)
    y = m - 10*q(m)
endfunction
function y=BL(m)
    y = 20 + 3*q(m) + 11*r(m)
endfunction
function y=BH(m)
    y = BL(m) + 40
endfunction
function y = normalized(x)
    y = 2*%pi*x/fs
endfunction
function y = bilinear(x)
    y = tan(x/2)
endfunction
function y=freq_transform(w)
    y = (b*w)/(w0square - w^2)
endfunction
function y=kth_root(k)
    y = %i*cos(%pi*(0.5/n + k/n)+%i*bk)
endfunction

//un-normalized discrete time specs for bandstop
//equiripple passband => chebyschev approximation
m = 55;
discrete_time_ws1 = BL(m);//90
discrete_time_wp1 = discrete_time_ws1 - 5;//85
discrete_time_ws2 = BH(m);//130
discrete_time_wp2 = discrete_time_ws2 + 5;//135

//normalized digital specs for bandstop
fs = 425;
digital_wp1 = normalized(discrete_time_wp1);//1.2566371
digital_ws1 = normalized(discrete_time_ws1);//1.3305569
digital_ws2 = normalized(discrete_time_ws2);//1.9219155
digital_wp2 = normalized(discrete_time_wp2);//1.9958353

//analog specs for bandstop
analog_wp1 = bilinear(digital_wp1);//0.7265425
analog_ws1 = bilinear(digital_ws1);//0.7845976
analog_ws2 = bilinear(digital_ws2);//1.4312732
analog_wp2 = bilinear(digital_wp2);//1.5502977

//parameters for frequency transformation
w0square = analog_wp2*analog_wp1;//1.1263572 = (1.0612998)^2
b = analog_wp2 - analog_wp1;//0.8237552

//analog specs for lowpass
ws1L = freq_transform(analog_ws1);//1.2653920
wp1L = freq_transform(analog_wp1);//1.0000000
wp2L = freq_transform(analog_wp2);//-1.0000000
ws2L = freq_transform(analog_ws2);//-1.2785044

wp = 1;
ws = min(abs(ws1L), abs(ws2L));//1.2653920

//parameters for chebyschev approximation
delta = 0.15;
d1 = 1 / (1 - delta)^2 - 1;//0.3840830
d2 = 1 / delta^2 - 1;//43.444444
e = sqrt(d1);//0.6197443
n = ceil(acosh(sqrt(d2/d1))/(acosh(ws/wp)));//=ceil(4.2829034) = 5
bk = asinh(1/e)/n;//0.2512306

//denominator of the analog lowpass transfer function
s = %s;
den = poly(1, "s", "coeff");

//a figure
poles = scf(2);
clf(2);
title('LHP Poles')
xlabel('$\mathbb{R}$');
ylabel('$j\mathbb{R}$');
xgrid(color('grey'));
a=get("current_axes");
a.data_bounds=[-1.2,-1.2;1.2,1.2];

//lhp poles
for k = n:2*n-1
    den = den * (s - kth_root(k));
    plot(real(kth_root(k)), imag(kth_root(k)), 'rx');
end
den = real(den);
angle = 0:0.001:2*%pi;
plot(sinh(bk)*cos(angle), cosh(bk)*sin(angle));

//analog lowpass transfer function
h_analog_lowpass = horner(den, 0) / den;
//analog bandstop transfer function
h_analog_bandstop = horner(h_analog_lowpass, (b*s) / (s^2 + w0square));
//discrete time bandstop transfer function
z = %z;
h_discrete_time_bandstop = horner(h_analog_bandstop, (z - 1)/(z + 1))

//a figure
mag_res = scf(0);
clf(0);
title('Magnitude Response')
xlabel('Un-normalized frequency (kHz)');
ylabel('Magnitude');
xgrid(color('grey'));

//magnitude
for omega = 0:0.001:1;
    eval(1000*omega+1) = horner(h_discrete_time_bandstop, exp(%i*%pi*omega));
end

xrange = (0:0.001:1)*fs/2
plot(xrange, abs(eval));

plot(xrange, .85, 'black')

for f = [digital_wp1, digital_ws1, digital_ws2, digital_wp2]
    plot(f*fs/(2*%pi), horner(h_discrete_time_bandstop, exp(%i*f)), 'rx');
end

//a figure
phi_res = scf(1);
clf(1);
title('Phase Response')
xlabel('Un-normalized frequency (kHz)');
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

plot(xrange, phi);
