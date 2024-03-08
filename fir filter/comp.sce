load('C:\Users\kanak\Documents\self\studies\IITB_EEDDCSP\sem6\EE338\filter_design\fir\test.dat');
a = 0:0.001:1
b = 0:0.0001:1
fs = [600, 425]
pi = %pi
for i = 1:2
    scf(i);clf();
    xtitle('Magnitude Response Comparison', 'Un-normalized frequency (in kHz)', 'Magnitude');
    xgrid(color('grey'))
    plot(a*fs(i)/2, iirdata(2*i-1, :), 'b');
    plot(b*fs(i)/2, firdata(2*i-1, :), 'r');
    legend([['Butterworth ', 'Chebyschev '](i)+'IIR filter';'Kaiser window FIR filter']);
end
for i = 3:4
    scf(i);clf();
    xtitle('Phase Response Comparison', 'Un-normalized frequency (in kHz)', 'Phase (in radians)');
    xgrid(color('grey'))
    plot(a*fs(i-2)/2, iirdata(i - modulo(i, 2), :)*[1, 4](i-2), 'b');
    plot(b*fs(i-2)/2, firdata(i - modulo(i, 2), :), 'r');
    legend([['Butterworth ', '4x scaled Chebyschev '](i-2)+'IIR filter';'Kaiser window FIR filter']);
end
