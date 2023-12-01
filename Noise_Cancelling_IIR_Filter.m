%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% IIR filter as a noise cancelling filter %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

[y, Fs] = audioread('Don_Giovanni_1.wav');

Ts = 1/Fs;                                                               %%Time axis step             
N = length(y);                                                           %%Length of the signal in samples
t = 0:Ts:(N-1)*Ts;                                                       %%Time axis
deltaF = Fs/N;                                                           %%Frequency axis step
f = 0:deltaF:(N-1)*deltaF;                                               %%Frequency axis

%%We can listen to the original signal:
%sound(0.05*y, Fs);

%% Temporal plot
figure(1)
subplot(2,1,1);
plot(t, y);
title('Signal in the time domain');
xlabel('Time (s)');
ylabel('Amplitude')
grid on;

%% Frequential plot

Y = fft(y);
Y = abs(Y);

subplot(2,1,2);
plot(f, 20*log10(Y));
title('Signal in the frequency domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;

%% We detect 2 noising frequencies : 1397.01 Hz and 5658.01 Hz

[M,I] = findpeaks(Y(1:N/2),f(1:N/2),'Threshold',10^4);                   %%We use findpeaks to detect the peaks of the noising frequencies

%% Normalized cutoff frequencies:
w1 = 2*pi*I(1)/Fs;      %%For frequency 1
w2 = 2*pi*I(2)/Fs;      %%For frequency 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filtering Frequency 1
%%Compute the zeros and the poles to construct the first notch filter
% Zeros
Z1f1 = exp(1i*w1);
Z2f1 = conj(Z1f1);

% Poles (for a notch filter, we have observed that the notch filters are at
% the same angles of the zeros). If the pole value is close to the zero
% value, the filter will be narrower. We will not exceed 1 for the
% stability of the filter. Here, we just multiply the value of the zero
% with a coefficient to get the value of the pole. The filter works with a
% multiplication coefficient of 0.9, but it is more precise with 0.99
P2 = 0.99*Z1f1;
P2f1 = conj(P2);

Z1 = [Z1f1, Z2f1];
Z1 = Z1';
P1 = [P2, P2f1];
P1 = P1';

[B1, A1] = zp2tf(Z1,P1,1);                                               %%Finding the coefficients of the filter

y_f1_filtered = filter(B1, A1, y);                                       %%Filtering the signal

%fvtool(B1,A1);                                                          %%Filter visualization on fvtool             

%%We can listen to the filtered signal
%sound(y_f1_filtered, Fs);

%%Plot to verify the deleting of the first frequency
Y_f1_filtered = fft(y_f1_filtered);
Y_f1_filtered = abs(Y_f1_filtered);

figure(2)
plot(f, 20*log10(Y_f1_filtered));
title('Signal filtered a first time in the frequency domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filtering Frequency 2
%%Compute the zeros and the poles to construct the second notch filter
%%Zeros
Z1f2 = exp(1i*w2);
Z2f2 = conj(Z1f2);

% Poles (for a notch filter, we have observed that the notch filters are at
% the same angles of the zeros). If the pole value is close to the zero
% value, the filter will be narrower. We will not exceed 1 for the
% stability of the filter. Here, we just multiply the value of the zero
% with a coefficient to get the value of the pole. The filter works with a
% coefficient of 0.9, but is more precise with 0.99
P1f2 = 0.99*Z1f2;
P2f2 = conj(P1f2);

Z2 = [Z1f2, Z2f2];
Z2 = Z2';
P2 = [P1f2, P2f2];
P2 = P2';

[B2, A2] = zp2tf(Z2,P2,1);                                                  %%Finds the coefficients of the filter

y_f2_filtered = filter(B2, A2, y_f1_filtered);                              %%Filtering the signal

%fvtool(B2,A2);                                                             %%Filter visualization on fvtool             

%%We can listen to the sound before the filtering:
sound(y, Fs);
%%And compare it to the sound after the filtering:
sound(y_f2_filtered, Fs);

%%Plot to verify the deleting of the second frequency 
Y_f2_filtered = fft(y_f2_filtered);
Y_f2_filtered = abs(Y_f2_filtered);

figure(4)
plot(f, 20*log10(Y_f2_filtered));
title('Signal filtered a second time in the frequency domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;

%%Final plot to compare the signal before and after the filtering
figure(5)
subplot(2,1,1);                                                             %%Time domain
plot(t, y);
hold on;
plot(t, y_f2_filtered);
title('Signal in the time domain before and after the filtering');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Before the filtering', 'After the filtering');
subplot(2,1,2);                                                             %%Frequency domain
plot(f, 20*log10(Y));
hold on;
plot(f, 20*log10(Y_f2_filtered));
title('Signal in the frequency domain before and after the filtering');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;
legend('Before the filtering', 'After the filtering');