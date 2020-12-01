load EEGdata_assignment3.mat;

% create time array for the two EEG signals
t1 = 0:1/Fs:(length(EEG1)-1)/Fs;
t2 = 0:1/Fs:(length(EEG2)-1)/Fs;

% Compute the spectra of the two signals
[MEEG1,phEEG1,fEEG1] = fourier_dt(EEG1,Fs,'full');
[MEEG2,phEEG2,fEEG2] = fourier_dt(EEG2,Fs,'full');

% plot the limited magnitude spectra of the two ECG signals
figure
subplot(2,1,1)
plot(fEEG1,MEEG1)
xlim([0 30])
xticks([0:10 15:5:30])
ylabel('|X(f)| (mV)')
title('EEG1')
subplot(2,1,2)
plot(fEEG2,MEEG2)
xlim([0 30])
xticks([0:10 15:5:30])
ylabel('|X(f)| (mV)')
title('EEG2')
xlabel('f (Hz)')

% get the indices of particular frequencies
f0i = find(fEEG1 >= 0, 1);
f3i = find(fEEG1 < 3, 1, 'last');
f8i = find(fEEG1 < 8, 1, 'last'); 
f13i = find(fEEG1 < 13, 1, 'last');
f25i = find(fEEG1 < 25, 1, 'last');
f100i = find(fEEG1 <= 100, 1, 'last');

f0i2 = find(fEEG2 >= 0, 1);
f3i2 = find(fEEG2 < 3, 1, 'last'); %use the find function to determine the most representitive index 
f8i2 = find(fEEG2 < 8, 1, 'last'); 
f13i2 = find(fEEG2 < 13, 1, 'last');
f25i2 = find(fEEG2 < 25, 1, 'last');
f100i2 = find(fEEG2 <= 100, 1, 'last');

% calculate band power of particular EEG frequency bands
Pd1 = sum(MEEG1(f0i:f3i).^2);
Pt1 = sum(MEEG1(f3i+1:f8i).^2); %comput all the sums of the squares to compute the band power
Pa1 = sum(MEEG1(f8i+1:f13i).^2);
Pb1 = sum(MEEG1(f13i+1:f25i).^2);
Pg1 = sum(MEEG1(f25i+1:f100i).^2);

Pd2 = sum(MEEG2(f0i2:f3i2).^2);
Pt2 = sum(MEEG2(f3i2+1:f8i2).^2);
Pa2 = sum(MEEG2(f8i2+1:f13i2).^2);
Pb2 = sum(MEEG2(f13i2+1:f25i2).^2);
Pg2 = sum(MEEG2(f25i2+1:f100i2).^2);

P1 = [Pd1 Pt1 Pa1 Pb1 Pg1];
P2 = [Pd2 Pt2 Pa2 Pb2 Pg2];

figure
bands = categorical({'delta', 'theta', 'alpha', 'beta', 'gamma'});
bands = reordercats(bands,{'delta', 'theta', 'alpha', 'beta', 'gamma'});
bar(bands, P1)
xlabel('Frequency Bands')
ylabel('Band Power (mV^2)')
title('EEG1')

figure
bands = categorical({'delta', 'theta', 'alpha', 'beta', 'gamma'});
bands = reordercats(bands,{'delta', 'theta', 'alpha', 'beta', 'gamma'});
bar(bands, P2)
xlabel('Frequency Bands')
ylabel('Band Power (mV^2)') 
title('EEG2')

% normalized by bandwidth
% example bandwidth = 3 - 0 = 3
Pd1n =  Pd1 / 3;
Pt1n = Pt1 / 5; 
Pa1n = Pa1 / 5;
Pb1n = Pb1 / 12;
Pg1n = Pg1 / 75;

Pd2n = Pd2 / 3;
Pt2n = Pt2 / 5;
Pa2n = Pa2 / 5;
Pb2n = Pb2 / 12;
Pg2n = Pg2 / 75;

P1n = [Pd1n Pt1n Pa1n Pb1n Pg1n];
P2n = [Pd2n Pt2n Pa2n Pb2n Pg2n];

figure
bands = categorical({'delta', 'theta', 'alpha', 'beta', 'gamma'}); %categorize the bands
bands = reordercats(bands,{'delta', 'theta', 'alpha', 'beta', 'gamma'}); %random cat of the bands
bar(bands, P1n)
xlabel('Frequency Bands')
ylabel('Normalized Band Power (mV^2/Hz)') 
title('EEG1')

figure
bands = categorical({'delta', 'theta', 'alpha', 'beta', 'gamma'});
bands = reordercats(bands,{'delta', 'theta', 'alpha', 'beta', 'gamma'});
bar(bands, P2n)
xlabel('Frequency Bands')
ylabel('Normalized Band Power (mV^2/Hz)') 
title('EEG2')

% BONUS: plot the spectrograms of the two EEG signals
winlen = 1e3;  % length of the windowed segments
overlap = 500; % number of samples overlapping for each window position
NFFT = 20e3;   % number of points in the FFT (the signal is zero-padded to this length)

figure
subplot(2,1,1)
[s_EEG1,f_EEG1,t_EEG1] = spectrogram(EEG1,winlen,overlap,NFFT,Fs);
imagesc(t_EEG1,f_EEG1,abs(s_EEG1)/winlen)
axis xy
ylim([0 30])
title('EEG1')
ylabel('f (Hz)')
cb1 = colorbar;
cb1.Label.String = '|X(f)| (mV)';
subplot(2,1,2)
[s_EEG2,f_EEG2,t_EEG2] = spectrogram(EEG2,winlen,overlap,NFFT,Fs);
imagesc(t_EEG2,f_EEG2,abs(s_EEG2)/winlen)
axis xy
ylim([0 30])
title('EEG2')
ylabel('f (Hz)')
xlabel('t (s)')
cb2 = colorbar;
cb2.Label.String = '|X(f)| (mV)';
