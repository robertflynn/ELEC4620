clc; clear all; close all;

% Program polyf.m - A program to demonstrate polyphase filtering
%
% The Problem - We wish to separate two tones at 100 and 200 Hz from
% tones at 1000 and 1100 Hz and the downsample the output by a factor
% of 20
%
% Author          - Brian Lovell (ft. MC Rob 2.0)  3/4/96

Ds = 128;                       % Downsample rate (same as upsample rate)
Fcd = 150;                      % Cutoff Frequency downsampler
Ftd = 100;                      % Transition Bandwidth downsampler
Fs = 32000;                     % Sampling Frequency
N =  80000;                     % Number of Samples
Nfd = 2000;                     % Length of Filter downsampler
Npd = ceil(Fcd/(Fs/Nfd)) + 3    % Number of passband samples downsampler
Loops = 100;                    % number of times to repeat filter

t = [1:N]/Fs;

% Generate Data
x = zeros(1,N);
x = x + sin(2*pi*50*t);       % tone at 50 Hz spacing
x = x + sin(2*pi*100*t);      % tone at 100 Hz spacing
x = x + sin(2*pi*1000*t);     % tone at 1000 Hz spacing
x = x + sin(2*pi*1100*t);     % tone at 1100 Hz spacing

% Plot the frequency response of the original data
X = Fs*([0:N-1]/N);
Y = (abs(fft(x)));
zoom on;
plot(X',Y');
title(sprintf('Original data sampled at %d Hz', Fs));
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');
pause

%% Create the filter
% Desired response vector
vect = zeros(1,Nfd);
vect(1:Npd)=ones(1,Npd);
vect(Nfd-Npd+2:Nfd)=ones(1,Npd-1);
Y = vect;
X = Fs*[0:length(vect)-1]/length(vect);

% Window coeffs.
wind = kaiser(Nfd,7.30626)';

% Convolve window and desired response
filt = real(fftshift(ifft(vect)));
filt = wind.*real(filt);

% Plot the frequency response of the filter
Y=20*log10(abs(fft(filt)));
X=Fs*[0:Nfd-1]/(Nfd);
plot(X',Y'); title('Filter Response'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude (dB)');
pause

%% Normal downsampling (filter then downsample)
% Time calcualtion
tic; % reset time counter
for n=1:Loops % do filtering multiple times so cpu time is significant
    y = conv(filt,x);
end

% Record time taken
tsingle=toc/Loops;
disp(' ');
disp(['Time for single rate filtering = ' num2str(tsingle,12)]);

y = y(length(filt):length(y)-length(filt)); % Remove edge effects

% Plot the filtered response
wy = y.*blackman(length(y))'; % Use window to reduce sidelobes
Y=20*log10(abs(fft(wy)));
X=Fs*[0:length(Y)-1]/length(Y);
plot(X',Y');
title(sprintf('Data filtered at %d Hz', Fs));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
pause

% Downsample and plot the above response
y = y(Ds:Ds:length(y)); 

wy = y.*blackman(length(y))';
Y=20*log10(abs(fft(wy)));
X=Fs*[0:length(Y)-1]/length(Y)/Ds;

subplot(2,1,1),plot(X',Y')
title(sprintf('After downsampling by a factor of %d', Ds));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');


%% Polyphase downsampling
% Create the polyphase filter matrix
filtp = zeros(1,Ds*ceil(length(filt)/Ds)); % zero pad to multiple of Ds
filtp(1:Nfd) = filt; % Copy the calculated filter coeffs.
polyfilt = reshape(filtp,Ds,ceil(length(filtp)/Ds)); % Enter the Matrix
polyfilt =flipud(polyfilt); % flip matrix top to bottom to get rows in right order

xp = reshape(x,Ds,ceil(length(x)/Ds)); % input data (sampled signal)
out = zeros(1,-1+(length(x)+length(filtp))/Ds); % output matrix

tic; % Reset timer counter
for n=1:Loops % do this multiple times so cpu time is significant
    % Carry out the filtering for each decomposition
    for i=1:Ds
	    out = out + conv(polyfilt(i,:),xp(i,:));
    end
end

% Calculate time taken
tpoly = toc/Loops;
disp(' ');
disp(['Time for polyphase filtering = ' num2str(tpoly,12)]);
disp(' ');
disp(['Ratio of CPU time (ideally should be M=' num2str(Ds) ') = ' num2str((tsingle/tpoly),3)]);

out = out(size(polyfilt,2):length(out)-size(polyfilt,2)); % Remove edge effects
wy = out.*blackman(length(out))'; % Use window to reduce sidelobes

% Plot the polyphase filter response
Y=20*log10(abs(fft(wy)));
X=Fs*[0:length(Y)-1]/length(Y)/Ds;
subplot(2,1,2),plot(X',Y');
title('Polyphase filtered');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
pause