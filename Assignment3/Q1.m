% Program polyf.m - A program to demonstrate polyphase filtering
%
% The Problem - We wish to separate two tones at 100 and 200 Hz from
% tones at 1000 and 1100 Hz and then downsample the output by a factor
% of 20
%
% Author          - Brian Lovell (ft. Rob) 3/4/96

Fs = 32000;	% Sampling Frequency
N =  128000;% Number of Samples
Nf = 2200;	% Length of Filter
Np = 14;	% Number of passband samples
Loops = 100; % number of times to repeat filter

subplot(1,1,1)  % reset axes

t = [1:N]/Fs;

% Generate Data
x = zeros(1,N);
x = x + sin(2*pi*50*t);
x = x + sin(2*pi*100*t);
x = x + sin(2*pi*1000*t);
x = x + sin(2*pi*1100*t);


X = Fs*([0:N-1]/N);
Y = (abs(fft(x)));
zoom on;
plot(X',Y');
title('Original data sampled at 32000 Hz');
xlabel('Frequency (Hz)');xlim([0 Fs]);
set(gca, 'XTickLabel', num2str(transpose(get(gca, 'XTick'))));
XX = X; YY = Y;
pause

%% Setting up the Kaiser window and filter
wind = kaiser(Nf,7.30626)';
plot(wind);
title('Kaiser Window');
pause
vect = zeros(1,Nf);
vect(1:Np)=ones(1,Np);
vect(Nf-Np+2:Nf)=ones(1,Np-1);
Y = vect;
X = Fs*[0:length(vect)-1]/length(vect); % Normalise Frequency axis
plot(X',Y');title('Desired Response'); xlim([0 Fs]);
set(gca, 'XTickLabel', num2str(transpose(get(gca, 'XTick'))));
xlabel('Frequency (Hz)');pause
filt = real(fftshift(ifft(vect)));
subplot(1,2,1);
plot(filt); title('Unwindowed Filter'); pause
filt = wind.*real(filt);
subplot(1,2,2);
plot(filt); title('Windowed Filter');pause
Y=20*log10(abs(fft(filt)));
X=Fs*[0:Nf-1]/Nf;
figure;
plot(X',Y'); title('Filter Response'); xlim([0 Fs]);xlabel('Frequency (Hz)');
%line([250 250],get(gca,'YLim'),'Color','r','LineStyle','--');
atten = -75*ones(1,Nf);
hold on;
plot(X',atten');
set(gca, 'XTickLabel', num2str(transpose(get(gca, 'XTick')))); grid;
pause

%% Recording number of operations for single rate FIR filter
tic; % reset time counter
for i=1:Loops % do filtering 100 times so cpu time is significant
    y = conv(filt,x);
end
y = y(length(filt):length(y)-length(filt)); % Remove edge effects

tsingle=toc/Loops;  % The toc function records elapsed time since tic call
disp(' ');
disp(['Time for single rate filtering = ' num2str(tsingle,12)]);

wy = y.*blackman(length(y))'; % Use window to reduce sidelobes
Y=20*log10(abs(fft(wy)));
X=Fs*[0:length(Y)-1]/length(Y);
plot(X',Y');
title('Data filtered at 20000 Hz');
xlabel('Frequency (Hz)');
pause

%% Downsample
downsample = 128;       % The factor to downsample by
y = y(downsample:downsample:length(y)); % Downsample y

wy = y.*blackman(length(y))';
Y=20*log10(abs(fft(wy)));
X=Fs*[0:length(Y)-1]/length(Y)/downsample;  % Normalise frequency axis
XXX = X; YYY = Y;

subplot(2,1,1),plot(X',Y')
title('After downsampling by a factor of 128');
xlabel('Frequency (Hz)');
pause;

%% Convert to Polyphase filter
%filtp = zeros(1,Nf); % zero pad to multiple of 128
%filtp(1:Nf) = filt;
Lf = ceil(length(filt)/downsample); % Must round up to avoid throw aways
filtp = zeros(1,downsample*Lf);
filtp(1:Nf) = filt;

polyfilt = reshape(filtp,downsample,Lf);
polyfilt =flipud(polyfilt); % flip matrix top to bottom to get rows in right order

Lx = ceil(length(x)/downsample);
xp = reshape(x,downsample,Lx);
out = zeros(1,-1+(length(x)+length(filtp))/downsample);
tic;
for i=1:Loops % do this 500 times so cpu time is significant
    for i=1:downsample
	    out = out + conv(polyfilt(i,:),xp(i,:));
    end
end
tpoly = toc/Loops;
disp(' ');
disp(['Time for polyphase filtering = ' num2str(tpoly,12)]);
disp(' ');
disp(['Ratio of CPU time (ideally should be M=128) = ' num2str((tsingle/tpoly),3)]);

out = out(size(polyfilt,2):length(out)-size(polyfilt,2)); % Remove edge effects
wy = out.*blackman(length(out))';
Y=20*log10(abs(fft(wy)));
X=Fs*[0:length(Y)-1]/length(Y)/downsample;
subplot(2,1,2),plot(X',Y');
title('Polyphase filtered');
pause;

subplot(2,1,1);
plot(XX',YY');
title('Original data sampled at 32000 Hz');
subplot(2,1,2);
plot(XXX', YYY');
title('Downsampled data');
