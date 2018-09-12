   x = sqrt(2)*randn(1000,1);
   Numlags = 50;
   [xc,lags] = xcorr(x,Numlags,'coeff');
   stem(lags(51:end),xc(51:end))
   % power spectrum
   Fs = 1; % sampling frequency
   [Pxx,F] = periodogram(x,[],length(x),Fs);
   figure;
   plot(F,10*log10(Pxx))