function c=spectha2M(a,nrows)
% c=spectha2M(a,nrows)
% Only nrows/2 are used, throw neg. frequ. away.
% Returns a time series of 2D spectra.
% Time series is Hamming-windowed and overlapped 50% (as opp to spectha.m).
% True power density spectrogram with all the appropriate normalizations.
% Amplitudes are linear in µPa^2/Hz. For dB, take 10*log10.
% Written by Christine Erbe, 1994

window=hamming(nrows);
n=size(a,1);
ncols=fix(n/(nrows/2))-1;
% c=zeros(nrows/2,ncols);
c = [];
for i=1:ncols
    seg = a((i-1)/2*nrows+1:(i+1)/2*nrows);
    seg = seg .* window;
    temp=abs(fft(seg))/44100;
    temp1=2*((temp(2:nrows/2+1)).^2);
    c(:,i)=temp1*44100/nrows;
end



