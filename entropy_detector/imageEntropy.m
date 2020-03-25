function e = imageEntropy(im, n)
% compute the entropy of the image in a sliding window of size n
% written by Christine Erbe, 2013
% im = input image 
% n = window size
% leave out the - sign in the entropy equation

if nargin < 2, n = 7; end

imE = im .* log(im);

imE(isnan(imE)) = 0; % so that 0*log(0) is 0 rather than NaN

% c = ones(n); % flat sliding window

height = n(1); width = n(2); sigma = n(1)/2; % Gaussian, ie smooth, sliding window
[n1,n2] = meshgrid(-(width - 1)/2 : 1 : (width-1)/2, -(height - 1)/2 : 1 : (height-1)/2);
c = exp(-(n1.^2 + n2.^2)/ 2 / sigma^2);
c = c / sum(sum(c));
% figure(4); clf; imagesc(c); colorbar; sum(sum(c))

e = conv2(imE,c);
a = conv2(im,c);

% e = conv2(imE,c,'same');
% a = conv2(im,c,'same');

% e = e ./ a;  % this used to normalise each entropy kernel box with the energy in
% each kernel box; removed bec we want to normalise by the energy of the
% total file; but that is just a one-number addition and falls away when
% the entropy threshold is taken as the mean + std of the whole file

return

%% test
load PSD1
PSD1 = PSD1 - min(PSD1(:)); PSD1 = PSD1+max(PSD1(:))/100; PSD1 = PSD1/sum(PSD1(:));
e = imageEntropy(PSD1,7);

figure;
subplot(211); imagesc(PSD1);
subplot(212); imagesc(e);

e = imageEntropy(PSD1,[7,3]);

figure;
subplot(211); imagesc(PSD1);
subplot(212); imagesc(e);

