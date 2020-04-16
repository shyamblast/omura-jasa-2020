%% OmuraEntropy2XSTD.m

% written by Christine Erbe, 17 July 2018

% The code is a 2D version (2D = operating on spectrograms) of the 1D
% version (1D = operating on spectra) that was published earlier:
% Erbe, C., and King, A. R. (2008). Automatic detection of marine mammals
% using information entropy. The Journal of the Acoustical Society of America,
% 124(5), 2833-2840, doi: 10.1121/1.2982368.

% Modifications:
% 17 Oct 2018: OmuraEntropy2XSTD computes 2 thresholds (merely for comparison)
% and combines the detections into 1 output file.

% If you use this software, we request that you please cite:
%   Shyam Madhusudhana, Anita Murray, and Christine Erbe. (2020). "Automatic
%   detectors for low-frequency vocalizations of Omura’s whales, Balaenoptera
%   omurai: A performance comparison." The Journal of the Acoustical Society
%   of America. 147(4). DOI: 10.1121/10.000108

% Notes:
% For performance assessment (e.g., ROC curves), need to vary Xstd on line 57.
% Let it range from 0 to 5 perhaps.

% Requirements:
% Function spectha2M.m written by Christine Erbe, 1994 or any other
% function that computes a spectrum.


%% set directories and parameters
clear all; close all;

wdir = '/Matlab/mProjects/Omuras whales/'; % working directory
wavdir = '/Volumes/TOSHIBA EXT/Testing Data/WAV Files/'; % test data
cd(wdir); % change to working directory
path(path,wdir);

% Omura call min f, max f, min t, max t [Hz, s]
minlowf = 17; maxlowf = 21; minhighf = 29; maxhighf = 70; mint = 2; maxt = 19;
minlowf = 15; maxlowf = 26; minhighf = 26; maxhighf = 62; mint = 2; maxt = 15;

bw = maxhighf - minlowf;

% limit fmax in PSD
fmax = 90; % Hz

% entropy kernel for smoothing
heightpx = 7; widthpx = 5; % kernel box; make them odd [pixels]
kernelsizey = bw / heightpx; kernelsizex = maxt / widthpx * 2;
% if NFFT = fs, then 1 px = 1 Hz, and then 9 kernels fit into a near-range Omura's call.
% Given FFT had 50% overlap, 8 kernels fit into max duration of an Omura's call.

minslope = 6.6/(15.3*2+widthpx); % min 90% bw / max dur * 50% overlap
maxslope = 41.8/(2.1*2+widthpx); % delta-f/delta-t

% entropy threshold: thresh = mean(Hvec) + Xstd * std(Hvec);
XstdA = 4 ; XstdB = 5;
% If set to 2, then it finds weak calls in strong ambient but misses strong calls in no ambient because
% the lumps are larger than the max pixel size.
% If set to 3, then it does better for strong calls in no ambient, but
% misses weak calls in strong ambient.
% Problem also is if a seismic pulse touches an Omura call.
% Note there are hardly any false detections w either Xstd = 2 or 3.
% So maybe combine results from 2 runs w diff Xstd.

%% recording system specifications

cal = -142; % calibration dB
callin = 10.^(abs(cal)/20); % linear pressure calibration factor


%% batch processing
cd(wavdir); nwavdir = dir(wavdir); % only run from 3 : end as 1,2,end are Mac folders

Tall = []; % summary table
Tfinfo = []; % wav file info

for kfolder = 3 : length(nwavdir) - 1 % loop over the set folders
    cd([nwavdir(kfolder).folder,'/',nwavdir(kfolder).name]) % loop over all wav files within 1 folder
    d = dir('*.wav');
    
    for k = 1 : length(d)
        
        % read audio file
        wavfname = d(k).name;
        [wf, fs] = audioread(wavfname);
        d(k).fs = fs;
        wf = wf - mean(wf); % use this to remove any potential DC offset
        p = wf * callin; % calibrated pressure time series
        clear wf; % remove wf to free memory
        d(k).dur = length(p) / fs;
        p = p(:);
        % figure(2); clf; plot(p)
        
        % compute spectrogram
        timewindow = fs;
        PSD = spectha2M(p,timewindow); % "2" = 50% overlap
        n = length(p); duration = n/fs; % in seconds
        t = ceil(size(PSD,2)/2) * timewindow / fs;
        times = .5:.5:t;
        
        % remove unwanted f
        PSD1 = PSD(1:fmax,:);
        
        %         figure(1); clf;
        %         h = imagesc(times,1:fmax,10*log10(PSD1));
        %         maxc = 10*log10(max(max(PSD1)));
        %         set(gca,'clim',[30 maxc]); set(gca,'tickdir','out');
        %         axis xy; % set(gca,'ylim',[0 300]);
        %         h1 = colorbar('peer',gca); colormap('jet'); set(get(h1,'ylabel'),'String','dB re 1 \muPa^2/Hz');
        %         xlabel('Time [s]'); ylabel('Frequency [Hz]'); title('Original PSD')
        
        % normalise to lie between 0 and 1
        PSD1 = 10*log10(PSD1);
        PSD1 = PSD1 - min(min(PSD1));
        PSD1 = PSD1 / max(PSD1(:));
        
        % convolve with entropy kernel
        height = heightpx; % 1 pixel = 1 Hz if NFFT = fs
        width = widthpx; % 1 pixel = 0.5 s if FFT had 50% overlap
        
        imE = PSD1 .* log(PSD1);
        % figure(2); imagesc(imE); axis xy; colorbar; colormap('jet'); title('Entropy')
        
        H = imageEntropy(PSD1,[height,width]);
        H = H((height-1)/2+1 : end - (height-1)/2, (width-1)/2+1 : end - (width-1)/2); % cut edges
        
        %     figure(5); clf;
        %     subplot(211); imagesc(times,1:fmax,PSD1); axis xy; colorbar;  colormap('jet');
        %     xlabel('Time [s]'); ylabel('Frequency [Hz]');
        %     title('Spectrogram, mean spectrum removed');
        %     subplot(212); imagesc(times,1:fmax,H); axis xy; colorbar;  colormap('jet');
        %     xlabel('Time [s]'); ylabel('Frequency [Hz]');
        %     title('-Entropy of the spectrogram');
        
        % h=findobj(gcf,'type','axes'); linkaxes(h,'x'); zoom xon
        
        % treshold entropy to build mask for original spectrum
        H1 = H((height-1)/2+1 : end - (height-1)/2, (width-1)/2+1 : end - (width-1)/2); % cut edges once more to avoid spill from half kernel
        Hvec = reshape(H1,size(H1,1)*size(H1,2),1);
        
        threshA = mean(Hvec) + XstdA * std(Hvec);
        threshB = mean(Hvec) + XstdB * std(Hvec);
        
        PSD2A = PSD1; PSD2A(H < threshA) = 0;
        PSD2B = PSD1; PSD2B(H < threshB) = 0;
        PSD2A(1:(height-1)/2+1, 1: end) = 0; PSD2A(end - (height-1)/2+1 : end, 1: end) = 0;
        PSD2A(1:end, 1: (width-1)/2+1) = 0; PSD2A(1: end, end - (width-1)/2+1: end) = 0;
        PSD2B(1:(height-1)/2+1, 1: end) = 0; PSD2B(end - (height-1)/2+1 : end, 1: end) = 0;
        PSD2B(1:end, 1: (width-1)/2+1) = 0; PSD2B(1: end, end - (width-1)/2+1: end) = 0;
        
        %         figure(8); clf; subplot(2,1,1)
        %         imagesc(times,1:fmax,PSD2A); axis xy; colorbar; colormap('jet'); ylabel('Frequency [Hz]');
        %         title(['Spectrogram, masked where -entropy < thresh. Xstd = ', num2str(XstdA)]);
        %         subplot(2,1,2); imagesc(times,1:fmax,PSD2B); axis xy; colorbar; colormap('jet');
        %         xlabel('Time [s]'); ylabel('Frequency [Hz]');
        %         title(['Spectrogram, masked where -entropy < thresh. Xstd = ', num2str(XstdB)]);
        
        PSD2A(H >= threshA) = 1; PSD2B(H >= threshB) = 1;
        PSD2A(1:(height-1)/2+1, 1: end) = 0; PSD2A(end - (height-1)/2+1 : end, 1: end) = 0;
        PSD2A(1:end, 1: (width-1)/2+1) = 0; PSD2A(1: end, end - (width-1)/2+1: end) = 0;
        PSD2B(1:(height-1)/2+1, 1: end) = 0; PSD2B(end - (height-1)/2+1 : end, 1: end) = 0;
        PSD2B(1:end, 1: (width-1)/2+1) = 0; PSD2B(1: end, end - (width-1)/2+1: end) = 0;
        % figure(9); clf; imagesc(times,1:fs/2,PSD2A); axis xy; colorbar
        % xlabel('Time [s]'); ylabel('Frequency [Hz]');
        % title('Spectrogram, forced to [0,1]');
        
        minlumpsize = (minhighf - maxlowf) * mint * 2;
        maxlumpsize = (maxhighf - minlowf) * maxt * 2;
        
        % find connecting pixels
        % im2 has zeros where nthg was detected, ones for the 1st detection, twos for the 2nd detection etc.
        [im2A,nA] = bwlabel(PSD2A); % nA blobs found, labelled 1:nA in im2A
        newcounter = 0; bbA = []; % now remove blobs that are too small or too large
        for i = 1 : nA  % loop over the blobs found
            [r, c] = find(im2A==i);
            if length(r) < minlumpsize, im2A(find(im2A==i)) = 0; % remove blobs that are too small
            elseif length(r) > maxlumpsize, im2A(find(im2A==i)) = 0; % remove blolbs that are too big
            else newcounter = newcounter + 1; % renumber the blobs
                im2A(find(im2A==i)) = newcounter;
                % compute bounding box for each blob
                bbA(1,newcounter) = min(r); % minf
                bbA(2,newcounter) = max(r); % maxf
                bbA(3,newcounter) = min(c); % startt
                bbA(4,newcounter) = max(c); % endt
            end
        end
        
        [im2B,nB] = bwlabel(PSD2B);
        newcounter = 0; bbB = [];
        for i = 1 : nB
            [r, c] = find(im2B==i);
            if length(r) < minlumpsize, im2B(find(im2B==i)) = 0;
            elseif length(r) > maxlumpsize, im2B(find(im2B==i)) = 0;
            else newcounter = newcounter + 1;
                im2B(find(im2B==i)) = newcounter;
                bbB(1,newcounter) = min(r); bbB(2,newcounter) = max(r);
                bbB(3,newcounter) = min(c); bbB(4,newcounter) = max(c);
            end
        end
        
        figure(10); clf;
        subplot(2,1,1); imagesc(times,1:fmax,im2A); axis xy; colorbar; ylabel('Frequency [Hz]');
        title(['Detected signals, each hue is different signal. Xstd = ', num2str(XstdA)]);
        subplot(2,1,2); imagesc(times,1:fmax,im2B); axis xy; colorbar; ylabel('Frequency [Hz]');
        title(['Detected signals, each hue is different signal. Xstd = ', num2str(XstdB)]);
        xlabel('Time [s]');
        % each different hue of colour is a different lump = signal component
        
        % find bb within Omura f- and t-range
        % bb has rows for: minf, maxf, startt, endt
        bbA = [bbA; bbA(4,:) - bbA(3,:)]; % add a 5th row that is duration of blob
        bbA = [bbA; (bbA(2,:) - bbA(1,:))./bbA(5,:)]; % add a 6th row that is the diagonal dimension
        
        fcriterionA = find(bbA(1,:) >= minlowf & bbA(2,:) <= maxhighf & bbA(5,:) >= mint * 2 & ...
            bbA(5,:) <= (maxt + width) * 2 & minslope < bbA(6,:) & bbA(6,:) < maxslope); % these are the blobs that meet the f-criterion
        
        bbB = [bbB; bbB(4,:) - bbB(3,:)]; % add a 5th row that is duration of blob
        bbB = [bbB; (bbB(2,:) - bbB(1,:))./bbB(5,:)]; % add a 6th row that is the diagonal dimension
        
        fcriterionB = find(bbB(1,:) >= minlowf & bbB(2,:) <= maxhighf & bbB(5,:) >= mint * 2 & ...
            bbB(5,:) <= (maxt + width) * 2 & minslope < bbB(6,:) & bbB(6,:) < maxslope); % these are the blobs that meet the f-criterion
        
        % plot spectrum after masking
        im3A = ismember(im2A,fcriterionA); % im3 only has ones and zeros for detections and not
        signalsA = PSD1; signalsA(im3A == 0) = 0;
        im3B = ismember(im2B,fcriterionB);
        signalsB = PSD1; signalsB(im3B == 0) = 0;
        
        %         figure(11); clf;
        %         subplot(2,1,1); imagesc(times,1:fmax,signalsA); axis xy; colorbar
        %         xlabel('Time [s]'); ylabel('Frequency [Hz]');
        %         title(['Detected signals. Xstd = ', num2str(XstdA)]);
        %         subplot(2,1,2); imagesc(times,1:fmax,signalsB); axis xy; colorbar
        %         xlabel('Time [s]'); ylabel('Frequency [Hz]');
        %         title(['Detected signals. Xstd = ', num2str(XstdB)]);
        
        % write to file starttime of detection in file
        % bb(row 3) has the start pixels of each detection
        % use only the fcriterion columns
        % times has the times corresponding to the pixels; remember FFT had 50% overlap
        % so far, code ran w 2 thresholds at once, and this outputs results for each threshold;
        % there is no averaging over thresholds or so
        
        detectstarttimes = [times(bbA(3,fcriterionA)), times(bbB(3,fcriterionB))];
        detectendtimes = [times(bbA(4,fcriterionA)), times(bbB(4,fcriterionB))];
        detectduration = detectendtimes - detectstarttimes;
        [sorted, isort] = sort(detectstarttimes);
        d(k).startdetections = detectstarttimes(isort)';  % s
        d(k).enddetections = detectendtimes(isort)';  % s
        minf = [bbA(1,fcriterionA), bbB(1,fcriterionB)];  % Hz
        maxf = [bbA(2,fcriterionA), bbB(2,fcriterionB)];  % Hz
        d(k).minf = minf(isort)';
        d(k).maxf = maxf(isort)';
        Xstd = [XstdA * ones(length(fcriterionA),1); XstdB * ones(length(fcriterionB),1)];
        Xstd = Xstd(isort);
        
        finame = {}; threshs = []; folders = {};
        for jj = 1 : length(detectstarttimes)
            finame{jj} = d(k).name;
            threshs(jj) = Xstd(jj)';
            folders{jj} = d(k).folder;
        end
        d(k).detectfnames = finame';
        d(k).threshs = threshs';
        d(k).detfolders = folders';
        
    end % over wav files in any k folder
    
    % write to xls
    fnames = []; startt = []; endt = []; minf = []; maxf = []; threshs = [];
    fsizes = []; FS = []; dur = []; folders = []; kfolders = []; knames = [];
    
    for k = 1 : length(d)
        folders = [folders; d(k).detfolders];
        fnames = [fnames; d(k).detectfnames];
        startt = [startt; d(k).startdetections];
        endt = [endt; d(k).enddetections];
        minf = [minf; d(k).minf];
        maxf = [maxf; d(k).maxf];
        threshs = [threshs; d(k).threshs];
        kfolders = [kfolders; d(k).folder]; % each folder
        knames = [knames; d(k).name]; % all the wav in each folder, no matter the detections
        fsizes = [fsizes; d(k).bytes]; % size of each wav
        FS = [FS; d(k).fs]; % fs
        dur = [dur; d(k).dur]; % dur [s]
    end
    T = table(folders,fnames,startt, endt, minf, maxf, threshs, 'VariableNames', ...
        {'Folder', 'Fname','Startt','Endt','Minf','Maxf','Xstd'});
    Tall = [Tall; T];
    
    Tfileinfo = table(kfolders, knames,fsizes, FS, dur,'VariableNames',...
        {'Folder','Fname','Sizebytes','FSHz','Durs'});
    Tfinfo = [Tfinfo; Tfileinfo];
    
    kfolder
end % over all kfolders

cd(wdir)
xlsname = ['threshs_',num2str(XstdA),'_',num2str(XstdB),'.csv'];
writetable(Tall,xlsname);

xlsname = 'wavfileinfo.csv'; writetable(Tfinfo,xlsname);

unique(Tfinfo.FSHz)
sum(Tfinfo.Sizebytes)
min(Tfinfo.Durs)/60, mean(Tfinfo.Durs)/60, max(Tfinfo.Durs)/60

