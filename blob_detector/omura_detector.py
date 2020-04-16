"""
A Python implementation of the blob detector algorithm described in
  Shyam Madhusudhana, Anita Murray, and Christine Erbe. (2020). "Automatic
  detectors for low-frequency vocalizations of Omura’s whales, Balaenoptera
  omurai: A performance comparison." The Journal of the Acoustical Society
  of America. 147(4). DOI: 10.1121/10.000108
for the detection of Omura's whale vocalizations.
If you use this software, we request that you please cite the above article.

Author: Shyam Madhusudhana <shyamm@cornell.edu>
"""

import numpy as np
from scipy.signal import butter, sosfiltfilt, spectrogram
from scipy.signal.windows import hann
import librosa
import argparse
from timeit import default_timer as timer
from detectsound import BlobExtractor
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def detect_omuras(audio_file, output_file=None):
    """Detect Omura's whale vocalizations in the recording 'audio_file'.
    If 'output_file' is None, the detections will be displayed on-screen.
    Otherwise, 'output_file' must specify a path where detections will be
    written in Raven Selection Table format."""

    t_start = timer()   # Measure data load time

    # Load audio from file. Up/downsample (as necessary) to 200 Hz.
    data, fs = librosa.load(audio_file, sr=200)

    time_load = timer() - t_start

    # Build a Butterworth filter to suppress energies outside the
    # 5-68 Hz band.
    wn = np.array([5, 68]) / (fs / 2)  # convert to angular frequency
    filter_sos = butter(64, wn, btype='bandpass', output='sos')

    t_start = timer()   # Measure data preparation time

    # Filter the audio
    data = sosfiltfilt(filter_sos, data).astype(np.float32)

    # Compute spectrogram with 1.0 s Hann window (and NFFT) and 70% overlap
    f, t, spec = spectrogram(data, fs=fs, window=hann(200),
                             nperseg=200, noverlap=140,
                             nfft=200, detrend=False, mode='psd')

    # Find out the indices of where to clip the frequency axis
    valid_f_idx_start = f.searchsorted(5, side='left')
    valid_f_idx_end = f.searchsorted(68, side='right') - 1

    # Clip frequency extents
    spec = spec[valid_f_idx_start:(valid_f_idx_end + 1), :]
    f = f[valid_f_idx_start:(valid_f_idx_end + 1)]

    # Convert to decibels
    spec = 10 * np.log10(spec + 1e-15)

    time_prep = timer() - t_start

    t_start = timer()   # Measure detector time

    # Create a detector instance with parameters tuned for detection
    # of Omura's whale vocalizations.
    be = BlobExtractor(len(f),
                       centroid_bounds=[17, 39],    # Translates to 22-44 Hz
                       min_snr=0.,                  # All non-negative SNRs
                       cost_std=[1.5, 2.5, 2.0],
                       min_frames=5,                # Translates to 1.5 s
                       first_scale=1.3,
                       t_blur_sigma=3)

    # Run the detector on the spectrogram
    detections = be.extract_blobs(spec)

    time_detect = timer() - t_start

    print('Source        : {:s}'.format(audio_file))
    print('Audio duration: {:.6f} s'.format(len(data) / fs))
    print('Time taken: ')
    print('    Audio load: {:.6f} s'.format(time_load))
    print('    Audio prep: {:.6f} s'.format(time_prep))
    print('    Detection : {:.6f} s'.format(time_detect))

    if output_file is not None:
        # Write out detections into a Raven Selection Table file.
        with open(output_file, 'w') as of:

            # Write header
            of.write('Selection\tBegin Time (s)\tEnd Time (s)\t' +
                     'Low Freq (Hz)\tHigh Freq (Hz)\tScore\n')

            discarded = 0
            for b_idx, blob in enumerate(detections):

                # Discard very narrow blobs (median bandwidth of 5 Hz)
                if np.median(f[blob.edges[:, 1]] - f[blob.edges[:, 0]]) < 5:
                    discarded += 1
                    continue

                of.write(
                    '{:d}\t{:.4f}\t{:.4f}\t{:.3f}\t{:.3f}\t{:.4f}\n'.format(
                        b_idx + 1 - discarded,
                        t[blob.first_frame], t[blob.last_frame],
                        f[blob.edges[:, 0].min()], f[blob.edges[:, 1].max()],
                        np.median(blob.snrs)))

        print('{:d} detections written to {:s}'.format(
            len(detections) - discarded, output_file))

    else:
        # Display detections
        fig, ax = plt.subplots()

        # Show the spectrogram
        ax.imshow(spec, cmap=plt.cm.gray_r,
                  extent=(t[0] - (t[1] - t[0]) / 2, t[-1] + (t[1] - t[0]) / 2,
                          f[0] - (f[1] - f[0]) / 2, f[-1] + (f[1] - f[0]) / 2),
                  interpolation='none', origin='lower', aspect='auto')

        # Overlay detections
        patches = []
        snrs = []
        discarded = 0
        for blob in detections:

            # Discard very narrow blobs (median bandwidth of 5 Hz)
            if np.median(f[blob.edges[:, 1]] - f[blob.edges[:, 0]]) < 5:
                discarded += 1
                continue

            # Build list of polygon points from reported blob edges
            temp = np.arange(blob.first_frame, blob.last_frame + 1)
            time_axis = t[np.concatenate([temp, np.flip(temp), [temp[0]]])]
            freq_axis = f[
                np.concatenate([blob.edges[:, 0],
                                np.flip(blob.edges[:, 1]),
                                blob.edges[0:1, 0]])]
            patches.append(Polygon(np.stack([time_axis, freq_axis]).T, True))
            snrs.append(np.median(blob.snrs))

        p = PatchCollection(patches, cmap=plt.cm.Greens, alpha=0.4, ec='k')
        p.set_array(np.asarray(snrs))
        ax.add_collection(p)

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Frequency (Hz)')
        ax.set_title('{:s}: {:d} detections'.format(
            repr(audio_file), len(detections) - discarded))
        plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='omura_detector', allow_abbrev=False,
        description='Detect Omura\'s whale vocalizations in recordings.')
    parser.add_argument(
        'src', metavar='<AUDIO SOURCE>',
        help='Path to the audio file to process.')
    parser.add_argument(
        'out', metavar='<OUTPUT FILE>', nargs='?',
        help='(Optional) Path to an output file into which detection ' +
             'results (Raven selection tables) will be written. If not ' +
             'specified, detection results will be displayed on screen.')

    args = parser.parse_args()

    detect_omuras(args.src, args.out)
