# PHAIN implementation using LTFAT
Computationally faster version of PHAIN [1].
## Goal
Reimplement PHAIN [1,2] such that fast routines of the LTFAT time-frequency toolbox [3] are utilized, resulting in a lower computational time while maintaining the reconstruction quality.

## Work done
- Functions for the Gabor transform were substituted by optimized routines from the LTFAT toolbox. Computation of the Gabor transform is about 35 % faster this way, resulting in approximately **10 % total speedup** of PHAIN.
Beware that the original code treats the signal as non-periodic, while it is considered periodic by LTFAT. In effect, if a gap would be present near the signal border, the reconstruction would be different; otherwise, the results are identical (see demo below).
- Projection to the reliable set has also been substituted by a more efficient code improving its computational speed by approximately 70 %. However, this improvement does not effect the total speedup due to its already fast computation. 
- paramsolver.x0?

## Experiment
The experiment was run to prove the goals, i.e. speedup not affecting the quality. The original and the optimized implementation was run on the same signal with the same gaps in the signal. Elepsed time was measured by the tic/toc commands, and the reconstruction quality was measured using the SNR only in the gaps.

The tests were run on signals from 
...
...

Tests were run in Matlab 2024b on PC with Intel Core i7-6829HQ CPU @2.7GHz, 16 GB RAM and Windows 10. The Matlab codes use the Signal Processing Toolbox and LTFAT [3].

## References
[1] Tanaka, Tomoro, Kohei Yatabe, and Yasuhiro Oikawa, “PHAIN: Audio inpainting via phase-aware optimization with instantaneous frequency,” IEEE/ACM Transactions on Audio, Speech, and Language Processing, Sep 2024.

[2] TomoroTanaka, “TomoroTanaka/PHAIN: v1.0.1”. Zenodo, Jun. 30, 2024. doi:10.5281/zenodo.12597334.

[3] Zdeněk Průša, Peter L. Søndergaard, Nicki Holighaus, Christoph Wiesmeyr, Peter Balazs, “The Large Time-Frequency Analysis Toolbox 2.0,” Sound, Music, and Motion, Lecture Notes in Computer Science 2014, pp 419-442.
