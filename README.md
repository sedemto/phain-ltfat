# PHAIN implementation using LTFAT
Computationally faster version of PHAIN [1].

## Goal
Reimplement PHAIN [1,2] such that fast routines of the LTFAT time-frequency toolbox [3] are utilized, resulting in a lower computational time while maintaining the reconstruction quality.

## Work done
- Functions for the Gabor transform were substituted by optimized routines from the LTFAT toolbox. Computation of the Gabor transform is about 45 % faster this way, resulting in approximately **10 % total speedup** of PHAIN.
Beware that the original code treats the signal as non-periodic, while it is considered periodic by LTFAT. In effect, if a gap would be present near the signal border, the reconstruction would be different; otherwise, the results are identical (see demos below).
- Projection to the reliable set has also been substituted by a more efficient code improving its computational speed by approximately 50 %. However, this improvement does not effect the total speedup due to its already fast computation.
- In the outer loop, phase gets repeatedly updated. During an inspection of the code, we found out that when the signal fed into the CP algorithm (the inner loop) gets updated as well, it leads to better stabilisation of the Chambolle-Pock algorithm (see figure below). For some signals, it can also lead to better results. Instead of leaving it fixed as in the original code, we added an option to change it using the parameter `updateInputCP`.

## Experiment
The experiment was run to prove the goals, i.e. speedup not affecting the quality. The original and the optimized implementation was run on the same signal with the same gaps in the signal. Elapsed time was measured by the tic/toc commands, and the reconstruction quality was measured using the SNR only in the gaps.


The tests were run on signals from DPAI dataset available at [DPAI](https://github.com/fmiotello/dpai). Multiple tests were done, all codes available in the `demos` folder :
- to prove that the reconstruction quality is the same for both implementations the `demo_reconstruction.m` can be run. It produces the following image:
  <img width="1920" height="973" alt="comparisonLTFATvsOriginalpng" src="https://github.com/user-attachments/assets/93469f57-eb27-445b-a819-16ea215d6e02" />

- to compare the difference in speed between the DGT from LTFAT and the DGT from the original code run `demo_DGT.m`. The table below is acquired. 

| Test Number | DGT Original code [s] | DGT LTFAT code [s] | DGT Improvement [%] | iDGT Original code [s] | iDGT LTFAT code [s] | iDGT Improvement [%] | Both Original code [s] | Both LTFAT code [s] | Both Improvement [%] |
|:--------------:|:------------------:|:---------------:|:---------------------:|:------------------:|:---------------:|:---------------------:|:------------------:|:---------------:|:---------------------:|
| 1            | 4.2116           | 2.2901       | 45.6249            | 5.1495           | 2.2625       | 56.0634            | 9.3007           | 4.4351       | 52.3139            |
| 2            | 4.0899           | 2.3402       | 42.7820            | 4.9792           | 2.4354       | 51.0895            | 9.7830           | 4.8119       | 50.8131            |
| 3            | 4.7095           | 2.3324       | 50.4758            | 5.0017           | 2.3381       | 53.2544            | 9.2260           | 4.7211       | 48.8278            |
| 4            | 4.0667           | 2.2608       | 44.4080            | 4.9769           | 2.3798       | 52.1835            | 9.4594           | 4.8343       | 48.8945            |
| 5            | 4.2012           | 2.3322       | 44.4888            | 4.9826           | 2.3158       | 53.5213            | 10.1460          | 4.9532       | 51.1808            |
| 6            | 4.2371           | 2.4842       | 41.3701            | 5.1303           | 2.6454       | 48.4365            | 9.2817           | 4.5700       | 50.7630            |
| 7            | 3.9911           | 2.3136       | 42.0322            | 4.7256           | 2.2643       | 52.0847            | 9.0098           | 5.9609       | 33.8402            |
| 8            | 4.1466           | 2.2976       | 44.5906            | 4.6344           | 2.4144       | 51.6367            | 8.7427           | 4.6189       | 47.0975            |
| 9            | 3.9556           | 2.2343       | 43.5145            | 4.6164           | 2.4551       | 46.8176            | 8.9491           | 4.6126       | 48.4828            |
| 10           | 4.0874           | 2.2765       | 44.3031            | 4.7099           | 2.2437       | 52.3615            | 8.9211           | 4.4772       | 49.8134            |
| **Average**  | 4.1697           | 2.3162   |**44.3590**       | 4.8906       | 2.3581   | **51.7449**      | 9.2955       | 4.7995   | **48.2878**      |

- to compare the differences between the projections a small demo called `demo_proj.m` was used.
- SNR of the reconstructed signals with the parameter `updateInputCP = false` (blue) and `updateInputCP = true` (orange):
<img width="1920" height="973" alt="diffWithWithoutUpdate" src="https://github.com/user-attachments/assets/0c148e33-5230-4809-bcf6-a035377c6256" />

### Notes
The `main.m` function in `U-PHAIN` folder can be used to process datasets in bulk with the default option being the DPAI dataset.

Even though LTFAT has functions for the calculation of instantaneous frequency, the original functions "calcInstFreq", "instPhaseCorrection" and "invInstPhaseCorrection" from [1] were utilized due to faster computational time.

Tests were run in Matlab 2024b on PC with Intel Core i7-6829HQ CPU @2.7GHz, 16 GB RAM and Windows 10. The Matlab codes use the Signal Processing Toolbox and LTFAT [3].

## References
[1] Tanaka, Tomoro, Kohei Yatabe, and Yasuhiro Oikawa, “PHAIN: Audio inpainting via phase-aware optimization with instantaneous frequency,” IEEE/ACM Transactions on Audio, Speech, and Language Processing, Sep 2024.

[2] TomoroTanaka, “TomoroTanaka/PHAIN: v1.0.1”. Zenodo, Jun. 30, 2024. doi:10.5281/zenodo.12597334.

[3] Zdeněk Průša, Peter L. Søndergaard, Nicki Holighaus, Christoph Wiesmeyr, Peter Balazs, “The Large Time-Frequency Analysis Toolbox 2.0,” Sound, Music, and Motion, Lecture Notes in Computer Science 2014, pp 419-442.
