# PHAIN implementation using LTFAT
Computationally faster version of PHAIN [1].

## Goal
Reimplement PHAIN [1,2] such that fast routines of the LTFAT time-frequency toolbox [3] are utilized, resulting in a lower computational time while maintaining the reconstruction quality.

## Work done
- Functions for the Gabor transform were substituted by optimized routines from the [LTFAT toolbox](https://github.com/ltfat/ltfat/releases). Computation of the Gabor transform is about 22 % faster this way. Beware that the original code treats the signal as non-periodic, while it is considered periodic by the LTFAT. In effect, if a gap would be present close to the signal border, the reconstruction would be different; otherwise, the results are identical (see demos below).
- Functions for the phase correction R and its adjoint R* were rewritten. Even though the instantaneous frequency is fixed inside the Chambolle-Pock algorithm (CPA), the phase rotations were originally calculated in each iteration, which was unnecessary. Now, the phase rotations are calculated in the outer loop, which speeds up the computation of R and R* by approximately 27 %.
- Time-directional variation operator D and its adjoint operator D* are calculated using Matlab `diff` function, which results in approximately 43 % speed up compared with the ones in the original algortihm.
- Projection to the reliable set has also been substituted by a more efficient code improving its computational speed by approximately 33 %. However, this improvement does not effect the total speedup due to its negligible share on the total time.
- All of the changes result in approximately **22 % total speedup** of PHAIN.
- In the outer loop, the instanteneous frequency gets repeatedly updated. During an inspection of the original code, we found out that when the signal fed into the CP algorithm (the inner loop) gets updated as well, opposite to the original code, it leads to a stabilisation of the Chambolle-Pock algorithm (see the figure below). For some signals, it can also lead to better results. Instead of leaving it fixed as in the original code, we added a switch `updateInputCP`.

## Experiment
The experiment was run to prove the goals, i.e. achieving a speedup not affecting the quality. The original and the optimized implementation were run on the same signal with the same gaps in the signal. Elapsed time was measured by the tic/toc commands, and the reconstruction quality was measured using the SNR (solely in the gaps).


The tests were run on signals from the [DPAI dataset](https://github.com/fmiotello/dpai). Multiple tests were performed, all codes available in the `demos` folder :
- To prove that the reconstruction quality is the same for both implementations the `demo_reconstruction.m` can be run. It also computes the SNR of both reconstructions, measures their execution time, and produces the following image:
  <img width="1920" height="973" alt="comparisonLTFATvsOriginalpng" src="https://github.com/user-attachments/assets/93469f57-eb27-445b-a819-16ea215d6e02" />
The reconstructions are practically identical and they result in SNRs 3.7397 and 3.7397, respectively.

- To compare the difference in speed between the DGT from LTFAT and the DGT from the original code run `demo_DGT.m`. The table below is acquired, which shows that on average, the LTFAT implementaion is about 22 % faster.

| Test Number | DGT Original code [s] | DGT LTFAT code [s] | DGT Improvement [%] | iDGT Original code [s] | iDGT LTFAT code [s] | iDGT Improvement [%] | Both Original code [s] | Both LTFAT code [s] | Both Improvement [%] |
|:--------------:|:------------------:|:---------------:|:---------------------:|:------------------:|:---------------:|:---------------------:|:------------------:|:---------------:|:---------------------:|
| 1 | 1.1229 | 0.9772 | 12.9740 | 1.4162 | 0.9792 | 30.8580 | 2.5186 | 1.9592 | 22.2130 |
| 2 | 1.0943 | 0.9786 | 10.5740 | 1.4228 | 0.9784 | 31.2350 | 2.5351 | 1.9738 | 22.1440 |
| 3 | 1.0967 | 0.9796 | 10.6750 | 1.4305 | 0.9795 | 31.5300 | 2.5370 | 1.9651 | 22.5410 |
| 4 | 1.0949 | 0.9792 | 10.5680 | 1.4299 | 0.9774 | 31.6450 | 2.5319 | 1.9712 | 22.1450 |
| 5 | 1.0993 | 0.9830 | 10.5840 | 1.4276 | 0.9766 | 31.5870 | 2.5330 | 1.9740 | 22.0680 |
| 6 | 1.0865 | 0.9755 | 10.2140 | 1.4254 | 0.9781 | 31.3820 | 2.5447 | 1.9725 | 22.4860 |
| 7 | 1.0937 | 0.9795 | 10.4400 | 1.4303 | 0.9778 | 31.6340 | 2.5387 | 1.9619 | 22.7200 |
| 8 | 1.0937 | 0.9756 | 10.8060 | 1.4220 | 0.9803 | 31.0630 | 2.5380 | 1.9665 | 22.5170 |
| 9 | 1.0925 | 0.9871 | 9.6460 | 1.4295 | 0.9788 | 31.5330 | 2.5359 | 1.9827 | 21.8150 |
| 10 | 1.1189 | 0.9781 | 12.5870 | 1.4315 | 0.9795 | 31.5750 | 2.5727 | 1.9844 | 22.8650 |
| **Average**  | 1.0993 | 0.9793 |**10.9068** | 1.4266 | 0.9786 | **31.4042** | 2.5386 | 1.9711 | **22.3514** |

- To compare the differences between the phase corrections a demo called `demo_phaseCor.m` was used.
- Similarly, to compare the differences between the time-directional variations a demo called `demo_timeVariation` was used. 
- To compare the differences between the projections a small demo called `demo_proj.m` was used.
- The SNR in the course of iterations; comparison between reconstruction with `updateInputCP = false` (blue) and `updateInputCP = true` (orange) :
<img width="1920" height="973" alt="diffWithWithoutUpdate" src="https://github.com/user-attachments/assets/0c148e33-5230-4809-bcf6-a035377c6256" />

### Notes
The `main.m` function in `U-PHAIN` folder can be used to process datasets in bulk with the default option being the DPAI dataset.

Even though LTFAT offers functions for the calculation of instantaneous frequency, the original function "calcInstFreq" from [1] was utilized due to its faster computational time.

The percentages change depending on the particular computer for the experiments.

Tests were run in Matlab 2025a on PC with AMD Ryzen 9 9900X, 4.4 GHz, 64 GB RAM and Windows 11. The Matlab codes use the Signal Processing Toolbox and the [LTFAT toolbox](https://github.com/ltfat/ltfat/releases) [3].

## References
[1] Tomoro Tanaka, Kohei Yatabe, and Yasuhiro Oikawa, “PHAIN: Audio inpainting via phase-aware optimization with instantaneous frequency,” IEEE/ACM Transactions on Audio, Speech, and Language Processing, Sep 2024.

[2] Tomoro Tanaka, “TomoroTanaka/PHAIN: v1.0.1”. Zenodo, Jun. 30, 2024. doi:10.5281/zenodo.12597334.

[3] Zdeněk Průša, Peter L. Søndergaard, Nicki Holighaus, Christoph Wiesmeyr, Peter Balazs, “The Large Time-Frequency Analysis Toolbox 2.0,” Sound, Music, and Motion, Lecture Notes in Computer Science 2014, pp 419-442.
