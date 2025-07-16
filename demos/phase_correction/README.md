# Instantaneous phase correction of complex spectrogram

This folder contains a set of MATLAB codes written as a supplementary material of the following tutorial paper explaining phase-related property of spectrograms:

- Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa, "**Representation of complex spectrogram via phase conversion**," *Acoustical Science and Technology*, vol.40, no.3, May 2019. (Open Access)

## Features

The codes provide some simple methods to perform the following things:

1. Calculating spectrogram of a real-valued signal;

2. Designing a synthesis window for perfectly reconstructing the signal;

3. Calculating the inverse transform reconstructing the signal from its spectrogram;

4. Calculating instantaneous frequency at each time-frequency bin;

5. Obtaining derivative of a window for computing instantaneous frequency; and

6. Calculating instantaneous-phase-corrected (iPC) spectrogram explained in the above paper.

## Demo

To demonstrate the usage of included functions, the following files are provided:

- **Executable demo codes**
    - demo1_DGTusage.m
    - demo2_windowUsage.m
    - demo3_IFandIPC.m
    - demo_FDGTusage.m

## Functions

The following MATLAB functions are included in the folder "DGT+iPC":

- **Window function utilities**
    - generalizedCosWin.m
    - calcCanonicalDualWindow.m
    - calcCanonicalTightWindow.m
    - numericalDiffWin.m
- **Discrete Gabor transform**
    - DGT.m
    - invDGT.m
    - zeroPaddingForDGT.m
- **Instantaneous phase correction**
    - calcInstFreq.m
    - instPhaseCorrection.m
    - invInstPhaseCorrection.m

In addition to these functions, a faster implementation of DGT and its inverse is also included in the subfolder for practical use:

- **Faster DGT implementation**
    - FDGT.m
    - invFDGT.m
    - precomputationForFDGT.m

It should be easier to open and read the demo files first as they illustrate how to use these functions. For detailed explanation and usage of each function, please open the file or type `help functionName` into the command line (for example, `help DGT`). For technical explanation, please read the above tutorial paper.

## Execution in Code Ocean

In Code Ocean, all figures are saved by the special function `codeOceanFigSave.m`. In your own MATLAB, these figures are not saved but appear in figure windows as usual.

## References

In addition to the associated tutorial paper above, some applications of the instantaneous-phase-corrected (iPC) spectrogram can be found in the papers below. We are happy if our codes and papers are helpful for you and if you cite some of them :)

- K. Yatabe and Y. Oikawa, "**Phase corrected total variation for audio signals**, *IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP2018)*, Apr. 2018.
- Y. Masuyama, K. Yatabe and Y. Oikawa, "**Low-rankness of complex-valued spectrogram and its application to phase-aware audio processing**," *IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP2019)*, May 2019.
- Y. Masuyama, K. Yatabe and Y. Oikawa, "**Phase-aware harmonic/percussive source separation via convex optimization**," *IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP2019)*, May 2019.