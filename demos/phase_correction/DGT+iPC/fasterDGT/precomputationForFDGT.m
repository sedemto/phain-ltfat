function [sigIdx,sumIdx,sumArray,ifftArray,rotIdx] = precomputationForFDGT(sigLen,winLen,shiftLen,fftLen)
winIdx   = 0:winLen-1;          % row vector
winIdx   = winIdx(:);           % column vector
fftIdx   = 0:fftLen-1;          % row vector
fftIdx   = fftIdx(:);           % column vector
idxShift = 0:shiftLen:sigLen-1; % row vector

sigIdx = mod(winIdx + idxShift - floor(winLen/2), sigLen); % 2D index matrix
sigIdx = sigIdx + 1; % MATLAB-type index (starting from 1)

rotIdx = mod(fftIdx - idxShift, fftLen) + fftLen*(0:length(idxShift)-1);
rotIdx = rotIdx + 1; % MATLAB-type index (starting from 1)

timeLen = size(sigIdx,2);
lap = ceil(winLen/shiftLen); % number of overlapping windows
can = lap:timeLen;           % candidates of integer ">= winLen/shiftLen"
rowLen = min(can(mod(timeLen-1,can)>=lap-1)); % best integer in candidates

ifftArray = zeros(fftLen,timeLen); % memory space for inverse FFT
sumArray = zeros(sigLen,rowLen);   % memory space for overlap add

rowIdx = mod(0:timeLen-1,rowLen);
rowIdx = rowIdx + 1; % MATLAB-type index (starting from 1)
rowIdx = repmat(rowIdx,winLen,1);

sumIdx = sub2ind(size(sumArray),sigIdx(:),rowIdx(:)); % linear index for overlap add
end