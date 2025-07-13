% MIT License
% 
% Copyright (c) 2023 TomoroTanaka
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

%   [1] Tanaka, Tomoro, Kohei Yatabe, and Yasuhiro Oikawa,"PHAIN: Audio
%       inpainting via phase-aware optimization with instantaneous frequency,"
%       IEEE/ACM Transactions on Audio, Speech, and Language Processing, Sep 2024.
%
% this is a copy of the Matlab function shortenForDGT.m from [1]

function [firstIdx, L] = shortenForDGT(w, a, s, f)

    center = ceil((s+f)/2);
    diff = 1 + floor((center - 1)/a)*a + ceil(a/2);
    offset = center - diff;
    offset = mod(offset, a);

    firstWin = ceil((s - ceil(w/2))/a) + 1;
    firstWinPos = 1 + (firstWin - 1)*a;
    firstWinPos = firstWinPos + offset;

    if firstWinPos - a + ceil(w/2) >= s + 1
        firstWin = firstWin - 1;
        firstWinPos = firstWinPos - a;
    end
    
    firstIdx = firstWinPos - ceil(floor(w/2)/a)*a;

    lastWin = firstWin + floor((f + floor(w/2) - firstWinPos)/a);
    lastWinPos = firstWinPos + (lastWin - firstWin)*a;

    lastIdx = lastWinPos + ceil(floor(w/2)/a)*a;
    
    L = lastIdx - firstIdx + 1;

end