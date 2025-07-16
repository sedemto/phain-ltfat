function varargout = generalizedCosWin(N, coeffs)
% generalizedCosWin: Collection of generalized cosine windows with derivatives.
%   
%   Usage:
%      win = generalizedCosWin(N);
%      win = generalizedCosWin(N,coeffs);
%      [win,diffWin] = generalizedCosWin(___);
%      generalizedCosWin(___);
%   
%   Input parameters:
%      N       : Window length.
%      coeffs  : Characters specifying the type of window (row vector).
%               (or may be a raw coefficient vector, see examples below)
%   
%   Output parameters:
%      win     : Specified window function.
%      diffWin : Derivative of the specified window.
%   
%   "[win,diffWin] = generalizedCosWin(N)" returns the Hann window and
%   its derivative. Other generalized cosine windows can be obtained by
%   specifying its name as "generalizedCosWin(N,'windowName')". Acceptable
%   window names are the followings:
%   
%      'rect' - returns the rectangular window.
%         The rectangular window is the simplest window whose non-zero
%         values are unity. Its derivative is zero.
%   
%      'hann' or 'hanning' - returns the Hann window.
%         The Hann window is a 2-term cosine window whose coefficients are
%         chosen so that it is continuous at the edges.
%   
%      'hamming' - returns the Hamming window.
%         The Hamming window is a 2-term cosine window whose coefficients
%         are chosen so that the first sidelobe is attenuated [2].
%   
%      'blackman' - returns the Blackman window.
%         The Blackman window is a 3-term cosine window whose coefficients
%         are obtained by truncating coefficients of the exact Blackman
%         window in order to remove discontinuity at the edges [2].
%   
%      'exactblackman' - returns the exact Blackman window.
%         The exact Blackman window is a 3-term cosine window whose
%         coefficients are chosen so that the third and fourth sidelobes
%         are attenuated [2].
%   
%      'blackmanharris' or 'blackmanharris4term' - returns the (4-term) Blackman-Harris window.
%         The Blackman-Harris window is a 4-term cosine window whose
%         coefficients are chosen to minimize the peak sidelobe level [2].
%         Improved versions of this window were later proposed by Nuttall [3].
%   
%      'blackmanharris3term' - returns the 3-term Blackman-Harris window.
%         The 3-term Blackman-Harris window is a 3-term cosine window whose
%         coefficients are chosen to minimize the peak sidelobe level [2].
%   
%      'nuttall3term', 'nuttall3termC1', 'nuttall3termC3', 'nuttall4term',
%      'nuttall4termC1', 'nuttall4termC3', 'nuttall4termC5' - return the Nuttall windows.
%         The Nuttall windows are windows whose coefficients are chosen to
%         minimize the peak sidelobe level [3]. The following variations
%         are acceptable in this code.
%   
%            'nuttall3term'  : the 3-term Nuttall window.
%            'nuttall3termC1': the 3-term Nuttall window with continuous first derivative.
%            'nuttall3termC3': the 3-term Nuttall window with continuous third derivative.
%            'nuttall4term'  : the 4-term Nuttall window.
%            'nuttall4termC1': the 4-term Nuttall window with continuous first derivative.
%            'nuttall4termC3': the 4-term Nuttall window with continuous third derivative.
%            'nuttall4termC5': the 4-term Nuttall window with continuous fifth derivative.
%   
%      'flattop' - returns the flattop window.
%         The flattop window is a 5-term cosine window whose mainlobe is
%         flattened for reducing scalloping losses.
%   
%      'kawahara5term' - returns Kawahara's 5-term window.
%         Kawahara's 5-term window is a 5-term cosine window with continuous
%         fifth derivative, where its coefficients are chosen to minimize
%         the peak sidelobe level [4].
%   
%      'kawahara6term' - returns Kawahara's 6-term window.
%         Kawahara's 6-term window is a 6-term cosine window with continuous
%         seventh derivative, where its coefficients are chosen to minimize
%         the peak sidelobe level [4].
%   
%      'albrecht7term' - returns Albrecht's 7-term minimum sidelobe window.
%         Albrecht's 7-term minimum sidelobe window is a 7-term cosine
%         window whose coefficients are chosen to minimize the peak sidelobe
%         level (with no edge condition) [5].
%   
%      'albrecht11term' - returns Albrecht's 11-term minimum sidelobe window.
%         Albrecht's 11-term minimum sidelobe window is an 11-term cosine
%         window whose coefficients are chosen to minimize the peak sidelobe
%         level (with no edge condition) [5].
%   
%   One can directly specify the coefficients of the cosine series as
%   "generalizedCosWin(N,[0.5,0.5])" (which results in the Hann window).
%   This code is associated with "DGT.m" implemented as a supporting
%   material of [1], where all windows were implemented so that "DGT.m"
%   can make the window zero-phased when "zeroPhaseFlag" is "true".
%   
%   This function plots the specified window and its derivative when the
%   output argument is omitted.
%   
%   Example 1: Generating the 128-point 4-term Nuttall window with
%      continuous first derivatives by its name.
%   
%      [win,diffWin] = generalizedCosWin(128, 'nuttall4termc1');
%   
%   Example 2: Generating the 128-point 4-term Nuttall window with
%      continuous first derivative by specifying its coefficients.
%   
%      coeffs = [0.355768, 0.487396, 0.144232, 0.012604];
%      [win,diffWin] = generalizedCosWin(128, coeffs);
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)
%   
%   [2] Fredric J. Harris, "On the use of windows for harmonic analysis
%       with the discrete Fourier transform," Proceedings of the IEEE,
%       vol.66, no.1, pp.51-83, Jan. 1978.
%   
%   [3] Albert H. Nuttall, "Some windows with very good sidelobe behavior,"
%       IEEE Transactions on Acoustics, Speech, and Signal Processing,
%       vol.29, no.1, pp.84-91, Feb. 1981.
%   
%   [4] Hideki Kawahara, Ken-Ichi Sakakibara, Masanori Morise, Hideki
%       Hideki Banno, Tomoki Toda and Toshio Irino, "A new cosine series
%       antialiasing function and its application to aliasing-free glottal
%       source models for speech and singing synthesis," Proceedings of
%       INTERSPEECH 2017, pp.1358-1362, Aug. 2017.
%   
%   [5] Hans-Helge Albrecht, "A family of cosine-sum windows for
%       high-resolution measurements," Proceedings of IEEE International
%       Conference on Acoustices, Speech, and Signal Processing (ICASSP), 
%       pp.3081-3084, May 2001.

%   Author: Tsubasa Kusano (2019)


narginchk(1,2);
nargoutchk(0,2);

if (length(N)~=1)
    error('N must be a scalar.')
end
if N <= 0
    error('N must be a positive number.')
end
if nargin == 1
    coeffs = 'hann';
    warning('Window type is not specified. The Hann window is output as default.')
end
if ~isnumeric(coeffs)
    coeffs = getcoeffs(lower(coeffs));
end

t = (2*pi*((-N/2:N/2-1)+mod(N,2)/2)/N).';
win = cos((0:length(coeffs)-1).*t) * coeffs(:);

if nargout == 1
    varargout{1} = win;
else  
    diffcoeffs = (0:length(coeffs)-1).' .* coeffs(:);
    diffWin = -sin((0:length(coeffs)-1).*t) * diffcoeffs;
    
    if nargout == 0
        disp_window(win, diffWin);
    else
        varargout{1} = win;
        varargout{2} = diffWin;
    end
end
end



function coeffs = getcoeffs(str)
switch str
    case 'rect'
        coeffs = 1;
    case {'hann', 'hanning'}
        coeffs = [0.5, 0.5];
    case 'hamming'
        coeffs =  [25, 21]./46;
    case 'blackman'
        coeffs = [0.42, 0.5, 0.08];
    case 'exactblackman'
        coeffs = [7938, 9240, 1430]./18608;
    case {'blackmanharris', 'blackmanharris4term'}
        coeffs = [0.35875, 0.48829, 0.14128, 0.01168];
    case 'blackmanharris3term'
        coeffs = [0.42323, 0.49755, 0.07922];
    case 'nuttall'
        error(['Choose one from the following Nuttall-type windows: \n\n', ...
               '   nuttall3term   : the 3-term Nuttall window \n',...
               '   nuttall3termc1 : the 3-term Nuttall window with continuous first derivative \n'...
               '   nuttall3termc3 : the 3-term Nuttall window with continuous third derivative \n',...
               '   nuttall4term   : the 4-term Nuttall window \n',...
               '   nuttall4termc1 : the 4-term Nuttall window with continuous first derivative \n',...
               '   nuttall4termc3 : the 4-term Nuttall window with continuous third derivative \n',...
               '   nuttall4termc5 : the 4-term Nuttall window with continuous fifth derivative \n',...
               ],[])
    case 'nuttall3term'
        coeffs = [0.4243801, 0.4973406, 0.0782793];
    case 'nuttall3termc1'
        coeffs = [0.40897, 0.5, 0.09103];
    case 'nuttall3termc3'
        coeffs = [3, 4, 1]/8;
    case 'nuttall4term'
        coeffs = [0.3635819, 0.4891775, 0.1365995, 0.0106411];
    case 'nuttall4termc1'
        coeffs = [0.355768, 0.487396, 0.144232, 0.012604];
    case 'nuttall4termc3'
        coeffs = [0.338946, 0.481973, 0.161054, 0.018027];
    case 'nuttall4termc5'
        coeffs = [10, 15, 6, 1]/32;
    case 'flattop'
        coeffs = [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368];
    case 'kawahara5term'
        coeffs = [0.2940462892, 0.4539870314, 0.2022629686, 0.0460129686, 0.0036907422];
    case {'kawahara6term', '6term'}
        coeffs = [0.2624710164, 0.4265335164, 0.2250165621, 0.0726831633, 0.0125124215, 0.0007833203];
    case {'albrecht7term', '7term'}
        coeffs = [2.712203605850388e-1, 4.334446123274422e-1, 2.180041228929303e-1, ...
                  6.578534329560609e-2, 1.076186730534183e-2, 7.700127105808265e-4, ...
                  1.368088305992921e-5];
    case {'albrecht11term', '11term'}
        coeffs = [2.151527506679809e-1, 3.731348357785249e-1, 2.424243358446660e-1, ...
                  1.166907592689211e-1, 4.077422105878731e-2, 1.000904500852923e-2, ...
                  1.639806917362033e-3, 1.651660820997142e-4, 8.884663168541479e-6, ...
                  1.938617116029048e-7, 8.482485599330470e-10];
    otherwise
        error(['Window type', ' ''', str, ''' ', 'is not included in this code. \n', ...
               'Please choose one of the following windows: \n\n', ...
               ' - ''rect'' \n', ...
               ' - ''hann'' (or ''hanning'') \n', ...
               ' - ''hamming'' \n', ...
               ' - ''blackman'' \n', ...
               ' - ''exactblackman'' \n', ...
               ' - ''blackmanharris'' (or ''blackmanharris4term'') \n', ...
               ' - ''blackmanharris3term'' \n', ...
               ' - ''nuttall3term'' \n', ...
               ' - ''nuttall3termc1'' \n', ...
               ' - ''nuttall3termc3'' \n', ...
               ' - ''nuttall4term'' \n', ...
               ' - ''nuttall4termc1'' \n', ...
               ' - ''nuttall4termc3'' \n', ...
               ' - ''nuttall4termc5'' \n', ...
               ' - ''flattop'' \n', ...
               ' - ''kawahara5term'' \n', ...
               ' - ''kawahara6term'' (or ''6term'') \n', ...
               ' - ''albrecht7term'' (or ''7term'') \n', ...
               ' - ''albrecht11term'' (or ''11term'') \n'],[])
end
end



function disp_window(win, diffWin)
min_nfft = 2^12;
nfft = max(length(win), min_nfft/2)*2;

h = figure('Units', 'normalized', 'Position', [.2, .2, .6, .6]);
subplot(2,2,1)
plot(0:length(win)-1, win, 'LineWidth', 1)
xlabel('Samples')
title('Window')
xlim([0, length(win)-1])

[specLeakage, bandwidth, hsl, hsl_f, lsl] = comp_windowcharacteristics(win,nfft);

[spec, freq] = freqz(win,1,nfft);
subplot(2,2,2)
plot(freq, 20*log10(abs(spec)), 'LineWidth', 1)
MSpec = max(20*log10(abs(spec)));
if isempty(specLeakage)
    str = {['  Spectral leakage: -%'], ['Bandwidth (-3 dB): -']};
    annotation('textbox',[0.75 0.8 0.1 0.1],'String',str, 'FitBoxToText','on')
else
    hold on, plot([0, hsl_f, hsl_f], [MSpec MSpec hsl+MSpec], 'r')
    
    plot([0, hsl_f+0.2], [1,1] * hsl + MSpec, 'r')
    text(hsl_f, MSpec+hsl, [num2str(hsl, '%2.2f'), 'dB'], 'VerticalAlignment', 'bottom')
    str = {['  Spectral leakage: ', num2str(specLeakage, '%2.4f'), '%'], ['Bandwidth (-3 dB): ', num2str(bandwidth, '%2.4f')]};
    annotation('textbox',[0.75 0.8 0.1 0.1],'String',str, 'FitBoxToText','on')
end
xlim([0 pi])
ylim([20*(ceil((lsl+MSpec)/20)-1), 20*ceil(MSpec/20)])
xlabel('Normalized frequency [rad/sample]')
ylabel('Amplitude [dB]')
xticks((0:4)/4 * pi)
xticklabels({'0','1/4 \pi','2/4 \pi','3/4 \pi','\pi'})

subplot(2,2,3)
plot(0:length(win)-1,diffWin, 'LineWidth', 1)
xlabel('Samples')
title('Differential window')
xlim([0, length(win)-1])

subplot(2,2,4)
diffspec = freqz(diffWin,1,nfft);
plot(freq, 20*log10(abs(diffspec)), 'LineWidth', 1)
xlim([0 pi])
if any(diffspec) % check rectangle or otherwise
    ylim([20*(ceil((lsl+MSpec)/20)-1), 20*ceil(max(20*log10(abs(diffspec)))/20)])
end
xlabel('Normalized frequency [rad/sample]')
ylabel('Amplitude [dB]')
xticks((0:4)/4 * pi)
xticklabels({'0','1/4 \pi','2/4 \pi','3/4 \pi','\pi'})

h.Units = 'pixels';
end



function [specLeakage, bandwidth, hsl, hsl_f, lsl] = comp_windowcharacteristics(w,nfft)
[h, f] = freqz(w, 1, nfft);
logSpec = 20*log10(abs(h));

% computation of bandwidth
r = logSpec-logSpec(1) + 3;
ind_minus3dB = find(r>=0, 1, 'last' );
bandwidth = (f(2)-f(1)) * (ind_minus3dB-1 + r(ind_minus3dB)./(r(ind_minus3dB) - r(ind_minus3dB+1)));
bandwidth = 2*bandwidth;

% computation of highest sidelobe level
d = diff(logSpec);
ind_peaks = find(diff(sign(d))==-2) + 1;

% Peak Amplitudes
logSpec_peaks = logSpec(ind_peaks);

% Reject mainlobe peak if present, keep highest peak
[hsl, ind_hsl] = max(logSpec_peaks(logSpec_peaks<logSpec(1)));
hsl = hsl - logSpec(1);
hsl_f = f(ind_peaks(ind_hsl));
lsl = min(logSpec_peaks(logSpec_peaks<logSpec(1)));
lsl = lsl - logSpec(1);

% computation of spectral leakage
if ~isempty(ind_minus3dB)
  firstzero = find(d(ind_minus3dB:end)>0);
  if isempty(firstzero)
      specLeakage = [];
  else
    firstzero = firstzero(1)+ind_minus3dB-1;
    specLeakage = 1 - (sum(abs(h(1:firstzero).^2))/sum(abs(h.^2)));
    specLeakage = 1e2*specLeakage; % percent
  end
end
end