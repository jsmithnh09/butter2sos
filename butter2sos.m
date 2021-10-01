function [sos, k] = butter2sos(N, Fc, Fs, type, dir)
% Even or odd butterworth HPF/LPF generation. Since zeros are known, poles 
% should be generated and ordered properly without conditionally checking 
% for odd-man-out singularities. The motivation for the function is for 
% real-time filter design in C.
%
% SOS = BUTTER2SOS(N, Fc, Fs, [Type, Dir])
%   N (scalar) is an even order filter to generate.
%   Fc (scalar) is the corner -3 dB frequency.
%   Fs (scalar) is the sampling rate of the filter.
%   TYPE (char) indicates "l(owpass)", "h(ighpass)", or "a(llpass)".
%   DIR (char) indicates "(u)p" or "(d)own" stage ordering for LPF/HPFs.
%     By default, the ordering is in an up configuration.
%
% NOTES:
% With odd order, pop the last pole and place it's coefficients in
% the first stage. zeros are at +1 or -1 based on the LPF/HPF requirement.
%
% If odd order, the butterworth analog prototype pre-warping would have an
% additional pole located at (-1, 0).
%
% [1] J.G. Proakis and D.G. Manolakis, Digital Signal Processing, Prentice
%     Hall, 2007, chapter 10, section 3.
% [2] A.V. Oppenheim and R.W. Schafer, Digital Signal Processing, Prentice
%     Hall, 1975, chapter 5, sections 1 through 3.
% [3] J.O. Smith III, Digital Filters with Audio Applications, BookSurge
%     Publishing, 2007, appendix I.4.


if (nargin < 4) || (isempty(type))
  type = 'h';
end
if (nargin < 5) || (isempty(dir))
  dir = 'up';
end

z = [];                                         % butterworth are all poles, no zeroes.
phi = pi*(1:2:N-1)/(2*N) + pi/2;                % angular positions of pole on unit circle
Wn = 2*Fc/Fs;                                   % corner frequency fractional to Nyquist.
warp = 4*tan(pi/2*Wn);
p = cos(phi) + 1i.*sin(phi);                    % eulers method (faster C-code?)
p = [p; conj(p)];                               % adding conjugates
p = p(:);                                       % row-vector enforced.
if (mod(N,2))
  p = [p; -1];                                  % if odd order, dangling pole (S-plane, no conjugates.)
end
k = real(prod(-p));                             % higher orders will increase quantization error around 1.
Fs = 2;                                         % corner frequency now normalized to w = Nq.

if (type(1) == 'h')
  [z, p, k] = static_hpfwarp(z, p, k, warp);    % static functions are order independent thankfully.
else
  [z, p, k] = static_lpfwarp(z, p, k, warp);    % LPF or APF default
end
[~, p, k] = static_bilinear(z, p, k, Fs);       % S to Z-plane conversion, (zeros are already solved for butterworth.)
a = repmat([1,  0, 0], ceil(N/2), 1);           % pre-allocating the denominator.


lastpole = [];
if (mod(N,2))                                   % pop the dangling pole for the coefficient quantization loop.
  lastpole = p(end);
  p = p(1:end-1);
end

p = p(imag(p)>0);                               % isolate the poles from their conjugates.
p = flip(p);                                    % assuming in the "up" stage direction...
if (~isempty(lastpole))
  a(1,:) = [1, -lastpole, 0];                   % at zero and the warped dangling pole positions. Can safely assume real.
  offset = 1;                                   % we put the dangling pole in the first section.
else
  offset = 0;                                   % no offset needed for an additional pole.
end
for iP = 1:length(p)
  a(iP+offset,2) = -2*real(p(iP));              % complex conjugate expansion, (x-z1)(z-z1*)
  a(iP+offset,3) = real(p(iP) .* conj(p(iP)));  % roots/poly functions may be slow in comparison to a direct evaluation.
end

if (type(1) == 'h')                             % given that the zeros are all +/-1, we can save a poly call here.
  b = repmat([1, -2, 1], ceil(N/2), 1);         % zeros at +1
  if (~isempty(lastpole))
    b(1,:) = [1, -1, 0];
  end
elseif (type(1) == 'l')
  b = repmat([1,  2, 1], ceil(N/2), 1);         % zeros at -1
  if (~isempty(lastpole))
    b(1,:) = [1, +1, 0];
  end
elseif (type(1) == 'a')
  b = fliplr(a);                                % if allpass, reverse the pole polynomial for the zeros to achieve unity gain but maintain phase.
  if (~isempty(lastpole))
    b(1,:) = [a(1,2), a(1,1), 0];               % dangling pole flip in first stage.
  end
  sos = [b, a];                                 % no gain incorporation, otherwise the unity gain is lost. normalize by a2(?)
  return;
end
sos = [b, a];                                   % form the overall SOS structure and normalize the LPF/HPF types.
if (strcmp(dir(1), 'd'))                        % currently in an "up" configuration, can flip for down ordering.
  sos = flip(sos,1);
end
if (nargout == 1)
  sos(1,1:3) = sos(1,1:3)*k;                    % encorporate the gain.
end
return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zl, pl, kl] = static_lpfwarp(z, p, k, w)
% lowpass is already "fitting" prototype, need to scale to the frequency
% axis defined by the warp frequency "w".
ord = length(p) - length(z);
zl = w .* z;
pl = w .* p;
kl = k .* w^(ord); % total gain change from shifting zeros and poles cancelled here.
return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zh, ph, kh] = static_hpfwarp(z, p, k, w)
% for highpass, invert singularity positions since we need the 
% HPF counterpart from the original LPF analog filter prototype.
ord = length(p) - length(z);
if (~isempty(z)) % since butterworth is all-pole, no zeros...
  zh = w ./ z;
else
  zh = [];
end
if (ord > 0) % zeros moved back to origin if any at Inf.
  zh = [zh; zeros(ord, 1)];
end
ph = w ./ p;
kh = k * real(prod(-z)/prod(-p)); % gain change from inversion.
return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zd, pd, kd] = static_bilinear(z, p, k, Fs)
% bilinear ZPK, H(z) = H(s) where s = 2*Fs*(z+1)/(z-1)
ord = length(p) - length(z);
Fs = 2*Fs; % 2/T equivalent in bilinear xform.
pd = (1+p/Fs)./(1-p/Fs);
zd = (1+z/Fs)./(1-z/Fs);
kd = k * real(prod(Fs-z)./prod(Fs-p));
if (ord > 0) % zeros at infinity moved to Nyq, (-1, 0.)
  zd = [zd; -ones(ord, 1)];
end
return;
end
