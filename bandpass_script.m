%% PARAMETERS
clear;
N = 4;
W1 = 300;
W2 = 2500;
Fs = 44100;


%% butterworth pole seeding
phi = pi*(1:2:N-1)/(2*N) + pi/2;
p = cos(phi) + 1i.*sin(phi);
p = [p; conj(p)];
p = p(:);
if (mod(N, 2))
  p = [p; -1];
end
z = [];                   % no zeros in butterworth.
k = 1;
p = p(:);
Wn1 = 2*W1/Fs;            % normalize to nyquist
Wn2 = 2*W2/Fs;

wp1 = 4*tan(pi/2 * Wn1);  % warped points
wp2 = 4*tan(pi/2 * Wn2);

bw = wp2 - wp1;           % normalized bandwidth
Wn = sqrt(wp1 * wp2);     % center frequency 
Fs = 2;                   % normalized Wq = 1, resetting Fs.

%% bandpass transformation.
order = length(p) - length(z);
zl = z * bw/2;
pl = p * bw/2;

pbp = zeros(2*length(pl), 1);
zbp = [z; zeros(order, 1)];

for i = 1:length(pl)
  ind1 = 2*i - 1;
  ind2 = 2*i;
  pbp(ind1) = pl(i) + sqrt(pl(i)^2 - Wn^2);
  pbp(ind2) = pl(i) - sqrt(pl(i)^2 - Wn^2);
end
kbp = k * bw^(order);

%% bilinear
ft = 2*Fs; % 2/T equivalent in bilinear xform.
pd = (1+pbp/ft)./(1-pbp/ft);
zd = (1+zbp/ft)./(1-zbp/ft);
if (order > 0) % zeros at infinity moved to Nyq, (-1, 0.)
  zd = [zd; -ones(order, 1)];
end
kd = kbp * real(prod(ft-zbp)./prod(ft-pbp));

%% pole sorting
% pole ordering needs to be sorted in proximity to the unit circle,
% (descending distance.)
%
% if odd order, we KNOW the last two poles are real. we can limit the
% sorting to 1:p-2.
%
% Since the p +- sqrt(p^2 - Wn^2) generates real pole pairs at the end of
% the vector in the case of an odd order bandpass, we can limit the range
% of the sort.

temp = [];
if (mod(N, 2))
  range = length(pd)-2;
else
  range = length(pd);
end

for iter = 1:range
  for comp = 1:range-1
    if (real(pd(comp)) < real(pd(comp+1))) % swap point.
      temp = pd(comp);
      pd(comp) = pd(comp+1);
      pd(comp+1) = temp;
    end
  end
end

% swap the real poles in the case of an odd order bandpass.
if (mod(N, 2)) && (pd(range+1) < pd(range+2))
  temp = pd(range+1);
  pd(range+1) = pd(range+2);
  pd(range+2) = temp;
end

% flips the imaginary component if the sign is different for complex
% conjugates.
for iter = 1:range/2
  if (imag(pd(2*iter-1) == 0))
    continue;
  end
  if (imag(pd(2*iter)) > 0)
    temp = pd(2*iter);
    pd(2*iter) = pd(2*iter-1);
    pd(2*iter-1) = temp;
  end
end

%% SOS formation
sos = repmat([1 0 0 1 0 0], N, 1);
for iStg = 1:N
  if ((zd(2*iStg-1) == 1) && (zd(2*iStg) == 1)) % +1/+1 pair.
    sos(iStg,1:3) = [1 -2 1];
  elseif ((zd(2*iStg-1) == -1) && (zd(2*iStg) == -1)) % -1/-1 pair.
    sos(iStg,1:3) = [1, +2, 1];
  elseif (sign(zd(2*iStg)) ~= sign(zd(2*iStg-1))) % -1/+1 zero pair.
    sos(iStg,1:3) = [1, 0, -1];
  end
  if (imag(pd(2*iStg-1)) ~= 0) % real pole pairing.
    sos(iStg,4:6) = [1, -2*real(pd(2*iStg-1)), pd(2*iStg-1)*conj(pd(2*iStg-1))];
  else
    z1 = pd(2*iStg-1);
    z2 = pd(2*iStg);
    sos(iStg,4:6) = [1, -z1+-z2, -z1 * -z2]; % second order poly expansion.
  end
end
sos(N,1:3) = sos(N,1:3) * kd;
  
