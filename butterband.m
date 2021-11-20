function sos = butterband(N, w1, w2, Fs, type)
    % BUTTERBAND designs bandpass/bandstop Butterworth filters.
    %   SOS = BUTTERBAND(N, W1, W2, FS, TYPE)
    %       N (scalar) is the order of the filter. This will match the number of biquad stages in the SOS matrix.
    %       W1 (scalar) is the lower corner frequency closer to DC.
    %       W2 (scalar) is the upper corner frequency towards Nyquist.
    %       Fs (scalar) is the discrete sampling rate of the digital filter.
    %       TYPE (char) indicates bandpass or bandstop.
    %
    % Copyright Jordan R. Smith
    % See Also: BUTTER2SOS
    
    if (nargin == 0)
        help('butterband.m');
        clear('sos'); % prevents overwriting "ans".
        return;
    end

    % seed the poles on the S-plane.
    phi = pi*(1:2:N-1)/(2*N) + pi/2;
    p = cos(phi) + 1i.*sin(phi); % euler
    p = [p; conj(p)];
    p = p(:);
    z = [];
    k = 1;
    if (mod(N,2))
        p = [p; -1];
    end

    % normalize the corner frequencies.
    Wn1 = 2*w1/Fs;
    Wn2 = 2*w2/Fs;

    wp1 = 4*tan(pi/2 * Wn1); % warped points
    wp2 = 4*tan(pi/2 * Wn2);

    bw = wp2 - wp1;
    Wn = sqrt(wp1 * wp2);
    Fs = 2; % normalized Wq = 1, resetting Fs.

    % prototype filter transformation.
    if (regexp(type, 'pass'))
        [zb, pb, kb] = static_xform_bp(z, p, k, bw, Wn);
    else
        [zb, pb, kb] = static_xform_bs(z, p, k, bw, Wn);
    end

    % bilinear transformation
    [zd, pd, kd] = static_bilinear(zb, pb, kb, Fs);
    pd = flipud(pd);

    % sort the poles/zeros together.
    pd = static_polesort(pd, N);

    % same stage count as the order...
    sos = zeros(N, 6);
    sos(:, [1, 4]) = 1; % a0/b0 are unity.

    % quantize the poles and zeros
    for iStg = 1:N
        z0 = zd(2*iStg-1);
        z1 = zd(2*iStg);
        p0 = pd(2*iStg-1);
        p1 = pd(2*iStg);
        
        if ((z0 == 1) && (z1 == 1)) % +1/+1 pair.
            sos(iStg, 1:3) = [1, -2, 1];
        elseif ((z0 == -1) && (z1 == -1)) % -1/-1 pair.
            sos(iStg, 1:3) = [1, +2, 1];
        elseif (imag(z0) ~= 0) % complex zeros.
            sos(iStg, 1:3) = [1, -2*real(z0), real(z0)*conj(z0)];
        elseif (sign(z1) ~= sign(z0)) % +1/-1 pair.
            sos(iStg, 1:3) = [1, 0, -1];
        end

        if (imag(p0) ~= 0) % complex conjugate poles.
            sos(iStg, 4:6) = [1, (-p0)+(-p1), (-p0)*(-p1)]; % second-order poly expansion.
        end
    end
    sos(1,1:3) = sos(1,1:3) * kd; % incorporate the gain in the first stage.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zbp, pbp, kbp] = static_xform_bp(z, p, k, bw, Wn)
% XFORM_BP transforms prototype to bandpass filter.
    order = length(p) - length(z);
    z1 = z * bw/2;
    p1 = p * bw/2;

    pbp = zeros(2*length(p1), 1);
    zbp = [z1; zeros(order, 1)]; % appending zeros at origin to manage order diff.

    for iSing = 1:length(p1)
        ind1 = 2*iSing-1;
        ind2 = 2*iSing;
        pbp(ind1) = p1(iSing) + sqrt(p1(iSing)^2 - Wn^2);
        pbp(ind2) = p1(iSing) - sqrt(p1(iSing)^2 - Wn^2);
    end
    kbp = k * bw^(order);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zbs, pbs, kbs] = static_xform_bs(z, p, k, bw, Wn)
% XFORM_BS transforms prototype to bandstop.
    
    z1 = (bw/2) ./ z; % similar inversion as the HPF counterpart.
    p1 = (bw/2) ./ p;

    pbs = zeros(2*length(p1), 1);
    zbs = zeros(2*length(p1), 1);

    for iSing = 1:length(p1)
        ind1 = 2*iSing - 1;
        ind2 = 2*iSing;
        pbs(ind1) = p1(iSing) + sqrt(p1(iSing)^2 - Wn^2);
        pbs(ind2) = p1(iSing) - sqrt(p1(iSing)^2 - Wn^2);

        zbs(ind1) = 0 + 1i*Wn;
        zbs(ind2) = 0 - 1i*Wn;
    end
    kbs = k * real(prod(-z) ./ prod(-p));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = static_polesort(p, N)
% STATIC_POLESORT sorts poles based on proximity to the unit circle, (absolute distance.)
% In another filter-type design case, we'd need to sort the poles based on proximity to other
% zeros as well, but since butterworth is all-pole, we can simplify this process a bit.

temp = [];
if (mod(N, 2))
    range = length(p)-2; % last two conjugates are real in odd order case.
else
    range = length(p); % all complex; even order.
end

% sort based on decreasing distance from the unit circle.
for iter = 1:range
    for comp = 1:range-iter
        if (abs(p(comp) < abs(p(comp+1))))
            temp = p(comp);
            p(comp) = p(comp+1);
            p(comp+1) = temp;
        end
    end
end

% swap the real poles in the case of an odd order bandpass
if ((mod(N, 2)) && (p(range+1) < p(range+2)))
    temp = p(range+1);
    p(range+1) = p(range+2);
    p(range+2) = temp;
end

% flip the imaginary sign if the sign is different for complex conjugates.
for iter = 1:range/2
    if (imag(p(2*iter-1)) == 0)
        continue;
    end
    if (imag(p(2*iter)) > 0)
        temp = p(2*iter);
        p(2*iter) = p(2*iter-1);
        p(2*iter-1) = temp;
    end
end

end % static_polesort







        



