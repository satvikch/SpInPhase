function [R_fil, B_fil, D_fil] = nlinsar(z1, z2, ...
                                         R_est, B_est, D_est, ...
                                         hW, hD, ...
                                         Lmin)
%NLINSAR estimator for pairs of SLC SAR images
%   R = NLINSAR(Z1, Z2) estimates reflectivity, coherence and
%   actual interferometric phase from a pair of coregistred SLC SAR
%   data as described in:
%
%       Charles-Alban Deledalle, Loïc Denis and Florence Tupin,
%       NL-InSAR : Non-Local Interferogram Estimation,
%       IEEE Trans. on Geoscience and Remote Sensing (to appear)
%
%   Please refer to these papers for a more detailed description of
%   the arguments. Note that this function ables also to treat
%   large images by preserving memory thanks to a block processing
%   on 256x256 subimages.%
%
%       ARGUMENT DESCRIPTION:
%               Z1       - First SLC SAR noisy image
%               Z2       - Second SAR noisy image
%               R_EST    - Image of pre-estimated reflectivities
%                          (default: image of values equal to 1)
%               B_EST    - Image of pre-estimated phases
%                          (default: image of values equal to 0)
%               D_EST    - Image of pre-estimated coherences
%                          (default: image of values equal to 0)
%               HW       - Half sizes of the search window width
%                          (default 10)
%               HD       - Half sizes of the  window width
%                          (default 3)
%               LMIN     - Minimum equivalent number of looks
%
%       OUTPUT DESCRIPTION:
%               R_FIL    - Image of estimated reflectivities
%               B_FIL    - Image of estimated actual
%                          interferommetric phase
%               D_FIL    - Image of estimated coherences
%
%       AUTHOR:
%               Charles Deledalle
%               email: deledalle@telecom-paristech.fr
%
%       VERSION: 23 August 2010

    if nargin < 4
        R_est = ones(size(z1));
        D_est = zeros(size(z1));
        B_est = zeros(size(z1));
    end
    if nargin < 5
        hW  = 10;
    end
    if nargin < 6
        hD  = 3;
    end
    if nargin < 7
        Lmin  = 10;
    end
    t = 2 * hW + 1;
    w = 2 * hD + 1;

    switch w
      case 3
        h_theo = 0.58; %alpha = 0.92
      case 7
        h_theo = 0.2440; %alpha = 0.92
      otherwise
        h_theo = quantile_insar(w, 0.92);
    end

    T_theo  = 0.2;
    h = h_theo .* w.^2;
    T = T_theo ./ h  * pi / 4 .* w.^2;

    width = size(z1, 1);
    height = size(z1, 2);
    sw = 256;
    sh = 256;
    overlap = hW + hD;

    z_est = R_est .* exp(sqrt(-1) * B_est);
    z1(abs(z1) <= 0) = ...
        min2(abs(z1(abs(z1) > 0)));
    z2(abs(z2) <= 0) = ...
        min2(abs(z2(abs(z2) > 0)));
    z1(isnan(abs(z1))) = ...
        min2(abs(z1(abs(z1) > 0)));
    z2(isnan(abs(z2))) = ...
        min2(abs(z2(abs(z2) > 0)));

    for i = 0:(ceil(width / sw) - 1)
        for j = 0:(ceil(height / sh) - 1)

            sx = 1 + i * sw;
            ex = sw + i * sw;
            sy = 1 + j * sh;
            ey = sh + j * sh;
            margesx = overlap;
            margeex = overlap;
            margesy = overlap;
            margeey = overlap;
            if ex > width
                ex = width;
            end
            if ey > height
                ey = height;
            end
            if sx - margesx < 1
                margesx = 0;
            end
            if ex + margeex > width
                margeex = 0;
            end
            if sy - margesy < 1
                margesy = 0;
            end
            if ey + margeey > height
                margeey = 0;
            end

            xrange = (sx - margesx):(ex + margeex);
            yrange = (sy - margesy):(ey + margeey);

            sub_z1 = z1(xrange, yrange);
            sub_z2 = z2(xrange, yrange);
            sub_z_est = z_est(xrange, yrange);
            sub_D_est = D_est(xrange, yrange);

            [sub_z_fil sub_D_fil] = ...
                nlInSAR (sub_z1, sub_z2, ...
                         sub_z_est, sub_D_est, ...
                         hW, hD, ...
                         1, ...
                         h, T, ...
                         Lmin);

            xrange = (1 + margesx):(ex - sx + 1 + margesx);
            yrange = (1 + margesy):(ey - sy + 1 + margesy);

            z_fil(sx:ex, sy:ey) = sub_z_fil(xrange, yrange);
            D_fil(sx:ex, sy:ey) = sub_D_fil(xrange, yrange);
            z_fil(abs(z_fil) <= 0) = ...
                min2(abs(z_fil(abs(z_fil)>0)));
            z_fil(isnan(abs(z_fil))) = ...
                min2(abs(z_fil(abs(z_fil)>0)));
        end
    end
    R_fil = abs(z_fil);
    B_fil = -angle(z_fil);

function res = min2(mat)

    res = min(min(mat));

function r = quantile_insar(ws, alphas)
%QUANTILE_INSAR Centered Quantile of the single-look InSAR
%patch-based log similarity likelihood
%   R = QUANTILE_INSAR(D, W, ALPHA) estimates the centered
%   alpha-quantile of the single-look InSAR patch-based log
%   similarity likelihood.
%
%       ARGUMENT DESCRIPTION:
%               WS       - 1-D array of the widths of the
%                          similarity windows (sizes WxW)
%               ALPHAS   - 1-D array of the alpha values of the
%                          desired alpha-quantiles
%
%       OUTPUT DESCRIPTION:
%               R        - 2-D array of the desired
%                          alpha-quantiles. The First dimension
%                          corresponds to the different alpha
%                          values, the second dimension corresponds
%                          to the different window widths.

for kw = 1:size(ws, 2);
    w = ws(kw);
    [ima_nse ima_nse_p] = insarrnd(ones(w * 256), ...
                                   ones(w * 256), ...
                                   zeros(w * 256));

    k = 1;
    for i = 1:(2*w):(size(ima_nse,1) - 2*w)
        for j = 1:w:(size(ima_nse,2) - w)
            sub_nse_1 = ima_nse(i:(i + w - 1), j:(j + w - 1));
            sub_nse_1p = ima_nse_p(i:(i + w - 1), j:(j + w - 1));
            sub_nse_2 = ima_nse((i + w):(i + 2 * w - 1), j:(j + w - 1));
            sub_nse_2p = ima_nse_p((i + w):(i + 2 * w - 1), j:(j + w - 1));

            a1 = abs(sub_nse_1);
            a1p = abs(sub_nse_1p);
            a2 = abs(sub_nse_2);
            a2p = abs(sub_nse_2p);
            ddp = angle(sub_nse_1 .* conj(sub_nse_1p)) - ...
                  angle(sub_nse_2 .* conj(sub_nse_2p));
            A = (a1.^2 + a1p.^2 + a2.^2 + a2p.^2).^2;
            B = 4 * ((a1 .* a1p).^2 + (a2 .* a2p).^2 ...
                     + 2 * a1 .* a1p .* a2 .* a2p .* cos(ddp));
            C = a1 .* a1p .* a2 .* a2p;
            lsl = -log((C ./ B).^(3/2) .* ...
                       ((A + B) ./ A .* sqrt(B ./ (A - B)) ...
                        - asin(sqrt(B./A))));
            v(k) = mean2(lsl);
            k = k + 1;
        end
    end

    for q = 1:size(alphas, 2)
        r(q, kw) = quantile(v, alphas(q)) - mean(v);
    end

end