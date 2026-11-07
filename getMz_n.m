function [Mz_n, pMri, Mo] = getMz_n(pMri, pRelax, n, Mo, nWidth)
% getMz_n  Longitudinal magnetization just before an RF pulse, for spins that entered fully relaxed (M0)
% and have since received n earlier RF pulses.
%
%   Mz_n = getMz_n(pMri, pRelax, n)              n integer (or Inf): exactly n earlier pulses
%   Mz_n = getMz_n(pMri, pRelax, x, Mo, xWidth)  cell average over the pulse-unit interval
%                                                [x - xWidth/2, x + xWidth/2] (see below)
%
% Unlike getMz_ss (util), which averages over a whole slice column of plug flow, getMz_n is for
% isochromats at known positions in the slab (simVesselSpins).
%
% INPUTS
%   pMri   - MRI acquisition struct. Required fields: TR, FA.
%   pRelax - Relaxation struct with field T1 [s]. Used to set (override) pMri.E1 and pMri.Q1.
%   n      - Number of earlier RF pulses (integer, scalar or array; Inf -> Mz_ss).
%            With nWidth: the centre x of an interval in "pulse units" x = depth/(v*TR), where a spin
%            at x has received floor(x) earlier pulses.
%   Mo     - (optional) Equilibrium magnetization. Default 1.
%   nWidth - (optional) Width of the interval in pulse units (scalar or same size as n), e.g. the
%            isochromat's thickness along the slice divided by (v_through-plane * TR).
%
% OUTPUTS
%   Mz_n   - Longitudinal magnetization just before the next pulse. Same size as n.
%   pMri   - Same struct with E1 and Q1 updated.
%   Mo     - Equilibrium magnetization used.
%
% Formula (exactly n earlier pulses):  Mz_n = Mz_ss + (Mo - Mz_ss) * Q1^n
%   with Mz_ss = Mo*(1-E1)/(1-Q1), Q1 = E1*cosd(FA), E1 = exp(-TR/T1).
% Cell average (nWidth given): the mean of Mz_floor(x) over x in [xa, xb],
%   Mz = Mz_ss + (Mo - Mz_ss) * (F(xb) - F(xa)) / (xb - xa),
%   F(x) = integral_0^x Q1^floor(t) dt = (1 - Q1^K)/(1 - Q1) + (x - K)*Q1^K,  K = floor(x).
% The pulse count is an integer: passing a non-integer n without nWidth is an error, because
% Q1^n with fractional n credits a spin with pulses it has not received (Q1^x < Q1^floor(x)).
% See docs/inflowPulseCount.md.
%
% Special cases:
%   n = 0    -> Mz_n = Mo     (not yet pulsed, fully relaxed)
%   n = Inf  -> Mz_n = Mz_ss  (stationary spins: steady state)
%
% INFLOW MODEL VERSION: exact (cell average of the integer pulse count, piecewise linear).
%   A spin at x = depth/(v*TR) has received floor(x) earlier pulses, never a fraction of one; the cell
%   average weights each integer-count population by its fractional width in the cell, so the result
%   is piecewise linear in x. Over a whole slice column of plug flow it equals util getMz_ss (exact),
%   which in turn equals Bianciardi et al. 2016 eq. 2.4 at integer vCrit/v. Other versions in older
%   code (not used here): "continuous count" (Q1^x with fractional x; too low, the bug fixed on this
%   branch), "staircase" (ceil(vCrit/v) in eq. 2.4) and "continuous-u" (eq. 2.4 at non-integer
%   vCrit/v; smooth, up to a few % high between integers; inflowEnhancementBianciardi).
%
% REFERENCES
%   Bianciardi M, Toschi N, Polimeni JR, et al. (2016) The pulsatility volume index: an indicator of
%     cerebrovascular compliance based on fast magnetic resonance imaging of cardiac and respiratory
%     pulsatility. Phil Trans R Soc A 374:20150184. doi:10.1098/rsta.2015.0184
%     (three-regime formulation; eq. 2.4, derived for integer vCrit/v)
%   Gullberg GT, Simons MA, Wehrli FW (1988) A mathematical model for signal from spins flowing during
%     the application of spin echo pulse sequences. Magn Reson Imaging 6:437-461.
%     doi:10.1016/0730-725X(88)90481-X (fractional-volume method: each spin population weighted by
%     its fractional volume)
%   Gao J-H, Holland SK, Gore JC (1988) Nuclear magnetic resonance signal from flowing nuclei in rapid
%     imaging using gradient echoes. Med Phys 15:809-814. doi:10.1118/1.596197
%     (that method for spoiled FLASH: integer-count populations plus the partly filled residual one,
%     sec. III.A, eq. 26)
%   Comparison of the four versions and the choice of the exact one (2026-11-06; private repo, the
%   summary above stands on its own):
%   https://github.com/Proulx-S/microscopyValidation/blob/9057997d65a96cad63419fed24daa1ba885ed6c6/docs/inflowBetweenIntegerPulseCounts.md

if exist('pRelax','var') && ~isempty(pRelax) && isfield(pRelax,'T1') && ~isempty(pRelax.T1)
    pMri.E1 = exp(-pMri.TR / pRelax.T1);
    pMri.Q1 = pMri.E1 * cosd(pMri.FA);
elseif ~isfield(pMri,'E1') || isempty(pMri.E1)
    error(['E1 not found in pMri' newline 'provide pRelax input']);
end
if ~exist('Mo','var') || isempty(Mo)
    Mo = 1;
end
Q1    = pMri.Q1;
Mz_ss = Mo * (1 - pMri.E1) / (1 - Q1);

if ~exist('nWidth','var') || isempty(nWidth)
    % exactly n earlier pulses
    fin = isfinite(n);
    if any(n(fin) ~= round(n(fin)))
        error('getMz_n:nonIntegerN', ['n must be an integer number of RF pulses (or Inf).' newline ...
            'For a position in pulse units x = depth/(v*TR), use floor(x), or pass nWidth for a cell average.']);
    end
    Mz_n = Mz_ss + (Mo - Mz_ss) .* Q1.^n;
else
    % cell average of Q1^floor(x) over [x - w/2, x + w/2]
    x  = n;
    w  = nWidth .* ones(size(x));
    xa = max(0, x - w/2);
    xb = x + w/2;
    f  = nan(size(x));
    st = ~isfinite(x);                 % stationary (v = 0): steady state
    f(st) = 0;
    pt = isfinite(x) & (xb - xa) <= 0; % zero width: a point
    f(pt) = Q1.^floor(x(pt));
    cl = isfinite(x) & (xb - xa) > 0;  % finite cell
    f(cl) = (Fint(xb(cl), Q1) - Fint(xa(cl), Q1)) ./ (xb(cl) - xa(cl));
    Mz_n = Mz_ss + (Mo - Mz_ss) .* f;
end
end


function F = Fint(x, Q1)
% integral_0^x Q1^floor(t) dt
K = floor(x);
F = (1 - Q1.^K)./(1 - Q1) + (x - K).*Q1.^K;
end
