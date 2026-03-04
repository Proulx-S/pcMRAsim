function [f,Mzss_v0,fMax,vCrit,v_cm_s] = inflowEnhancementBianciardi(v_cm_s, inflow)
%INFLOWENHANCEMENTBIANCIARDI  Inflow enhancement per Bianciardi et al. (2016) formulation.
%
% Predicts relative inflow signal (slice thickness, TR, flip angle, T1) using the
% Bloch-equation solution for non-stationary spins in three regimes [1]:
%
%   Regime 1 (v = 0):           stationary spins, saturated longitudinal magnetization
%   Regime 2 (0 < v < vCrit):     flow-related enhancement; Mz depends on velocity
%   Regime 3 (v >= vCrit):        complete spin replacement per TR; Mz = M0
%
% With q = exp(-TR/T1)*cos(FA) and vCrit = slice_thickness_cm / TR:
%
%   Mz(v=0)   = M0 * (1 - exp(-TR/T1)) / (1 - q)
%   Mz(0<v<vCrit) = Mz(v=0) + (M0 - Mz(v=0)) * (1 - q^(vCrit/v)) / ((vCrit/v)*(1 - q))
%   Mz(v>=vCrit) = M0
%   
%   !! All for M0=1 !!
%   so Mz           ranges from 0 to 1
%      f          ranges from 1 to fMax
%      Mzss_v0*fMax = 1
%
% At FA = 90°, q = 0 and regime-2 factor reduces to v/vCrit (linear in v).
%
% [1] Bianciardi M et al. The pulsatility volume index: an indicator of cerebrovascular
%     compliance based on fast magnetic resonance imaging of cardiac and respiratory
%     pulsatility. Phil Trans R Soc A 374:20150184, 2016.
%
% Inputs:
%   v_cm_s              velocity [cm/s] (scalar or array)
%   T1_s                blood T1 [s], e.g. 1.2–1.6
%   sliceThickness_mm   slice thickness [mm]
%   TR_s                repetition time [s]
%   FA_deg              flip angle [deg]
%
% Outputs:
%   f            inflow enhancement factor -> Regime 1 (stationary): f(v=0) = 1; Regime 2: f(v=0)<fMax; Regime 3: f(v>=vCrit) = fMax
%   Mzss_v0          longitudinal magnetization for stationary spins (regime 1)
%   fMax           max relative signal (regime 3 vs regime 1)
%   vCrit             critical velocity

% Example:
%   v = linspace(0, 6, 50); vCrit = 3.6;
%   [f,Mzss_v0,fMax,vCrit] = inflowEnhancementBianciardi(v, 2.5, 1.2, 0.033, 90);


T1_s              = inflow.T1;
sliceThickness_mm = inflow.sliceThickness;
TR_s              = inflow.TR;
FA_deg            = inflow.FA;




E1    = exp(-TR_s / T1_s);
alpha = FA_deg * pi/180;
q     = E1 * cos(alpha);   % paper: q = e^(-TR/T1)*cos(FA)

% Longitudinal magnetization for stationary spins (regime 1), Eq. (2.3)
Mzss_v0 = (1 - E1) / (1 - q);   % Mz(v=0)/M0

% Max relative signal (regime 3 vs regime 1): M0 / Mz(v=0)
fMax = 1 / Mzss_v0;            % = (1 - q) / (1 - E1)

% Critical velocity: vCrit = ST/TR [cm/s], with slice thickness in cm
slice_cm = sliceThickness_mm / 10;
vCrit       = slice_cm / TR_s;


if isempty(v_cm_s)
    v_cm_s = linspace(0, vCrit+2*vCrit/32, 34);
end

% v_cm_s = max(0, v_cm_s);

% Regime 2 factor: (1 - q^(vCrit/v)) / ((vCrit/v)*(1 - q));  u = vCrit/v
u = nan(size(v_cm_s));
u(v_cm_s > 0) = vCrit ./ v_cm_s(v_cm_s > 0);

if abs(1 - q) < 1e-12
    % q ≈ 1 (e.g. very low FA): limit → 1
    f = double(v_cm_s >= vCrit);
    f(v_cm_s > 0 & v_cm_s < vCrit) = v_cm_s(v_cm_s > 0 & v_cm_s < vCrit) / vCrit;
else
    % Paper Eq. (2.4): factor (1 - q^u) / (u*(1 - q))
    f = (1 - q.^u) ./ (u * (1 - q));
    f = min(1, max(0, f));
    f(v_cm_s >= vCrit) = 1;
    f(v_cm_s <= 0)  = 0;
end

% Relative signal: f = Mz / Mz(v=0) = 1 + (fMax - 1)*f
f = 1 + f .* (fMax - 1);

end
