> **Frozen record, copied verbatim on 2026-11-03. Do not edit it; correct the maintained docs instead.**
> This is the investigation that found two inflow errata: pcMRAsim `simVesselSpins` / `getMz_n`
> ([`docs/inflowPulseCount.md`](../inflowPulseCount.md)) and util `getMz_ss` (util
> `docs/inflowColumnAverage.md`, branch `dev-erratum-inflowColumn`). It was written on
> 2026-11-02 by a tako session of the vesTrack project (Proulx-S/vesTrack, private), first named
> vesTrack-example-vesSynth-inflow and later vesTrack-vesSynth-model, at
> `vesTrack/scratch/inflowCheck/inflowColumnAverage.md`. Paths below are on the takoyaki server.
> The probe scripts it names are not copied. The committed tests reproduce its numbers: util
> `tests/test_getMz_ss.m` and pcMRAsim `tests/test_inflowCount.m`. Its "Recommendation" section
> addresses the vesTrack port, not this repo.

# Inflow: does pcMRAsim account for the pulse-count spread along the slice column?

2026-11-02, vesTrack-example-vesSynth-inflow, for Seb. Read-and-think task, no repo edits.
Code read: pcMRAsim `origin/dev-qvfMRI` (13c5890), `origin/Miller2026`, `origin/main`; util `main`
(`getMz_ss.m`, `getMxy_ss.m`) and its `Miller2026`/`dev` variants. Probes in this folder:
`columnAverage.py` (+ `.out`), `realCode_qvf.m` (+ `.out`; runs exported copies of the real dev-qvfMRI
functions in `code_qvf/`).

Target: 7 T, 2D spoiled GRE (`fl2d1`), slice 1.2 mm, TR_exc 8.4 ms, FA 40°, TE 3.68 ms, blood
T1 2.58 s → Q1 = E1·cos FA = 0.764, Mz_ss = 0.0137·M0, vCrit = 14.3 cm/s.

## Short answer

- **Bianciardi (`getMz_ss`, Miller2026 / xinAndSava2026 / main): yes, implicitly and correctly.** The
  formula *is* the column average over the slice. Exact for integer vCrit/v; between integers the
  `ceil(u)` in util `getMz_ss` gives up to **−10 %** error in the inflow excess (Mz − Mz_ss). The
  paper's continuous u is closer (≤ +2 %).
- **dev-qvfMRI (`getMz_n` on a 3D grid): yes, explicitly, but with the wrong pulse count.** It
  uses a continuous `n = d/(v·TR)` in `Q1^n`. The true count is an integer, `floor(d/(v·TR))`.
  This gives **~12 % too little inflow excess at all slow velocities**. Blood never reaches M0
  above vCrit: Mz = 0.89 at vCrit, 0.94 at 2·vCrit, instead of 1. With the default 9 planes
  through the slice, discretization adds −11 % to +17 % on top.

So Seb's intuition is right for the Bianciardi formulation. The 3D formulation, the one we planned
to port, models the spread explicitly but gets the per-spin count wrong.

## Algebra

Take a perpendicular vessel, plug flow at velocity v, and a rectangular slice [0, Δz] entered
fresh (M0) at z = 0. RF pulses come every TR. Just before a pulse, a spin at depth z entered
z/v ago, so it has had **m = floor(z/(v·TR)) prior pulses**. That holds at every pulse, because the
pattern is stationary for steady flow. Its Mz just before this pulse (`getMz_n` with n = m) is

  Mz_m = Mz_ss + (M0 − Mz_ss)·Q1^m,  and its signal is Mz_m·sin FA·e^(−TE/T2*).

The voxel sums the column, and the signal is linear in Mz with the same phase for all spins of a
column (same v). So **the signal-weighted average equals the Mz average** times sin FA·E2. With
U = vCrit/v = Δz/(v·TR) and x = z/(v·TR) uniform on [0, U]:

  ⟨Q1^m⟩ = (1/U)∫₀^U Q1^floor(x) dx = [ (1 − Q1^K)/(1 − Q1) + (U − K)·Q1^K ] / U,  K = floor(U).

- **Integer U = u:** ⟨Q1^m⟩ = (1/u)·Σ_{m=0}^{u−1} Q1^m = (1 − Q1^u)/(u·(1 − Q1)). This is exactly
  Bianciardi's Eq. 2.4. The paper says so in words: the last-entered velocity segment "experienced
  1 RF pulse", the first-entered one vC/v pulses, assuming integer vC/v. Counting the readout pulse,
  1..u pulses means 0..u−1 prior pulses. **No off-by-one.**
- **U < 1 (v > vCrit):** every spin in the slice entered less than one TR ago, so Mz = M0.
- **Continuous n (dev-qvfMRI as coded):** (1/U)∫₀^U Q1^x dx = (1 − Q1^U)/(U·ln(1/Q1)). That
  differs from the discrete sum by (1 − Q1)/ln(1/Q1) = **0.876** here. Q1^x lies below the true
  step function Q1^floor(x), because a spin is credited fractional pulses it has not received.
- **`ceil(u)` (util `getMz_ss` since e566edf):** this applies the integer formula with u rounded up.
  It is exact at integer U and too saturated between integers. The error is worst just below
  vCrit/k, where u jumps by 1.

A brute-force pulse-by-pulse recursion (`truth` in `columnAverage.py`, 20 000 stratified spins)
reproduces the closed form to 4 decimals. So the closed form is the reference below.

## Numbers (perpendicular, central column, plug v)

Mz/M0 before the readout pulse:

| v/vCrit | v cm/s | exact | 3D as coded, fine grid | 3D as coded, 9 planes | 9 planes, floor n | Bianciardi `ceil` (util) | Bianciardi continuous u |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.71 | 0.221 | 0.196 | 0.256 | 0.208 | 0.221 | 0.221 |
| 0.10 | 1.43 | 0.403 | 0.355 | 0.408 | 0.409 | 0.403 | 0.403 |
| 0.20 | 2.86 | 0.632 | 0.555 | 0.597 | 0.636 | 0.632 | 0.632 |
| 0.30 | 4.29 | 0.752 | 0.664 | 0.697 | 0.750 | 0.702 | 0.756 |
| 0.45 | 6.43 | 0.854 | 0.756 | 0.781 | 0.851 | 0.785 | 0.860 |
| 0.70 | 10.0 | 0.930 | 0.832 | 0.850 | 0.922 | 0.883 | 0.948 |
| 0.90 | 12.9 | 0.977 | 0.866 | 0.880 | 0.974 | 0.883 | 0.986 |
| 1.00 | 14.3 | 1.000 | 0.878 | 0.891 | 1.000 | 1.000 | 1.000 |
| 2.00 | 28.6 | 1.000 | 0.936 | 0.943 | 1.000 | 1.000 | 1.000 |

The relative error of the inflow excess (Mz − Mz_ss) against exact, over v/vCrit 0.05–0.9:

- 3D as coded, fine grid: **−11 to −12 %**.
- 3D as coded, 9 planes: **+17 % to −10 %**. At slow v the first plane (d = 0) counts as fully fresh.
- `ceil`: **0 to −10 %**.
- Continuous u: **0 to +2 %**.

The real dev-qvfMRI `simVesselSpins` reproduces the "9 planes as coded" column to 4 decimals
(`realCode_qvf.out`).

Parabolic flow, lumen-mean Mz (each radius with its own column average):

| vMean/vCrit | exact | 3D as coded, fine grid | 3D as coded, 9 planes | Bianciardi `ceil` |
|---|---|---|---|---|
| 0.10 | 0.373 | 0.329 | 0.384 | 0.358 |
| 0.25 | 0.618 | 0.546 | 0.586 | 0.585 |
| 0.50 | 0.780 | 0.691 | 0.720 | 0.734 |
| 1.00 | 0.890 | 0.802 | 0.822 | 0.867 |

How the error of floor-n at plane centres shrinks with the plane count (max over v/vCrit
0.05–1): 9 planes 12 %, 21 planes 4.5 %, 41 planes 1.7 %, 81 planes 0.8 %. Point sampling converges
slowly because it samples a step function. The **cell-averaged** form (below) is exact, to
1e−13 %, with any number of planes, including 1.

## The other details checked

- **Pulse counts in the 3D grid (default `isochromatPerVoxPerDirMin` = 2).** `setGrid` makes the
  count odd (3) on the smallest voxel side: 0.4 mm / 3 = 0.133 mm. The slice then gets
  round(1.2/0.133) = **9 planes**. This was confirmed by running `setGrid`: 3 × 3 × 9, first plane at
  −0.533 mm.
- **Slab entry.** `gridD = gridSLC − gridSLC(:,:,1)` measures depth from the first isochromat plane,
  not from the slice face. That is a left-edge Riemann sum, offset by half a plane (0.067 mm). The
  bottom plane counts as never pulsed, which explains the +17 % at slow flow. The entry is also
  `fovSLC`, not `pMri.sliceThickness`. The runSim defaults disagree: `fovSLC` = 2 mm, slice 1 mm.
  Only `getMz_ss` uses `sliceThickness`.
- **`n_hat(3)` (tilt).** `gridD/n_hat(3)` gives the path length along the axis. `d/(v·TR)` is then
  `z/(v·n3·TR)`, i.e. the through-plane velocity. That is right for an infinite slab with entry
  through the bottom face.
  - **Bug: n3 < 0 with v > 0 gives negative n.** The real code returns Mz up to 2.2·M0 down the
    column. costFun avoids this only because it builds nz ≥ 0.
  - n3 = 0 gives 0/0.
- **NaN edge case.** A lumen isochromat exactly on r = ID/2 has v = 0. If it sits in the first plane
  (d = 0), n = 0/0 = NaN, and the voxel sum becomes NaN. Not hit in my runs.
- **Side bug: the surround uses the blood T1.** When `S.surround` is empty, the surround is
  computed with `getMz_n(pMri, pMri.relax.blood, inf)` and then GM T2*. The real code gives Mz =
  0.01375, which is the blood Mz_ss; with GM T1 it should be 0.01822.
- **Plug vs parabolic.** Both formulations use each spin's own v. Miller2026 does it per in-plane
  spin, dev-qvfMRI per isochromat. So each radius gets its own column average, as it should. Under
  `plug`, Miller2026 uses vMean only.
- **Bianciardi variants.** The pcMRAsim `main` `inflowEnhancementBianciardi` and the paper use
  continuous u. util `main`, `Miller2026` and `dev` use `ceil(u)`. util `dev` also moves where vCrit
  is set.

## Literature (full-text pass, 2026-11-02, in Seb's Chrome)

The pulse-count convention matters for the off-by-one question. In the literature, n counts the
pulses received **including the readout pulse**, and the signal comes from Mz just before the n-th
pulse. That equals `getMz_n(n−1)`, i.e. m = n − 1 prior pulses, which is the floor convention used
above.

- **Bianciardi et al. 2016** [1]. Read: the whole PMC article text, plus the image of Eq. 2.4.
  - Eq. 2.4 is exactly (1 − q^(vC/v))/((vC/v)(1 − q)) inside Mz(v=0) + (M0 − Mz(v=0))·…, times
    G·sin FA·e^(−TE/T2*).
  - The text says the slice holds velocity segments. The last-entered segment "experienced 1 RF
    pulse"; the first-entered one experienced vC/v pulses. It also says it "assumed integer values
    of vC/v". So it is the column average, with the readout counted in n.
  - Regime 3 (v ≥ vC) is full replacement each TR: Mxy = M0·sin FA·E2.
  - Figure 1 uses only integer vC/v (2–100).
  - It cites the Haacke textbook [13] for the Bloch steady-state solutions; see Brown et al. 2014
    [8] below, where Eq. 24.9 is the same formula. It states no flow-profile assumption and no
    slice-profile assumption.
- **Kim & Parker 2012** [6], book chapter "Time-of-flight angiography", read through the Stanford
  library. Read: pp. 39–42 (TOF signal theory); not the later 3D-TOF / MOTSA / artifact sections.
  - Fig. 2.1b divides the slice into dz = v·TR segments that see 1, 2, 3 … pulses.
  - Eq. 2.4: Mz(n⁻) = Mz_ss + (e^(−TR/T1) cos θ)^(n−1)·(M0 − Mz_ss), n ≥ 1.
  - For v < Vc, the inflow enhancement "can be modeled by subdivision of imaging slices into
    multiple compartments and summation of the geometric series".
  - The chapter also notes that the enhancement drops for oblique flow: only the perpendicular
    velocity component counts. That supports using v_z = v·n3.
- **Ashenagar et al. 2025** [2]. Read: the whole PMC text and Appendix A, plus the image of Eq. 1.
  - Spin ensembles are point particles whose position is the integral of the velocity.
  - n is the integer count of RF events an ensemble's trajectory crosses; Δt are the actual
    inter-pulse intervals.
  - The signal is "averaged together within each slice" after each repetition. So this is the
    history-aware column average.
  - The slice profile is a step function, named as a limitation.
  - I checked numerically that their Eq. 1 with all Δt = TR equals sin θ·`getMz_n(n−1)` to 6
    decimals (`ashenagarEq1.py`). Appendix A says the same expression reduces to Kim & Parker's
    equation at steady state.
  - Relevant to vesTrack time series: for pulsatile flow, their simulations and data show a
    **frequency-dependent attenuation**. Faster velocity oscillations give relatively smaller
    inflow signal, i.e. the signal is "a weakly low-pass filtered version" of the flow.
    Instantaneous-v models cannot produce this.
  - Their Fig. 7d shows signals saturating at integer multiples of vCrit across slices.
- **Gao & Liu 2012** [5]. Read: the whole review through ScienceDirect.
  - Its Eq. 1 is a cruder schematic: S ∝ M_sat·L + (M_fresh − M_sat)·V·TR for V·TR < L. Only the
    newest segment is fresh and all older ones sit at M_sat, which makes S linear in V.
  - With M_sat = Mz_ss that gives an inflow excess of 1/U. That is the q = 0 (90°) limit of
    Bianciardi. At 40° here it underestimates the excess by up to 1/(1 − Q1) ≈ 4×.
  - The review points to Gao 1988 for the full spoiled-FLASH treatment ("RF pulse history, slice
    thickness, flow velocity…"), validated with laminar phantom flow.
- **Gullberg et al. 1988** [3] (spin echo). Read: the scanned PDF, pp. 437–439 (introduction,
  fractional-volume method); not the derivations for the 10 pulse sequences.
  - "The total signal is the sum of the echo amplitudes from each spin population, each weighted
    by the fractional volume occupied by that population."
  - Laminar flow uses paraboloid volume segments (Eqs. 3–6). That equals giving each radius its own
    column average and integrating over the cross-section, as pcMRAsim does per spin.
- **Brown, Cheng, Haacke, Thompson & Venkatesan 2014** [8], ch. 24 "MR angiography and flow
  quantification". This is the 2nd edition of the Haacke textbook that Bianciardi cites as [13].
  Seb's PDF of the chapter; read §24.1 (pp. 702–711) in full, not §24.2–24.4.
  - It is the source of Bianciardi's formula. Eq. 24.5: Mz(m⁻) = Mze + q^(m−1)(M0 − Mze), m ≥ 1,
    "as long as the spin is affected by exactly (m − 1) pulses".
  - §24.1.3 (2D, flow perpendicular to the slice) cuts the slice into n = TH/(v·TR) "velocity
    segments" that have seen m = 1 … n pulses. Eq. 24.9: S(p ≥ n) = Λ (sin θ/n) Σ_{m=1..n}
    [Mze + (M0 − Mze) q^(m−1)] = Λ sin θ [Mze + (M0 − Mze)(1 − qⁿ)/(n(1 − q))]. That is the column
    average, written out with the summation.
  - Footnote 4: the segment's Mz just before the pulse "depend[s] on m − 1 pulses". No off-by-one.
  - On integers, the book "assume[s] that an exact integer n of velocity segments fit into the
    slice". For "arbitrary, but large values of TH/vTR, a subdivision according to the nearest
    integer will lead to a reasonable approximation". So **rounding**, not util's `ceil`. And only
    for large n; at small n (v near vCrit) neither rounding is good, which matches the −10 % I found
    for `ceil` just below vCrit. The exact non-integer form above removes the issue.
  - §24.1.5 gives an even cruder continuous approximation: a "pseudo-T1", 1/T1,eff = 1/T1 + v/TH.
  - Laminar flow is only mentioned qualitatively: the centre sees the fewest pulses and is
    brightest.
- **Gao et al. 1996** [7], MRM 36:314. Seb's PDF, in place of Gao 1988; read in full (6 pp.).
  - It applies the Gao 1988 / Gao & Gore 1994 model to single-slice fMRI with through-plane flow
    and a laminar profile.
  - Eq. 2 (spoiled GRE, laminar): I_b(V) = Σ_{n=1..m} M_y^(n)(TE)·V·TR·πa² + Σ_{i=1..N−m}
    M_y^(m+i)(TE)·Δ_(m,i) for 0 < V ≤ L/(2TR). Here M_y^(n) is the transverse magnetization "after
    experiencing the nth RF pulse", and Δ_(m,i) is the volume of each partially filled population.
    For V ≥ L/(2TR) it is Σ_i M_y^(i)(TE)·Δ_(0,i).
  - So it is a volume-weighted sum over spin populations with integer pulse counts, the readout
    included in n. The paraboloid segments give per-radius column averages.
  - The threshold is L/(2TR) in mean velocity, i.e. peak velocity = L/TR.
  - For laminar flow even V above that threshold does not fully refresh the slice: "there are always
    residual portions of each flow population … in the image slice", from the slow flow near the
    wall. pcMRAsim reproduces this, because each radius gets its own v.
- **Gao et al. 1988** [4]. Not read: the Wiley/AAPM links bounce to the abstract. Gao 1996 cites it
  for the details of Eq. 2.

Bottom line: every source models slow through-plane inflow as a sum over slice segments, each with
an integer pulse count, weighted by volume. That is the column average. Closed forms (Bianciardi,
Kim & Parker) assume steady flow, integer vC/v or summation over whole segments, a perpendicular
vessel and a rectangular slice profile. Ashenagar generalizes this to time-varying flow by tracking
each spin's history. None uses a fractional pulse count, so dev-qvfMRI's continuous n has no
support in the literature.

## Recommendation for the vesTrack port

1. **Port the dev-qvfMRI 3D structure, but fix the per-isochromat pulse count.**
   - Replace `Q1^n` (continuous n) by the exact **cell average of Q1^floor over each isochromat's
     slab of thickness**: ⟨Q1^m⟩_cell = [F(x_b) − F(x_a)]/(x_b − x_a), with
     F(x) = (1 − Q1^K)/(1 − Q1) + (x − K)·Q1^K, K = floor(x), and x = d/(v_z·TR) at the
     isochromat's lower and upper faces.
   - Measure d from the slice face (`−sliceThickness/2`), not from the first plane.
   - The result is exact for a perpendicular vessel with any plane count. It is ~10 lines and costs
     nothing.
   - The plane count is then set by **geometry only**: how r varies along z for a tilted vessel.
     The default 9 planes is fine for perpendicular vessels; use ≥ 9 for tilt.
2. Use v_z = v·|n3| and pick the entry face from sign(v·n3). Assert n3 ≠ 0. Treat v = 0 as Mz_ss
   before dividing.
3. Use one slab thickness: tie `fovSLC` to `pMri.sliceThickness` (2D single slice: `matSLC` = 1).
4. Compute the surround with the GM T1.
5. If a 2D shortcut is wanted (perpendicular vessel only), use the **exact** formula, not `ceil(u)`.
   Seb's e566edf change should be revisited.
6. **For time series, the bigger approximation is quasi-static flow.** Ashenagar 2025 [2] shows
   its signature: the inflow signal low-pass filters the velocity. Every formulation here uses
   the instantaneous v. The pulse count really depends on the velocity history over the transit
   time Δz/v_z: 84 ms at 1.4 cm/s, 170 ms at 0.7 cm/s, a sizeable fraction of a cardiac cycle.
   - If vesTrack synthesizes pulsatile slow flow, consider the Ashenagar-style history model [2]:
     entry time from ∫v dt, then integer counts.
   - Check this against the quasi-static model before deciding.
   - run-1.json has frame TR 0.840 s = 100 × 8.4 ms, so no RF gap between frames; the
     steady-pulsing assumption holds.

## Not verified

- Slice-profile effects. The real excitation is not rectangular at 40°, and none of the models
  include it.
- Spoiling imperfections and T2* differences within the column.
- Gao et al. 1988 itself (Gao 1996 read in its place). Brown et al. 2014 beyond §24.1.
- Gullberg 1988 beyond pp. 437–439.
- Kim & Parker 2012 beyond pp. 39–42.
- Tilted-vessel numbers: I checked the algebra, not a simulation.
- The NaN boundary case: from reading only, not triggered.
- vesTrack's own synthesize code: none exists yet; I grepped for these functions and found none
  outside scratch.

## References

1. Bianciardi M et al. (2016) The pulsatility volume index… *Phil Trans R Soc A* 374:20150184.
   https://doi.org/10.1098/rsta.2015.0184 ·
   [Consensus](https://consensus.app/papers/details/96e214338860540c8786fe95c054b86d/?utm_source=claude_desktop)
2. Ashenagar B et al. (2025) Modeling dynamic inflow effects in fMRI to quantify cerebrospinal fluid
   flow. *Imaging Neuroscience*. https://doi.org/10.1162/imag.a.9 ·
   [Consensus](https://consensus.app/papers/details/111fd3d1d77454099995b971edf34dc5/?utm_source=claude_desktop)
3. Gullberg GT et al. (1988) A mathematical model for signal from spins flowing during the
   application of spin echo pulse sequences. *Magn Reson Imaging*.
   https://doi.org/10.1016/0730-725x(88)90481-x ·
   [Consensus](https://consensus.app/papers/details/a46d2fc8624157b0a1e167c642e11adb/?utm_source=claude_desktop)
4. Gao JH et al. (1988) NMR signal from flowing nuclei in rapid imaging using gradient echoes.
   *Med Phys*. https://doi.org/10.1118/1.596197 ·
   [Consensus](https://consensus.app/papers/details/1d7983289f2b5ff8a0cffae462e138ca/?utm_source=claude_desktop)
5. Gao JH, Liu HL (2012) Inflow effects on functional MRI. *NeuroImage*.
   https://doi.org/10.1016/j.neuroimage.2011.09.088 ·
   [Consensus](https://consensus.app/papers/details/907fc97d823f54b59ba9cbae0bec6188/?utm_source=claude_desktop)
6. Kim SE, Parker DL (2012) Time-of-flight angiography. In: Carr JC, Carroll TJ (eds), *Magnetic
   Resonance Angiography: Principles and Applications*, pp. 39–50. Springer.
   https://doi.org/10.1007/978-1-4419-1686-0_2
7. Gao JH, Miller I, Lai S, Xiong J, Fox PT (1996) Quantitative assessment of blood inflow effects
   in functional MRI signals. *Magn Reson Med* 36:314–319. https://doi.org/10.1002/mrm.1910360219
8. Brown RW, Cheng YCN, Haacke EM, Thompson MR, Venkatesan R (2014) MR angiography and flow
   quantification. Ch. 24 in *Magnetic Resonance Imaging: Physical Principles and Sequence Design*,
   2nd ed. Wiley. https://doi.org/10.1002/9781118633953.ch24
