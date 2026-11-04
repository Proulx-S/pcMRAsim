# Inflow in the 3D isochromat model: counting RF pulses per isochromat

> **Erratum (2026-11-03), on the separate branch `errata` only.** `dev-qvfMRI` and the other
> existing branches are deliberately left unchanged, so they remain the reproducible default; check
> out `errata` to get the fix. On `dev-qvfMRI` (13c5890),
> `simVesselSpins` gave each lumen isochromat $Q_1^{\,n}$ with a continuous pulse count
> $n = d/(v\,TR)$, with $d$ measured from the first isochromat plane. A spin receives a whole number
> of pulses, so this is wrong in three ways:
>
> - the bulk is too saturated: the excess is scaled by 0.876 at 7 T, 40°;
> - the slab never refreshes above $v_c$: $0.89\,M_0$ at $v_c$;
> - the result depends on the plane count: −11 % to +39 % of the inflow excess at the default
>   9 planes.
>
> The fix gives each isochromat the exact average of the integer count over its cell (§3). It also
> fixes five related defects (§4). Affected branches and callers are in §7, and how this was found
> in §8.

`simVesselSpins` samples the vessel on a 3D isochromat grid that spans the excited slab. Each
lumen isochromat gets its own longitudinal magnetization from the number of RF pulses it has
received since it entered the slab (`getMz_n`). The voxel signal is the sum over its
isochromats. This note gives the model, the corrections made to it, and the validation in
[`tests/test_inflowCount.m`](../tests/test_inflowCount.m).

## 1. Per-spin magnetization

The model is a spoiled gradient-echo sequence with flip angle $\alpha$ and repetition time $TR$,
applied to a rectangular slab of thickness $TH$, with blood $T_1$. Write

$$
E_1 = e^{-TR/T_1}, \qquad Q_1 = E_1\cos\alpha, \qquad M_{ss} = M_0\,\frac{1-E_1}{1-Q_1}.
$$

A spin that entered fully relaxed and has received exactly $m$ earlier pulses has, just before
the next pulse,

$$
M_z(m) = M_{ss} + (M_0 - M_{ss})\,Q_1^{\,m},
$$

and contributes $M_z(m)\sin\alpha\,e^{-TE/T_2^*}$ to the signal of that pulse. Kim & Parker (2012,
eq. 2.4) and Brown et al. (2014, eq. 24.5) write the same expression as
$M_z(n^-) = M_{z,ss} + (e^{-TR/T_1}\cos\theta)^{n-1}(M_0 - M_{z,ss})$ for the $n$-th pulse,
i.e. with $n-1$ earlier pulses. So does Ashenagar et al. (2025, eq. 1) with all inter-pulse
intervals equal to $TR$. $m$ is a whole number of pulses.

## 2. The pulse count of an isochromat

Let blood move along the vessel axis $\hat n$ at speed $v$. Its through-plane speed is
$v_z = v\,|\hat n_z|$. Blood enters through the slab face that is upstream for the sign of
$v\,\hat n_z$. Pulses fire at $t_k = k\,TR$. Just before a pulse, a spin at through-plane depth
$d_z$ below the entry face entered $d_z/v_z$ ago, so it has received

$$
m = \left\lfloor x \right\rfloor, \qquad x = \frac{d_z}{v_z\,TR}
$$

earlier pulses. Under steady flow this pattern is the same at every pulse.

Only the through-plane component of the velocity matters. For a tilted vessel, the path length
along the axis, $d_z/|\hat n_z|$, divided by $v$, gives the same transit time $d_z/v_z$. This is
geometry, and it is how the literature treats it:

- Ashenagar et al. (2025, Discussion) note that the inflow method "measures flow orthogonal to the
  imaging slice, and can only quantify flow in that direction".
- Kim & Parker (2012, book chapter) note that the enhancement drops for oblique flow.

No source read states the transit time $TH/(v\cos\theta)$ in these words. The test checks it
directly: a 30° tilt matches a perpendicular vessel at equal $v\,\hat n_z$ (§5).

## 3. One isochromat stands for a cell: exact cell average

An isochromat at depth $d_z$ represents the cell $[d_z - \Delta/2,\ d_z + \Delta/2]$, where
$\Delta$ is the plane spacing. Within the cell, $x$ runs over $[x_a, x_b]$ with
$x_{a,b} = (d_z \mp \Delta/2)/(v_z TR)$, and the pulse count steps from $\lfloor x_a\rfloor$ to
$\lfloor x_b\rfloor$. The isochromat is given the exact mean over its cell:

$$
\left\langle Q_1^{\lfloor x\rfloor}\right\rangle_{[x_a,x_b]} = \frac{F(x_b) - F(x_a)}{x_b - x_a},
\qquad
F(x) = \int_0^x Q_1^{\lfloor t\rfloor}\,dt = \frac{1-Q_1^{K}}{1-Q_1} + (x-K)\,Q_1^{K},\quad K=\lfloor x\rfloor,
$$

$$
M_z^{\text{cell}} = M_{ss} + (M_0-M_{ss})\,\frac{F(x_b)-F(x_a)}{x_b-x_a}.
$$

This is `getMz_n(pMri, pRelax, x, [], xWidth)`. The cells tile the slab, so the column sum is the
exact integral over the slab for any number of planes. For a perpendicular plug-flow column it
equals the closed-form slice average. That average is Bianciardi et al. (2016, eq. 2.4) and Brown
et al. (2014, eq. 24.9) at integer $v_c/v$, and util `getMz_ss` at any $v$. The plane count now
only matters for the geometry, i.e. how the radius, and hence $v$, varies through the slab for a
tilted vessel.

### Provenance: what comes from the literature and what is derived here

| Ingredient | Peer-reviewed journal articles | Textbooks |
|---|---|---|
| $M_z$ after $m$ pulses, $M_{ss} + (M_0-M_{ss})Q_1^m$ | Ashenagar et al. (2025, eq. 1, with all inter-pulse intervals equal to $TR$); Bianciardi et al. (2016, eqs. 2.3–2.4) | Brown et al. (2014, eq. 24.5: valid "as long as the spin is affected by exactly $(m-1)$ pulses"); Kim & Parker (2012, eq. 2.4) |
| The pulse count is a **whole number** that grows with depth into the slab | Bianciardi et al. (2016, after eq. 2.4: the last-entered velocity segment "experienced 1 RF pulse, whereas the first entered segment experienced $v_C/v$ RF pulses"); Ashenagar et al. (2025: "fluid within the first slice can be divided into partitions each having different amounts of excitation", and each ensemble's excitations are tracked to define $n$); Söderström et al. (2026, §2.5: "the number of RF pulses applied to a specific element") | Brown et al. (2014, §24.1.3); Kim & Parker (2012, Fig. 2.1) |
| Signal = volume-weighted average over populations, partly filled ones included | Gullberg et al. (1988, abstract: "contributions from a finite number of spin populations, with each population signal weighted by the fractional volume of that spin population"); Gao et al. (1996, eq. 2: $\Delta_{(m,i)}$, "the volume occupied by the $(m+i)$-labeled spin population in the selected slice"); Ashenagar et al. (2025: signals "averaged together within each slice"); Söderström et al. (2026, §2.5: "the mean transverse magnetization of all elements localized within the slice") | Brown et al. (2014, eq. 24.9) |
| Point-sampled sub-slices limit how finely the inflow is resolved | Markl et al. (2003, Methods): with the slice cut into $N_S$ sub-slices, "the number of subslices determined the minimal spin replacement per TR"; 20 sub-slices gave 5 % per TR | |
| Cell average $F(x_b)-F(x_a)$ over an isochromat's thickness | **Derived here** from the rows above. It is the exact integral of the integer count over the cell, so it removes the sub-slice resolution limit that the previous point sampling had. No source states it. It is validated against a pulse-by-pulse simulation that uses no closed form (§5) | |

**No source read uses a fractional pulse count** $Q_1^{x}$ with continuous $x$, which is what the
previous code did. The only continuous closed forms in the literature are Bianciardi et al.'s
eq. 2.4 evaluated at non-integer $v_C/v$ (which they exclude: "we assumed integer values of
$v_C/v$") and Brown et al.'s "pseudo-$T_1$" approximation (§24.1.5). Both are whole-column
formulas, not per-spin counts.

## 4. What was wrong, and what changed

| | Before | Now |
|---|---|---|
| Pulse count | continuous $n = d/(v\,TR)$ in $Q_1^{\,n}$ | integer $\lfloor x\rfloor$, averaged exactly over each cell |
| Depth origin | the first isochromat plane | the slab face |
| Plane-count dependence | yes (9 planes by default: −11 % to +39 % of the excess) | none (spread $< 10^{-15}$ from 3 to 21 planes) |
| Axis toward $-z$ with $v>0$ | negative $n$: $M_z$ up to $2.2\,M_0$ | entry face chosen from the sign of $v\,\hat n_z$ |
| In-plane vessel ($\hat n_z = 0$) | $0/0$ | explicit error |
| Isochromat with $v = 0$ at $d = 0$ (lumen edge, first plane) | $0/0$ → NaN; it reaches the voxel sum when `WT > 0` (with `WT = 0` the tissue value overwrote it, see [#3](https://github.com/Proulx-S/pcMRAsim/issues/3)) | $M_{ss}$ |
| Surround (`S.surround` empty) | $M_{ss}$ with the **blood** $T_1$ | $M_{ss}$ with the GM $T_1$ |
| Slab thickness | `fovSLC`, silently independent of `pMri.sliceThickness` (runSim defaults: 2 mm vs 1 mm) | error if they differ when `matSLC = 1`; runSim default `fovSLC = 1` |
| `getMz_n(…, n)` with non-integer $n$ | accepted | error; pass `nWidth` for a cell average |

Why the continuous count is wrong:

- **Too saturated in the bulk.** $Q_1^{\,x} \le Q_1^{\lfloor x\rfloor}$, so a spin is credited
  with a fraction of a pulse it has not received. Integrated over a column this scales the inflow
  excess by $(1-Q_1)/\ln(1/Q_1)$, which is 0.876 at $Q_1 = 0.764$.
- **No full refresh above $v_c$.** For $v \ge v_c$, every spin in the slab entered less than one
  $TR$ ago, so $M_z = M_0$. The continuous count gives $0.89\,M_0$ at $v_c$ and $0.94\,M_0$ at
  $2v_c$.

## 5. Validation

[`tests/test_inflowCount.m`](../tests/test_inflowCount.m) uses the target protocol: 7 T 2D spoiled
GRE, TR 8.4 ms, 40°, TE 3.68 ms, slab 1.2 mm, voxel 0.4 × 0.4 × 1.2 mm, blood $T_1$ 2.58 s,
$v_c = 14.3$ cm/s, default grid with 9 planes. The reference is a pulse-by-pulse simulation:
200 003 spins with stratified entry times cross the slab, each spin's received pulses are counted
from its entry time, and its $M_z$ is evolved explicitly. The test compares the reference with the
central isochromat column of a perpendicular vessel ($r = 0$, so $v = 2\,v_{mean}$).

| $v/v_c$ | reference $\overline{M_z}/M_0$ | new | previous | previous: error of the excess |
|---|---|---|---|---|
| 0.03 | 0.1389 | 0.1389 | 0.1872 | +38.6 % |
| 0.05 | 0.2214 | 0.2214 | 0.2557 | +16.5 % |
| 0.20 | 0.6315 | 0.6315 | 0.5967 | −5.6 % |
| 0.45 | 0.8539 | 0.8539 | 0.7805 | −8.7 % |
| 0.90 | 0.9767 | 0.9767 | 0.8802 | −10.0 % |
| 1.00 | 1.0000 | 1.0000 | 0.8912 | −11.0 % |
| 2.00 | 1.0000 | 1.0000 | 0.9433 | −5.7 % |

New vs reference: $\le 6\times10^{-7}\,M_0$, the reference's sampling error. The test also checks:

- the result does not change from 3 to 21 planes ($8\times10^{-16}$);
- an axis toward $-z$ gives the mirrored column and $M_z \le M_0$;
- a vessel tilted by 30° matches a perpendicular one at equal $v\,\hat n_z$ ($10^{-6}$);
- the surround is at the GM steady state (0.01822, not the blood 0.01375);
- `getMz_n` is unchanged for integer $n$ and rejects non-integer $n$;
- there is no NaN.

`runSim` runs end to end with its defaults. `checkcode` is clean. Run:

```matlab
addpath(pcMRAsimDir); addpath(fullfile(pcMRAsimDir,'tests')); addpath(utilDir); test_inflowCount
```

## 6. Assumptions and limits (unchanged)

- **Rectangular slab profile** and ideal spoiling. Ashenagar et al. (2025) list the step-function
  profile as a limitation of their model too.
- **Steady flow.** Each isochromat's count uses the current velocity. For pulsatile flow the count
  depends on the velocity history over the transit time; Ashenagar et al. (2025) show that the
  inflow signal then behaves like a low-pass-filtered version of the velocity.
- **Fresh blood upstream** and a single excited slab equal to `fovSLC`. With `matSLC > 1` the
  whole FOV is treated as one slab, not as separately excited 2D slices.
- **Cell average along the slice only.** The radius, and so $v$, is evaluated at the isochromat
  centre. In-plane partial volume comes from the in-plane grid.
- **Lumen edge with `WT = 0`.** An isochromat exactly at $r = ID/2$ matches both the lumen and the
  tissue mask, and the tissue value is kept. This is unchanged here and filed as
  [#3](https://github.com/Proulx-S/pcMRAsim/issues/3).

## 7. Affected branches and callers

Checked on 2026-11-03 on every `origin` branch:

- **The continuous count** (`Q1.^n` in `getMz_n`, depth measured from the first plane in
  `simVesselSpins`) exists only on **`dev-qvfMRI`**. The other branches have no `getMz_n`; they use
  util `getMz_ss` per spin (2D column model), or `inflowEnhancementBianciardi` on `main`.
- **Projects:** no project on takoyaki loads pcMRAsim `dev-qvfMRI` from a doIt.
  multiVencISMRM2026 carries project-local copies of `getMz_n.m` and `simVesselSpins.m` with the
  same continuous count; this branch does not reach them.

Related errors in the 2D column model, on other branches, not fixed here:

- util `getMz_ss`, which every branch except `main` calls from `simVesselSpins`, used
  `ceil(vCrit/vel)`. That is fixed on util `dev-erratum-inflowColumn`.
- pcMRAsim `main`'s `inflowEnhancementBianciardi` uses the continuous $u = v_c/v$ in Bianciardi
  et al.'s eq. 2.4. The investigation (§8) found that to be off by up to +2 % of the excess.
  util `getMz_ss` on `dev-erratum-inflowColumn` is the exact replacement.

Found in the same code and filed separately:

- the lumen/tissue mask overlap, [#3](https://github.com/Proulx-S/pcMRAsim/issues/3);
- `costFun_inflow` setting fields that the simulation never reads,
  [#4](https://github.com/Proulx-S/pcMRAsim/issues/4).

## 8. How this was found

The error turned up on 2026-11-02, while the vesTrack project was porting this 3D model to
synthesize vessel time series. The question was whether the existing inflow models account for
the spread of pulse counts along the slice. The investigation:

- derived the column average from first principles;
- ran the real `dev-qvfMRI` functions against a pulse-by-pulse recursion;
- compared them with the 2D column variants in util;
- read the sources below in full text.

Its write-up is kept unedited as a frozen record:
[`docs/errata/inflowCheck-2026-11-02.md`](errata/inflowCheck-2026-11-02.md).

## References

How each source was read:

- **Full text, 2026-11-03:**
  - Bianciardi et al. (2016), PMC4822444, equations read from their images;
  - Ashenagar et al. (2025), PMC12319909;
  - Söderström et al. (2026), PMC13418971;
  - Brown et al. (2014) §24.1 and Gao et al. (1996), from the PDFs.
- **Full-text excerpt only, 2026-11-03:** Markl et al. (2003), the Methods passage quoted above.
- **Full text, 2026-11-02** (the investigation in §8):
  - Gullberg et al. (1988), pp. 437–439 only; its abstract was re-read on 2026-11-03;
  - Kim & Parker (2012), pp. 39–42 only.
- **Not read:** Gao, Holland & Gore (1988, *Med. Phys.* 15:809), which Gao et al. (1996) cite for
  the details of their eq. 2. It is not used as a source here.

Every DOI was resolved on 2026-11-03, and its title and first author checked.

Ashenagar, B., Gomez, D. E. P., & Lewis, L. D. (2025). Modeling dynamic inflow effects in fMRI to
quantify cerebrospinal fluid flow. *Imaging Neuroscience*, 3, IMAG.a.9.
https://doi.org/10.1162/imag.a.9

Bianciardi, M., Toschi, N., Polimeni, J. R., Evans, K. C., Bhat, H., Keil, B., Rosen, B. R.,
Boas, D. A., & Wald, L. L. (2016). The pulsatility volume index: an indicator of cerebrovascular
compliance based on fast magnetic resonance imaging of cardiac and respiratory pulsatility.
*Philosophical Transactions of the Royal Society A*, 374(2067), 20150184.
https://doi.org/10.1098/rsta.2015.0184

Brown, R. W., Cheng, Y.-C. N., Haacke, E. M., Thompson, M. R., & Venkatesan, R. (2014). MR
angiography and flow quantification. In *Magnetic Resonance Imaging: Physical Principles and
Sequence Design* (2nd ed., pp. 701–737). Wiley. https://doi.org/10.1002/9781118633953.ch24
(textbook)

Gao, J.-H., Miller, I., Lai, S., Xiong, J., & Fox, P. T. (1996). Quantitative assessment of blood
inflow effects in functional MRI signals. *Magnetic Resonance in Medicine*, 36(2), 314–319.
https://doi.org/10.1002/mrm.1910360219

Gullberg, G. T., Simons, M. A., & Wehrli, F. W. (1988). A mathematical model for signal from spins
flowing during the application of spin echo pulse sequences. *Magnetic Resonance Imaging*, 6(4),
437–461. https://doi.org/10.1016/0730-725X(88)90481-X

Kim, S.-E., & Parker, D. L. (2012). Time-of-flight angiography. In J. C. Carr & T. J. Carroll
(Eds.), *Magnetic Resonance Angiography: Principles and Applications* (pp. 39–50). Springer.
https://doi.org/10.1007/978-1-4419-1686-0_2 (book chapter)

Markl, M., Alley, M. T., Elkins, C. J., & Pelc, N. J. (2003). Flow effects in balanced steady
state free precession imaging. *Magnetic Resonance in Medicine*, 50(5), 892–903.
https://doi.org/10.1002/mrm.10631

Söderström, P., Björnfot, C., Andersson, B. M., Malm, J., Eklund, A., & Wåhlin, A. (2026).
Quantifying cardiac, respiratory, and low frequency components of CSF motion from fMRI inflow
effects. *Magnetic Resonance in Medicine*, 96(4), 1916–1928. https://doi.org/10.1002/mrm.70438
