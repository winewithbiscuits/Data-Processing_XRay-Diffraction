# Local Background & Noise Estimation for Spot Finding

References: 

- Barty, A., Kirian, R.A., Maia, F.R.N.C., Hantke, M., Yoon, C.H., White, T.A., Chapman, H.N. (2014)
- Cheetah: software for high-throughput reduction and analysis of serial femtosecond X-ray diffraction data, J. Appl. Cryst. 47, 1118–1131 (origin of peakfinder8); Kabsch, W. (2010) Acta Cryst. D66, 133–144 (XDS INIT/COLSPOT/DEFPIX);
- Hadian-Jazi, M. et al. (2021) J. Appl. Cryst. 54 ("Data reduction for serial crystallography using a robust peak finder" — the PF8-vs-robust comparison paper);
- XDS documentation (xds.mr.mpg.de)

## Restating the chicken-and-egg problem precisely

You want two numbers per pixel (or per local region): a background level λ̂_b and a noise scale σ, so you can threshold (N_i − λ̂_b) > n·σ. But the most natural way to estimate λ̂_b — average the pixels around the candidate — is corrupted by the fact that some of those surrounding pixels might themselves belong to a real Bragg peak. A true signal pixel included in the background window pulls the mean up and the variance up, which simultaneously (a) overestimates the background level, causing you to under-subtract real signal, and (b) overestimates σ, raising the effective threshold and causing you to miss weak neighboring peaks. This is a contaminated-population estimation problem, and every method discussed below is a different strategy for estimating "the mean and spread of the background-only subpopulation" when you can't yet label which pixels are background versus signal — that's precisely the labeling task you're trying to solve in the first place.

## Sigma-clipping — the generic algorithm

The formal iterative procedure, applicable whether your "local region" is a square window, a detector module, or a resolution ring:

S₀ = { all pixel values in the region }

repeat for t = 0, 1, 2, ...:

    μ_t = mean(S_t)
    σ_t = std(S_t)
    S_(t+1) = { x ∈ S_t : x ≤ μ_t + k·σ_t }        (k typically 2–4)
until S_(t+1) = S_t  (no further pixels rejected), or a fixed max-iteration count
output: λ̂_b = μ_final,  σ = σ_final

**Why clip only the high side.** The physical picture from the Poisson derivation two turns back tells you exactly why: background+noise is (to good approximation) a single-sided-contaminated population — a true Bragg peak can only add photons on top of the background rate, it can never subtract them. So any genuine outlier in the fluctuation distribution that arises from real diffraction signal must lie above the mean, never below. Low-value outliers (dead pixels, masked gaps, detector module edges, beamstop shadow) are a different failure mode — handled by a static bad-pixel/trusted-region mask applied before sigma-clipping even starts, not by symmetric clipping. Symmetric (two-sided) clipping would be the wrong tool here and would additionally introduce a small systematic underestimate of σ, since truncating both tails of a distribution mechanically shrinks its sample variance even when there's no real contamination to remove.

**Why this converges, and converges fast.** In any reasonably sized region, true background pixels vastly outnumber peak-contaminated pixels (a Bragg spot occupies a tiny fraction of a resolution ring's or a window's total pixel count). Each pass removes the worst offenders and pulls μ_t, σ_t back toward the true background statistics, which in turn tightens the threshold μ_t+k·σ_t and catches slightly-less-extreme contaminants on the next pass. In practice this stabilizes in 2–4 iterations — you're not searching a large space, you're peeling off a small, well-separated high-value tail from an otherwise clean population. Formally this is a special case of an M-estimator / iteratively-reweighted robust location-scale estimator (closely related to a trimmed mean or, in the limit, to the Expectation-Maximization-style adaptive MeanShift procedure — this is literally the comparison the peakfinder8 evaluation literature draws: the algorithm is described as similar to adaptive MeanShift, an expectation-maximization method that iteratively updates the model prior to outlier detection). 

https://journals.iucr.org/j/issues/2021/05/00/te5078/index.html

**Window/ring-size trade-off — the bias-variance argument, made concrete.**
From the earlier Poisson derivation, the standard error of a sample mean over N pixels is σ_pixel/√N. So:

- Too small a region (few tens of pixels): σ(λ̂_b) = σ_pixel/√N is itself large — your background estimate is noisy just from having too few samples, degrading the very threshold you're trying to set precisely (and this noise in λ̂_b propagates directly into the reflection-intensity variance formula from before: Var(Î_signal) includes a λ̂_b/N_bg term that grows as N_bg shrinks).
- Too large a region: bias, not variance, dominates — background is not actually constant across a big region (it varies with distance from the beam center / resolution), so you're now averaging over a population whose true mean itself varies across the window, producing a systematically wrong local estimate regardless of how many pixels you have.

This is the exact reason the two tools below choose different region shapes: they're both solving this same bias-variance trade-off, but under different constraints on how much real data they have available to pool.

 ## XDS's approach — an empirically-measured, per-pixel 2D table (not a parametric model)

XDS has an advantage CrystFEL's still-shot problem doesn't: with rotation data, you have multiple real images of the same detector before you even start finding spots, so you don't need to assume any particular functional form (radial or otherwise) for the background — you can measure it directly, pixel by pixel, from real frames.

INIT step. XDS builds the background lookup table BKGINIT.cbf (BKGINIT.pck in older versions) directly by estimating the initial background at each pixel from a small number of data images specified by the BACKGROUND_RANGE parameter, obtained by adding the X-ray background from each of those images. Concretely, BACKGROUND_RANGE= 2 6 tells XDS to build the initial background table from images 2 through 6. Because this table is a genuine per-pixel empirical measurement rather than a fitted radial/parametric curve, it automatically captures whatever detector module gain differences, radial scatter falloff, beamstop shadow geometry, or panel-gap structure actually exist on that specific detector — nothing about circular symmetry needs to be assumed. This is the sense in which XDS operates "per module": not because it explicitly bins by module, but because its background table has full per-pixel resolution, so module boundaries, gain differences, and any non-radial structure fall out for free

One caveat worth being precise about (correcting the looser gloss from a couple turns ago): the INIT-step background construction as documented is a straightforward accumulation across the BACKGROUND_RANGE images, not an explicit multi-pass sigma-clipped rejection loop the way peakfinder8's ring statistics are. It relies instead on the fact that, within a small handful of rotation images, any single pixel is background on most of those images and only briefly touched by a transient Bragg reflection sweeping through — so contamination is naturally diluted by construction rather than explicitly clipped out. Robustness against genuinely bad data is handled separately: badly spoiled images (e.g., from insufficient electromagnetic-pulse shielding) included in this range can corrupt the whole background table, which is why users are expected to inspect it

DEFPIX step. This refines the table by identifying and masking untrusted regions — DEFPIX recognizes regions in the initial background table that are obscured by intruding hardware (e.g. the beamstop) and marks the shaded pixels as untrusted. It does this via a value-range test: pixels are classified as unreliable if the initial-background control image falls outside VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS, defaulting to 6000–30000, with unshaded, trusted pixels typically showing control-image values around 10000 and shaded/obscured regions falling below roughly 7000. This is the equivalent of the static bad-pixel mask mentioned in §2 — it's applied as its own dedicated step, separately from the noise-based clipping used to find actual spots.

COLSPOT step — the actual thresholding. Given the (now-masked) background table, COLSPOT classifies each pixel as background or "strong" using two decision constants: spots are defined as sets of strong pixels that are adjacent in three dimensions, with the classification of a strong pixel controlled by the decision constants STRONG_PIXEL and BACKGROUND_PIXEL (in current XDS, the older STRONG_PIXEL parameter has been renamed SIGNAL_PIXEL). A real run shows these as concrete n·σ-style multipliers — a representative log shows BACKGROUND_PIXEL=6.00 and SIGNAL_PIXEL=3.00, i.e., pixels within 6σ of background are treated as background for statistics purposes, and pixels exceeding roughly 3σ above background are flagged strong. Critically, the connected-component step operates in 3D — x, y, and frame number — exactly matching the discussion from two turns ago about spots spanning several rotation images: strong pixels adjacent across x, y, and z (frame) are grouped into a single spot, then accepted only if the group contains a minimum number of pixels (MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT) and its centroid falls close enough to its strongest pixel (SPOT_MAXIMUM-CENTROID) — this connectivity+size requirement is exactly the false-positive-rate control mechanism derived mathematically in §6 of the previous turn. IUCr + 3

https://journals.iucr.org/d/issues/2010/02/00/dz5179/

This pipeline has a documented failure mode that directly illustrates the "assumptions break" theme: sharp intensity edges such as ice rings can cause an excessive number of pixels to be erroneously classified as strong, contaminating the spot list — because an ice ring is a real, sharp, resolution-localized intensity jump, not random noise, so it defeats the STRONG_PIXEL/BACKGROUND_PIXEL threshold at exactly the pixels around that ring regardless of how well the background table elsewhere is estimated

## CrystFEL peakfinder8 — pooled statistics on a resolution ring, per shot

A single XFEL still gives you exactly one frame to work with — there's no multi-image averaging trick available the way XDS uses BACKGROUND_RANGE. peakfinder8 solves the small-sample problem from §2 not by pooling across frames, but by pooling spatially, around a ring, exploiting the physical fact (established two turns ago) that background from air/parasitic/solvent scatter is, to first order, a function of scattering angle alone.

The ring-statistics model. Peakfinder8 begins by modeling the pixel intensities on one resolution ring around the center of the diffraction pattern using a single average value μ_X and a scale σ_X describing the Gaussian noise on that ring. A ring at even a modest radius on a modern megapixel detector contains thousands of pixels — which is exactly what makes single-shot, single-frame background estimation statistically viable at all: the √N denominator in the standard-error formula from §2 becomes large enough for a stable μ_X, σ_X even though you only get one look at the data.

https://journals.iucr.org/j/issues/2021/05/00/te5078/index.html

Outlier rejection, iteratively, on the ring population. This is a direct instance of the §2 sigma-clipping algorithm, applied per ring rather than per rectangular window: peakfinder8 uses non-robust statistics (ordinary mean/std) combined with careful, iterative algorithmic outlier removal to keep bright, peak-contaminated pixels from biasing μ_X, σ_X upward — confirmed independently by a later paper describing the same underlying method explicitly as a sigma-clipping algorithm that produces the background average and deviation, used to identify pixels likely belonging to Bragg peaks, the same approach peakfinder8 uses. 

https://journals.iucr.org/j/issues/2025/01/00/jo5109/

**Why radial pooling specifically, and its limits.** This algorithm was brought into CrystFEL from Cheetah because it performs better than alternatives when the background varies radially but retains approximate circular symmetry. The practical payoff is described directly: peakfinder8 uses the radial average and standard deviation to set a dynamically varying threshold, which lets it pick up weak peaks at low resolution without being swamped by false positives from the strong solvent/"water" ring — precisely the position/resolution-dependent-background problem the whole discussion started from, solved here by making the threshold itself a function of ring radius rather than a single global number. The explicit trade-off is also documented: radial averaging is simple and computationally cheap, but it assumes the background is genuinely radially symmetric, an assumption that breaks down when the sample environment (e.g. a liquid jet, membrane, or cell) introduces anisotropic scattering features — this is why CrystFEL geometry/mask files and per-panel bad-pixel handling still matter on top of the ring model; radial symmetry is a first-order approximation, not an exact physical law, exactly as flagged when this was first introduced two turns ago.

https://journals.iucr.org/d/issues/2019/02/00/ba5291/

https://www.desy.de/~barty/cheetah/Cheetah/Best_practice.html

https://arxiv.org/html/2605.19199

Local peak validation on top of the ring threshold. Passing the ring-based n·σ test isn't sufficient on its own to call a pixel a peak — a related implementation of the same family of algorithm adds explicit local-neighborhood checks: within a small local patch (typically 3×3 or 5×5 pixels), the algorithm verifies that the candidate pixel is the local maximum and that enough neighboring pixels also satisfy the signal-to-noise condition before accepting it as a peak — the ring-level μ_X/σ_X sets where the bar is, but the local patch check is what actually enforces the "connected cluster, not a lone spike" requirement discussed in §6 of the previous turn. The real implementation exposes this directly as tunable parameters — the C-level function signature includes ADCthresh, hitfinderMinSNR, hitfinderMinPixCount, hitfinderMaxPixCount, and hitfinderLocalBGRadius, mapping cleanly onto: an absolute floor threshold, the ring-based SNR cut, the minimum/maximum connected-pixel-count bounds (the same false-positive-control logic as XDS's MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT), and a local radius parameter refining the background estimate right around the candidate on top of the coarser ring statistic.

https://journals.iucr.org/j/issues/2025/01/00/jo5109/

https://www.desy.de/~barty/cheetah/Cheetah/SFX_hitfinding.html

## Why the two designs diverge — and what they share

Both are, underneath, solving exactly the estimation problem in §1 with exactly the statistical tool in §2 (robust local mean+σ via outlier-suppressed averaging) — they differ only in which axis they pool over to get enough samples for a stable estimate: XDS pools over time (multiple real frames of the same physical pixel), CrystFEL pools over space (many different pixels sharing, to first order, the same expected background), because time-pooling simply isn't available when every shot is a different, randomly-oriented crystal.
