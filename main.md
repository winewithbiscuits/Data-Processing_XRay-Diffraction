# Summary of Data Processing

Poisson-Distributed Signal.md

XDS was built for rotation/oscillation data (one physical crystal, continuous ω-sweep, single detector), and CrystFEL was built for serial crystallography (thousands of still, randomly-oriented crystal snapshots at XFELs or serial synchrotron setups, no rotation, partial reflections by default). 

The physical problem, briefly

A crystal is a 3D periodic array of unit cells. X-rays scattered by that periodicity interfere constructively only at points on a reciprocal lattice — the Ewald sphere construction: reciprocal lattice points that intersect a sphere of radius 1/λ produce a diffraction spot. As the crystal rotates (rotation method) or 
as many differently-oriented crystals are fired through (serial method), different reciprocal lattice points cross the sphere, producing a stream of spots on the detector.
The data-processing pipeline exists to invert this: from a set of 2D spot positions/intensities on a detector, recover (a) the unit cell and orientation matrix (indexing), 
(b) accurate integrated intensities per reflection (integration), 
(c) a self-consistent, correctly scaled, merged dataset of unique reflection intensities |F(hkl)|² suitable for phasing/refinement (scaling & merging).


# Spot finding (peak/spot search)

Goal: locate pixel clusters on each detector image that are Bragg peaks, not background, cosmic rays, ice rings, or detector artifacts, and record centroid position + approximate intensity for each.

Core algorithmic logic (shared across tools):

1. Background/noise modeling per pixel. Raw pixel counts = background (diffuse scatter, air scatter, detector dark current, fluorescence) + Poisson-distributed signal. You can't threshold on raw counts because background varies smoothly across the detector (stronger near the beam center, varies with resolution shell). So first estimate a local background level and local noise (σ), typically via:

- A local mean/variance in a sliding window or resolution-shell-based statistics, iteratively excluding outlier (high) pixels so peaks don't bias the background estimate (sigma-clipping).
- XDS does this per module/resolution shell with iterative rejection.
- CrystFEL's peakfinder8 computes background statistics per resolution ring (annulus around the beam center) since background is radially symmetric to first order — faster and more robust for the high peak-density, high-noise regime of XFEL shots.

2. Thresholding. Flag pixels where (pixel_value − local_background) > n·σ (typically n ~ 3–6, tunable). This gives a binary "candidate signal" mask.
3. Connected-component labeling. Group adjacent flagged pixels into connected blobs (4- or 8-connectivity flood fill). Each blob is a candidate spot.
4. Filtering candidate spots by shape/size heuristics: minimum pixel count (rejects hot pixels / single-pixel noise spikes), maximum pixel count (rejects saturated blooms, ice rings, or overlapping-spot blobs), roundness/eccentricity checks, and resolution-dependent size expectations (spots get smaller/more compact at high resolution due to the Ewald sphere geometry — XDS uses this).
5. Centroid calculation. For each accepted blob, compute an intensity-weighted centroid (x̄, ȳ) = Σ(I_i · x_i)/ΣI_i over the pixels in the blob — this sub-pixel position is what indexing will use. Some tools also fit a 2D Gaussian or use moment analysis to get more precise centroid + spot shape (σx, σy, correlation) for later use in integration profile modeling.
6. (For rotation data) A spot's centroid also has a z-coordinate — the frame number/rotation angle at which it's brightest — found either from a single strongest frame or by 3D connected-component labeling across frames (spots can span multiple adjacent oscillation images).

## XDS specifics (COLSPOT step)

XDS's COLSPOT uses local background subtraction per pixel based on statistics computed in INIT (which builds a background/gain image from a subset of frames), applies the threshold SIGNAL_PIXEL / SPOT_MAXIMUM-CENTROID, and does 3D (x,y,frame) connected component analysis to merge spot fragments seen on adjacent images into a single strong reflection with an interpolated centroid in all three dimensions. Output: a spot list (SPOT.XDS) of (x, y, z, intensity) — this is what indexing consumes.

## CrystFEL specifics (indexamajig peak search)

Because a still XFEL shot has no rotation (each image is one 2D slice through reciprocal space, no z-refinement across frames possible for a single hit), CrystFEL's peak finders (zaef, peakfinder8, peakfinder9) are purely 2D per-frame:

- zaef — classic threshold + local-mean rejection, simple and fast but less robust to strong background gradients.
- peakfinder8 — radial (resolution-ring) background statistics, robust to the strong low-resolution background typical of XFEL data; industry-standard default for a long time.
- peakfinder9 — a more recent refinement using a local signal-to-noise + gradient-based approach with better handling of overlapping/close peaks, generally higher accuracy at the cost of more tuning parameters.

Since each XFEL shot might have hundreds of crystals worth of noise/background from jet streaming, gas dynamic virtual nozzle scatter, etc., background modeling per shot is a much harder and more central problem in CrystFEL than in XDS.

# Indexing

Goal: from the observed spot centroids (a partial, noisy sampling of the reciprocal lattice), determine (a) the unit cell parameters, (b) the crystal orientation matrix A (sometimes decomposed as U·B, orientation times reciprocal metric tensor), such that predicted reciprocal lattice vectors h = A·(h,k,l) match observed spot positions.

Core algorithmic logic

Step 1 — convert 2D spot positions to 3D reciprocal space vectors. Using the known detector geometry (distance, beam center, pixel size, tilt), beam wavelength, and (for rotation data) the rotation angle at which each spot was observed, back-project each spot centroid onto the Ewald sphere to get its reciprocal space coordinate r = (S − S₀)/λ, where S is the scattered beam unit vector and S₀ the incident beam unit vector. This gives a cloud of 3D points that should coincide with points of an (unknown) reciprocal lattice, just missing the (h,k,l) integer labels.

Step 2 — find the lattice basis from the point cloud. This is the actual "indexing" problem: given an unlabeled 3D point cloud that's (approximately) a lattice, recover 3 primitive basis vectors a* , b* , c* that generate it. Two dominant algorithmic families:

(a) Real-space / difference-vector methods (classic: DPS — Duisenberg, or the "1997 Kabsch/Steller/Rossmann" family; used by LABELIT, MOSFLM-style indexing, and as one option in CrystFEL):

- Compute all pairwise difference vectors between observed reciprocal points. Differences between points that both belong to the true lattice are themselves lattice vectors (or lattice vector combinations), so peaks recur at multiples of the true basis vectors.
- Build a histogram (in effect a Patterson-like function in reciprocal-difference space) and look for the shortest, most frequently recurring, linearly-independent set of 3 vectors — that's a candidate primitive cell.
- Reduce this to a Niggli-reduced cell (canonical smallest-volume basis) to remove basis-choice ambiguity, then compare against space-group symmetry constraints if known.

(b) Fourier-transform-based methods (XDS's IDXREF uses this — the "1993 Kabsch" algorithm, closely related to what's often called the FFT/DPS method used also by DIALS' fft3d indexer):

- For each observed reciprocal vector, and for a large set of trial 1D directions in reciprocal space, project the reciprocal vectors onto that direction and compute a 1D Fourier transform of the resulting projected-point density.
- A trial direction that's parallel to a true lattice row will show strong periodicity (sharp peaks) in this projection — because true lattice points that lie on lines parallel to a true axis are evenly spaced.
- Search over trial directions (a coarse grid over the unit sphere, then refined) to find the 3 directions/periods giving the strongest, mutually independent periodicity signals → these become candidate basis vectors.
- Reduce to Niggli cell, test cell candidates against Bravais lattice constraints (XDS reports candidate lattices ranked by "quality of fit" / distortion index from ideal cell symmetry), and let the user or software pick the correct/expected Bravais lattice.

(c) Global/exhaustive & "brute force" methods (relevant for CrystFEL's harder problem of indexing every single still shot fast and often from very few spots):

- xgandalf — extends the FFT gradient-descent idea, refines using a fast approximate nearest-neighbor lattice-fit gradient descent from many random starting orientations; good balance of speed and robustness, now a common CrystFEL default.
- pinkindexer — designed for pink-beam (broad bandwidth) data; solves indexing without a single fixed wavelength assumption via a Radon-transform-like search over orientation space.
- taketwo — builds up an orientation matrix incrementally from pairs/triplets of spots consistent with the known cell, using graph-based consistency checking; strong when the unit cell is already known (common in high-throughput serial pipelines where you index the same crystal type repeatedly).
- mosflm/dirax/xds indexers wrapped inside CrystFEL — CrystFEL can call out to several external indexing engines and try them in sequence/parallel per pattern until one succeeds, since no single indexer works on 100% of noisy still shots. indexamajig typically runs a cascade like xgandalf,mosflm,xds and takes the first internally-consistent solution.

Step 3 — refine the orientation/cell. Once a trial basis is found, predict where all reflections up to some resolution should fall given that cell + orientation + detector geometry, match predictions to observed spots, and least-squares refine the 9 (or fewer, with symmetry constraints) unit-cell/orientation parameters (plus, iteratively, detector geometry parameters — distance, beam center, tilts) to minimize the position residuals (Δx, Δy) and, for rotation data, Δφ (angular residual). This is a nonlinear least-squares problem, normally Levenberg-Marquardt or straightforward Gauss-Newton, iterated a few cycles.

Step 4 — symmetry/Bravais lattice determination. The primitive triclinic cell found by pure indexing math doesn't know about crystal symmetry. Tools compare the reduced cell's metric to the 44 conventional Niggli-reduced cell types (Burzlaff/de Wolff/International Tables classification) to suggest the likely Bravais lattice (and hence possible space groups), letting the user confirm or override. XDS reports a table of candidate lattices ranked by a "quality" (distortion from ideal symmetric metric) — you pick the one matching your known/expected space group. CrystFEL does the equivalent via cell_tolerance matching against a target cell you supply (--pdb=), since with hundreds of thousands of stills you cannot manually vet each one.

# Integration

Goal: given indexed orientation + cell + detector geometry, for every predicted reflection (h,k,l), extract a single reliable intensity I(hkl) and its uncertainty σ(I) by summing/fitting signal in the pixels around the predicted spot position, correctly separating signal from local background.

## Two competing approaches

(a) Summation integration. Define a fixed-shape box/ellipse around the predicted spot center, sum pixel counts inside it, subtract a background estimated from an annulus/ring around that box (background is fit, often as a locally-planar or locally-constant function, from surrounding non-peak pixels). Simple, robust for strong reflections, but statistically suboptimal — every pixel (even in the low-signal tails of a weak, noisy reflection) is weighted equally, which is not the maximum-likelihood-optimal weighting.

(b) Profile-fitting integration (the standard in serious modern pipelines — this is what makes XDS/DIALS/etc. much better than naive summation). The key insight (from Kabsch's and Rossmann's work, and originally Ford's 1974 profile-fitting method): reflections close together in the same region of the detector, at similar resolution and rotation angle, share nearly the same intrinsic spot shape (the "reflection profile"), because that shape is dominated by beam divergence, crystal mosaicity, and detector PSF — properties that vary smoothly, not randomly, across the image. So:

1. Build a 3D reference profile (a normalized 2D-spatial × 1D-rotation intensity distribution, or a "learned average shape") from strong, well-measured, non-overlapping reflections in a local region of detector + rotation range.
2. Scale that reference profile (single multiplicative scale factor = the reflection's integrated intensity) to best fit the observed pixel counts of a given (possibly weak) reflection, via a maximum-likelihood or least-squares fit — essentially projecting the observed pixel-count vector onto the known profile-shape vector, which is statistically the minimum-variance unbiased estimator (this is literally why it beats summation for weak reflections: the profile carries prior information about where the signal should be concentrated within the box, downweighting noisy tail pixels).
3. This profile is typically stored as a lookup varying smoothly across detector position and rotation angle ("learned profiles" partitioned into batches), updated as integration proceeds through the image sequence.

XDS's INTEGRATE step implements exactly this 3D profile fitting (spatial × φ), which is why XDS output typically outperforms simple box summation especially in the weak/high-resolution shells. It also does both — reports summation-integration intensity too, for comparison/QC.

## Handling partiality

For rotation data, a reflection's full intensity is often spread across 2–3 consecutive oscillation images because the reciprocal lattice point takes a finite rotation range to fully cross the Ewald sphere's finite-width "reflecting condition." XDS profile-fits and sums the partial contributions across the relevant frames into one full reflection intensity (using the Lorentz factor and rocking-curve-width formalism to combine partial measurements correctly).

For still shots (CrystFEL/serial crystallography) there is no rotation to sweep through — a still image intersects the Ewald sphere at one instant, so essentially every reflection is a partial measurement of the "true" fully-integrated structure factor, with a partiality factor 0 < p ≤ 1 depending on how close the reciprocal lattice point's excitation error is to zero (how exactly it satisfies the diffraction condition) and on the crystal's mosaicity/domain size. indexamajig estimates per-reflection partiality using a chosen partiality model (e.g. "unity," or physical models like Monte Carlo / dynamical models available in more recent CrystFEL versions), and integration here uses fixed small boxes around predicted spot positions (since orientation is known from indexing) with the same profile-fit-vs-summation choice available (--integration=rings or prof2d etc.). This partiality problem is precisely why serial-crystallography merging needs an extra step (post-refinement) that rotation data doesn't strictly need.

# Scaling and merging

Goal: combine many individual, partial, differently-scaled measurements of the same unique reflection (from different images/frames/crystals) into one best-estimate merged intensity per unique hkl, correcting for systematic experimental effects (beam intensity drift, absorption, radiation damage, differing partiality, detector nonuniformity) — this is what actually produces the dataset you hand to phasing/refinement software.

## Core algorithmic logic (rotation data — XDS's CORRECT/XSCALE)

1. Apply per-reflection corrections first: Lorentz-polarization factor (geometric correction for scattering geometry and beam polarization), absorption correction (empirical, based on symmetry-equivalent reflection comparison across different crystal orientations — this is essentially what XDS's "empirical absorption correction" via spherical harmonics does), detector-response/geometric corrections.
2. Scale model. Assume each measurement's true intensity I(hkl) is related to observed intensity I_obs by I_obs = K(t) · I(hkl) · exp(−2B(t)·sin²θ/λ²) roughly, where K(t) is a smoothly time/frame-varying scale factor (accounts for beam intensity drift, illuminated-volume changes) and B(t) is a smoothly-varying relative B-factor (accounts for global effects that fall off with resolution, notably radiation damage — as damage accumulates, high-resolution intensities fade faster than low-resolution ones). These scale/B parameters are modeled as smooth (often B-spline or low-order polynomial) functions of batch/frame number, not independently free per frame, or the fit is underdetermined.
3. Iterative least-squares refinement. Initialize all scale factors to 1, compute a provisional merged average per unique hkl (averaging over symmetry equivalents too, using the space group), then refine K(t), B(t) to minimize Σ_hkl Σ_i w_i·(I_i(hkl) − K(t_i)·exp(...)·⟨I(hkl)⟩)² , recompute the merged average, iterate to convergence. This is the classic scaling algorithm going back to Hamilton/Rollett/Sutton and refined by Fox & Holmes, implemented by XDS's XSCALE and equivalently by AIMLESS in the CCP4 world (Evans' algorithm, similar idea with Bayesian/robust refinements).
4. Outlier rejection. After a provisional scale model, flag and downweight/reject measurements deviating too many σ from their symmetry-equivalent group mean (protects against misindexed reflections, detector artifacts, ice ring contamination).
5. Final merge. Weighted mean of all corrected, scaled, symmetry-equivalent measurements per unique hkl gives the merged I(hkl) and its σ, typically weighted inverse-variance (w = 1/σ²), sometimes with additional robust down-weighting.

## Serial crystallography specifics — why "merging" ≠ "scaling" here (CrystFEL's partialator)

Because each still crystal is a different physical crystal with random, uncorrelated orientation and unknown absolute scale (unlike a single rotating crystal where frame-to-frame drift is small and smooth), you cannot merge naively by simple weighted averaging of raw partial intensities — partiality must be modeled and refined per-crystal first. This is called post-refinement:

1. For each indexed crystal/pattern, refine per-crystal parameters: orientation (small corrections), unit cell, and crystal-specific partiality-model parameters (mosaic domain size, reflecting-range/"profile radius," and an overall relative scale factor for that crystal) — refined by requiring that, once partiality-corrected, each crystal's reflections agree with a running estimate of the merged "true" structure factor amplitudes.
2. This is inherently a chicken-and-egg / EM-style iterative problem: you need the merged reference intensities to refine per-crystal partiality, but you need per-crystal-corrected intensities to build the merged reference. partialator iterates: merge → refine per-pattern parameters against current merge → re-merge → repeat, similar in spirit to the rotation-data scale-refinement loop but with many more free per-crystal parameters (potentially thousands of crystals × several parameters each) and an explicit physical partiality model rather than a smooth scale curve.
3. Because a single still crystal only samples each reflection partially/once, you fundamentally need many crystals (hundreds to thousands) per unique reflection to get a statistically solid merged intensity — this is why serial datasets need such high multiplicity/redundancy compared to rotation datasets.
4. partialator supports multiple partiality models (unity = ignore partiality/treat as full, useful as a crude baseline; physical models that use the refined profile radius/reflecting range) and refinement targets (least-squares residual to merged reference, or maximum-likelihood formulations).

# Data-quality metrics used to judge every stage above

These are the numbers you actually look at to decide if indexing/integration/merging worked and where to cut resolution:

- R_merge / R_sym — Σ|I_i − ⟨I⟩| / ΣI_i over symmetry equivalents; classic but biased upward by high multiplicity (more measurements → more deviation terms summed → inflated R even for good data), so largely superseded by:
- R_meas / R_pim — multiplicity-corrected versions of R_merge (R_meas corrects for the redundancy-dependent bias; R_pim isolates the precision of the merged mean specifically).
- CC1/2 — Pearson correlation coefficient between two randomly-split half-datasets' merged intensities, per resolution shell; the modern gold-standard for deciding the "true" resolution cutoff (Karplus & Diederichs, 2012) because unlike R-factors it doesn't blow up simply from high multiplicity/noise the same way and has a clear statistical meaning (correlation with the "true" signal).
- *CC ** — derived from CC1/2, estimates correlation of your merged half-dataset with the true underlying (infinite-multiplicity) dataset — used to judge whether adding more data would actually still help.
- I/σ(I) — per-shell signal-to-noise, classic simple cutoff criterion (commonly resolution is cut around I/σ ≈ 1–2, though CC1/2-based cutoffs are now preferred and often push usable resolution further than the old I/σ>2 rule suggested).
- Completeness & multiplicity per resolution shell — fraction of unique possible reflections actually measured, and average number of independent observations per unique reflection; both indexing failures and geometric blind-spots (e.g., along the rotation axis in rotation data) show up as completeness gaps.
- Rsplit — CrystFEL's serial-crystallography analogue of R_merge computed between two random half-datasets, used because R_merge's usual per-crystal-scale assumptions don't map cleanly onto still-image partial data.

# Putting it together — the two concrete pipelines

### XDS (rotation data), step by step, keyword-named:

XYCORR (detector geometric distortion correction) → INIT (background/gain/bad-pixel map from a few images) → COLSPOT (spot finding, produces SPOT.XDS) → IDXREF (indexing: FFT-based basis search, Niggli reduction, Bravais lattice ranking, geometry refinement) → DEFPIX (refine bad-pixel/background mask using now-known reflection positions to exclude spot pixels from background stats) → INTEGRATE (3D profile-fitting integration over all images) → CORRECT (apply Lorentz-polarization/absorption corrections, do provisional scaling, produce XDS_ASCII.HKL) → then XSCALE (multi-dataset scaling and merging into a final merged HKL file) → downstream POINTLESS/AIMLESS (CCP4) or XPREP etc. for space-group confirmation and final merging/statistics, then phasing (SHELX, Phaser) or refinement (REFMAC, phenix.refine).

### CrystFEL (serial/still data), step by step:

indexamajig does spot-finding + indexing + integration all in one pass per image (calls peak finder → tries cascade of indexing engines → integrates predicted reflections with chosen partiality/profile model) → outputs a stream file of per-crystal reflection lists with orientations → partialator performs post-refinement (per-crystal scale + partiality parameters) and merges into final reflection list → check_hkl/compare_hkl/get_hkl and stats.py-style tools compute CC1/2, Rsplit, completeness per shell → merged intensities go on to phasing/refinement same as above.


