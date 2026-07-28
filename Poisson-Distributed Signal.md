# **Poisson-Distributed Signal in Diffraction Images**

References:

1. Feller, W. An Introduction to Probability Theory and Its Applications, Vol. 1 (for the formal Poisson process derivation);
2. Rossmann, M.G. & Erickson, J.W. (1983) J. Appl. Cryst. 16, 629–636 (counting statistics in oscillation crystallography);
3. Diederichs, K. (2010) Acta Cryst. D66, 733–740 ("Quantifying instrument errors in macromolecular X-ray data sets");
4. French, S. & Wilson, K. (1978) Acta Cryst. A34, 517–525 ("On the treatment of negative intensity observations" — the classic paper on what Poisson noise does to weak/background-subtracted intensities);
5. Leslie, A.G.W. (1999) Acta Cryst. D55, 1696–1702 (background modeling, integration).

## Why photon detection is a Poisson process — the physical origin

X-ray photons arrive at a detector pixel as discrete, quantized events. Whether a given pixel registers N photons during an exposure of duration T depends on two independent stochastic layers:

Photon emission from the source (synchrotron bending magnet/undulator, or XFEL self-amplified spontaneous emission) is itself a random process — individual photons are emitted at random times with some average rate, not at fixed regular intervals.
Scattering/absorption interactions (a photon being Compton/Bragg-scattered into a given solid angle, or absorbed in the sensor) are each independent probabilistic events, uncorrelated with whether any other photon was scattered.

### What is a Poisson process

A process satisfying three axioms is, a Poisson process:

(a) Independence — the number of events in disjoint time (or space) intervals are statistically independent.

(b) Stationarity / proportionality — the probability of exactly one event in a small interval [t, t+dt] is λ·dt + o(dt), where λ is a constant rate (constant over the interval considered).

(c) Rarity of coincidence — the probability of two or more events in the same infinitesimal dt is o(dt) (negligible compared to dt itself).

**Formal derivation of the Poisson distribution**

Derive it as the limit of a binomial process. Split the exposure time T into n tiny sub-intervals of width Δt = T/n, small enough that each sub-interval has probability p = λΔt = λT/n of containing exactly one photon-detection event, and probability ≈0 of containing more than one (by axiom (c) above).
The number of sub-intervals that register a hit, out of n independent trials each with success probability p, is Binomial(n, p):

P(N = k) = C(n,k) · p^k · (1−p)^(n−k)

Now take the limit n → ∞, p → 0, holding np = λT ≡ μ fixed (μ = the expected total count):

P(N=k) = [n! / (k!(n−k)!)] · (μ/n)^k · (1 − μ/n)^(n−k)

= (μ^k / k!) · [n!/((n−k)! n^k)] · (1 − μ/n)^n · (1 − μ/n)^(−k)

As n→∞: n!/((n−k)! n^k) → 1, (1 − μ/n)^n → e^(−μ), and (1 − μ/n)^(−k) → 1.

This leaves:
P(N=k) = μ^k · e^(−μ) / k!         k = 0, 1, 2, ...

This is the Poisson distribution with parameter μ (the expected count). μ plays the role of the true underlying signal level you're trying to measure; k is what you actually observe on a given exposure.

**The defining statistical property: mean = variance**

Compute the mean and variance directly from the distribution (via the moment generating function M(t) = E[e^(tN)] = exp(μ(e^t − 1)), or directly by summation):

E[N] = μ

Var[N] = μ

Mean equals variance. This single fact is the entire basis of "shot noise" statistics in X-ray detection: if you observe a pixel with N counts, your best estimate of the true rate is μ̂ = N, and your best estimate of the uncertainty on that measurement is:

σ(N) = √μ ≈ √N


This is why you'll see "σ = √I" used everywhere in crystallographic data processing without further derivation — it is not a modeling choice, it's the exact variance of a Poisson-distributed count, valid down to N=0 (unlike a Gaussian model, which would allow negative counts and doesn't have variance tied to the mean).

Signal-to-noise ratio for a raw photon count is therefore:

SNR = μ/σ = μ/√μ = √μ

— SNR grows only as the square root of exposure/photon count. Doubling SNR requires quadrupling exposure time/flux, not doubling it. This √N scaling is the fundamental reason serial crystallography and weak/high-resolution data collection are so photon-hungry, and it's the quantitative reason longer exposures or brighter sources give diminishing returns per unit time spent.

**Sums of independent Poisson variables are Poisson**

A pixel's total observed count is the sum of background (diffuse scatter, air scatter, fluorescence, dark current) and signal (Bragg diffraction), each itself Poisson-distributed with its own rate:

N_pixel = N_background + N_signal,     N_background ~ Poisson(λ_b),  N_signal ~ Poisson(λ_s)

Claim: the sum of two independent Poisson random variables is itself Poisson, with rate = sum of rates.

Proof via MGFs: M_N(t) = M_b(t)·M_s(t) = exp(λ_b(e^t−1))·exp(λ_s(e^t−1)) = exp((λ_b+λ_s)(e^t−1)), which is exactly the MGF of Poisson(λ_b+λ_s). So:

N_pixel ~ Poisson(λ_b + λ_s),   Var(N_pixel) = λ_b + λ_s

This closure property is why background-and-signal decomposition is mathematically clean: the total pixel variance is just the sum of the two rates, no cross terms, no covariance correction needed (independence is assumed and is physically reasonable — background scatter and Bragg scatter are different physical processes).

**Consequence for background subtraction.**

You never observe λ_b and λ_s separately — you estimate λ_b (from surrounding non-peak pixels, see §6) and subtract it from the total observed N_pixel to estimate signal:

Î_signal = N_pixel − λ̂_b

But the variance of this estimate is not just λ̂_b — it's the sum of the variance of the raw pixel count and the variance of your background estimate itself (since λ̂_b is itself a noisy statistic computed from a finite number of background pixels):

Var(Î_signal) = Var(N_pixel) + Var(λ̂_b) = (λ_s + λ_b) + λ_b/N_bg_pixels

where the last term comes from λ̂_b being a sample mean over N_bg_pixels background pixels, each Poisson(λ_b), so Var(mean) = λ_b/N_bg_pixels. **This is why spot-finding/integration algorithms want a large, well-sampled background region** — a background estimated from too few pixels adds its own noise on top of the unavoidable Poisson shot noise of the peak itself, degrading σ(I) even for otherwise strong spots.

## The Gaussian approximation — when it's valid, and when it isn't

By the Central Limit Theorem (or directly from the Poisson MGF via a saddle-point/Stirling expansion), for large μ the Poisson distribution converges to a Gaussian:

Poisson(μ) → Normal(μ, μ)     as μ → ∞

In practice this approximation is reasonable once μ ≳ 15–20 (skewness ~1/√μ becomes small, the distribution looks symmetric and bell-shaped). **This is why weighted-least-squares treatments (§2.2 of the earlier profile-fitting derivation) are valid approximations for strong reflections** — high photon counts, Gaussian-like Poisson, ordinary χ² minimization works fine.

**It breaks down for weak signals** (μ roughly < 10, common in high-resolution shells, weak reflections, or serial-crystallography single-shot integration where exposures are extremely short): the true distribution is discrete, strictly non-negative, and noticeably right-skewed. 

Two concrete failure modes of the Gaussian approximation in this regime:

**Negative "intensities."** Background-subtracted Î = N_pixel − λ̂_b can legitimately go negative purely from Poisson fluctuation (N_pixel came in below the true background mean by chance), even though the true underlying signal intensity is physically ≥ 0. A naive Gaussian/least-squares pipeline will happily report negative I(hkl) for weak reflections. French & Wilson (1978) is the standard treatment: they derive a Bayesian correction (using a prior that true intensity J ≥ 0) that maps the observed (possibly negative) I_obs and its σ to a properly non-negative posterior estimate of the structure factor amplitude — this is applied routinely in the final stages of processing weak/low-resolution data before it's usable for phasing (it's baked into TRUNCATE/ctruncate in the CCP4 pipeline, and equivalent French-Wilson-style treatments exist in phenix).

**Poisson-correct maximum likelihood, not Gaussian least-squares,** matters more for weak reflections — this is exactly why §2.4 of the earlier integration derivation (the IRLS Poisson-ML profile-fit estimator) exists as a refinement on top of the simpler weighted-least-squares estimator: at low counts the Gaussian-weighted estimator is measurably biased/suboptimal relative to the true Poisson likelihood.

## Application to spot-finding thresholding

Recall the spot-finding logic from earlier: flag pixel i as candidate signal if (N_i − λ̂_b) > n·σ. Now we can be precise about σ and about what "n" should be.

σ construction. Using the variance decomposition from §4:

σ_i = √( N_i + λ̂_b + λ̂_b/N_bg )   ≈   √( λ̂_b(x,y) · (1 + 1/N_bg) )   for background-dominated pixels

(the local λ̂_b(x,y) is itself position-dependent — stronger near the beam center, falls off with scattering angle — which is why background/σ maps are computed per-pixel or per-resolution-ring, not as one global number; this is the direct statistical justification for XDS's per-module background modeling and CrystFEL's peakfinder8 resolution-ring approach mentioned earlier).

**Choosing the threshold multiplier n.** Under the Gaussian approximation, the false-positive probability for a single pixel exceeding an n·σ threshold by chance is given by the standard normal tail:

P(false positive per pixel) = 1 − Φ(n)

n=3 → P ≈ 1.3×10⁻³ per pixel; n=5 → P ≈ 2.9×10⁻⁷; n=6 → P ≈ 9.9×10⁻¹⁰.

Why this matters quantitatively — multiple-testing.

A modern detector has ~ 10⁶–10⁸ pixels per image. If you used n=3σ across a 4-megapixel detector, the expected number of pure-noise pixels exceeding threshold per image is:

E[false positives] = N_pixels × P(false positive) ≈ 4×10⁶ × 1.3×10⁻³ ≈ 5,200 spurious pixels/image

— clearly unworkable. This is precisely why real spot-finders don't rely on the per-pixel threshold alone: they additionally require a minimum connected-pixel count (a single noise spike is almost never several contiguous pixels all independently exceeding threshold — the joint probability of k adjacent pixels all crossing n·σ by pure chance falls off roughly as P(single)^k, assuming rough independence between neighboring pixels' noise), which is why "connected component with ≥ M pixels" (step 3–4 of the spot-finding algorithm from earlier) is not an arbitrary heuristic but a direct, derivable consequence of wanting to control the family-wise false-positive rate across millions of pixels using per-pixel Poisson/Gaussian tail statistics. In practice this lets pipelines use a relatively permissive per-pixel n (~3–4σ) combined with a connectivity/size requirement, rather than needing an impractically strict per-pixel n alone.

Also note: because the true distribution is discretely Poisson (not exactly Gaussian) especially at low background rates (modern hybrid pixel detectors can have background rates of a few counts/pixel or less in a typical exposure — see §7), the Gaussian tail probabilities above are only approximate. Rigorous implementations use the exact Poisson CDF (or its complement) for the tail probability at low λ̂_b, since the Gaussian approximation understates the tail probability at low counts (Poisson is right-skewed — heavier right tail than a Gaussian with matched mean/variance), meaning a "3σ Gaussian" threshold is actually somewhat looser than 3σ-worth of true statistical significance would suggest when background is low.

## Detector-specific complications — Poisson isn't the whole story

**Photon-counting hybrid pixel detectors** (Pilatus, Eiger, Jungfrau — the dominant modern detector class): each pixel has its own discriminator/counter circuit that increments by exactly 1 for each photon above an energy threshold. The recorded number is a literal photon count, no gain conversion needed, negligible dark current, and essentially zero electronic read noise. This makes the pure-Poisson model in §1–6 essentially exact for these detectors — one of the major reasons they superseded CCDs/image plates for weak-signal work: the noise floor is limited almost purely by unavoidable photon shot noise, not by added instrumental noise. The one Poisson-violating correction that matters at high flux is dead-time / pile-up: if two photons arrive within the detector's minimum resolvable time window (paralyzable or non-paralyzable dead-time models), they can be miscounted as one, causing a systematic (not random) undercount at high count rates — this is a deterministic nonlinear correction applied to the raw counts, separate from the statistical Poisson noise model, and matters mainly for very strong low-resolution spots or intense XFEL pulses.

Integrating detectors (CCDs, image plates, older technology, still relevant in some legacy pipelines and CrystFEL default test cases): the sensor accumulates charge proportional to absorbed photon energy, then that charge is converted to a digital "ADU" (analog-digital unit) value via a gain factor g (electrons or ADU per photon), and the readout electronics add their own approximately-Gaussian, signal-independent read noise σ_read. The correct combined variance model (in ADU units) is:

σ²_ADU = g · S_ADU + σ²_read

where S_ADU is the background-subtracted signal in ADU (g·S_ADU recovers the Poisson-variance-in-photon-units term scaled back into ADU²), and σ_read is a fixed additive term independent of signal level. There's also a dark current contribution — thermally generated electrons accumulate over the exposure exactly like an extra Poisson process with its own rate λ_dark(T, exposure time), typically subtracted via a "dark frame" or "pedestal" calibration image taken with no beam, and its own shot noise contributes an extra λ_dark term into the total variance, exactly analogous to λ_b in §4. This is why CCD-generation processing pipelines needed careful pedestal/dark subtraction and gain calibration as explicit preprocessing steps (this is what XDS's INIT step and its "gain map"/"bad pixel mask" construction is doing under the hood) — steps that are largely unnecessary (or much simpler) for modern photon-counting detectors.

## Why this all matters downstream

Every later stage in the pipeline inherits this statistical foundation directly:

- Integration (§2 from before): the weighting w_i = 1/σ_i² in both the Gaussian weighted-least-squares profile-fit estimator and the Poisson-ML IRLS scheme comes directly from the σ² = μ relationship derived here — it's not a free/tunable parameter, it's fixed by the physics of photon counting (plus, for CCD data, the extra read-noise/gain terms from §7).
- Reported σ(I) for every merged reflection ultimately traces back to propagating these per-pixel Poisson variances through background subtraction, profile-fitting, partiality-correction, and scaling — which is why careful crystallographic software reports σ(I) as a first-class, physically meaningful quantity (used directly in CC1/2, I/σ resolution cutoffs, and refinement weighting), rather than an ad hoc error bar.
- French-Wilson correction (§5) is applied precisely because this Poisson noise model predicts — correctly — that weak reflections will sometimes show negative background-subtracted intensity, and downstream structure-factor amplitude estimation needs a statistically principled (not just clipped-to-zero) way to handle that.
