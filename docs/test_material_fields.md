# Test-material fields used by the ERM alignment

This page describes the data fields that set the **particle geometry, size
distribution and composition of a test material** in the alignment of
effect concentrations to ecologically relevant metrics (ERMs). It is kept
identical in two repositories:

- `ToMEx2.0_EcoToxRisk` (R package **PSSDplusplus**: `align_data()`,
  `align_mixture()`), and
- `ToMEx_AquaticOrganisms` (the ToMEx Shiny app and the database tidying).

Sections 2 and 3 describe **new, optional fields**: measured size
distributions and co-dosed mixtures. A record that leaves them empty is
aligned exactly as before, so adding the fields to the database changes no
existing threshold. Section 1 describes existing fields whose meaning is easy
to misread.

Units: lengths, widths and heights in µm unless the field name says `mm`;
concentrations in particles/mL (water) or particles/kg (sediment) and mg/L or
mg/kg for mass.

## 1. Existing geometry fields (how they are used)

| Field | Meaning | Used for |
|---|---|---|
| `shape_f` | `"Sphere"`, `"Fiber"`, `"Fragment"` or `"Not Reported"` (capitalised; never match lower case) | particle geometry |
| `polydispersity` | `"monodisperse"` (one size) or `"polydisperse"` (a range between a minimum and maximum length) | which formulas apply |
| `size.length.um.used.for.conversions` | particle length; for a polydisperse material the midpoint of the range | monodisperse geometry; SSD size labels |
| `size.length.{min,max}.um.used.for.conversions` | reported minimum and maximum lengths (measured, else nominal) | polydisperse range |
| `size.width.{min,max}.mm.{measured,nominal}` | **reported** minimum/maximum widths | fiber diameters (below) |
| `size.width.um.used.for.conversions` | **reported** single width (measured, else nominal; mean of min and max if only those are given) | monodisperse width; fiber diameter |
| `size.width.{min,max}.um.used.for.conversions` | **not measurements**: the tidying sets them to R.ave × length (spheres: length) | nothing in the alignment (it recomputes widths) |
| `density.g.cm3` | polymer density: the record's own value (template, or ToMEx 1.0), else the polymer table `polymer_densities.csv` of the database repository (literature values, or the median of the ToMEx 2.0 v2.0.1 template densities); `density.source` says which ("record", "polymer table") | particle mass |
| `dose.particles.mL.master`, `dose.mg.L.master` (+ sediment `dose.particles.kg.sediment.master`, `dose.mg.kg.sediment.master`) | the test's dose as particles and as mass | dose |
| `dose.particles.mL.master.converted.reported` (+ sediment equivalent) | `"reported"`: the study gave a particle count; `"converted"`: the database derived it from the dosed mass | particle counts (below) |

### Geometry rules (PSSDplusplus ≥ 0.4.1, app)

- **Sphere:** width = height = length.
- **Fragment / Not Reported:** width = R.ave × length, height = H:W × width,
  where R.ave and H:W are the environment's width:length ratio and
  height:width ratio (Kooi et al. 2021).
- **Fiber:** one diameter, from the first of: reported minimum (or maximum)
  width, measured > nominal; reported single width; 15 µm (existing default).
  It is never wider than the fiber is long: the diameter is capped at the
  length (for polydisperse fibers, at the minimum length for the smallest
  particle and at the bioaccessible maximum length for the largest particle).
  Height = H:W × width.
- Polydisperse materials: the smallest and largest particles get these
  geometries at the minimum length and at the maximum bioaccessible length
  (the length truncated at the mouth size for food dilution, at the
  translocation limit for tissue translocation).

### Particle counts of tests dosed by mass

When the particle count was derived from the mass (`"converted"`), the
alignment recomputes it from the dosed mass and the mean particle mass of
its *own* geometry and size distribution, so that the particles it aligns
carry exactly the dosed mass:

`count = dose (mg/L) / mean particle mass (mg) / 1000` (particles/mL),
`mean particle mass = density × mean particle volume`; the mean volume is the
particle's volume (monodisperse) or the power-law mean over the whole
reported size range (polydisperse; or the measured distribution of
section 2). Outputs: `particle_count_basis` (`"as reported"`,
`"from mass (alignment geometry)"`; with
`align_data(particle_count_source = "database")`, used only for comparisons
with PSSDplusplus ≤ 0.4.0: `"from mass (database conversion)"`) and
`mean.mass.per.particle.mg`. Reported counts are always used as reported.

## 2. Measured size distribution of a polydisperse test material (optional)

By default a polydisperse test material is assumed to follow the
**environmental** power-law length distribution (exponent α of Kooi et al.
2021 for the environment) between its reported minimum and maximum lengths,
because no study in ToMEx reports its material's size distribution (SI of
Coffin et al. 2026). A record can instead carry its own measured
distribution:

| Field | Type, units | Meaning |
|---|---|---|
| `size.dist.family` | text: `"powerlaw"` or `"lognormal"` | family of the **number-based** length distribution. Empty = environmental default |
| `size.dist.alpha` | number (dimensionless) | power-law exponent of the number-based length distribution, p(L) ∝ L^−α (`"powerlaw"`) |
| `size.dist.median.um` | number, µm | median length (`"lognormal"`) |
| `size.dist.gsd` | number (> 1) | geometric standard deviation (`"lognormal"`) |
| `size.D10.um`, `size.D50.um`, `size.D90.um` | number, µm | number-based 10th/50th/90th percentile lengths; used for `"lognormal"` when the median or GSD is missing (D50 = median; GSD from D10 and D90, or from D50 and one of them) |

Rules for data entry:

- **Number-based only.** Percentiles and medians must refer to particle
  numbers. Laser-diffraction instruments usually report volume-based D10/D50/D90:
  do not enter those here (convert first or leave the fields empty).
- **Length.** The distribution is of particle length (longest dimension),
  like the other size fields.
- The reported minimum and maximum lengths still bound the material: the
  distribution is truncated to [min, max].
- Only polydisperse records use these fields. A family without usable
  parameters falls back to the default, with a warning; an unknown family is
  an error.

How it is used: the bioaccessible fraction (CF_bio of the material, the
number fraction of its particles no longer than the mouth size or
translocation limit) is taken from the distribution's CDF, and the mean
particle volume and surface area of the bioaccessible particles are
integrated numerically over the distribution with the shape geometry of
section 1. Outputs: `poly_size_method` (`"monodisperse"`,
`"assumed_env_powerlaw"`, `"measured_powerlaw"`, `"measured_lognormal"`) and,
in PSSDplusplus, `poly_volume_factor_ingest` (mean ingestible particle volume
÷ smallest particle volume; large values flag alignments that rest heavily on
the assumed distribution).

## 3. Co-dosed mixtures of particle types (optional)

A test that doses several particle types together (two shapes, two polymers,
two size ranges) is entered as **one row per component**:

| Field | Type | Meaning |
|---|---|---|
| `mixture_id` | text | identifier shared by all components of the mixture (one test treatment). Empty = an ordinary single-material record |
| `component.fraction` | number, 0–1 | the component's share of the mixture; normalised to sum to 1 within a mixture |
| `component.fraction.basis` | text: `"number"` (default) or `"mass"` | whether the share is of the particle number or of the dosed mass; the same for all components of a mixture |

Each component row carries that component's own shape, polymer, density and
sizes, and the **whole test's** dose (the same dose on every component row),
plus the usual record fields (species, endpoint, effect metric, ...).

How it is combined (`PSSDplusplus::align_mixture()`; app:
`combine_mixture_components()` at the end of each alignment):

1. Each component is aligned as if it were the whole test material.
2. Its aligned concentrations are scaled by its share of the dose. The share
   is its **number share** when the particle count was reported (the count is
   the whole test's count) and its **mass share** when the count was derived
   from the dosed mass (count = whole mass ÷ the component's mean particle
   mass). Number and mass shares are converted into each other with the
   components' mean particle masses.
3. The scaled concentrations of the **bioaccessible** components are added:
   under both ERMs the ecologically relevant dose (ingestible particle volume
   for food dilution, translocatable surface area for tissue translocation)
   is additive, and a component that is not ingestible (or not
   translocatable) adds nothing.

The combined record keeps the labels of its **dominant component** (largest
share), including `shape_f` and `poly_f`, so it remains selectable by shape
and polymer, and adds `n_components`, `mixture_shapes` and
`mixture_polymers` (the components' shapes and polymers, `"; "`-separated).
`ingestible`/`translocatable` become `"ingestible (mixture)"` /
`"translocatable (mixture)"` when any component is, and
`size.length.um.used.for.conversions` becomes the smallest component length.
The component rows (whole-test doses) and their shares
(`.component_share`) are attached as attribute `"components"`.

Example (a 70:30 mixture by number of PS spheres and PET fibers at
1,000 particles/mL):

| mixture_id | shape_f | poly_f | size.length.um.used.for.conversions | dose.particles.mL.master | component.fraction | component.fraction.basis |
|---|---|---|---|---|---|---|
| T12-high | Sphere | Polystyrene | 10 | 1000 | 0.7 | number |
| T12-high | Fiber | Polyethylene Terephthalate | 200 | 1000 | 0.3 | number |

## 4. Where this is implemented

| | PSSDplusplus (`ToMEx2.0_EcoToxRisk/package`) | ToMEx app (`ToMEx_AquaticOrganisms`) |
|---|---|---|
| geometry, mass-dosed counts | `R/functions.R::align_data()` (`particle_count_source`) | `app.R` alignment blocks |
| size distributions | `R/poly_size_dist.R` | `test_material_functions.R` (verbatim copy) |
| mixtures | `R/align_mixture.R` (`align_mixture()`, `combine_mixture_components()`) | `test_material_functions.R` (verbatim copy of `combine_mixture_components()`) |
| database fields | — | `ToMEx2.0_Onboarding/ToMEx2.0_Data_Tidying.R` (`add_test_material_fields()`: passes the fields through from the templates, else NA; validates values) |
| tests | `tests/testthat/test-polydisperse-geometry.R`, `test-mass-dosed-counts.R`, `test-align_mixture.R` | `tests/test_alignment_geometry.R`, `tests/test_alignment_erms.R` |

Template columns use the field names of sections 2 and 3 exactly and are
appended at the **end** of the data template (the tidying moves the screening
columns by position).

The analysis scripts of Coffin et al. (`ToMEx2.0_EcoToxRisk/scripts`) mirror
`align_data()` without these optional inputs: ToMEx records none of them yet.

## References

- Kooi, M., Primpke, S., Mintenig, S. M., Lorenz, C., Gerdts, G., Koelmans, A. A.
  (2021). Characterizing the multidimensionality of microplastics across
  environmental compartments. *Water Research* 202, 117429.
  doi:10.1016/j.watres.2021.117429
- Koelmans, A. A., et al. (2020). Solving the nonalignment of methods and
  approaches used in microplastic research to consistently characterize risk.
  *Environmental Science & Technology*. doi:10.1021/acs.est.0c02982
- Coffin, S., et al. (2026). A Probabilistic Risk Framework for Microplastics
  Integrating Uncertainty Across Toxicological and Environmental Variability:
  Development and Application to Marine and Freshwater Ecosystems. *Journal of
  Hazardous Materials*. Supplementary Information.
