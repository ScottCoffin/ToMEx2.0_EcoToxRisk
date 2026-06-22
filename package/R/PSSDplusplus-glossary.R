#' Glossary of Terms and Acronyms for PSSDplusplus
#'
#' Quick reference for the field-specific terminology and acronyms used
#' throughout the PSSDplusplus package, the bundled ToMEx 2.0 dataset, and
#' the associated manuscript (Coffin et al. 2026,
#' \doi{10.1016/j.jhazmat.2025.141021}).
#'
#' @section Statistical and modelling methods:
#' \describe{
#'   \item{LHS (Latin Hypercube Sampling)}{
#'     An efficient stratified-random parameter-space sampling method used to
#'     generate the alignment parameter matrix. Unlike simple Monte Carlo
#'     sampling, LHS ensures all regions of the parameter space are
#'     represented, reducing the number of simulations needed for stable
#'     estimates. Implemented via \code{\link[sensobol]{sobol_matrices}}.
#'   }
#'   \item{MC / Monte Carlo}{
#'     A simulation technique that propagates uncertainty by repeatedly
#'     sampling from parameter distributions. Used here to propagate
#'     uncertainty in ERM alignments across thousands of simulated datasets.
#'   }
#'   \item{Sobol' sequence}{
#'     A quasi-random low-discrepancy sequence used as the basis for LHS
#'     to fill the parameter space uniformly. Produces more even coverage
#'     than pseudo-random sequences.
#'   }
#'   \item{Bootstrapping}{
#'     A resampling method for estimating uncertainty (e.g., confidence
#'     intervals) around SSD-derived thresholds.
#'   }
#'   \item{MLE (Maximum Likelihood Estimation)}{
#'     A statistical method for fitting distribution parameters. Used to
#'     derive power-law exponents for environmental particle distributions.
#'   }
#'   \item{Power law exponent (alpha)}{
#'     Describes the frequency-size relationship of particles in a
#'     polydisperse mixture. Higher alpha = steeper decline in abundance
#'     with increasing particle size. Environment- and compartment-specific
#'     values are sampled from distributions reported in Kooi et al. (2021).
#'   }
#'   \item{RCI (Relative Confidence Interval)}{
#'     A standardised uncertainty measure: (95th percentile - 5th percentile)
#'     / median. Allows comparison of uncertainty widths across methods and
#'     environments.
#'   }
#' }
#'
#' @section Ecological risk assessment framework:
#' \describe{
#'   \item{ERA (Ecological Risk Assessment)}{
#'     A framework for evaluating the probability and magnitude of adverse
#'     ecological effects from a contaminant. The PSSD++ approach is one ERA
#'     methodology.
#'   }
#'   \item{SSD (Species Sensitivity Distribution)}{
#'     A statistical model representing variation in sensitivity to a
#'     stressor across species. The HC5 or HC10 of the SSD is used as the
#'     threshold (PNEC) below which 95\% or 90\% of species are protected.
#'   }
#'   \item{PSSD (Probabilistic SSD)}{
#'     A variant of the SSD that preserves species-level variability by
#'     building a probability distribution for each species, then combining
#'     them (Gottschalk & Nowack, 2013).
#'   }
#'   \item{PSSD+ }{
#'     Extension of PSSD that incorporates probabilistic assessment factors
#'     (uncertainty factors) to harmonise acute/sub-lethal endpoints to
#'     chronic NOECs (Wigger et al., 2020).
#'   }
#'   \item{PSSD++ }{
#'     The method implemented in this package. Extends PSSD+ by also
#'     propagating uncertainty in ERM-based particle-trait alignments via
#'     Monte Carlo simulation.
#'   }
#'   \item{HC5 / HC10 (Hazard Concentration for x percent of species)}{
#'     The concentration below which x percent of species are protected.
#'     HC5 is used for tiers 1-3; HC10 for tier 4 (reflecting higher
#'     confidence at organismal/population endpoints).
#'   }
#'   \item{PNEC (Predicted No Effect Concentration)}{
#'     The threshold concentration below which no adverse ecological effects
#'     are expected. Derived from the HC5/HC10 of the SSD or PSSD.
#'   }
#'   \item{Tier}{
#'     One of four data-filtering levels (1-4) that trade off protectiveness
#'     vs. predictive confidence. Tiers 1-2 include all biological endpoints
#'     and are more protective; tiers 3-4 restrict to organism- and
#'     population-level endpoints and are more predictive.
#'   }
#'   \item{AF (Assessment Factor)}{
#'     A multiplicative safety factor applied to toxicity values to
#'     extrapolate from tested endpoints (e.g., acute, EC50) to chronic
#'     no-effect concentrations (NOECs). Two types are used: \code{af.time}
#'     for acute-to-chronic conversion and \code{af.noec} for
#'     LOEC/ECx-to-NOEC conversion.
#'   }
#' }
#'
#' @section Effect metrics (dose descriptors):
#' \describe{
#'   \item{NOEC (No Observed Effect Concentration)}{
#'     The highest tested concentration at which no statistically significant
#'     adverse effect was observed. The target endpoint for all PSSD++ tiers.
#'   }
#'   \item{LOEC (Lowest Observed Effect Concentration)}{
#'     The lowest tested concentration at which a statistically significant
#'     adverse effect was observed.
#'   }
#'   \item{ECx (Effect Concentration at x percent)}{
#'     The concentration causing a specified percentage of the maximum
#'     observed effect. EC50 = 50 percent effect.
#'   }
#'   \item{LCx (Lethal Concentration at x percent)}{
#'     The concentration causing x percent mortality.
#'   }
#'   \item{HONEC (Highest Observed No Effect Concentration)}{
#'     Excluded from analyses because it is not a rigorously determined
#'     no-effect level (no statistical test establishing the boundary).
#'   }
#' }
#'
#' @section Microplastic-specific terms:
#' \describe{
#'   \item{ERM (Ecologically Relevant Metric)}{
#'     A particle- and species-specific measure of MP exposure effects.
#'     Unlike conventional chemical concentrations (mg/L), ERMs account for
#'     physical MP traits relevant to the specific toxicity mechanism. Two
#'     ERMs are used: \strong{food dilution} and
#'     \strong{tissue translocation-mediated effects}. Defined in
#'     Koelmans et al. (2020).
#'   }
#'   \item{Food dilution ERM}{
#'     MPs are modelled as caloric dilutants in the organism's diet. The
#'     alignment converts lab doses to environmentally relevant particle
#'     volumes per litre (particles/L) based on ingestible particle sizes
#'     (determined by species body length via an allometric model).
#'   }
#'   \item{Tissue translocation ERM}{
#'     MPs are modelled as particles that cross the gut epithelium into
#'     tissues. The alignment uses surface area and a logistic regression
#'     model predicting the particle length at which there is 50\%
#'     probability of translocation.
#'   }
#'   \item{Bioaccessibility}{
#'     The fraction of environmental MPs that can be absorbed by an organism
#'     based on particle size, shape, and species physiology. Two
#'     bioaccessibility limits are used: (1) ingestion (body-length-based
#'     allometric model) and (2) tissue translocation (logistic regression).
#'   }
#'   \item{Translocation}{
#'     Movement of particles from the digestive tract into tissues or the
#'     circulatory system. Particle length is the primary predictor: smaller
#'     particles have a higher probability of translocation.
#'   }
#'   \item{Eco-corona}{
#'     A layer of biomolecules (proteins, lipids, polysaccharides) that
#'     adsorbs onto MP surfaces in environmental waters. Alters particle
#'     behavior, bioavailability, and toxicity.
#'   }
#'   \item{ToMEx 2.0}{
#'     Toxicity of Microplastics Explorer version 2.0. The largest
#'     published MP ecotoxicity database (~13,000 data points, ~300
#'     studies). Bundled in this package as \code{\link{tomex2}}.
#'   }
#' }
#'
#' @section Polymer abbreviations:
#' \describe{
#'   \item{PA}{Polyamide (e.g., nylon)}
#'   \item{PE}{Polyethylene}
#'   \item{PET}{Polyethylene Terephthalate}
#'   \item{PP}{Polypropylene}
#'   \item{PS}{Polystyrene}
#'   \item{PTFE}{Polytetrafluoroethylene (Teflon)}
#'   \item{PUR}{Polyurethane}
#'   \item{PVC}{Polyvinyl Chloride}
#'   \item{NIAS}{Non-Intentionally Added Substances: impurities or
#'     by-products present in plastic formulations.}
#' }
#'
#' @references
#' Coffin, S. et al. (2026). A probabilistic risk framework for microplastics
#' integrating uncertainty across toxicological and environmental variability.
#' \emph{Journal of Hazardous Materials} 503, 141021.
#' \doi{10.1016/j.jhazmat.2025.141021}
#'
#' Gottschalk, F. & Nowack, B. (2013). A probabilistic method for species
#' sensitivity distributions. \emph{Integr. Environ. Assess. Manag.} 9, 79-86.
#'
#' Koelmans, A.A. et al. (2020). Solving the nonalignment of methods and
#' approaches used in microplastic research. \emph{Environ. Sci. Technol.}
#' 54, 12307-12315.
#'
#' Kooi, M. et al. (2021). Characterizing the multidimensionality of
#' microplastics across environmental compartments. \emph{Water Res.} 202,
#' 117429.
#'
#' Wigger, H. et al. (2020). Systematic consideration of parameter
#' uncertainty and variability in probabilistic species sensitivity
#' distributions. \emph{Integr. Environ. Assess. Manag.} 16, 211-222.
#'
#' @name PSSDplusplus-glossary
#' @aliases PSSDplusplus_glossary glossary
NULL
