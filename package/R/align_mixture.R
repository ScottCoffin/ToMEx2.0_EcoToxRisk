#' Align tests that dose a mixture of particle types
#'
#' Aligns a test whose material is a mixture of components (for example two
#' shapes or two polymers dosed together) component by component, then adds
#' the components' aligned effect concentrations. Under both ERMs the
#' ecologically relevant dose is additive across components: ingestible
#' particle volume for food dilution and translocatable particle surface area
#' for tissue translocation. A component that is not bioaccessible contributes
#' nothing.
#'
#' [align_data()] treats every record as a single material (one shape, one
#' polymer and density, one size range). ToMEx does not yet record
#' component-resolved test materials, so this function is for data extracted
#' that way (one row per component); it changes nothing for existing records.
#' The data fields are described in `docs/test_material_fields.md` of the
#' repository.
#'
#' Each component is aligned with the **whole test's** dose, then its aligned
#' concentrations are scaled by its share of that dose (aligned concentrations
#' are proportional to the dose). The share depends on how the component's
#' particle count was obtained:
#' * a count reported by the study (`particle_count_basis` `"as reported"`) is
#'   the whole test's count, so the share is the component's **number** share;
#' * a count derived from the dosed mass (`"from mass ..."`, see
#'   `particle_count_source` in [align_data()]) is the whole test's mass divided
#'   by the component's mean particle mass, so the share is the component's
#'   **mass** share.
#' Number and mass shares are converted into each other with the components'
#' mean particle masses (`mean.mass.per.particle.mg`, over each component's
#' whole reported size range).
#'
#' @param components A data frame with one row per mixture component. Each row
#'   carries the record's usual columns (as for [align_data()]), with the
#'   component's own `shape_f`, sizes, `density.g.cm3` and `poly_f`, and the
#'   **whole mixture's** dose (`dose.particles.mL.master` and/or
#'   `dose.mg.L.master`; the sediment equivalents for sediment tests),
#'   repeated on every component row.
#' @param mixture_id Name of the column identifying the mixture (test record)
#'   each component belongs to.
#' @param fraction Name of the column giving each component's share of the
#'   mixture. Shares are normalised to sum to 1 within a mixture (with a
#'   warning if they did not).
#' @param fraction_basis `"number"` (share of particles) or `"mass"` (share of
#'   the dosed mass). Used for rows without a `component.fraction.basis`
#'   column value.
#' @param ... Further arguments passed to [align_data()].
#'
#' @return One row per mixture: the record columns of its dominant component
#'   (largest share; its `shape_f`, `poly_f` and other labels, so the record
#'   stays selectable by shape and polymer), with `n_components`,
#'   `mixture_shapes` and `mixture_polymers` (the components' shapes and
#'   polymers, `"; "`-separated), the summed ERM-aligned concentrations
#'   (`particles.mL.food.dilution`, `particles.mL.ox.stress`,
#'   `particles.kg.food.dilution`, `particles.kg.ox.stress` and the matching
#'   `EC_env_*` columns), the summed component particle doses,
#'   `ingestible`/`translocatable` (`"ingestible (mixture)"` if any component
#'   is, else `"not ingestible"`), and `size.length.um.used.for.conversions`
#'   set to the smallest component length, so the record passes the SSD size
#'   filters when any component does. The aligned component rows (whole-test
#'   doses, before scaling) are attached as attribute `"components"`, with
#'   their shares in `.component_share`.
#'
#' @seealso [align_data()]
#' @examples
#' \dontrun{
#' comp <- PSSDplusplus::tomex2_mini[c(1, 1), ]
#' comp$mixture_id <- "test-1"
#' comp$component.fraction <- c(0.7, 0.3)
#' comp$shape_f <- c("Sphere", "Fiber")
#' align_mixture(comp)
#' }
#' @export
align_mixture <- function(components,
                          mixture_id = "mixture_id",
                          fraction = "component.fraction",
                          fraction_basis = c("number", "mass"),
                          ...) {
  fraction_basis <- match.arg(fraction_basis)
  for (col in c(mixture_id, fraction)) {
    if (!col %in% names(components)) {
      stop("`components` has no column ", dQuote(col, FALSE), call. = FALSE)
    }
  }
  comp <- as.data.frame(dplyr::ungroup(components), stringsAsFactors = FALSE)
  comp$shape_f <- as.character(comp$shape_f)
  if ("poly_f" %in% names(comp)) comp$poly_f <- as.character(comp$poly_f)
  id <- as.character(comp[[mixture_id]])
  f <- suppressWarnings(as.numeric(comp[[fraction]]))
  if (anyNA(id) || any(id == "") || anyNA(f) || any(f < 0)) {
    stop("mixture ids and component fractions must be present and fractions >= 0",
      call. = FALSE
    )
  }
  total <- stats::ave(f, id, FUN = sum)
  if (any(abs(total - 1) > 1e-6)) {
    warning("component fractions do not sum to 1 within every mixture; normalising",
      call. = FALSE
    )
  }
  if (!"component.fraction.basis" %in% names(comp)) {
    comp$component.fraction.basis <- fraction_basis
  }
  basis <- as.character(comp$component.fraction.basis)
  basis[is.na(basis) | basis == ""] <- fraction_basis
  comp$component.fraction.basis <- basis
  comp$mixture_id <- id
  comp$component.fraction <- f

  aligned <- suppressWarnings(align_data(comp, verbose = FALSE, ...))
  combine_mixture_components(aligned,
    ec_pattern = "^EC_env_|^particles[.](mL|kg)[.](food[.]dilution|ox[.]stress)$"
  )
}

#' Combine the aligned components of co-dosed mixtures into one record each
#'
#' Shared with the ToMEx Shiny app (SCCWRP/ToMEx_AquaticOrganisms,
#' `test_material_functions.R`, verbatim copy): keep the two identical.
#' Input: aligned rows, one per component, each aligned with the WHOLE test's
#' dose, carrying `mixture_id`, `component.fraction` and (optionally)
#' `component.fraction.basis` ("number" or "mass"; default "number"). Rows
#' without a `mixture_id` are returned unchanged.
#'
#' A component's aligned concentrations are scaled by its share of the whole
#' test's dose and added over the bioaccessible components. The share is the
#' number share when the component's particle count was reported (the whole
#' test's count), and the mass share when it was derived from the dosed mass
#' (`particle_count_basis` "from mass ...": whole mass / the component's mean
#' particle mass). Number and mass shares are converted with
#' `mean.mass.per.particle.mg`. Bioaccessibility of a column: names containing
#' "ingest" or "food" use `ingestible`, "trans" or "ox.stress" use
#' `translocatable`; other columns follow `switch_col` (the app's
#' ingestion/translocation switch).
#' @noRd
combine_mixture_components <- function(df, ec_pattern = "^EC_env_",
                                       switch_col = "ingestion.translocation.switch") {
  if (!"mixture_id" %in% names(df) || all(is.na(df$mixture_id) | df$mixture_id == "")) {
    return(df)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  is_mix <- !is.na(df$mixture_id) & df$mixture_id != ""
  singles <- df[!is_mix, , drop = FALSE]
  mix <- df[is_mix, , drop = FALSE]
  f <- suppressWarnings(as.numeric(mix$component.fraction))
  basis <- if ("component.fraction.basis" %in% names(mix)) {
    tolower(trimws(as.character(mix$component.fraction.basis)))
  } else {
    rep("number", nrow(mix))
  }
  basis[is.na(basis) | basis == ""] <- "number"
  if (any(is.na(f)) || any(f < 0) || any(!basis %in% c("number", "mass"))) {
    stop("mixture components need component.fraction >= 0 and ",
      "component.fraction.basis 'number' or 'mass'",
      call. = FALSE
    )
  }
  id <- as.character(mix$mixture_id)
  if (any(stats::ave(seq_along(basis), id, FUN = function(i) length(unique(basis[i]))) > 1)) {
    stop("all components of a mixture need the same component.fraction.basis", call. = FALSE)
  }
  f <- f / stats::ave(f, id, FUN = sum)
  mean_mass <- if ("mean.mass.per.particle.mg" %in% names(mix)) {
    as.numeric(mix$mean.mass.per.particle.mg)
  } else {
    rep(NA_real_, nrow(mix))
  }
  from_mass <- if ("particle_count_basis" %in% names(mix)) {
    grepl("^from mass", as.character(mix$particle_count_basis))
  } else {
    rep(FALSE, nrow(mix))
  }
  # a share in the other basis is needed where the count basis and the
  # fraction basis differ
  need_conv <- (basis == "number" & from_mass) | (basis == "mass" & !from_mass)
  conv_ids <- unique(id[need_conv])
  if (any(id %in% conv_ids & !(is.finite(mean_mass) & mean_mass > 0))) {
    stop("converting between number and mass shares needs every component's ",
      "mean.mass.per.particle.mg (density and sizes)",
      call. = FALSE
    )
  }
  n_share <- ifelse(basis == "number", f, NA_real_)
  m_share <- ifelse(basis == "mass", f, NA_real_)
  if (length(conv_ids)) {
    k <- id %in% conv_ids
    nm <- ifelse(basis == "number", f * mean_mass, NA_real_) # number -> mass
    mn <- ifelse(basis == "mass", f / mean_mass, NA_real_) # mass -> number
    m_share[k & basis == "number"] <- (nm / stats::ave(ifelse(is.na(nm), 0, nm), id, FUN = sum))[k & basis == "number"]
    n_share[k & basis == "mass"] <- (mn / stats::ave(ifelse(is.na(mn), 0, mn), id, FUN = sum))[k & basis == "mass"]
  }
  scale <- ifelse(from_mass, m_share, n_share)
  sw <- if (switch_col %in% names(mix)) as.character(mix[[switch_col]]) else rep("ingestion", nrow(mix))
  ingest_ok <- !as.character(mix$ingestible) %in% c("not ingestible", "unknown")
  trans_ok <- !as.character(mix$translocatable) %in% "not translocatable"
  switch_ok <- dplyr::case_when(
    sw == "translocation" ~ trans_ok,
    sw == "ingestion" ~ ingest_ok,
    TRUE ~ TRUE
  )
  ec_cols <- grep(ec_pattern, names(mix), value = TRUE)
  dose_cols <- intersect(c(
    "dose.particles.mL.master", "dose.particles.kg.master",
    "dose.particles.kg.sediment.master"
  ), names(mix))
  ids <- unique(id)
  out <- do.call(rbind, lapply(ids, function(k) {
    r <- which(id == k)
    lead <- r[which.max(f[r])] # dominant component carries the record's labels
    row <- mix[lead, , drop = FALSE]
    for (col in ec_cols) {
      ok <- if (grepl("ingest|food", col)) {
        ingest_ok
      } else if (grepl("trans|ox[.]stress", col)) {
        trans_ok
      } else {
        switch_ok
      }
      row[[col]] <- sum(ifelse(ok[r], scale[r] * mix[[col]][r], 0))
    }
    for (col in dose_cols) row[[col]] <- sum(scale[r] * mix[[col]][r])
    row$n_components <- length(r)
    row$mixture_shapes <- paste(sort(unique(as.character(mix$shape_f[r]))), collapse = "; ")
    if ("poly_f" %in% names(mix)) {
      row$mixture_polymers <- paste(sort(unique(as.character(mix$poly_f[r]))), collapse = "; ")
    }
    row$ingestible <- if (any(ingest_ok[r])) "ingestible (mixture)" else "not ingestible"
    row$translocatable <- if (any(trans_ok[r])) "translocatable (mixture)" else "not translocatable"
    row$size.length.um.used.for.conversions <- min(mix$size.length.um.used.for.conversions[r], na.rm = TRUE)
    row
  }))
  mix$.component_share <- scale
  if (nrow(singles)) {
    for (col in setdiff(names(out), names(singles))) singles[[col]] <- NA
    out <- rbind(singles, out[, names(singles), drop = FALSE])
  }
  attr(out, "components") <- mix
  out
}
