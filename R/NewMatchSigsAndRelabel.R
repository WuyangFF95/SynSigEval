

#' A wrapper function around \code{\link[ICAMSxtra]{match_two_sig_sets()}}
#' and substitutes deprecate function \code{\link[ICAMSxtra]{MatchSigsAndRelabel()}},
#'
#'
#' @param ex.sigs Newly extracted signatures to be compared to \code{gt.sigs}.
#'
#' @param gt.sigs "Ground truth" signatures.
#'
#' @param exposure If \code{NULL}, then match
#'   \code{ex.sigs} against all signatures in \code{gt.sigs}.
#'   Otherwise this should be ground-truth exposures used generate the
#'   synthetic spectra from which \code{ex.sigs} were extracted.
#'   In this case we do not
#'   match to ground-truth signatures that were not in the ground
#'   truth exposure.
#'
#' @param similarity.cutoff A ground-truth signature must have been
#'   the best match of an extracted signature with a cosine
#'   similarity \eqn{$ge$} this value to be considered a true positive.
#'   Otherwise we consider the ground-truth signature to be a false
#'   negative.
#'
#' @return A list with the elements
#'
#' * \code{averCosSim} The average of the cosine similarities
#'    between each signature in
#'    \code{ex.sigs} and its closest match in \code{gt.sigs}
#'    and the closest match
#'    between each signature in \code{gt.sigs}
#'    and its closest match in \code{ex.sigs}.
#'    This may not be what you want. Often one wants
#'    the average of the cosine similarities of the true
#'    positives to their matching ground-truth signatures.
#'
#' * \code{match1} The match from extracted signatures to ground truth
#'         signatures. Associated with each extracted signature is
#'         a ground truth signature with best cosine similarity.
#'         TODO: Ties?
#'
#' * \code{match2} The match from ground truth signatures to extracted
#'         signatures. Associated with each extracted signature is
#'         a ground truth signature with best cosine similarity.
#'         TODO: Ties?
#'
#' * \code{extracted.with.no.best.match}
#' * \code{ground.truth.with.no.best.match}
#' * \code{ex.sigs} Echo input argument
#' * \code{gt.sigs} Echo input argument
#' * \code{gt.mean.cos.sim}
#'
#' @export
NewMatchSigsAndRelabel <- function(ex.sigs,
                                   gt.sigs,
                                   exposure          = NULL,
                                   similarity.cutoff = 0.9) {


  if (is.null(colnames(ex.sigs))) {
    colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
  }

  if (!is.null(exposure)) {
    # Remove signatures that are not present in
    # the exposure from which the synthetic data were
    # generated
    exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
    # Make sure we do not have any signatures in exposures that
    # are not in gt.sigs.
    stopifnot(
      setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
    gt.sigs <- gt.sigs[  , exposed.sig.names]
  }

  gt.sigs.all.sim <-
    suppressMessages(philentropy::distance(t(gt.sigs), method = "cosine"))
  if (is.null(dim(gt.sigs.all.sim))) {
    if (gt.sigs.all.sim >= similarity.cutoff) {
      warning("The two ground truth signatures have cosine similarity >= ",
              similarity.cutoff)
    }
  } else {
    # browser()
    gt.sigs.all.sim[lower.tri(gt.sigs.all.sim, diag = TRUE)] <- 0
    if (any(gt.sigs.all.sim >= similarity.cutoff)) {
      rownames(gt.sigs.all.sim) <- colnames(gt.sigs)
      colnames(gt.sigs.all.sim) <- colnames(gt.sigs)
      warning(
        "Some ground truth signatures have cosine similarity >= ",
        similarity.cutoff)
      print(gt.sigs.all.sim)
    }
  }

  sim <- ICAMSxtra::match_two_sig_sets(gt.sigs, ex.sigs,
                                       cutoff = 0.9)

  # According to Hungarian algorithm,
  # the match in table minimizes the sum of cosine distance
  # and thus maximizes the sum of cosine similarity.
  sim$table[, 3] <- 1 - sim$table[, 3]
  colnames(sim$table) <- c("gt", "ex", "cossim")

  # List true positive signatures
  sim$ground.truth.with.best.match <- unique(sim$table$gt)
  # List false negative signatures
  # Ground-truth sigs without a best match from any extracted sig
  sim$ground.truth.with.no.best.match <-
    setdiff(colnames(gt.sigs), sim$table$gt)
  # List false positive signatures.
  # Extracted sigs without a best match from any ground-truth sig
  sim$extracted.with.no.best.match <-
    setdiff(colnames(ex.sigs), sim$table$ex)

  # Calculate average cosine similarity.
  # Only best matches listed in sim$table are considered.
  sim$averCosSim <- mean(sim$table$cossim)

  # Added reduced ground-truth signatures
  # to list sim
  sim$gt.sigs <- gt.sigs

  # Added extracted signatures
  # to list sim
  # Renaming is temporarily disabled.
  sim$ex.sigs <- ex.sigs

  invisible(sim)
}
