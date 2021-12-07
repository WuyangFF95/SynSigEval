

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
#'
#' @examples
#' gt.sigs <- matrix(c(1, 3,   4, 1, 2, 4), ncol = 2)
#' ex.sigs <- matrix(c(1, 3.1, 4, 5, 1, 1, 1, 2.8, 4), ncol = 3)
#' colnames(gt.sigs) <- c("gt.1", "gt.2")
#' colnames(ex.sigs) <- c("ex.1", "ex.2", "ex.3")
#' tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
#' tout
#'

NewMatchSigsAndRelabel <- function(ex.sigs,
                                   gt.sigs,
                                   exposure          = NULL,
                                   similarity.cutoff = 0.9) {

  # Sanity check --------------------------------------------------------------

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



  # Fetch true positive, false positive and false negative signatures ---------


  # According to Hungarian algorithm,
  # the match in table minimizes the sum of cosine distance
  # and thus maximizes the sum of cosine similarity.
  #
  # sim$table stores matches and cosine similarities
  # between two signature sets: gt and ex.
  #
  # This match suffices:
  # 1. The match is unique bi-directionally. i.e. Only one gt sig matches
  #   to one ex sig, and vice versa.
  #
  # 2. The cosine similarity of each match is higher than the threshold,
  #   by default, 0.9.
  #
  # 3. The sum of pairwise cosine similarity is maximized.
  sim <- ICAMSxtra::TP_FP_FN_avg_sim(gt.sigs, ex.sigs, similarity.cutoff)


  # List true positive signatures
  # as ground-truth signature with a best match.
  TP_gt_names <- sim$table$gt
  # List extracted signatures
  # with a best match to a ground-truth signature.
  extracted.with.best.match <- sim$table$ex
  sim$extracted.with.best.match <- extracted.with.best.match
  # List false negative signatures
  # Ground-truth sigs without a best match from any extracted sig
  sim$ground.truth.with.no.best.match <-
    setdiff(colnames(gt.sigs), sim$table$gt)
  # List false positive signatures.
  # Extracted sigs without a best match from any ground-truth sig
  extracted.with.no.best.match <-
    setdiff(colnames(ex.sigs), sim$table$ex)
  sim$extracted.with.no.best.match <- extracted.with.no.best.match


  # Calculate average cosine similarity.
  # Only best matches listed in sim$table are considered.
  sim$averCosSim <- mean(sim$table$cosSim)

  # Output ground-truth signatures with non-zero exposures
  sim$gt.sigs <- gt.sigs



  # Re-order and re-name extracted signatures with best match -----------------

  # Extracted signatures with best match are ranked first,
  # and their internal ranking is determined by the ID of
  # ground-truth signature matching to.
  new_order_TP_ex <-
    gtools::mixedorder(ground.truth.with.best.match)

  # Names of labels are old extracted signature names
  # with new order.
  #
  # Values of labels_TP_ex are new extracted signature
  # names with new order.
  labels_TP_ex <- character(length(extracted.with.best.match))
  # Names of labels are old extracted signature names
  # with new order.
  names(labels_TP_ex) <- extracted.with.best.match[new_order_TP_ex]
  for (ii in 1:length(labels_TP_ex)) {
    # Let the labels follow the new-order.
    jj <- new_order_TP_ex[ii]
    labels_TP_ex[ii] <- paste0(extracted.with.best.match[jj],
                           " (",
                           ground.truth.with.best.match[jj],
                           " ",
                           sim$table[jj, 3],
                           ")")
  }



  # Re-order and re-name extracted signatures without best match --------------


  # Starting from false positive extracted signatures
  # without a best bi-lateral match listed in sim$table,
  #
  # We first match unilaterally from these false positive signatures
  # to their most similar ground-truth signature.
  #
  # We next re-sort these false positive signatures according to the
  # ID of these ground-truth signature in unilateral match.
  #
  # E.g.,
  # sigs c(ex.1, ex.2, ex.3)
  # are most similar to
  # sigs c(ID10, ID2, ID5),
  # then these false positive sigs are re-sorted as
  # c(ex.2, ex3, ex1).
  gt_sigs_FP <- character(0)
  cosSim_FP <- numeric(0)
  names(gt_sigs_FP) <- extracted.with.no.best.match
  names(cosSim_FP) <- extracted.with.no.best.match

  for(ex_sig_name in extracted.with.no.best.match) {
    gt_sigs_FP[ex_sig_name] <-
      rownames(sim$sim.matrix)[which.max(sim$sim.matrix[, ex_sig_name])]
    cosSim_FP[ex_sig_name] <- max(sim$sim.matrix[, ex_sig_name])
  }

  new_order_FP <- gtools::mixedorder(extracted.with.no.best.match)

  # Names of labels_FP are old extracted signature names
  # with new order.
  #
  # Values of labels_FP are new extracted signature
  # names with new order.
  labels_FP <- character(length(extracted.with.no.best.match))
  # Names of labels are old extracted signature names
  # with new order.
  names(labels_FP) <- extracted.with.no.best.match[new_order_FP]
  for (ii in 1:length(labels_FP)) {
    # Let the labels follow the new-order.
    jj <- new_order_FP[ii]
    labels_FP[ii] <- paste0(extracted.with.no.best.match[jj],
                           "(False-Positive) (",
                           gt_sigs_FP[jj],
                           " ",
                           cosSim_FP[jj],
                           ")")
  }



  # Rename and sort extracted signatures altogether ---------------------------
  #
  # Extracted signatures with best match are ranked first
  ex.sigs.sorted <- ex.sigs[, c(names(labels_TP_ex), names(labels_FP))]
  ex.sigs.sorted.renamed <- ex.sigs.sorted
  colnames(ex.sigs.sorted.renamed) <- c(labels_TP_ex, labels_FP)

  sim$ex.sigs <- ex.sigs.sorted.renamed



  # Calculate BEST cosine similarity for each gt signature --------------------

  # gt sigs in sim$table are gt sigs with best match.
  # For these sigs, we directly record the cosSim values in sim$table.
  cosSim_TP_gt <- sim$table$cosSim
  names(cosSim_TP_gt) <- sim$table$gt

  # For gt sigs which do not have a best match (and thus FN sigs),
  # their best cosine similarity values are calculated in cosSim_FN
  cosSim_FN <- numeric(length(ground.truth.with.no.best.match))
  names(cosSim_FN) <- ground.truth.with.no.best.match
  for(gt_sig_name in ground.truth.with.no.best.match) {
    cosSim_FN[gt_sig_name] <- max(sim$sim.matrix[gt_sig_name,])
  }

  sim$cosSim <- c(cosSim_TP_gt, cosSim_FN)




  invisible(sim)
}
