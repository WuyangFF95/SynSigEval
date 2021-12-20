#' Append most similar ground-truth sig and pairwise cosine similarity
#' to the name of each extracted sig in matrix of extracted sigs.
#'
#' @param sigAnalysis A list returned by function \code{\link{ReadAndAnalyzeSigs}},
#' at least including: \enumerate{
#'   \item \code{ex.sigs}: matrix of extracted sigs
#'   \item \code{sim.matrix}: full matrix between extracted sigs and
#'     ground-truth sigs
#'   \item \code{table}: matrix include pairs of true positive ground-truth
#'   sigs and true positive extracted sigs.
#'  }
#'
#' @return Matrix of extracted sigs, yet the names of signatures changed.
#'   New name: <old_name> (<name_of_most_similar_ground_truth_sig> cosine_similarity)
#'   e.g., Sig.A --> Sig.A (SBS1 0.998)
RelabelExSigs <- function(sigAnalysis) {

  # Old names of extracted sigs
  old_names <- colnames(sigAnalysis$ex.sigs)

  # Use old name to label new names
  new_names <- character(nrow(sigAnalysis$sim.matrix))
  names(new_names) <- old_names

  # Generate new names
  for (old_name in old_names) {
    # cosine sim > 0.9 to at least one gt sig
    if (old_name %in% sigAnalysis$table[, 1]) {
      index <- which(sigAnalysis$table[, 1] == old_name)
      matched_gt_sig <- sigAnalysis$table[index, 2]
      cossim <- sigAnalysis$table[index, 3]
      # Create new name
      new_names[old_name] <- paste0(
        old_name, " (", matched_gt_sig, " ",
        round(cossim, 3), ")"
      )
    } else {
      # For false positive sig,
      # find its most similar ground-truth sig.
      index <- match(old_name, old_names)
      matched_gt_sig_index <-
        which.max(sigAnalysis$sim.matrix[index, ])
      matched_gt_sig <-
        colnames(sigAnalysis$sim.matrix)[matched_gt_sig_index]
      cossim <- max(sigAnalysis$sim.matrix[index, ])
      # Create new name
      new_names[old_name] <- paste0(
        old_name, " (False Positive) (", matched_gt_sig, " ",
        round(cossim, 3), ")"
      )
    }
  }

  # Return ex.sigs matrix, with new signature names
  ex.sigs.renamed <- sigAnalysis$ex.sigs
  colnames(ex.sigs.renamed) <- new_names
  return(ex.sigs.renamed)
}

