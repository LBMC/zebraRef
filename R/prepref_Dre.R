#' Build Interpolated Gene Expression References
#' 
#' Builds the interpolation of the reference datasets.
#' **These functions are internally called by \code{\link[RAPToR]{prepare_refdata}} from RAPToR.**
#' 
#' @param n.inter the resolution of the interpolation, as in \code{seq(start, end, length.out = n.inter)}.
#' 
#' @return A list with \code{interpGE} the interpolated gene expression matrix and 
#' \code{time.series} the time of the interpGE matrix columns.
#' 
#' @seealso \code{\link[RAPToR]{prepare_refdata}} \code{\link[RAPToR]{ge_im}}
#' 
#' @name Dre_prep
NULL

#' @rdname Dre_prep
#' @export
#' @importFrom RAPToR ge_im
#' @importFrom stats predict
#' 
.prepref_Dre_emb_larv <- function(n.inter){
  # utils::data("Dre_emb_larv", envir = environment())
  m <- RAPToR::ge_im(
    X = zebraRef::Dre_emb_larv$g,
    p = zebraRef::Dre_emb_larv$p,
    formula = zebraRef::Dre_emb_larv$geim_params$formula,
    method = zebraRef::Dre_emb_larv$geim_params$method,
    dim_red = zebraRef::Dre_emb_larv$geim_params$dim_red,
    nc = zebraRef::Dre_emb_larv$geim_params$nc
  )
  ndat <- data.frame(age = seq(min(zebraRef::Dre_emb_larv$p$age),
                               max(zebraRef::Dre_emb_larv$p$age),
                               l = n.inter))
  return(
    list(interpGE = predict(m, ndat), time.series = ndat$age)
  )
}