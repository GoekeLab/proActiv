
### PromoterAnnotation Class Definition ###

#' S4 class for promoter annotation data for a specific annotation version
#'
#' @slot reducedExonRanges A GRanges object. The reduced first exon ranges for
#'   each promoter without metadata.
#' @slot annotatedIntronRanges A GRanges object. The intron ranges annotated with the promoter
#'   information.
#' @slot promoterIdMapping A data.frame object. The id mapping between transcript ids, names,
#'   TSS ids, promoter ids and gene ids.
#' @slot promoterCoordinates A GRanges object. Promoter coordinates (TSS) with gene id and internal promoter state
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportClass PromoterAnnotation
#'
setClass(
  "PromoterAnnotation",
  slots = c(
    reducedExonRanges = "GRanges",
    annotatedIntronRanges = "GRanges",
    promoterIdMapping = "data.frame",
    promoterCoordinates = "GRanges"
  ),
  prototype = list(
    reducedExonRanges = GRanges(),
    annotatedIntronRanges = GRanges(),
    promoterIdMapping = data.frame(),
    promoterCoordinates = GRanges()
  )
)

#' Constructor for PromoterAnnotation class
#'
#' @param reducedExonRanges A GRanges object containing reduced exon ranges
#' @param annotatedIntronRanges A GRanges object containing annotated intron ranges
#' @param promoterIdMapping A data.frame containing mapping between transcript, tss, promoter and gene ids
#' @param promoterCoordinates A GRanges object containing promoter coordinates
#'
#' @name PromoterAnnotation
#' @rdname PromoterAnnotation-class
#'
#' @importFrom methods new
#'
#' @export
#'
PromoterAnnotation <-
  function(reducedExonRanges = GRanges(),
           annotatedIntronRanges = GRanges(),
           promoterIdMapping = data.frame(),
           promoterCoordinates = GRanges()) {
    new(
      "PromoterAnnotation",
      reducedExonRanges = reducedExonRanges,
      annotatedIntronRanges = annotatedIntronRanges,
      promoterIdMapping = promoterIdMapping,
      promoterCoordinates = promoterCoordinates
    )
  }

setValidity("PromoterAnnotation", function(object) {
  check <- TRUE
  if (is(object@reducedExonRanges, 'GRanges') == FALSE) {
    check <- FALSE
  }
  if (is(object@annotatedIntronRanges, 'GRanges') == FALSE) {
    check <- FALSE
  }
  if (is(object@promoterIdMapping, 'data.frame') == FALSE) {
    check <- FALSE
  }
  if (is(object@promoterCoordinates, 'GRanges') == FALSE) {
    check <- FALSE
  }
  return(check)
})

###############
### Getters ###

#' Getter for the reduced exon ranges
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod reducedExonRanges
#'
setGeneric("reducedExonRanges", function(x) standardGeneric("reducedExonRanges"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases reducedExonRanges,PromoterAnnotation-method
#'
setMethod("reducedExonRanges", "PromoterAnnotation", function(x) x@reducedExonRanges)

#' Getter for the annotated intron ranges
#'
#' @param x A PromoterAnnotation object
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod annotatedIntronRanges
#'
setGeneric("annotatedIntronRanges", function(x) standardGeneric("annotatedIntronRanges"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases annotatedIntronRanges,PromoterAnnotation-method
#'
setMethod("annotatedIntronRanges", "PromoterAnnotation", function(x) x@annotatedIntronRanges)

#' Getter for the promoter id mapping
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod promoterIdMapping
#'
setGeneric("promoterIdMapping", function(x) standardGeneric("promoterIdMapping"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases promoterIdMapping,PromoterAnnotation-method
#'
setMethod("promoterIdMapping", "PromoterAnnotation", function(x) x@promoterIdMapping)

#' Getter for the promoter coordinates
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod promoterCoordinates
#'
setGeneric("promoterCoordinates", function(x) standardGeneric("promoterCoordinates"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases promoterCoordinates,PromoterAnnotation-method
#'
setMethod("promoterCoordinates", "PromoterAnnotation", function(x) x@promoterCoordinates)

###############
### Setters ###

#' Setter for the reduced exon ranges
#'
#' @param value User specified value to set
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod reducedExonRanges<-
#'
#' @importFrom methods validObject
#'
setGeneric("reducedExonRanges<-", function(x, value) standardGeneric("reducedExonRanges<-"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases reducedExonRanges<-,PromoterAnnotation-method
#'
setMethod("reducedExonRanges<-", "PromoterAnnotation", function(x, value) {
  x@reducedExonRanges <- value
  validObject(x)
  x
})

#' Setter for the annotated intron ranges
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod annotatedIntronRanges<-
#'
#' @importFrom methods validObject
#'
setGeneric("annotatedIntronRanges<-", function(x, value) standardGeneric("annotatedIntronRanges<-"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases annotatedIntronRanges<-,PromoterAnnotation-method
#'
setMethod("annotatedIntronRanges<-", "PromoterAnnotation", function(x, value) {
  x@annotatedIntronRanges <- value
  validObject(x)
  x
})

#' Setter for the promoter id mapping
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod promoterIdMapping<-
#'
#' @importFrom methods validObject
#'
setGeneric("promoterIdMapping<-", function(x, value) standardGeneric("promoterIdMapping<-"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases promoterIdMapping<-,PromoterAnnotation-method
#'
setMethod("promoterIdMapping<-", "PromoterAnnotation", function(x, value) {
  x@promoterIdMapping <- value
  validObject(x)
  x
})

#' Setter for the promoter coordinates
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod promoterCoordinates<-
#'
#' @importFrom methods validObject
#'
setGeneric("promoterCoordinates<-", function(x, value) standardGeneric("promoterCoordinates<-"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases promoterCoordinates<-,PromoterAnnotation-method
#'
setMethod("promoterCoordinates<-", "PromoterAnnotation", function(x, value) {
  x@promoterCoordinates <- value
  validObject(x)
  x
})

