
### PromoterAnnotation Class Definition ###

#' S4 class for promoter annotation data for a specific annotation version
#'
#' @slot intronRanges A GRanges object. The intron ranges annotated with each intron's 
#'   associated transcript 
#' @slot promoterIdMapping A data.frame object. The id mapping between transcript ids, names,
#'   TSS ids, promoter ids and gene ids.
#' @slot promoterCoordinates A GRanges object. Promoter coordinates (TSS) with gene id,
#'   internal promoter state, coordinates of the end of the first reduced exon corresponding
#'   to each TSS, and associated intron ids
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportClass PromoterAnnotation
#'
setClass(
  "PromoterAnnotation",
  slots = c(
    intronRanges = "GRanges",
    promoterIdMapping = "data.frame",
    promoterCoordinates = "GRanges"
  ),
  prototype = list(
    intronRanges = GRanges(),
    promoterIdMapping = data.frame(),
    promoterCoordinates = GRanges()
  )
)

#' Constructor for PromoterAnnotation class
#'
#' @param intronRanges A GRanges object containing intron ranges
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
  function(intronRanges = GRanges(),
           promoterIdMapping = data.frame(),
           promoterCoordinates = GRanges()) {
    new(
      "PromoterAnnotation",
      intronRanges = intronRanges,
      promoterIdMapping = promoterIdMapping,
      promoterCoordinates = promoterCoordinates
    )
  }

setValidity("PromoterAnnotation", function(object) {
  check <- TRUE
  if (is(object@intronRanges, 'GRanges') == FALSE) {
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

#' Getter for the intron ranges
#'
#' @param x A PromoterAnnotation object
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod intronRanges
#'
setGeneric("intronRanges", function(x) standardGeneric("intronRanges"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases intronRanges,PromoterAnnotation-method
#'
setMethod("intronRanges", "PromoterAnnotation", function(x) x@intronRanges)

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

#' Setter for the intron ranges
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportMethod intronRanges<-
#'
#' @importFrom methods validObject
#'
setGeneric("intronRanges<-", function(x, value) standardGeneric("intronRanges<-"))

#'
#' @rdname PromoterAnnotation-class
#' @aliases intronRanges<-,PromoterAnnotation-method
#'
setMethod("intronRanges<-", "PromoterAnnotation", function(x, value) {
  x@intronRanges <- value
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

