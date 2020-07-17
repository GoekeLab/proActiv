
### PromoterAnnotation Class Definition ###

#' S4 class for promoter annotation data for a specific annotation version
#'
#' @slot intronRanges A GRanges object. The intron ranges annotated with the promoter
#'   information.
#' @slot promoterIdMapping A data.frame object. The id mapping between transcript ids, names,
#'   TSS ids, promoter ids and gene ids.
#' @slot promoterCoordinates A GRanges object. Promoter coordinates (TSS) with gene id and internal promoter state
#'
#' @name PromoterAnnotation-class
#' @rdname PromoterAnnotation-class
#' @exportClass PromoterAnnotation

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

###################
### Constructor ###

#' @param intronRanges A GRanges object containing annotated intron ranges
#' @param promoterIdMapping A data.frame containing mapping between transcript, tss, promoter and gene ids
#' @param promoterCoordinates A GRanges object containing promoter coordinates
#'
#' @name PromoterAnnotation
#' @rdname PromoterAnnotation-class
#'
#' @importFrom methods new
#'
#' @export
#' @return A promoter annotation object with three slots: intronRanges, promoterIdMapping 
#'   and promoter Coordinates
#'   
#' @examples 
#' 
#' promoterAnnotation <- PromoterAnnotation()
#' intronRanges(promoterAnnotation) <- intronRanges(promoterAnnotation.gencode.v19)
#' promoterIdMapping(promoterAnnotation) <- promoterIdMapping(promoterAnnotation.gencode.v19)
#' promoterCoordinates(promoterAnnotation) <- promoterCoordinates(promoterAnnotation.gencode.v19)
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

#' @param x A PromoterAnnotation object
#'
#' @describeIn PromoterAnnotation-class Getter for intronRanges
#' @exportMethod intronRanges

setGeneric("intronRanges", function(x) standardGeneric("intronRanges"))

#' @describeIn PromoterAnnotation-class Getter for intronRanges
#' @aliases intronRanges,PromoterAnnotation-method

setMethod("intronRanges", "PromoterAnnotation", function(x) x@intronRanges)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("promoterIdMapping", function(x) standardGeneric("promoterIdMapping"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("promoterIdMapping", "PromoterAnnotation", function(x) x@promoterIdMapping)

#' @describeIn PromoterAnnotation-class Getter for promoterCoordinates
#' @exportMethod promoterCoordinates

setGeneric("promoterCoordinates", function(x) standardGeneric("promoterCoordinates"))

#' @describeIn PromoterAnnotation-class Getter for promoterCoordinates
#' @aliases promoterCoordinates,PromoterAnnotation-method

setMethod("promoterCoordinates", "PromoterAnnotation", function(x) x@promoterCoordinates)

###############
### Setters ###

#' @param value intronRanges, promoterIdMapping or promoterCoordinates to be assigned
#' 
#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @exportMethod intronRanges<-
#' @importFrom methods validObject

setGeneric("intronRanges<-", function(x, value) standardGeneric("intronRanges<-"))

#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @aliases intronRanges<-,PromoterAnnotation-method

setMethod("intronRanges<-", "PromoterAnnotation", function(x, value) {
    x@intronRanges <- value
    validObject(x)
    x
})

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @exportMethod promoterIdMapping<-
#' @importFrom methods validObject

setGeneric("promoterIdMapping<-", function(x, value) standardGeneric("promoterIdMapping<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases promoterIdMapping<-,PromoterAnnotation-method

setMethod("promoterIdMapping<-", "PromoterAnnotation", function(x, value) {
    x@promoterIdMapping <- value
    validObject(x)
    x
})

#' @describeIn PromoterAnnotation-class Setter for promoterCoordinates
#' @exportMethod promoterCoordinates<-
#' @importFrom methods validObject

setGeneric("promoterCoordinates<-", function(x, value) standardGeneric("promoterCoordinates<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterCoordinates
#' @aliases promoterCoordinates<-,PromoterAnnotation-method

setMethod("promoterCoordinates<-", "PromoterAnnotation", function(x, value) {
    x@promoterCoordinates <- value
    validObject(x)
    x
})
