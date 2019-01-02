newCCSParams <- function() {
    new("CCSParams")
}

setMethod("getParam", "CCSParams", function(object, name){
    slot(object, name)
})
setMethod("setParam<-", "CCSParams", function(object, name, value){
    slot(object, name) <- value
    validObject(object)
    return(object)
})
setMethod("show", "CCSParams", function(object){
    cat("A", crayon::bold("CCSParams"), "object of class",
        crayon::bold(class(object)), "\n")
})
