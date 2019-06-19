context("CCSParamsFunctions")

# newCCSParams
params <- newCCSParams()
expect_true(class(params) =="CCSParams")

# getParam
expect_true(is.numeric(getParam(params, "nGene")))
expect_true(is.numeric(getParam(params, "nCell")))
expect_true(is.list(getParam(params, "cciInfo")))
expect_true(is.numeric(getParam(params, "lambda")))
expect_true(is.numeric(getParam(params, "seed")))

# setParam
setParam(params, "nGene") <- 20000
expect_true(identical(getParam(params, "nGene"), 20000))

setParam(params, "nCell") <- c(12, 43, 323)
expect_true(identical(getParam(params, "nCell"), c(12, 43, 323)))

setParam(params, "cciInfo") <- list(nPair=2000,
                                   CCI1=list(
                                       LPattern=c(1,0,0),
                                       RPattern=c(0,1,1),
                                       nGene=200,
                                       fc="E10"),
                                   CCI2=list(
                                       LPattern=c(0,0,1),
                                       RPattern=c(1,1,1),
                                       nGene=200,
                                       fc="E10"),
                                   CCI3=list(
                                       LPattern=c(1,1,1),
                                       RPattern=c(1,0,1),
                                       nGene=200,
                                       fc="E10")
                                   )
expect_true(identical(getParam(params, "cciInfo"),
    list(nPair=2000,
       CCI1=list(
           LPattern=c(1,0,0),
           RPattern=c(0,1,1),
           nGene=200,
           fc="E10"),
       CCI2=list(
           LPattern=c(0,0,1),
           RPattern=c(1,1,1),
           nGene=200,
           fc="E10"),
       CCI3=list(
           LPattern=c(1,1,1),
           RPattern=c(1,0,1),
           nGene=200,
           fc="E10")
       )))

setParam(params, "lambda") <- 0.1
expect_true(identical(getParam(params, "lambda"), 0.1))

setParam(params, "seed") <- 123
expect_true(identical(getParam(params, "seed"), 123))
