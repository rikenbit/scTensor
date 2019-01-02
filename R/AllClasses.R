setClass("CCSParams",
    slots = c(
        nGene = "numeric",
        nCell = "numeric",
        cciInfo = "list",
        lambda = "numeric",
        seed = "numeric"
        ),
    prototype = prototype(
        nGene=1000,
        nCell=c(50, 50, 50),
        cciInfo=list(
            nPair=500,
            CCI1=list(
                LPattern=c(1,0,0),
                RPattern=c(0,1,0),
                nGene=50,
                fc="E10"),
            CCI2=list(
                LPattern=c(0,1,0),
                RPattern=c(0,0,1),
                nGene=50,
                fc="E10"),
            CCI3=list(
                LPattern=c(0,0,1),
                RPattern=c(1,0,0),
                nGene=50,
                fc="E10")
            ),
        lambda = 1,
        seed = 1234
    )
)
