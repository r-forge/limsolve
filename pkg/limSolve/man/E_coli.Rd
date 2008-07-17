\name{E_coli}
\docType{data}
\alias{E_coli}
\title{An underdetermined linear inverse problem: the Escherichia Coli Core Metabolism Model}
\description{Input matrices and vectors for performing Flux Balance Analysis of the E.coli metabolism

(as from \url{http://gcrg.ucsd.edu/Downloads/Flux_Balance_Analysis}).

The original input file can be found in the package subdirectory \url{/inst/docs/E_coli.input}

There are 53 substances:

GLC,G6P,F6P,FDP,T3P2,T3P1,13PDG,3PG,2PG,PEP,PYR,ACCOA,CIT,ICIT,AKG,SUCCOA,
SUCC,FUM ,MAL,OA,ACTP,ETH,AC,LAC,FOR,D6PGL,D6PGC,RL5P,X5P,R5P,S7P,E4P,RIB,
GLX,NAD,NADH,NADP,NADPH,HEXT,Q,FAD,FADH,AMP,ADP,ATP,GL3P,CO2,PI,PPI,O2,COA,GL,QH2

and 13 externals:

Biomass,GLCxt,GLxt,RIBxt,ACxt,LACxt,FORxt,ETHxt,SUCCxt,PYRxt,PIxt,O2xt,CO2xt


There are 70 unknown reactions (named by the gene encoding for it):

GLK1,PGI1,PFKA,FBP,FBA,TPIA,GAPA,PGK,GPMA,ENO,PPSA,PYKA,ACEE,ZWF,PGL,GND,RPIA,RPE,
TKTA1,TKTA2,TALA,GLTA,ACNA,ICDA,SUCA,SUCC1,SDHA1,FRDA,FUMA,MDH,DLD1,ADHE2,PFLA,
PTA,ACKA,ACS,PCKA,PPC,MAEB,SFCA,ACEA,ACEB,PPA,GLPK,GPSA1,RBSK,NUOA,FDOH,GLPD,CYOA,
SDHA2,PNT1A,PNT2A,ATPA,GLCUP,GLCPTS,GLUP,RIBUP,ACUP,LACUP,FORUP,ETHUP,SUCCUP,
PYRUP,PIUP,O2TX,CO2TX,ATPM,ADK,Growth

The lsei model contains:
\item 54 equalities (Ax=B): the 53 mass balances (one for each substance) and one equation
that sets the ATP drain flux for constant maintenance requirements to a fixed value (5.87)
\item 70 unknowns (x), the reaction rates
\item 62 inequalities (Gx>h). The first 28 inequalities impose bounds on some reactions.
The last 34 inequalities impose that the reaction rates have to be
positive (for unidirectional reactions only).
\item 1 functions that has to be maximised, the biomass production (growth).

As there are more unknowns (70) than equations (54), there exist an infinite amount of solutions
(it is an underdetermined problem).
}
\usage{E_coli}
\format{A list with the matrices and vectors that constitute the mass balance problem:
\code{A}, \code{B}, \code{G} and \code{H} and

\code{Maximise}, with the function to maximise.

The columnames of \code{A} and \code{G} are the names of the unknow reaction rates;
The first 53 rownames of \code{A} give the names of the components (these rows consitute the mass balance equations).}


\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{


# 1. parsimonious (simplest) solution
pars <- lsei(E=E_coli$A,F=E_coli$B,G=E_coli$G,H=E_coli$H)$x

# 2. the optimal solution - solved with linear programming
#    as some unknowns can be negative -> the problem is reformulated
#    define X = x1 - x2, where x1>0 and x2>0
#    and solve for x1, x2.

     # double all columns (unknowns)
     AA   <- cbind(E_coli$A, -1*E_coli$A)
     BB   <- E_coli$B
     GG   <- cbind(E_coli$G, -1*E_coli$G)
     HH   <- E_coli$H
     Cost <- -1*c(E_coli$Maximise,-1*E_coli$Maximise)
     
     LP<-linp(E=AA,F=BB,G=GG,H=HH,Cost=Cost)
     Nx   <- ncol(E_coli$A)   # number unknowns

     # X = x1 - x2
     (Optimal <- LP$X[1:Nx] - LP$X[(Nx+1):(2*Nx)])
     

# 3.ranges of all unknowns, including the central value and all solutions
xr   <- xranges(E=E_coli$A,F=E_coli$B,G=E_coli$G,H=E_coli$H,
                central=TRUE, full=TRUE)

# plotting these values already shows the shape of the uncertainty region
plot(xr[1,],xr[2,])
# and the central value in red...
points(xr[1,3],xr[2,3],col="red",pch=16,cex=2)

# the central point is a valid solution:
X <- xr[,"central"]
max(abs(E_coli$A\%*\%X-E_coli$B))
min(E_coli$G\%*\%X-E_coli$H)

# 4. Sample solution space; the central value is a good starting point
#   for algorithms cda and rda - but these need many iterations
xs<-xsample(E=E_coli$A,F=E_coli$B,G=E_coli$G,H=E_coli$H,
            iter=50000,out=5000,type="rda",x0=X)$X
pairs(xs[,10:20],pch=".",cex=2)

# using mirror algorithm takes less iterations,
# but an iteration takes more time ; it is better to start in a corner...
# (i.e. no need to use X as starting value)
xs<-xsample(E=E_coli$A,F=E_coli$B,G=E_coli$G,H=E_coli$H,
            iter=10000,out=2500,jmp=50,type="mirror")$X
pairs(xs[,10:20],pch=".",cex=2)

}
\references{\item \url{http://gcrg.ucsd.edu/Downloads/Flux_Balance_Analysis }
\item Edwards,J.S., Covert, M., and Palsson, B.., " Metabolic Modeling of Microbes: the Flux Balance Approach", Environmental Microbiology, 4(3): pp. 133-140 (2002).

}
\keyword{datasets}

