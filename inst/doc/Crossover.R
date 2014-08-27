## ----OptionsAndLibraries, include=FALSE, message=FALSE------------------------------------------------------------------------
if (exists("opts_chunk")) {
  opts_chunk$set(concordance=TRUE)
  opts_chunk$set(tidy.opts=list(keep.blank.line=FALSE, width.cutoff=95))
  #opts_chunk$set(size="footnotesize")
  #opts_chunk$set(size="tiny")
  opts_chunk$set(size="scriptsize") 
  opts_chunk$set(cache=TRUE)
  opts_chunk$set(autodep=TRUE)
}

# See http://yihui.name/knitr/hooks for the following code.
CrossoverNamespace <- function(before, options, envir) {
  if (before) {
    ## code to be run before a chunk
    #attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3)
    attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3, warn.conflicts=FALSE)
  } else {
    ## code to be run after a chunk
    detach("namespace:Crossover")
  }
}

if (exists("opts_chunk")) {
  knit_hooks$set(withNameSpace = CrossoverNamespace)
}

library(Crossover, quietly=TRUE)
options(width=128)
options(digits=4)
startGUI <- function(...) {invisible(NULL)}
#options(prompt="> ", continue="+ ")
library(MASS)
library(multcomp)
library(ggplot2)
library(Matrix)



## ----StandardAdditive, echo=TRUE----------------------------------------------------------------------------------------------
# Design:
design <- rbind(c(3,2,1),
                c(2,1,3),
                c(1,2,3),
                c(3,2,1))
design
v <- 3 # number of treatments
# Link matrix:
H <- Crossover:::linkMatrix(model="Standard additive model", v)
H
# Row-Column-Design: (cf. John et al. 2004, Table II and page 2649f.)
rcDesign <- Crossover:::rcd(design, v=v, model=1)
rcDesign
# Design Matrix of Row-Column Design:
Xr <- Crossover:::rcdMatrix(rcDesign, v, model=1)
Xr
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----FullInteractions, echo=TRUE----------------------------------------------------------------------------------------------

H <- Crossover:::linkMatrix(model="Full set of interactions", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----SelfAdjacency, echo=TRUE-------------------------------------------------------------------------------------------------

H <- Crossover:::linkMatrix(model="Self-adjacency model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------------------------
# Link matrix:
H <- Crossover:::linkMatrix(model="Placebo model", v, placebos=1)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----NoIntoSelf, echo=TRUE----------------------------------------------------------------------------------------------------

H <- Crossover:::linkMatrix(model="No carry-over into self model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----TreatmentDecay, echo=TRUE------------------------------------------------------------------------------------------------

H <- Crossover:::linkMatrix(model="Treatment decay model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----Proportionality, echo=TRUE-----------------------------------------------------------------------------------------------

H <- Crossover:::linkMatrix(model="Proportionality model", v)
H
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----SecondOrder, echo=TRUE---------------------------------------------------------------------------------------------------
# Link matrix:
H <- Crossover:::linkMatrix(model="Second-order carry-over effects", v)
H
# Row-Column-Design: (cf. John et al. 2004, Table II and page 2649f.)
rcDesign <- Crossover:::rcd(design, v=v, model=8)
rcDesign
# Design Matrix of Row-Column Design:
Xr <- Crossover:::rcdMatrix(rcDesign, v, model=8)
Xr
# Design Matrix of Cross-Over Design:
X <- Xr %*% H
X


## ----SearchStrategy, echo=TRUE, eval=TRUE, cache=TRUE, dev='png', dpi=150-----------------------------------------------------

set.seed(42)
x <- searchCrossOverDesign(s=9, p=5, v=4, model=4)
plot(x)
plot(x, type=2)


## ----attachNameSpace, echo=TRUE, eval=FALSE, include=FALSE--------------------------------------------------------------------
#   # We will use a lot of internal commands since we will test and evaluate things the normal user will probably not be interested in. Therefore we load and attach the namespace.
#  # attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3)
#  # When we are finished we call 'detach("namespace:Crossover")'.

## ----TestOfDifferentApproaches, echo=TRUE, eval=TRUE, withNameSpace=TRUE------------------------------------------------------

attach(loadNamespace("Crossover"), name="namespace:Crossover", pos=3, warn.conflicts=FALSE)

s <- 6
p <- 3
v <- 3
model <- 1
data(williams)
design <- williams3t
  
rcDesign <- rcd(design, v, model)
# JRW, p 2650, first equation on that page, whithout number
Ar <- infMatrix(rcDesign, v, model)
Xr <- rcdMatrix(rcDesign, v, model)
# JRW, p 2650, second equation on that page, number 11
Ar2 <- t(Xr) %*% (diag(s*p)-Crossover:::getPZ(s,p)) %*% Xr
max(abs(Ar-Ar2))

fXr <- cbind(Xr, getZ(s,p))
Ar3 <- t(fXr) %*% fXr
ginv(Ar3)[1:12,1:12]-ginv(Ar2)

H <- linkMatrix(model=model,v=v)
fX <- cbind(Xr%*%H, getZ(s,p))
A1 <- t(fX) %*% fX
A2 <- t(H)%*%Ar%*%H

# While A1 and A2 differ:

ginv(A1)[1:6,1:6]
ginv(A2)
max(abs(ginv(A1)[1:6,1:6]-ginv(A2)))

# The variances for the estimable contrasts are the same:

C <- matrix(0,nrow=15,ncol=1)
C[1:2,1] <- c(-1,1)
tdiff1 <- t(C)%*%ginv(A1)%*%C
tdiff2 <- t(C[1:6,])%*%ginv(A2)%*%C[1:6,]
tdiff1 - tdiff2

C <- matrix(0,nrow=6,ncol=1)
C[1:2,1] <- c(-1,1)
tdiff1 <- t(C)%*%ginv(A1)[1:6,1:6]%*%C
tdiff2 <- t(C)%*%ginv(A2)%*%C
tdiff1 - tdiff2


## ----detachNameSpace, echo=TRUE, eval=FALSE, include=FALSE--------------------------------------------------------------------
#  detach("namespace:Crossover")

