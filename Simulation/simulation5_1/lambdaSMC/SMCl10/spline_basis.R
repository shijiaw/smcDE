
tobs = times
#knots = tobs
#knots = c(1,45,90)
#knots = c(times[1], seq(times[2], max(times),by=2))
knots = seq(times[1], max(times),by=knotsNpoints)
nknots   = length(knots)
norder   = 4
nbasis   = length(knots) + norder - 2
bsbasis = create.bspline.basis(range(knots),nbasis,norder,knots)


# basis values at sampling points
basismat   = eval.basis(tobs, bsbasis)

# square of basis matrix (Phi). It is used frequent for optimization,
# so we calculate in advance
basismat2  = t(basismat)%*%basismat

# values of the first derivative of basis functions at sampling points
Dbasismat  = eval.basis(tobs, bsbasis, 1)

# values of the second derivative of basis functions at sampling points
D2basismat = eval.basis(tobs, bsbasis, 2)

# number of quadrature points (including two knots) between two knots
nquad                        = 5;

# set up quadrature points and weights
#[sinbasis, quadpts, quadwts] 
quadre = quadset(nquad, bsbasis)
quadpts=quadre$quadvals[,1]
quadwts=quadre$quadvals[,2]

# values of the second derivative of basis functions at quadrature points
D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
#D2quadbasismat = eval.basis(quadpts, bsbasis, 2)

#Rmat  = t(D2quadbasismat)%*%(D2quadbasismat*((quadwts)%*%matrix(1,1,nbasis))) 


psi <- basismat  # bsplineS(norder=4,nderiv=0)
psi1 <- Dbasismat # bsplineS(norder=4,nderiv=1)
#psi2 <- D2basismat # bsplineS(norder=4,nderiv=2)


VV <- diag(quadwts)
D0VVD0 <- t(D0quadbasismat)%*%VV%*%D0quadbasismat
D1VVD1 <- t(D1quadbasismat)%*%VV%*%D1quadbasismat
bsmat2 <- t(basismat)%*%basismat



