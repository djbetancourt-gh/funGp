# ----------------------------------------------------------------------------------------------------------
ysurfmax <- function(){
  # >---------< LOO using D_H + D_JB + D_M for training
  sIn.cv <- rbind(sIn.H, sIn.JB, sIn.M)
  fIn.cv <- Map(rbind, fIn.H, fIn.JB, fIn.M)
  sOut.cv <- rbind(sOut.H[,1,drop=F], sOut.JB[,1,drop=F], sOut.M[,1,drop=F])

  # only scalar inputs
  m1 <- funGp(sIn = sIn.cv, sOut = sOut.cv)
  plotLOO(m1)

  # scalar and all functional
  m2 <- funGp(sIn = sIn.cv, fIn = fIn.cv, sOut = sOut.cv)
  plotLOO(m2)

  # scalar and (Td + Sg + Hs + Tp) functional
  m3 <- funGp(sIn = sIn.cv, fIn = fIn.cv[1:4], sOut = sOut.cv)
  plotLOO(m3)

  # optimized sctructure
  xm <- funGp_factory(sIn = sIn.cv, fIn = fIn.cv, sOut = sOut.cv, setup = list(n.pop = 10), quietly = F)
  plotX(xm)
  plotEvol(xm)
  # hola
}
# ----------------------------------------------------------------------------------------------------------
