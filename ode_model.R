

viral4dS <- function(t,y,p){
  ka <- p[[1L]]
  Imax <- exp(p[[2L]])
  beta <- 0.8
  ki <- 1
  pp <- 80
  gamma <- 13
  A50 <- 1
  kc <- 0.33 # was 0.4 before 10000 - > 1000
  c(beta*y[[2L]]*(1-(y[[1L]]/Imax)) - ka*y[[1L]]*(y[[5L]]**2)/((y[[5L]]**2)+(A50**2)),
    pp*y[[1L]] - gamma*y[[2L]],
    ki*(y[[1L]])/((y[[1L]]) + (1000)) - kc*y[[3L]], #was 1000
    kc*y[[3L]] - kc*y[[4L]],
    kc*y[[4L]])
}
#New IC function for model S (above)
# Ic_fnS <- function(y, ts){
#   IC <- NA
#   slop <- 0.24
#   m1 <- Matrix(c(0,0.8,80,-13), nc = 2)
#   e1 <- expm(4.2*m1)
#   mmy <- min(which(y>1))
#   #mmy <- which.max(y)
#   if(ts[mmy] > 6){ #5 + 1
#     aux3 <- 10**(log10(y[mmy]) + (ts[mmy]-6)*slop  )
#     IC <- aux3/e1[2,2]
#     if(ts[mmy] > 11){
#       aux3 <- 10**(log10(y[mmy]) + (5)*slop  ) #i.e. put a limit on interpolation
#       IC <- aux3/e1[2,2]
#     }
#     #print('late!')
#   }else{ 
#     aux3 <- y[1]#max(y)
#     IC <- aux3/e1[2,2]
#     #print('early!')
#     if(ts[mmy] < 5){ #adjust if first data pt has t<0
#       e1 <- expm((ts[mmy])*m1)
#       aux3 <- y[1]#max(y)
#       IC <- aux3/e1[2,2]
#     }
#   }
#   return(IC)
# }
