#Function for the ODEs

viral4dS <- function(t,y,p){
  #Parameters to vary (vector p)
  ka <- p[[1L]]
  Imax <- exp(p[[2L]])
  #Fixed params
  beta <- 0.8
  ki <- 1
  pp <- 80
  gamma <- 13
  A50 <- 1
  kc <- 0.33 
  #RHS for the ODEs. In this order: I, V, A1, A2, A3
  c(beta*y[[2L]]*(1-(y[[1L]]/Imax)) - ka*y[[1L]]*(y[[5L]]**2)/((y[[5L]]**2)+(A50**2)),
    pp*y[[1L]] - gamma*y[[2L]],
    ki*(y[[1L]])/((y[[1L]]) + (1000)) - kc*y[[3L]], #was 1000
    kc*y[[3L]] - kc*y[[4L]],
    kc*y[[4L]])
}