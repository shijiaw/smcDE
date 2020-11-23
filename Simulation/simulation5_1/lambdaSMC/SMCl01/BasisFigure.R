library(ggplot2)
library(magrittr)
library(reshape2)
library(stringr)
autoplot.basis <- function(basis, n=1000) {
  all.knots <- sort(c(attr(basis,"Boundary.knots") ,attr(basis, "knots"))) %>%
    unname
  bounds <- range(all.knots)
  knot.values <- predict(basis, all.knots) %>%
    set_colnames(str_c("S", seq_len(ncol(.))))
  newx <- seq(bounds[1], bounds[2], length.out = n+1)
  interp.values <- predict(basis, newx) %>%
    set_colnames(str_c("S", seq_len(ncol(.))))
  knot.df <- data.frame(x=all.knots, knot.values) %>%
    melt(id.vars="x", variable.name="Spline", value.name="y")
  interp.df <- data.frame(x=newx, interp.values) %>%
    melt(id.vars="x", variable.name="Spline", value.name="y")

  ggplot(interp.df) +
    aes(x=x, y=y, color=Spline, group=Spline) +
    geom_line() +
    geom_point(data=knot.df) +
    scale_color_discrete(guide=FALSE)+ theme_bw() +rremove("x.text")+rremove("ylab")+ xlab('Cubic B-spline basis')

}

library(splines)
x <- seq(0, 1, by=0.001)
spl <- bs(x,df=13)

gname = c("CubicBspline.eps",sep="")  
postscript(gname,width=6,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
autoplot(spl)
dev.off()


