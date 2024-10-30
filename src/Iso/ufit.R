if (!require("Iso")) {
  install.packages("Iso")
}
library(Iso)

y <- c(0,1,2,3,3,2)
f1 <- ufit(y,lmode=0.4) # The third entry of the default
                        # value of x = c(0.0,0.2,0.4,0.6,0.8,1.0).
f2 <- ufit(y,imode=3)   # Identical to f1.
f3 <- ufit(y,lmode=3,x=1:6)   # Effectively the same as f1 and f2.
                              # But is different  in appearance.
f4 <- ufit(y,imode=3,x=1:6)   # Identical to f3.

x <- c(0.00,0.34,0.67,1.00,1.34,1.67,2.00,2.50,3.00,3.50,4.00,4.50,
       5.00,5.50,6.00,8.00,12.00,16.00,24.00)
y <- c(0.0,61.9,183.3,173.7,250.6,238.1,292.6,293.8,268.0,285.9,258.8,
       297.4,217.3,226.4,170.1,74.2,59.8,4.1,6.1)
z <- ufit(y,x=x,type="b")
plot(x,y)
lines(z,col="red")
plot(z$h,do.points=FALSE,col.hor="blue",col.vert="blue",add=TRUE)
abline(v=z$mode,col="green",lty=2)
