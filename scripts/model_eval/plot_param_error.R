
f <- read.table('foo')
hist(f[[3]],20,freq=FALSE,main="param estimate error / sd")
curve(dnorm(x), add=TRUE,col="red")

