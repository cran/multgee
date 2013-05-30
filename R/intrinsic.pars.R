intrinsic.pars <-
function(response,id,repeated,ncategories,rscale="ordinal")
{
rscales <- c("ordinal", "nominal")
icheck <- as.integer(match(rscale, rscales, -1))
if (icheck < 1) {
        stop("'rscale' must be 'ordinal' or 'nominal'")
                 }
cdata <- datacounts(response,id,repeated,ncategories)
if(rscale=="ordinal"){
cmod <- gnm(counts~(factor(x)+factor(y))*factor(tp)+factor(tp):x:y,
                                             family=poisson,data=cdata)
ans <- as.vector(coef(cmod)[pickCoef(cmod,"x:y")])
                     } else {
ans <- rep(0,max(cdata$tp))
cdata$x <- factor(cdata$x)
cdata$y <- factor(cdata$y)
for(i in 1:max(cdata$tp))
{
cmod <- gnm(counts~x+y+MultHomog(x,y),
       family=poisson,data=cdata[cdata$tp==i,])
cscores <- coef(cmod)[pickCoef(cmod,"MultHomog")]
cscores <- c(tcrossprod(normscores(cscores)))
cmod <- gnm(counts~factor(x)+factor(y)+cscores,
        family=poisson,data=cdata[cdata$tp==i,])
ans[i] <- coef(cmod)["cscores"]
}
}
ans
}
