arg = commandArgs(trailingOnly=T)

lst = read.table(arg[1],as.is=T)[,1]
names = gsub(".wgt.RDat","",basename(lst))
N = length(lst)

model = c() # <name of all models>
colnames = c( "id" , "nsnps" , "hsq" , "hsq.se" , "hsq.pv" , paste(model,"r2",sep='.') , paste(model,"pv",sep='.') )

mat.snps = matrix(nrow=N,ncol=1)
mat.hsq = matrix(nrow=N,ncol=2)
mat.hsqpv = matrix(nrow=N,ncol=1)
mat.mod.rsq = matrix(nrow=N,ncol=length(model))
mat.mod.pv = matrix(nrow=N,ncol=length(model))

for ( i in 1:N ) {
load(lst[i])
m = match( model , colnames(cv.performance) )
mat.snps[i] = nrow(snps)
mat.hsq[i,] = hsq_afr
mat.hsqpv[i] = hsq_afr.pv
mat.mod.rsq[i,] = cv.performance[1,m]
mat.mod.pv[i,] = cv.performance[2,m]
print(cv.performance)
}


write.table( 
cbind( names , mat.snps , format(mat.hsq,digits=2) , format(mat.hsqpv,digits=2) , format(round(mat.mod.rsq,3),digits=3) , format(mat.mod.pv,digits=2) ),
file = "magepro.wgt.profile.txt",
quote=F,
col.names=colnames,
row.names=F,sep='\t')
best = apply(mat.mod.rsq,1,max,na.rm=T)

options(digits=3)
cat( "Average hsq:" , mean(mat.hsq[,1]) , '(' , sd(mat.hsq[,1])/sqrt(N) , ') \n' , file=stderr() )

cat( "\nAverage crossvalidation R2:\n" , file=stderr() )
write.table( cbind( format(apply(mat.mod.rsq,2,mean,na.rm=T),digits=3) , format(apply(mat.mod.rsq,2,sd,na.rm=T) / sqrt(N),digits=3) ), quote=F , col.names=c("R2","SE") , row.names=model , sep='\t' , file=stderr() )
cat( "BEST" , format(mean(best),digits=3) , '\n' , sep='\t' , file=stderr() )

cat( "\n% Model is best:\n" , file=stderr() )
for ( i  in 1:length(model) ) {
cat( model[i] , ": " , 100*mean( mat.mod.rsq[,i] == best , na.rm=T ) , "%" , '\n' , sep='' , file=stderr() )
}
