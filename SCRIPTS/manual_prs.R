
# How to make a GRS command in R

# beta snp in a vector:
Beta.snp <- c(0.12897038129696, -0.132389188, -0.123102197, 0.215671536475509, 0.15785808461558, 0.135819723142535, -0.218156009803171, 0.201306856705035, 0.107059072293408) 
# Dosage output in rows: 
output.snp <- Output_excel_bestand

#option 1
grs=data.frame()

for (i in 5:ncol(output.snp)){
  
  index=colnames(output.snp)[i]
  grs_i=t(Beta.snp)%*%as.vector(unlist(output.snp[,i]))
  
  row=data.frame(index=index, grs_i=grs_i)
  grs=rbind(grs, row)}


#option 2
grs=data.frame(index=colnames(output.snp)[5:ncol(output.snp)],
                grs= as.vector(t(Beta.snp)%*%as.matrix(output.snp[,5:ncol(output.snp)])))

#option 2: generic script

start=5 #first column with dosage
stop=ncol(output.snp) #last column with dosage

index=colnames(output.snp)[start:stop]
grs_i=as.vector(t(Beta.snp)%*%as.matrix(output.snp[,start:stop]))

grs=data.frame(index=index, grs_i=grs_i)

# Export your R output to csv file
write.csv(grs1, "GRS1.csv")
# Look at your desktop >> All files


# How to import a SPSS file into R?
# with the following package: haven

