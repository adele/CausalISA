library("pcalg")
library("dplyr")
library("corrplot")
library("DescTools")

#####pacotes##################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("graph")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")

###########################


#lendo banco
banco_FCI <- read.table("banco_FCI.txt", header=TRUE)
str(banco_FCI) #informação sobre o banco

#usando somente as variáveis fenotipicas para teste
#banco_FCI <- (banco_FCI[,1:10])
#banco_FCI <- banco_FCI[,-which(colnames(banco_FCI) %in% c("conjugal3","sexo","fumo"))]



#####  subgrafo para obesidade #######
banco_FCI <- select(banco_FCI, obesidade, rs2826583, rs6042994, faixa_etaria, pcr_cat )
#rs5868995, rs2826583, rs11678946
#str(banco_FCI)

###### subgrafo para pressao ######
#banco_FCI <- select(banco_FCI, has2, rs34500244, rs726164, conjugal3, sexo)
#rs34500244, rs10503700, rs12155372, rs36097237, rs9492001, rs881641, rs726164

#### subgrafo para diabetes #####
#banco_FCI <- select(banco_FCI, diabetes2, rs7638558, rs74062309, has2, conjugal3)
#rs4766012, rs11063692, rs74062309, rs34019131, rs60430780, rs7638558, rs9855667



#retirando as linhas com NA
banco_FCI <- banco_FCI[complete.cases(banco_FCI),]
dim(banco_FCI)

#criando a lista com a matriz das variaveis (banco)
suffStat <- list(dm = banco_FCI, adaptDF=TRUE)

#teste utilizado - variáveis discretas
disCItest (1, 2, S=NULL, suffStat)
indepTest <- disCItest

alpha <- 0.2
labels <- colnames(banco_FCI)


sink("verbose.txt") 
#rodando o algoritmo
resultado_FCI <- fci(suffStat, indepTest, alpha, labels,
    skel.method = "stable",
    type = "normal",
    fixedGaps = NULL, fixedEdges = NULL,
    NAdelete = FALSE, m.max = Inf, pdsep.max = Inf,
    rules = rep(TRUE, 10), doPdsep = TRUE, biCC = FALSE,
    conservative = FALSE, maj.rule = FALSE,
    numCores = 1, selectionBias = TRUE, 
    verbose = TRUE)
sink()

#plotar o grafo
plot(resultado_FCI)

##################################################################


### correlação entre obesidade e seus SNPs ######

D <- select(banco_FCI, c("obesidade", "rs2826583", "rs6042994", "rs5868995", "rs2826583", "rs11678946"))

#names(D) <- c("obesidade", "rs2826583", "rs6042994", "rs5868995", "rs2826583", "rs11678946")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

M <- cor(D,method = "pearson",use="pairwise.complete.obs")

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

D <- select(banco_FCI, c("obesidade", "rs2826583", "rs6042994", "rs5868995", "rs2826583", "rs11678946"))

corrplot(PairApply(D, CramerV),diag = F,is.corr = F, addCoef.col = "black")

#############################################################

#### correlação entre pressão e seus SNPs #####

D <- select(banco_FCI, c("has2", "rs34500244", "rs10503700", "rs12155372", "rs36097237", "rs9492001", "rs881641", "rs726164" ))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

M <- cor(D,method = "pearson",use="pairwise.complete.obs")

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

D <- select(banco_FCI, c("has2", "rs34500244", "rs10503700", "rs12155372", "rs36097237", "rs9492001", "rs881641", "rs726164" ))
corrplot(PairApply(D, CramerV),diag = F,is.corr = F, addCoef.col = "black")

###############################################################

#### correlação entre diabetes e seus SNPs ####

D <- select(banco_FCI, c("diabetes2", "rs4766012", "rs11063692", "rs74062309", "rs34019131", "rs60430780", "rs7638558", "rs9855667"  ))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

M <- cor(D,method = "pearson",use="pairwise.complete.obs")

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

D <- select(banco_FCI, c("diabetes2", "rs4766012", "rs11063692", "rs74062309", "rs34019131", "rs60430780", "rs7638558", "rs9855667"  ))
corrplot(PairApply(D, CramerV),diag = F,is.corr = F, addCoef.col = "black")

###############################################################

#### correlação dos SNPs ######
banco_FCI_c <- (banco_FCI[,11:33])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

M <- cor(banco_FCI_c,method = "pearson",use="pairwise.complete.obs")

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

####################################################################

### correlção entre dislipidemia e seus SNPs ###

D <- select(banco_FCI, c("dlp2","rs2644247", "rs347845", "rs4376525", "rs16878983", "rs4277749"))
#"rs2644247", "rs347845", "rs4376525", "rs16878983", "rs4277749" 

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

M <- cor(D,method = "pearson",use="pairwise.complete.obs")

corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

###############################################





# ## Simulate data
# n <- 100
# set.seed(123)
# x <- sample(0:2, n, TRUE) ## three levels
# y <- sample(0:3, n, TRUE) ## four levels
# z <- sample(0:1, n, TRUE) ## two levels
# dat <- cbind(x,y,z)
# dat[1,1]=NA
# 
# ## Analyze data
# gSquareDis(1,3, S=2, dat, nlev = c(3,4,2)) # but nlev is optional:
# gSquareDis(1,3, S=2, dat, verbose=TRUE, adaptDF=TRUE)
# ## with too little data, gives a warning (and p-value 1):
# gSquareDis(1,3, S=2, dat[1:60,], nlev = c(3,4,2))
# 
# suffStat <- list(dm = dat, adaptDF = TRUE)
# disCItest(1,3,2,suffStat)





