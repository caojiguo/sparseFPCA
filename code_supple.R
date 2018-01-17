dir  = "/Users/cao/Dropbox/Research/Yunlong/sparse_FPCA/JCGS2018Jan10/SimulationCodes"

# Step 0: install the required packages 
install.packages('fda')
install.packages('dplyr')
install.packages('MASS')
install.packages('DiceKriging')
library(fda)
library(dplyr)
library(MASS)
library(DiceKriging)
#####################
# Step 1: Simulation Setup
#####################

# Step 1.1: Load the true PCs
load(sprintf("%s/simulated_pcs_updated.rda",dir), verbose = 1)

# Step 1.2: Simulate the PC scores 
seed=1
set.seed(seed)
n=100; 
noise_sd=1
mu = rep(0,4)
Sigma=diag(c(30,20,10,3))
scores =mvrnorm(n = 100, mu, Sigma,empirical = TRUE)
coef_simfd = do.call(cbind,lapply(1:nrow(scores),function(i){
	coef = scores[i,]
	(matrix(rep(coef,nrow(coef(simulated_pcs))),nrow=nrow(coef(simulated_pcs)),byrow=TRUE)*coef(simulated_pcs))%>%rowSums+rnorm(nrow(coef(simulated_pcs)),sd=noise_sd)

}))
simfd = (fd(coef_simfd,simulated_pcs$basis))

# Step 1.3: Generate the observations 
xmat = eval.fd(seq(0,60,len=20),simfd) 
xmat = xmat%>%t 

#####################
# Step 2: Estimate the sparse FPCs 
#####################

# Step 2.1: Express each curve's observations using bspline functions
source(paste0(dir, '/fscad.R'))
basisobj = create.bspline.basis(c(0,60), norder=4, nbasis=10) 
simfds = Data2fd(seq(0,60,len=20),xmat%>%as.matrix%>%t, basisobj)

# Step 2.2: Set the initial values of alpha(t) using the conventional FPCs
D2Lfd <- int2Lfd(m=2)
D2fdPar<- fdPar(basisobj, D2Lfd, 1e2)
fpca_basis = pca.fd(simfds,4,D2fdPar,centerfns = F)
total_var=  mean(fpca_basis$scores[,1]^2)/fpca_basis$varprop[1]
(lambdas_pc_basis= colMeans(fpca_basis$scores^2))
npc = 4
pca_sim= pca.fd(simfds,npc,D2fdPar,centerfns = F) 
pca_sim$harmonics[[1]] = -pca_sim$harmonics[[1]]
(lambdas_pc= colMeans(pca_sim$scores^2))

# Step 2.3 Compute the sparse FPCs
wtemp = function(xbasis0) {
        xbasis0 %>% knots(, interior = FALSE) %>% unique %>% 
            data.frame(knot = .) %>% mutate(knotlead = lead(knot)) %>% 
            dplyr::filter(!is.na(knotlead)) %>% rowwise() %>% 
            do(temp = eval.penalty(xbasis0, int2Lfd(0), rng = c(.$knot, 
                .$knotlead)))
    }

W_m = wtemp(basisobj)$temp
W = inprod(basisobj, basisobj)
reg_mat  = inprod(simfds,basisobj)

npc=4
gamma_ridge = 10
gamma_R=1e2
lam1= 20
lambda_1 = rep(lam1,npc); 

Y_mats = inprod(simfds, pca_sim$harmonics)
coefs_fd= do.call(cbind,lapply(1:ncol(Y_mats), function(x){
 fscad_beta(Y_mats[,x],lambda=lambda_1[x],gamma_r=gamma_R,gamma=gamma_ridge,basisobj,reg_mat )

}))

difference_before = -1000
difference=1
threshold=1
j = 1
fitted_simfd_before = fd(matrix(0,coef(simfds)%>%nrow,coef(simfds)%>%ncol),basisobj)
afds_before = fd(matrix(0, basisobj$nbasis,npc),basisobj)
Y_mats_before  = Y_mats
Y_mats_before[]=0

while(threshold>1e-2){
print(threshold)
if (j>20) break
j=j+1
ests = fd(coefs_fd, basisobj);
if(any(inprod(ests)%>%as.numeric<0)) {
	index_flip = which(inprod(ests)%>%as.numeric<0)
	for (i in index_flip){
	coefs_fd[,i] = - coefs_fd[,i]	
	}
	ests = fd(coefs_fd, basisobj);

}
estsn= fd(coef(ests)*matrix(rep(sqrt(1/diag(inprod(ests, ests))),nrow(coef(ests))), nrow=nrow(coef(ests)),byrow=TRUE), ests$basis)
betascores = inprod(ests, fpca_basis$harmonics)
betascores= betascores*matrix(rep(lambdas_pc_basis,nrow(betascores)),nrow=nrow(betascores),byrow=TRUE)
pca_coefs= coef( fpca_basis$harmonics)

betatilda_coefs = do.call(cbind,lapply(1:nrow(betascores), function(x){
rowSums(matrix(rep(betascores[x,], nrow(pca_coefs)),nrow=nrow(pca_coefs),byrow=T)*pca_coefs)    
}))

betatilda = fd(betatilda_coefs, pca_sim$harmonics$basis);

D2Lfd <- int2Lfd(m=2)

betatilda_pca = pca.fd(betatilda,npc,centerfns=FALSE)

betatilda_pca_scores = apply(betatilda_pca$scores,1,function(x) x/sqrt(colSums(betatilda_pca$scores^2)))%>%t


pca_coefs_new = coef(betatilda_pca$harmonics)
afd_coefs = do.call(cbind,lapply(1:nrow(betatilda_pca_scores ), function(x){
rowSums(matrix(rep(betatilda_pca_scores [x,], nrow(pca_coefs_new)),nrow=nrow(pca_coefs_new),byrow=T)*pca_coefs_new) 
}))

afds = fd(afd_coefs, betatilda_pca$harmonics$basis);
if(any(inprod(afds)%>%as.numeric<0)) {
	index_flip = which(inprod(afds)%>%as.numeric<0)
	for (i in index_flip){
	afd_coefs[,i] = - afd_coefs[,i]	
	}
	afds = fd(afd_coefs, betatilda_pca$harmonics$basis);

}




Y_mats = inprod(simfds,afds)
threshold=1
print(j)
if (j>2){
Y_mats_before = Y_mats_before[,apply(cor(Y_mats_before,Y_mats), 1, which.max)]
threshold =abs(Y_mats_before - Y_mats)%>%max
} 

afds_before = afds
Y_mats_before = Y_mats

coefs_fd= do.call(cbind,lapply(1:ncol(Y_mats), function(x){
 fscad_beta(Y_mats[,x],lambda=lambda_1[x],gamma_r=gamma_R,gamma=gamma_ridge ,basisobj,reg_mat)

}))

ests = fd(coefs_fd, basisobj);

fitted_simfds = fd((coef(afds)%*%inprod(ests,simfds)),afds$basis)


difference = (inprod(simfds-fitted_simfds,simfds-fitted_simfds))%>%diag%>%mean + gamma_ridge *sum(diag(inprod(ests,ests)))+sum(DiceKriging::SCAD(x=diag(inprod(ests, ests)), lambda=lambda_1))+gamma_R*sum(diag(inprod(ests,ests,2,2)))

difference_before= difference
}

########################
# Step 3: Summary the results 
########################

# Step 3.1:  Compute the sparse FPC scores 
scores0 = inprod(simfds,estsn)
vars0 = apply(scores0,2,var)
scores= scores0[,vars0%>%order(decreasing = TRUE)]
estsn = estsn[vars0%>%order(decreasing = TRUE)]


# Step 3.2: Compute the adjusted variation explained 
varexp = c()
for (i in 1:ncol(scores)){
	if (i==1) varexp[i] = var(scores[,i])  else {
		varexp[i]  = var((lm(scores[,i]~scores[,1:(i-1)])%>%residuals))
	}
}


# Step 3.3: Plot the spase FPCs against the true FPCs
plot(simulated_pcs,ylim=c(-0.2,0.5),lty=1,main="sFPCs")
lines(estsn,lwd=2,lty=2)

# Step 3.4: Compute the AIC value 
coef_mat = lapply(1:nrow(scores0), function(x){
	rowSums(matrix(rep(scores0[x,], each=estsn$basis$nbasis), nrow=estsn$basis$nbasis)*coef(estsn))
})%>%do.call(cbind,.)
simfd_fit = fd(coef_mat,estsn$basis )
mse = mean((inprod(simfd_fit-simfds,simfd_fit-simfds )%>%diag/60)%>%sqrt)
k = sum(abs(estsn%>%coef)>1e-3)
n=length(simfds$fdnames$reps)
AIC= 2*k+ n*log(mse)


