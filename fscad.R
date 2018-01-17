
fscad_beta = function(Y,lambda,gamma_r,gamma,bspi, reg_mat)
{
X = reg_mat

W = inprod(bspi,bspi,1,1)
R = inprod(bspi,bspi)
r = inprod(bspi,bspi,2,2)

range_knots = knots(bspi,interior=FALSE)%>%range

W0=  function(beta_j){
 
 zero = NULL

    wi = function(i){
    beta_m_norm =  t(as.matrix(beta_j))%*%W_m[[i]]%*%as.matrix(beta_j)%>%as.numeric%>%sqrt
            # print(beta_m_norm)
    if(beta_m_norm < 10^-8)   {
        zero = c(i:i+4)
        res = matrix(0, nrow=bspi$nbasis, ncol=bspi$nbasis)
    } else{
    temp1 = DiceKriging::SCAD.derivative(sqrt((length(W_m))/((bspi$rangeval[2]-bspi$rangeval[1]))) * beta_m_norm, lambda = lambda)
    temp2 = beta_m_norm * sqrt(((bspi$rangeval[2]-bspi$rangeval[1]))/(length(W_m)))
    res = temp1/temp2*W_m[[i]]

    }
    return(list(res=res, zero=zero))
    }
            
res_list = lapply(1:length(W_m), wi)
res = lapply(res_list, function(x) x[[1]])
zero = do.call(c,lapply(res_list, function(x) ifelse(is.null(x[[2]]),0,x[[2]])))%>%unique

W_j0= Reduce('+', res)*0.5
return(list( W = W_j0, zero = zero))
}

range_knots = knots(bspi,interior=FALSE)%>%range

beta_now = runif(bspi$nbasis,-10,10)


beta_before = rep(0,bspi$nbasis)
value_before = 10^8 

non_zeros = 1:bspi$nbasis

threshhold =1
i=1
stoprule = 1e-3
# print(threshhold)
while(threshhold >stoprule|i==1){ # this shold be 
# print(threshhold)
W0_now= W0(beta_now)

W0_res = W0_now[[1]]

non_zeros = intersect(non_zeros,setdiff(1:bspi$nbasis,W0_now$zero))
R_non= R[non_zeros,non_zeros]
r_non= r[non_zeros,non_zeros]

W0_res_non = W0_res[non_zeros,non_zeros]
X_non = X[,non_zeros]



# t(X)%*%Y

beta_now_non_zero  = solve(t(X_non)%*%X_non   + nrow(X)*gamma*R_non+ nrow(X)*gamma_r*r_non+nrow(X)*W0_res_non, t(X_non)%*%Y,tol=1e-20)

beta_now = rep(0, length=bspi$nbasis)
beta_now[non_zeros] = beta_now_non_zero[,1]

threshhold = max(abs(beta_now - beta_before)) 

# cat(threshhold)
beta_before = beta_now
# if ((beta_now%>%abs%>%max) < stoprule) {
#     stoprule = 1e-2*(beta_now%>%abs%>%max)
#     print(stoprule)
# }
i = i+1

# print(non_zeros)
# beta_now_mat = matrix(beta_now_all, nrow=ncol(snps),byrow=TRUE)
}

return(beta_now)
}
