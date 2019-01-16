#this file is for core NMF functions including:
#NMF
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Update of W matrix in iteration
#'
#' `Update_W_Tang` conducts one simple iteration of W in the NMF process, however, it is designed as
#' a sub-function for NMF procedure, do not run it indenpendently
#'
#' @param Y  Y is the gene expression matrix that has #{signature genes} rows and #{samples} columns
#'
#' @param W_itr W matrix from previous step of iteration
#'
#' @param H_itr H matrix from previous step of iteration
#'
#' @param beta_W $\beta$ is the parameter controlling sparsity ridge penalty, see our paper for more information
#'
#' @param lam_W $\lambda$ is the parameter imposing A's structure into resulting W, see our paper for more info
#' @param A The reversed trichotomous guide matrix \tilde{A}
#'
#' @return 'Update_W_Tang' returns the udpated W matrix after iteration
#' @examples This function is not supposed to be run independently
#' @export

Update_W_Tang <- function(Y,W_itr,H_itr,beta_W,lam_W,A){
  #dimension
  m <- nrow(Y)
  n <- ncol(Y)
  r <- ncol(W_itr)

  #generate diagonal matrix
  I_rr <- diag(x=rep(1,r))

  #numerator
  num <- Y%*%t(H_itr)
  denom <- W_itr%*%H_itr%*%t(H_itr)+beta_W*W_itr%*%I_rr+lam_W*A
  denom = pmax(denom, 0.000000000000001*array(1,c(dim(denom)[1],dim(denom)[2])))

  #update
  W_new <- W_itr * num/denom
  W_new[which(W_new < 0)] <- 0
  ##return
  return(W_new)
}

#' Update of H matrix in iteration
#'
#' `Update_H_Tang` conducts one simple iteration of H in the NMF process, however, it is designed as
#' a sub-function for NMF procedure, do not run it indenpendently
#'
#' @param Y  Y is the gene expression matrix that has #{signature genes} rows and #{samples} columns
#'
#' @param W_itr W matrix from previous step of iteration
#'
#' @param H_itr H matrix from previous step of iteration
#'
#' @param eta $\lambda$ is the parameter for ensuring column sums of H are 1s
#'
#' @return 'Update_H_Tang' returns the udpated H matrix after iteration
#' @examples This function is not supposed to be run independently
#' @export
Update_H_Tang <- function(Y,W_itr,H_itr,eta){
  #dimension
  m <- nrow(Y)
  n <- ncol(Y)
  r <- ncol(W_itr)

  #generate diagonal matrix
  E_rr <- matrix(rep(1,r*r),ncol=r)
  E_rn <- matrix(rep(1,r*n),ncol=n)

  #update
  num <- t(W_itr)%*%Y  + eta*E_rn
  denom <- t(W_itr)%*%W_itr%*%H_itr + eta*E_rr%*%H_itr #+ beta_H*H_itr
  denom= pmax(denom, 0.000000000000001*array(1,c(dim(denom)[1],dim(denom)[2])))
  H_new <- H_itr*(num/denom)
  H_new[which(H_new<0)] <- 0

  #return
  return(H_new)

}

#' NMF core function
#'
#' `NMF_RNA_Tang` is the core function for NITUMID's factorization, it takes the gene expression matrix Y,
#' guide matrix A, initial W and H, as well as other parameters as input, and optimizes til converged/ or
#' reached maximum steps
#'
#' @param Y  Y is the gene expression matrix that has #{signature genes} rows and #{samples} columns
#'
#' @param ini.W Initial W matrix
#'
#' @param ini.H Initial H matrx
#'
#' @param beta_W $\beta$ is the parameter controlling sparsity ridge penalty, see our paper for more information. In
#' real calculation, upper level functions will automatically choose suitable $\beta$ no need to set by yourself
#'
#' @param lam_W $\lambda$ is the parameter imposing A's structure into resulting W, see our paper for more info. In
#' real calculation, upper level functions will automatically choose suitable $\lambda$ no need to set by yourself
#'
#' @param A The reversed trichotomous guide matrix \tilde{A}
#'
#' @param eta $\eta$ is the parameter for ensuring column sums of H are 1s, upper level function will choose suitable
#' $\eta$ automatically, no need to spcify yourself
#'
#' @param tol numeric scalar, the cutoff value for stop iterating, default $10^{-5}$
#'
#' @param prterr logical value indicating if it should print each step's error, default FALSE
#'
#' @param max.itr numerical scalar, maximum iteration number, default 2000
#'
#' @param alaways logical, if it should output results whatsover when max.itr is reached, default TRUE
#'
#' @return 'NMF_RNA_Tang' returns the resulted W and H, as well as track of errors along iteration
#'
#' @examples This function is not supposed to be run independently, it's supposed to be called by upper level functions
#' @export
NMF_RNA_Tang <- function(Y,ini.W,ini.H,lam_W, beta_W, A, eta, tol=10^-5,prterr=F,max.itr=2000,always=T){
  #obtain the dimensions
  m <- nrow(ini.W)
  r <- ncol(ini.W)

  #initialization
  itr_num <- 0
  W_itr <- ini.W
  H_itr <- ini.H
  #record err
  err_change =NULL; err_change_diff =NULL; err_change_W =NULL ; err_change_H =NULL

  while (itr_num < max.itr  ){
    W_itr0=W_itr; H_itr0=H_itr;
    H_itr=Update_H_Tang(Y = Y,W_itr = W_itr,H_itr = H_itr,eta = eta)
    W_itr <- Update_W_Tang(Y = Y, W_itr = W_itr,H_itr = H_itr,beta_W = beta_W, lam_W=lam_W,A=A )

    err <-   norm(Y-W_itr%*%H_itr,type = "f")
    err_diff= norm(W_itr%*%H_itr-W_itr0%*%H_itr0,type = "f")
    err_W <- norm(W_itr-W_itr0,type = "f")
    err_H = norm(H_itr-H_itr0,type = "f")
    err_change=c(err_change,log(err));
    err_change_diff=c(err_change_diff,log(err_diff));
    err_change_W=c(err_change_W,log(err_W));
    err_change_H=c(err_change_H,log(err_H));

    itr_num <- itr_num + 1
    if (prterr==T){print(err)}
    #if (err_diff < tol ) {return(list(W=W_itr,H=H_itr,err=err_change,step=itr_num,err_diff=err_change_diff,err_W=err_change_W, err_H=err_change_H))}

  }
  if (always==T){
    return(list(W=W_itr,H=H_itr,err=err_change,step=itr_num,err_diff=err_change_diff,err_W=err_change_W, err_H=err_change_H))
  }

  else {return(NA)}
}



##NMF Wrapped up Function
##This function includes data_initialization, data_scaling and normalizations

#' NMF core function
#'
#' `NMF_RNA_Tang_ind` is a upper level wrapper function for `NMF_RNA_Tang`. On top of core NMF, it also initialtes $W$ and $H$, and scales/normalizes
#' according to the input arguments: scale_Data. Besides, it also does random permutations to check results stability and helps choose the
#' best parameters. Although it has full sets of NITUMID's functions, the `NMF_RNA_Tang_ind` is still not supposed to be run independently, has hyper-function
#' will call it with different parameters choices.
#'
#' @param Y  Y is the gene expression matrix that has #{signature genes} rows and #{samples} columns
#'
#' @param ini.W Initial W matrix
#'
#' @param ini.H Initial H matrx
#'
#' @param beta_W $\beta$ is the parameter controlling sparsity ridge penalty, see our paper for more information. In
#' real calculation, upper level functions will automatically choose suitable $\beta$ no need to set by yourself
#'
#' @param lam_W $\lambda$ is the parameter imposing A's structure into resulting W, see our paper for more info. In
#' real calculation, upper level functions will automatically choose suitable $\lambda$ no need to set by yourself
#'
#' @param A The reversed trichotomous guide matrix \tilde{A}
#'
#' @param eta $\eta$ is the parameter for ensuring column sums of H are 1s, upper level function will choose suitable
#' $\eta$ automatically, no need to spcify yourself
#'
#' @param tol numeric scalar, the cutoff value for stop iterating, default $10^{-5}$
#'
#' @param prterr logical value indicating if it should print each step's error, default FALSE
#'
#' @param max.itr numerical scalar, maximum iteration number, default 2000
#'
#' @param alaways logical, if it should output results whatsover when max.itr is reached, default TRUE
#'
#' @param num_cell numeric scalar, number of component cell types, dfault 11. Noted that if you were to change component
#' cell types number, you should change your $A$ matrix accordingly.
#'
#' @param scale_Data the parameter controlling which scaling method to use for $Y$: 0 indicates no scaling, 1 indicates
#' scaling $||Y||_{F}=\sqrt{nm}$, where $n,m$ are dimensions of Y; 2 indicates scaling such that $||Y^{scaled}||_{F}=||log(Y^{raw}+1)||_{F}$
#'
#' @param perm_type logical variable, default TRUE. If TRUE, the function will do a permutation and measure consistency between results.
#'
#' @import gtools
#' @import rlist
#'
#' @return 'NMF_RNA_Tang' returns the resulted W and H, as well as track of errors along iteration
#'
#' @examples This function is not supposed to be run independently, it is supposed to be called by top level function, which will feed it with different pre-determined
#' parameters and choose thr best one
#'
#' @export
NMF_RNA_Tang_ind <- function(Y,lam_W, beta_W, A,eta, tol=10^-5,prterr=F,max.itr=2000,always=F,num_cell=11,scale_Data, perm_type=T){


  #make sure Y is a matrix bc sometimes we have only 1 sample
  Y <- as.matrix(Y)
  n1=dim(Y)[1]; n2=dim(Y)[2];



  #need to specify data scaling method, 0 means do not scale, 1 means scale by Frobeniues Norm, 2 means scaling by log-frobeniues norm
  if (scale_Data==0){
    Y_scale=as.matrix(Y) }
  else if (scale_Data==1){
    Y_scale=as.matrix(sqrt(n1*n2)*Y/norm(Y,'f'));}
  else {
    log_norm=norm(log(Y+1),'f')
    Y_scale <- as.matrix(log_norm*Y/norm(Y,'f'))
  }


  ##add a little noise to make results more robust
  if (perm_type==T){
    #browser()
    Y_scale=Y_scale+ matrix(rnorm(dim(Y_scale)[1]*dim(Y_scale)[2],0, 0.1*sd(Y_scale)),ncol=dim(Y_scale)[2])}

  Y_scale=pmax(Y_scale, array(0,c(dim(Y_scale)[1],dim(Y_scale)[2])))


  #browser()
  ini.W=abs(matrix(rnorm(dim(Y_scale)[1]*num_cell),dim(Y_scale)[1]))
  ini.H=abs(matrix(rnorm(num_cell*dim(Y_scale)[2]),num_cell))
  ini.W=1-A
  ini.H=array(1,c(dim(ini.H)[1],dim(ini.H)[2]))
  #NMF
  return(NMF_RNA_Tang(Y_scale,ini.W,ini.H,lam_W, beta_W, A, eta, tol=tol,prterr=F,max.itr=max.itr))

}




##If do not use test mode, then set test_mode=F and do not need to input a cell_structure
##Always remember A_input = 1-A_ref

#' NITUMID top level function
#'
#' `NITUMID` is the top level function that should be used by user. It takes gene expression, trichotomous guide matrix and other signature
#' genes' information as input, and output the deconvolution results. It has two different modes for pure cell's gene expression profile and
#' bulk gene expression profiles. Also, it has a test mode that allows you to input data as well as known cells fractions and test NITUMID's
#' performance.
#'
#' @param Y  Y is the gene expression matrix that has #{signature genes} rows and #{samples} columns
#'
#' @param if.bulk Logical, if Y is bulk gene expression data. If Y is purified cell gene expression, please set it as FALSE.
#'
#' @param A The reversed trichotomous guide matrix \tilde{A}
#'
#' @param tol numeric scalar, the cutoff value for stop iterating, default $10^{-5}$
#'
#' @param max.itr numerical scalar, maximum iteration number, default 2000
#'
#' @param num_cell numeric scalar, number of component cell types, dfault 11. Noted that if you were to change component
#' cell types number, you should change your $A$ matrix accordingly.
#'
#' @param test_mode Logical, default FALSE. This argument indicates if running NITUMID on test mode, which enables you to input
#' known true cells' fractions and measure NITUMID's outcomes' correlation with underlying true. For details, please see function
#' description as well as information regarding to arguments `real_H`, `cell_structure` and `correlation`
#'
#' @param real_H Optional, only needed when `test_mode`=TRUE, this should be a matrix with #{cell types} rows and #{samples} columns
#' indicating each cell types' fraction in each sample
#'
#'
#' @param row_index Default c(1,2,3,...,53). A numeric vector, indicating among all the siganture genes, which of them CAN be found in input $Y$, when you have
#' gene expression matrix $Y$, you could use `Signature_Match` function from this package to do the searching for you. A suggested way is to first create a vector
#' by seq(1,m), where m is the number of signature genes, and use `Signature_Match`, which will return a vector called `missing_row_index` then you just need to remove
#' those indexs from seq(1,m) by setdiff(seq(1,m),missing_row_index)
#'
#' @param cell_structure Optional, only needed when `test_mode`=TRUE, it is a dataframe with 11 rows and 2 variables: `origina_index` and
#' `destin_index`, you can modify index in `destin_index` to customize your cell types' matching relationship with NITUMID's 11 cell types.
#' See more information and samples from description of dataset `cell_structure_example` that comes with this package.
#'
#' @param correlation="p" Optional, a character of "p","s" or "kl",only needed when `test_mode`=TRUE. When you are in test mode with underlying
#' true cells' fractions input, you can evaluate their consistency with NITUMID's outcome by Pearson Correlation ("p"), Spearman Correlation ("s")
#' or K-L Divergence ("kl")
#'
#' @param row_mean a numeric vector speciying each signature gene's mean expression level across all cell types, defualt value is row_mean_v5 from
#' us
#'
#' @import gtools
#' @import rlist
#'
#' @return 'NITUMID' returns a list, depending on what specific mode you are using, it can return sevral different outcomes:
#' 1. For test_mode=F, it returns a list of 2 elements: the 'result' element itself is a list, containing output Ws and Hs under multiple parameters;
#'  the second element is the `consistency_table` for NITUMID's outcomes under different parameters, see our paper for more details, simply speaking, if the
#'  i-th value in `consistency_table` has the largest value, then you should extract 'result'[[i]] and use its W and H, e.g. result[[i]]$W
#' 2. For test_mode=T, it returns a list of 3 elements: on top of the 2 elements mentioned above, it has another element of `real_corr_table`, which is a
#' numeric vector contianing the correlations between NITUMID's outcome and real cells' fractions under different parameters. Again, you should choose the i-th value
#' (which you got from consistency_table) since that represents the final result.
#'
#' @examples TBA
#'
#' @export
NITUMID<-function(Y, A=A_melanoma_v5,row_index=seq(1:53), if.bulk,row_mean=row_mean_v5,tol=10^-5,num_cell=11,max.itr=2000,test_mode=F,real_H=NA,cell_structure=cell_structure_example,correlation="p"){


  Y <- as.matrix(Y)
  ##########################################no permutation case#############################################
  ##Iterate over beta=0.05,beta=100 and data-dependent beta
  Y <- Y[row_index,]
  Y_normal_row <- as.matrix(Y/row_mean[row_index])
  Y_normal_row_col <- sqrt( nrow(Y_normal_row) )*t(t(Y_normal_row)/sqrt(apply(Y_normal_row^2, 2,sum)))
  A <- A[row_index,]


##Loading some internal functions needed by NITUMID#######################################################
##in the test mode, if the labelled data has different cell types with us, we rearrange the results based on cell_structure
  transform_H <- function(input_H,cell_structure){
    dest_index <- sort(unique(cell_structure$destin_index))
    outcome_H <- matrix(rep(0,ncol(input_H)*max(dest_index)),ncol=ncol(input_H))
    for (i in 1:num_cell){
      if ( !is.na(cell_structure$destin_index[i]) ){
        outcome_H[cell_structure$destin_index[i],] <- outcome_H[cell_structure$destin_index[i],]+input_H[i,]
      }
    }
    return(outcome_H)
  }


#fxn for K-L divergence
kl_div_disc <- function(aa,bb){
  aa=aa/sum(aa); bb=bb/sum(bb); clength=length(aa); aa=pmax(aa,10^(-30)); bb=pmax(bb,10^(-30));
  return(sum(aa*log(aa/bb)))}

##################end of function##########################


##the case when if.bulk=TRUE
  if (if.bulk==T){
    consistency_table <- -1
    real_corr_table   <- -1
    #specify parameters
    lam_W =5000;
    beta_W = 0.05;
    eta=200;
    #non permutation
    deconvo_result_bulk <- NMF_RNA_Tang_ind(Y_normal_row_col,lam_W = lam_W, beta_W = beta_W, A=A,eta = eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data = 1,perm_type=0)

    ##starting permutation + noise case
    Y <- as.matrix(Y)
    ncolY=ncol(Y); rand_perm=sample(1:ncolY)
    Y2=as.matrix(Y[,rand_perm])
    #if (test_mode==T){
    #  real_H_permut <- real_H[,rand_perm]
    #}
    Y2_normal_row <- as.matrix(Y2/row_mean[row_index])
    Y2_normal_row_col <- as.matrix(sqrt( nrow(Y2_normal_row) )*t(t(Y2_normal_row)/sqrt(apply(Y2_normal_row^2, 2,sum))))
    deconvo_result_bulk_perm <- NMF_RNA_Tang_ind(Y2_normal_row_col,lam_W = lam_W, beta_W = beta_W, A=A, eta=eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data = 1)
    #calculate consistency
    consistency_table=mean( diag( cor(deconvo_result_bulk$H[,rand_perm],deconvo_result_bulk_perm$H)   ) )

    if (test_mode==T){

      if (correlation=="p"){
            browser()
            real_corr_table=mean(diag(cor(transform_H(deconvo_result_bulk$H, cell_structure),(real_H))))
      }
      if (correlation=="s"){
            real_corr_table=mean(diag(cor(transform_H(deconvo_result_bulk$H,cell_structure),(real_H),method = "spearman" ) ) )
      }

      if (correlation=="kl"){
            temp_kl <- rep(-1,ncol(Y))
            trans_H<-transform_H(deconvo_result_bulk$H, cell_structure)
            for (kl in 1:ncol(Y)){
              temp_kl[kl] <- kl_div_disc( aa = real_H[,kl], bb=trans_H[,kl])
            }
            real_corr_table=mean(temp_kl)
      }

    }

    report_consistency_table<-consistency_table

    if (test_mode==T){

      return(list(result=list(deconvo_result_bulk),consistency_table=report_consistency_table,real_corr_table=real_corr_table))
      }

    else {
      return(list(result=list(deconvo_result_bulk),consistency_table=report_consistency_table))

    }

  }

 #the case for pure cell gene expression samples
  else{
    consistency_table <- rep(-1,3)
    names(consistency_table) <- c("N1S1","N1S2","N2S1")
    real_corr_table   <- rep(-1,3)
    names(real_corr_table) <- c("N1S1","N1S2","N2S1")
    #specify parameters
    lam_W =5000;
    beta_W = 100;
    eta=200;
    #browser()
    deconvo_result_pure1 <- NMF_RNA_Tang_ind(Y,lam_W=lam_W, beta_W = beta_W, A=A, eta = eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data=1,perm_type=0)
    deconvo_result_pure2 <- NMF_RNA_Tang_ind(Y,lam_W=lam_W, beta_W = beta_W, A=A, eta = eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data=2,perm_type=0)
    deconvo_result_pure3 <- NMF_RNA_Tang_ind(Y_normal_row_col,lam_W = lam_W, beta_W = beta_W, A=A, eta= eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data=1,perm_type=0)


    result_non_permut_pure=list(deconvo_result_pure1,deconvo_result_pure2,deconvo_result_pure3)

    ##starting permutation + noise case
    Y <- as.matrix(Y)
    ncolY=ncol(Y); rand_perm=sample(1:ncolY)
    Y2=as.matrix(Y[,rand_perm])
    #if (test_mode==T){
    #  real_H_permut <- real_H[,rand_perm]
    #}
    Y2_normal_row <- as.matrix(Y2/row_mean[row_index])
    Y2_normal_row_col <- as.matrix(sqrt( nrow(Y2_normal_row) )*t(t(Y2_normal_row)/sqrt(apply(Y2_normal_row^2, 2,sum))))

    deconvo_result_pure1_perm <- NMF_RNA_Tang_ind(Y2,lam_W=lam_W, beta_W = beta_W, A=A, eta = eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data=1,perm_type=0)
    deconvo_result_pure2_perm <- NMF_RNA_Tang_ind(Y2,lam_W=lam_W, beta_W = beta_W, A=A, eta = eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data=2,perm_type=0)
    deconvo_result_pure3_perm <- NMF_RNA_Tang_ind(Y2_normal_row_col,lam_W = lam_W, beta_W = beta_W, A=A, eta= eta, tol=tol,max.itr=max.itr,num_cell=num_cell,scale_Data=1,perm_type=0)

    result_permut_pure=list(deconvo_result_pure1_perm,deconvo_result_pure2_perm,deconvo_result_pure3_perm)

    #calculate consistency
    for (ipm in 1:3){
    consistency_table[ipm]=mean( diag( cor(result_non_permut_pure[[ipm]]$H[,rand_perm],result_permut_pure[[ipm]]$H)   ) )
    }

    if (test_mode==T){


      if (correlation=="p"){
        for (ipm in 1:3){
        real_corr_table[ipm]=mean(diag(cor(transform_H(result_non_permut_pure[[ipm]]$H, cell_structure),(real_H))))
      }}
      if (correlation=="s"){
        for (ipm in 1:3){
        real_corr_table[ipm]=mean(diag(cor(transform_H(result_non_permut_pure[[ipm]]$H,cell_structure),(real_H),method = "spearman" ) ) )
      }}

      if (correlation=="kl"){
        for (ipm in 1:3){
        temp_kl <- rep(-1,ncol(Y))
        trans_H<-transform_H(result_non_permut_pure[[ipm]]$H, cell_structure)
        for (kl in 1:ncol(Y)){
          temp_kl[kl] <- kl_div_disc( aa = real_H[,kl], bb=trans_H[,kl])
        }
        real_corr_table[ipm]=mean(temp_kl)
      }}
     }

    report_consistency_table<-consistency_table

    if (test_mode==T){
      return(list(result=result_permut_pure,consistency_table=report_consistency_table,real_corr_table=real_corr_table))
    }

    else {
      return(list(result=result_non_permut_pure,consistency_table=report_consistency_table))

    }
  } }


#' Simple version of NITUMID
#'
#' `NITUMID_simple` This is a simplified version of NITUMID, which uses most of the default settings (generic signature genes, etc.),
#' you can only input a gene expression data frame whose rownames is gene symbol, and tells it if this data is from pure cells or bulk tumor (melanoma)
#'
#' @param Y  Y is the gene expression matrix that has #{signature genes} rows and #{samples} columns, Y's rowname must be gene symbols
#'
#' @param if.bulk Logical, if Y is bulk gene expression data. If Y is purified cell gene expression, please set it as FALSE.
#'
#' @return `NITUMID_simple` returns a list of two object, W matrix and H matrix; However, when the consistency table have multiple results
#' that share consistency, then it would return all the raw outcomes (consistency table and results) and you are supposed to select manually
#' @examples TBA
#'
#' @export
NITUMID_simple <- function(Y,if.bulk){
  match_outcome<-Signature_Match(gene_expression = Y,signature_genes = as.character(signature_marker_melanoma$gene_symbol) )
  NITUMID_out <- NITUMID(Y = Y[match_outcome$matched_index,],if.bulk = if.bulk,A = A_melanoma_v5,row_index = setdiff(seq(1,53),match_outcome$missing_row_index),row_mean=row_mean_v5)
  #browser()
  best_index <- which.max(NITUMID_out$consistency_table)

  consistency_round <- round(NITUMID_out$consistency_table,4)
  #browser()
  if (length(which(consistency_round==max(consistency_round)))>1){
    message("There are multiple 'best' consistency, please look at all results manually")
    return( composite_result = NITUMID_out )
  }
  return(list(W=NITUMID_out$result[[best_index]]$W,H=NITUMID_out$result[[best_index]]$H))
}


# #A function to test the outcome (purified cells)
# Test_Deconvo <- function(deconvo_result,A_ref,correct_H){
#
#   #test the structure concordance with A-ref
#   rot_ind2=array(0,c(nrow(deconvo_result$W),ncol(deconvo_result$W)))
#   for (ii in 1:nrow(deconvo_result$W)){
#     rot_ind2[ii,which(deconvo_result$W[ii,]==max(deconvo_result$W[ii,]))]<-1
#   }
#   #how diff
#   W_diff<-length(which(rot_ind2 != 1-A_ref))/(dim(A_ref)[1]*dim(1-A_ref)[2])
#
#   #how accurate H is
#   #by max
#   len <- ncol(deconvo_result$H)
#   indi <- 0
#   for (i in 1:len){
#     if (which(deconvo_result$H[,i]==max(deconvo_result$H[,i]))==which(correct_H[,i]==max(correct_H[,i]))   ){
#       indi <- indi +1
#     }
#   }
#   max_right_ratio <- indi/len
#
#   return(list(err=deconvo_result$err[deconvo_result$step],err_diff=deconvo_result$err_diff[deconvo_result$step],W=rot_ind2,W_diff_rate=W_diff,H_max_right=max_right_ratio))
#
# }




