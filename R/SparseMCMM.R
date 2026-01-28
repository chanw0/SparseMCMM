
SparseMCMM=function(Treatment,otu.com,outcome,n.split=10,
                    dirichlet.penalty=seq(0,1,0.1),
                    lm.penalty1=seq(0,1,0.1),lm.penalty2=seq(0,2,0.2),
                    low.bound1=NULL,up.bound1=NULL,low.bound2=NULL,up.bound2=NULL,
                    num.per=NULL,
                    parallel = TRUE,             # Enable or disable parallelization
                    ncores = future::availableCores() - 1, # Number of cores
                    seed = NULL )
{

  ### using a small value to replace 0
  pseudo=min(otu.com[otu.com>0])/2
  otu.com=t(apply(otu.com,1,function(x) {if(min(x)==0) return((x+pseudo)/sum(x+pseudo)) else return(x)}))

  ### Splitting strategy for post-selection estimation


  ### Configure parallel backend based on `parallel` parameter
  options(future.rng.onMisuse = "ignore")  # Automatically configures RNG for all parallel tasks
  if (parallel) {
    plan(multisession, workers = ncores)  # Use parallel backend with specified cores
    options(future.seed = seed)  # Set global RNG seed for reproducibility
  } else {
    plan(sequential)  # Disable parallelism; run computation sequentially
  }

  ### Outer Parallel Loop: Parallelize over splits
  results <- future_map(1:n.split, function(tt) {


    ### variable selection
    Training=c(sample(which(Treatment==1),ceiling(sum(Treatment==1)/2)),
               sample(which(Treatment==0),ceiling(sum(Treatment==0)/2)))


    Treatment.training=Treatment[Training]
    otu.com.training=otu.com[Training,]
    outcome.training=outcome[Training]


    # Estimate initial betas and alphas
    betas=beta.estimates(Treatment.training,otu.com.training,penalty.lambda=dirichlet.penalty,
                         low.bound=low.bound2,up.bound=up.bound2)

    beta.estimation=rep(0,length(betas))
    Index_beta=betas[-c(1:ncol(otu.com))]==0
    Index_0 <- if(sum(Index_beta) > 0) which(Index_beta) else NULL


    Alpha.Est0=alpha.estimates(Treatment.training,otu.com.training,outcome.training,
                               penalty.lambda1=lm.penalty1,penalty.lambda2=lm.penalty2,
                               low.bound=low.bound1,up.bound=up.bound1)

    Index_alpha=Alpha.Est0==0
    Index_alpha0 <- if(sum(Index_alpha) > 0) which(Index_alpha) else NULL
    alpha.estimation=rep(0,length(Alpha.Est0))

    ### post-selection estimation
    Treatment.est=Treatment[-Training]
    otu.com.est=otu.com[-Training,]
    outcome.est=outcome[-Training]


    est.beta=beta.estimates2(Treatment.est,otu.com.est,Index_0)
    if(!is.null(Index_0))  beta.estimation[-c(Index_0+ncol(otu.com))]=est.beta else  beta.estimation=est.beta


    est.alpha=alpha.estimates2(Treatment.est,otu.com.est,outcome.est,Index_alpha0)
    if(!is.null(Index_alpha0))  alpha.estimation[-Index_alpha0]=est.alpha else  alpha.estimation=est.alpha

    ### mediation effect
    CausalEffect=CausalE(otu.com,alpha.estimation=alpha.estimation,
                         beta.estimation=beta.estimation)


    # Base results for this split
    All.per <- list(c(tt, 0, CausalEffect[[1]]))
    Individual.per <- list(c(tt, 0, CausalEffect[[2]]))

    if(!is.null(num.per)) {

      # Parallelize over permutations within the split
      perm_results <- future_map(1:num.per, function(rr) {

        Treatment.per=sample(Treatment,length(Treatment))
        outcome.per=sample(outcome, length(outcome))

        Treatment.training=Treatment.per[Training]
        otu.com.training=otu.com[Training,]
        outcome.training=outcome.per[Training]

        betas=beta.estimates(Treatment.training,otu.com.training,penalty.lambda=dirichlet.penalty,
                             low.bound=low.bound2,up.bound=up.bound2)

        beta.estimation=rep(0,length(betas))
        Index_beta=betas[-c(1:ncol(otu.com))]==0
        Index_0 <- if(sum(Index_beta) > 0) which(Index_beta) else NULL

        Alpha.Est0=alpha.estimates(Treatment.training,otu.com.training,outcome.training,
                                   penalty.lambda1=lm.penalty1,penalty.lambda2=lm.penalty2,
                                   low.bound=low.bound1,up.bound=up.bound1)

        Index_alpha=Alpha.Est0==0
        Index_alpha0 <- if(sum(Index_alpha) > 0) which(Index_alpha) else NULL
        alpha.estimation=rep(0,length(Alpha.Est0))

        ### post-selection estimation
        Treatment.est=Treatment.per[-Training]
        otu.com.est=otu.com[-Training,]
        outcome.est=outcome.per[-Training]


        est.beta=beta.estimates2(Treatment.est,otu.com.est,Index_0)
        if(!is.null(Index_0))  beta.estimation[-c(Index_0+ncol(otu.com))]=est.beta else  beta.estimation=est.beta


        est.alpha=alpha.estimates2(Treatment.est,otu.com.est,outcome.est,Index_alpha0)
        if(!is.null(Index_alpha0))  alpha.estimation[-Index_alpha0]=est.alpha else  alpha.estimation=est.alpha


        ### mediation effect
        CausalEffect=CausalE(otu.com,alpha.estimation=alpha.estimation,
                             beta.estimation=beta.estimation)

        perm_res <- list(c(tt, rr, CausalEffect[[1]]), c(tt, rr, CausalEffect[[2]]))
        return(perm_res)

      })

      # Combine permutation results
      All.per <- c(All.per, lapply(perm_results, "[[", 1))
      Individual.per <- c(Individual.per, lapply(perm_results, "[[", 2))
    }

    # Return results for this split
    return(list(All.per = do.call(rbind, All.per),
                Individual.per = do.call(rbind, Individual.per)))
  })


  # Combine results across splits
  All.per <- do.call(rbind, lapply(results, `[[`, "All.per"))
  Individual.per <- do.call(rbind, lapply(results, `[[`, "Individual.per"))

  ### Effects

  colnames(All.per)=c("split","permu",colnames(All.per)[-c(1:2)])
  All.per=data.frame(All.per)

  if(n.split==1) MEs=unlist(All.per[All.per$permu==0,-c(1:2)]) else {

    MEs=rbind(colMeans(All.per[All.per$permu==0,-c(1:2)]),
              apply(All.per[All.per$permu==0,-c(1:2)], 2, sd)/sqrt(sum(All.per$permu==0)))
    rownames(MEs)=c("mean","se") }

  colnames(Individual.per)=c("split","permu",colnames(Individual.per)[-c(1:2)])
  Individual.per=data.frame(Individual.per)

  if(n.split==1) Com_MEs=unlist(Individual.per[Individual.per$permu==0,-c(1:2)]) else {

    Com_MEs_mean=colMeans(Individual.per[Individual.per$permu==0,-c(1:2)])
    Com_MEs_se=apply(Individual.per[Individual.per$permu==0,-c(1:2)], 2, sd)/sqrt(sum(Individual.per$permu==0))

    Com_MEs= rbind(Com_MEs_mean,
                   Com_MEs_mean-Com_MEs_se*1.96,
                   Com_MEs_mean+Com_MEs_se*1.96)

    rownames(Com_MEs)=c("Mean", "Lower_CI","Upper_CI")
  }


  Effect.estimates=list(MEs,Com_MEs)
  names(Effect.estimates)=c('Esitmated Causal Effects', 'Compontent-wise ME')

  if(!is.null(num.per)) {

    ### Testing
    OMD <- All.per %>%
      group_by(permu) %>%  # Group by two variables
      summarise(MeanME = mean(ME, na.rm = TRUE)) %>%  # Calculate mean
      select(-permu) %>% unlist()  %>% as.numeric() %>% na.omit()


    if (n.split > 1) {
      CMD <- Individual.per %>% as.data.frame() %>%
        group_by(permu) %>%  # Group by unique permutation values
        summarise(across(-c(split), ~ mean(.x, na.rm = TRUE))) %>%  # Calculate mean for all columns except split/permu
        select(-permu) %>%
        as.matrix()         # Convert to matrix if needed
    } else {
      CMD <- Individual.per[, -c(1:2)]
    }

    CMD=rowSums(CMD^2) %>% na.omit()

    p.values=c(sum(abs(OMD)>=abs(OMD[1]))/length(OMD),
               sum(CMD>=CMD[1])/length(CMD))

    names(p.values)=c("OME","CME")
    Effect.estimates$Test= p.values
  }


  return(Effect.estimates)

}


