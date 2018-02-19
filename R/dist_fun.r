# Function to fit discrete heavy tail functions
# Mostly using poweRlaw package
#
fit_dis_heavy_tail <- function (data_set,xmins,options.output)
{
  require(fitdistrplus)
  
  n <- length(unique(data_set))
  n_models <- 6
  # List of models
  model_list = list(list(model=vector("list", length=0), 
                         GOF=vector("list", length=0),
                         xmin_estimation=vector("list", length=0),
                         uncert_estimation=vector("list", length=0),
                         k=0,
                         LL=0,
                         n=n,
                         AICc=0,
                         delta_AICc=0,
                         AICc_weight=0,
                         model_name=character(0),
                         model_set=character(0))) 
  fit_ht <- rep(model_list,n_models)
  dim(fit_ht) <- c(n_models)

  # Declare models
  #
  # Discrete power law
  fit_ht[[1]]$model <- displ$new(data_set)
  fit_ht[[1]]$k <- 1
  # Discrete log-normal
  fit_ht[[2]]$model <- dislnorm$new(data_set)
  fit_ht[[2]]$k <- 2
  # Discrete exponential
  fit_ht[[3]]$model <- disexp$new(data_set)
  fit_ht[[3]]$k <- 1
  # Poisson 
  fit_ht[[4]]$model <- dispois$new(data_set)
  fit_ht[[4]]$k <- 1
  # power law with exponential cutoff
  fit_ht[[5]]$model <- "" 
  fit_ht[[5]]$k <- 2
  # Uniform 
  fit_ht[[6]]$model <- "" 
  fit_ht[[6]]$k <- 1
  
  # Estimate Xmin with complete data_set for power law model
  #
  if(xmins>1) {
    fit_ht[[1]]$xmin_estimation <- estimate_xmin(fit_ht[[1]]$model,
                                                 xmins = 1:xmins, 
                                                 pars = NULL, 
                                                 xmax = max(data_set))
    fit_ht[[1]]$model$setXmin(fit_ht[[1]]$xmin_estimation)
  } else {
    fit_ht[[1]]$model$setXmin(1)
    
  }
  
  
  
  model_names <- c("Power", "LogNorm","Exp","Poisson","PowerExp","Uniform")
  
  AICc_weight <- matrix( nrow = n_models, ncol = 1,
                         dimnames = list(model_names))
  
  delta_AICc <- AICc_weight
  GOF <- delta_AICc
  aic_min=Inf
  norm_aic_weight=0

  for (i in 1:(n_models))
  {
      # Set cut-off (x_min)
      fit_ht[[i]]$model$xmin <- fit_ht[[1]]$model$xmin
      
      # Correct n with xmin
      fit_ht[[i]]$n <-length(data_set[data_set>=fit_ht[[1]]$model$xmin])
      
      # Fit models
      #            PowerExp is different!!!!!!!!!
      #
      if(i<5){
        fit_ht[[i]]$model$setPars(estimate_pars(fit_ht[[i]]$model))
        # 
        # Get Loglikelihood
        fit_ht[[i]]$LL <- dist_ll(fit_ht[[i]]$model)
      } else if(i==5) {
        # Fit Power Exponential (not poweRlaw package)
        #
        fit_ht[[i]]$model<- discpowerexp.fit(data_set,fit_ht[[1]]$model$xmin)
        # 
        # Get Loglikelihood
        fit_ht[[i]]$LL <- fit_ht[[i]]$model$loglike 			
      } else if(i==6) {
        # Fit Power Exponential (not poweRlaw package)
        #
        fit_ht[[i]]$model<- max(data_set)
        teta<- max(data_set)
        # 
        # Get Loglikelihood
        fit_ht[[i]]$LL <- log(teta^(-length(data_set)))			
      }
      
      # 
      # Get Loglikelihood
      LL <- fit_ht[[i]]$LL
      k <- fit_ht[[i]]$k

      # Compute AICc
      #
      fit_ht[[i]]$AICc <- (2*k-2*LL)+2*k*(k+1)/(n-k-1)
      aic_min <- min(aic_min,fit_ht[[i]]$AICc)
      fit_ht[[i]]$model_name <- model_names[i]
  }
    
  for (i in 1:n_models){
    delta_AICc[i] <- fit_ht[[i]]$AICc - aic_min
    fit_ht[[i]]$delta_AICc <- delta_AICc[i]
    norm_aic_weight <- norm_aic_weight + exp(-0.5*fit_ht[[i]]$delta_AICc)
  }

  # Akaike weigths and dataframe with parameters
  #
  daf<-data.frame()
  
  for (i in 1:n_models){
    AICc_weight[i]  <- exp(-0.5*fit_ht[[i]]$delta_AICc)/norm_aic_weight 
    fit_ht[[i]]$AICc_weight <- AICc_weight[i]
    if(i<5){
      daf <- rbind(daf,data.frame(ModelNames=model_names[i],
                                  par1=fit_ht[[i]]$model$pars[1],
                                  par2=fit_ht[[i]]$model$pars[2],
                                  xmin=fit_ht[[1]]$model$getXmin(),
                                  n=fit_ht[[i]]$n,
                                  AICc=fit_ht[[i]]$AICc,Delta_AICc=delta_AICc[i],AICc_weight=AICc_weight[i]))
    } else if(i==5){
      daf <- rbind(daf,data.frame(ModelNames=model_names[i],
                                  par1=-fit_ht[[i]]$model$exponent,
                                  par2=fit_ht[[i]]$model$rate,
                                  xmin=fit_ht[[1]]$model$getXmin(),
                                  n=fit_ht[[i]]$n,
                                  AICc=fit_ht[[i]]$AICc,Delta_AICc=delta_AICc[i],AICc_weight=AICc_weight[i]))
    }else if(i==6){
      daf <- rbind(daf,data.frame(ModelNames=model_names[i],
                                  par1=fit_ht[[i]]$model,
                                  par2= NA,
                                  xmin=fit_ht[[1]]$model$getXmin(),
                                  n=fit_ht[[i]]$n,
                                  AICc=fit_ht[[i]]$AICc,Delta_AICc=delta_AICc[i],AICc_weight=AICc_weight[i]))
      
    }
  }

  # Plots
  #
  if (options.output$ploting){
    #setwd(options.output$resultsDir)
    
    require(RColorBrewer)
    colp <-brewer.pal(8,"Dark2")
    fnam <-paste0(options.output$data_set_name, "_xmin",fit_ht[[1]]$model$xmin,".png")
    png(filename=fnam, res=300,units = "mm", height=200, width=200,bg="white")
    po <-plot(fit_ht[[1]]$model,xlab="Degree",ylab="log[P(X > x)]",main=options.output$data_set_name)
    for (i in 1:n_models){
      if(i<5) {
        lines(fit_ht[[i]]$model, col=colp[i])
      } else if(i==5) {
        est1 <- fit_ht[[i]]$model
        est1$xmin <- fit_ht[[1]]$model$xmin
        x <- sort(unique(data_set))
        x <- x[x>=est1$xmin]
        shift <- max(po[po$x>=est1$xmin,]$y)
        y <- ppowerexp(x,est1$xmin,est1$exponent,est1$rate,lower.tail=F)*shift
        lines(x,y,col=colp[i])
      } else if(i==6) {
        est1 <- fit_ht[[i]]$model
        xmin <- fit_ht[[1]]$model$xmin
        x <- sort(unique(data_set))
        x <- x[x>=xmin]
        shift <- max(po[po$x>=xmin,]$y)
        y <- 1-((x - xmin+1)/(est1-xmin+1)*shift)
        lines(x,y,col=colp[i])
      }
  }
    
      
    legend("topright",model_names,bty="n",col=colp,lty=c(1,1,1,1),cex=1)
    dev.off()
  }

  return(list(fitted_models=fit_ht,da_fit=daf))  

}

fit_ht_dplyr_helper <-function(df,xmin=1){
  opt.output$data_set_name <- unique(df$Network)
  temp <-fit_dis_heavy_tail(df$Degree,xmin,opt.output)
  temp <-data.frame(temp$da_fit)
  temp$Network <- unique(df$Network)
  return(temp)
}


# Plot of frequencies of patch sizes with fitted continuous heavy tail functions
# x: data
# fit_ht_df : dataframe with fitted parameters
#
freq_plot_con_ht <- function(x,fit_ht,tit="") # PLOT ALL FROM X=1 ?????????????????
{
  require(dplyr)
  require(ggplot2)
  xx <-as.data.frame(table(x))
  xx$x <- as.numeric(as.character(xx$x))
  ff <- filter(fit_ht,ModelNames=="PowerExp")
  if(nrow(ff)>0){
    xmin <- ff$xmin
    xx$pexp <-dpowerexp(xx$x,1,ff$par1,ff$par2)
  } else {
    xx$pexp<-0
  }
  ff <- filter(fit_ht,ModelNames=="Power")
  xmin <- ff$xmin
  xx$pow  <-dpareto(xx$x,xmin,ff$par1)
  
  ff <- filter(fit_ht,ModelNames=="Exp")
  xx$exp  <-dexp(xx$x,ff$par1)	
  
  xx <-mutate(xx,pexp=ifelse(x<xmin,NA,pexp),pow=ifelse(x<xmin,NA,pow),
              exp=ifelse(x<xmin,NA,exp),
              Freq=Freq/sum(Freq))
  
  minFreq <- min(xx$Freq) - 0.5*min(xx$Freq)
  minFreq <- ifelse(minFreq<0,0,minFreq)
  g <- ggplot(xx, aes(y=Freq,x=x)) +  theme_bw() + geom_point(alpha=0.3) + 
    coord_cartesian(ylim=c(1,min(minFreq)))+
    scale_y_log10() +scale_x_log10() + ylab("Frequency") + xlab("Degree") +ggtitle(tit)
  
  mc <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  g <- g + geom_line(aes(y=pow,x=x,colour="P.law")) + 
#    geom_line(aes(y=pexp,x=x,colour="P.law with\nexp. cutoff"))+
    geom_line(aes(y=exp,x=x,colour="Exp."))+
    scale_colour_manual(values=mc,name="")  
  
  fil <- gsub(" ", "", tit, fixed = TRUE)
  fil <- paste0(fil,".png")
  if(tit=="")
    print(g)
  else
    ggsave(fil,plot=g,width=6,height=4,units="in",dpi=600)
  
}

# Calculates The Corrected Akaike Criterion
# k: number of parameters
# n: number of points
# mdl1,2: fitted nonlinear models
# 
Akaike_criterion<-function(k,n,mdl1,mdl2,names)
{
  Akaike <-data.frame(Model=names, AICc=0)
  Akaike$AICc[1] <- 2*k - 2* logLik(mdl1) + 2*k*(k+1)/(n-k-1)
  Akaike$AICc[2] <- 2*k - 2* logLik(mdl2) + 2*k*(k+1)/(n-k-1)
  Akaike$Delta <- Akaike$AICc - min(Akaike$AICc)
  return(Akaike)
}



# Function to plot CCDF of fitted continuous heavy tail functions using ggplo2
#
# x: data
# fit_ht : dataframe with fitted parameters
# tit: file name to save graph
# fit_ht1: second set of parameters to superimpose in the same graph

cdfplot_displ_exp <- function(x,fit_ht,tit="",xmax=0)
{
  require(poweRlaw)
  require(viridis)
  m <- displ$new(x)
  tP <- plot(m,draw=F)
  require(ggplot2)
  require(dplyr)
  
  ff <- filter(fit_ht,ModelNames=="Power")
  xmin <- ff$xmin
  
  tP1 <- cdfplot_conpl_exp_helper(x,tP,fit_ht,xmin)
  
  #tP2 <-filter(tP2, powl>= min(tP$Rank))
  #tP1 <-filter(tP1, powl>= min(tP$Rank))
  
  #mc <- c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2","#D55E00", "#CC79A7")
  # Brewer
  #mc <- c("#d7191c","#fdae61","#abd9e9","#2c7bb6")
  
  g <- ggplot(tP, aes(x=x,y=y)) +  theme_bw() + geom_point(alpha=0.3) + 
    coord_cartesian(ylim=c(1,min(tP$y)))+
    scale_y_log10() +scale_x_log10() + ylab("log[P(X > x)]") + xlab("Degree") #+ggtitle(tit)
  if(xmax>0) {
    g<-g + xlim(0,xmax+1) + scale_x_log10()
  }
  brk<-unique(tP1$model)
  g <- g + geom_line(data=tP1,aes(y=powl,x=psize,colour=model)) + 
    #scale_colour_discrete(name="",breaks=brk)
    #scale_colour_manual(values=mc,name="",breaks=brk)  
    #scale_colour_brewer(type="dark2",palette=7,name="",breaks=brk)
    scale_colour_viridis(discrete = TRUE)
  
  fil <- gsub(" ", "", tit, fixed = TRUE)
  fil <- paste0(fil,".png")
  if(tit=="")
    print(g)
  else
    ggsave(fil,plot=g,width=6,height=4,units="in",dpi=600)
  
}

cdfplot_conpl_exp_helper <- function(x,tP,fit_ht,xmin,mode="gt")
{
  x1 <- unique(x)
  
  # Select model and generate a data frame 
  #
  ff <- filter(fit_ht,ModelNames=="PowerExp")
  
  if(mode=="gt") {
    x1 <- x1[x1>=xmin]
    shift <- max(filter(tP,x>=xmin)$y)
  } else {
    x1 <- x1[x1<xmin]
    #xmin<-1
    #shift <-1
    # First select lower subset the change to the Xmin of this subset
    xmin <- ff$xmin
    shift <- max(filter(tP,x>=xmin)$y)
  }
  
  tP2 <- data_frame(psize=x1, powl=ppowerexp(x1,xmin,ff$par1,ff$par2,lower.tail=F)*shift,model="PowerExp")
  
  ff <- filter(fit_ht,ModelNames=="Power")
  m <- displ$new(x)
  m$setPars(ff$par1)
  m$setXmin(xmin)
  tP1 <- data_frame(psize=x1,powl=dist_cdf(m,x1,lower_tail=F)*shift,model="Power")
  
  ff <- filter(fit_ht,ModelNames=="Exp")
  m <- disexp$new(x)
  m$setPars(ff$par1)
  m$setXmin(xmin)
  tP3 <- data_frame(psize=x1,powl=dist_cdf(m,x1,lower_tail=F)*shift,model="Exp")
  
  ff <- filter(fit_ht,ModelNames=="LogNorm")
  m <- dislnorm$new(x)
  m$setPars(c(ff$par1,ff$par2))
  m$setXmin(xmin)
  tP4 <- data_frame(psize=x1,powl=dist_cdf(m,x1,lower_tail=F)*shift,model="LogNorm")

  ff <- filter(fit_ht,ModelNames=="Poisson")
  m <- dispois$new(x)
  m$setPars(c(ff$par1,ff$par2))
  m$setXmin(xmin)
  tP5 <- data_frame(psize=x1,powl=dist_cdf(m,x1,lower_tail=F)*shift,model="Poisson")
  tP1 <- bind_rows(tP1,tP2,tP3,tP4,tP5)
}


# Complementary cumulative distribution plot using base graphics and poweRlaw package
#
ccdf_base_plot_ht <- function(fit_ht,data_set,netName=""){
    # Plots
  #
  
    require(RColorBrewer)
    colp <-brewer.pal(8,"Dark2")
    if(netName!="") {
        fnam <-paste0(netName, "_xmin",fit_ht[1]$xmin,".png")
        png(filename=fnam, res=300,units = "mm", height=200, width=200,bg="white")
    }
    
    po <-plot(fit_ht[1]$model,xlab="Degree",ylab="log[P(X > x)]",main=options.output$data_set_name)
    for (i in 1:n_models){
      if(i!=n_models) {
        lines(fit_ht[[i]]$model, col=colp[i])
      } else {
        est1 <- fit_ht[[i]]$model
        x <- sort(unique(data_set))
        x <- x[x>=est1$xmin]
        shift <- max(po[po$x>=est1$xmin,]$y)
        y <- ppowerexp(x,est1$xmin,est1$exponent,est1$rate,lower.tail=F)*shift
        lines(x,y,col=colp[i])
      }
    }
    
      
    legend("topright",model_names,bty="n",col=colp,lty=c(1,1,1,1),cex=1)
    if(netName!="") dev.off()
}