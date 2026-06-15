#' Extract Canonical Frequency Bands
#'
#' Applies a linear bandpass filter to extract one or more canonical EEG
#' frequency bands from a univariate time series. Available bands are:
#' \itemize{
#'   \item \code{"delta"}: 0.5--4 Hz
#'   \item \code{"theta"}: 4--8 Hz
#'   \item \code{"alpha"}: 8--12 Hz
#'   \item \code{"beta"}: 12--30 Hz
#'   \item \code{"gamma"}: 30--50 Hz
#'   \item \code{"dta"}: 0.5--12 Hz
#'   \item \code{"bg"}: 12--50 Hz
#'   \item \code{"dtabg"}: 0.5--45 Hz
#' }
#'
#' @param x A numeric vector containing the time series to be filtered.
#' @param samprate A positive numeric scalar giving the sampling frequency (Hz).
#' @param fband A character vector specifying the frequency bands to extract.
#'   Possible values are \code{"delta"}, \code{"theta"}, \code{"alpha"},
#'   \code{"beta"}, \code{"gamma"}, \code{"dta"}, \code{"bg"}, and
#'   \code{"dtabg"}.
#' @param method A character scalar specifying the filtering method:
#'   \code{"butter"} for a Butterworth filter or \code{"fir"} for a finite
#'   impulse response filter.
#' @param twosided Logical scalar indicating whether a two-sided
#'   zero-phase filter (\code{TRUE}) or a one-sided causal filter
#'   (\code{FALSE}) should be used.
#'
#' @return A named list whose elements contain the filtered time series
#'   corresponding to the frequency bands specified in \code{fband}.
#' @examples
#' # Extracting frequency bands from a simulated time series
#' series <- rnorm(1000)
#' series.fband <- extract.fband(x = series,samprate = 100,fband = c("delta","theta"))
#' plot.ts(series.fband$delta)
#' plot.ts(series.fband$theta)
#' @export
extract.fband <- function(x,samprate,fband = c("delta","theta","alpha","beta","gamma","dta","bg","dtabg"),
                            method = "butter",twosided = TRUE){

  y <- vector("list",length(fband))
  names(y) <- fband

  for(i in 1:length(fband)){

    if(fband[i] == "delta"){

      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(0.5,4)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(0.5,4)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "theta"){

      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(4,8)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(4,8)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "alpha"){

      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(8,12)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(8,12)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "beta"){

      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(12,30)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(12,30)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "gamma"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(30,45)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(30,45)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "dta"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(0.5,12)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(0.5,12)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "bg"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(12,45)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(12,45)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    } else if(fband[i] == "dtabg"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(0.5,45)/(0.5 * samprate), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(0.5,45)/(0.5 * samprate), type="pass")
      }

      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }

    }

  }

  return(y)

}

#' Generate a Step-Wise Modulating Signal
#'
#' Generates a time-varying step-wise modulating signal from a sequence of
#' peak magnitudes and the intervals between successive peaks. Transitions
#' between adjacent peaks are smoothed according to the parameters
#' \code{adj} and \code{beta}.
#'
#' @param time A numeric vector specifying the time intervals between
#'   successive peaks.
#' @param peaks A numeric vector containing the magnitudes of the peaks.
#' @param adj A positive numeric scalar controlling the width of the
#'   transition region between adjacent peaks.
#' @param beta A positive numeric scalar controlling the abruptness of
#'   transitions between adjacent peaks.
#'
#' @return A numeric vector containing the generated modulating signal.
#'
#' @keywords internal
aux_modulating.signal <- function(time,peaks,adj = 6,beta = 0.75){

  adj.time <- time - (2*adj+1)

  y <- c()
  for(t in 1:(length(time)-1)){

    if(peaks[t] < peaks[t+1]){
      add <- (peaks[t+1]-peaks[t])*(1 - 1/(1+exp(beta*(-adj:adj)))) + peaks[t]
    } else {
      add <- -(peaks[t]-peaks[t+1])*(1 - 1/(1+exp(beta*(-adj:adj)))) + peaks[t]
    }

    y <- c(y,rep(peaks[t],adj.time[t]),add)

  }

  y <- c(y,rep(peaks[length(time)],adj.time[length(time)]))
  return(y)
}


#' Generate Band-Specific Oscillations (Carrier Signals)
#'
#' Simulates band-specific carrier signals using autoregressive processes of
#' order two (AR(2)) and subsequently applies a linear bandpass filter to
#' isolate the desired frequency bands. Available bands are:
#' \itemize{
#'   \item \code{"delta"}: 0.5--4 Hz
#'   \item \code{"theta"}: 4--8 Hz
#'   \item \code{"alpha"}: 8--12 Hz
#'   \item \code{"beta"}: 12--30 Hz
#'   \item \code{"gamma"}: 30--45 Hz
#' }
#'
#' @param nsec A positive numeric scalar specifying the duration of the
#'   simulated signal in seconds.
#' @param samprate A positive numeric scalar giving the sampling frequency
#'   (Hz).
#' @param fband A character vector specifying the frequency bands to
#'   generate. Possible values are \code{"delta"}, \code{"theta"},
#'   \code{"alpha"}, \code{"beta"}, and \code{"gamma"}.
#' @param minlength A positive numeric scalar specifying the minimum segment
#'   length used during amplitude normalization. Within each segment, the
#'   minimum (trough) and maximum (peak) values are rescaled to -1 and 1,
#'   respectively. If \code{NA}, no amplitude normalization is performed.
#'
#' @return A numeric matrix containing the simulated band-specific
#'   oscillations. Rows correspond to time points and columns correspond to
#'   the frequency bands specified in \code{fband}.
#' @examples
#' # Simulating frequency bands
#' set.seed(1)
#' series.fband <- sim.fband(nsec = 10,samprate = 100,fband = c("delta","theta"))
#' plot.ts(series.fband[,1])
#' plot.ts(series.fband[,2])
#' @export
sim.fband <- function(nsec,samprate,fband,minlength = NA){

  Tlength <- nsec*samprate
  Z <- array(NA,dim = c(Tlength,length(fband)))

  for(f in 1:length(fband)){

    if(fband[f] == "delta"){

      peaklocation <- 2
      sharpness <- 0.03

    } else if(fband[f] == "theta"){

      peaklocation <- 6
      sharpness <- 0.03

    } else if(fband[f] == "alpha"){

      peaklocation <- 10
      sharpness <- 0.03

    } else if(fband[f] == "beta"){

      peaklocation <- 22.5
      sharpness <- 0.05

    } else if(fband[f] == "gamma"){

      peaklocation <- 37.5
      sharpness <- 0.1

    }

    psi<- peaklocation/samprate
    par <- c(exp(-sharpness)*2*cos(2*pi*psi), -exp(-2*sharpness))

    Z_f <- scale(extract.fband(stats::arima.sim(model = list(ar = par),n = Tlength),
                                 samprate = samprate,fband = fband[f])[[1]])

    if(!is.na(minlength)){
      Z_pos <- diff(Z_f > 0)
      Z_ind0 <- (1:Tlength)[Z_pos==1]

      Z_ind <- c()
      for(s in 1:length(Z_ind0)){
        if(Z_ind0[s] > (length(Z_ind)+1)*minlength){
          Z_ind <- c(Z_ind,Z_ind0[s])
        }
      }
      Z_ind <- unique(c(0,Z_ind,Tlength))

      Z_f2 <- c()
      for(t in 1:(length(Z_ind)-1)){
        z <- Z_f[(Z_ind[t]+1):(Z_ind[t+1])]
        z_max <- max(abs(z))
        Z_f2 <- c(Z_f2,(z/z_max))
      }

      Z[,f] <- Z_f2[1:Tlength]

    } else {

      Z[,f] <- Z_f[1:Tlength]

    }

  }

  return(Z)

}



#' Generate Causal Band-Specific Oscillations
#'
#' Simulates band-specific oscillations with known causal interactions.
#' Dependence among signals is introduced through a vector autoregressive
#' (VAR) process of order one used to generate the peak magnitudes of the modulating
#' signals.
#'
#' Available frequency bands are:
#' \itemize{
#'   \item \code{"delta"}: 0.5--4 Hz
#'   \item \code{"theta"}: 4--8 Hz
#'   \item \code{"alpha"}: 8--12 Hz
#'   \item \code{"beta"}: 12--30 Hz
#'   \item \code{"gamma"}: 30--45 Hz
#' }
#'
#' @param nsec A positive numeric scalar specifying the duration of the
#'   simulated signal in seconds.
#' @param samprate A positive numeric scalar giving the sampling frequency
#'   in Hz.
#' @param fband A character vector specifying the frequency bands to
#'   generate. Possible values are \code{"delta"}, \code{"theta"},
#'   \code{"alpha"}, \code{"beta"}, and \code{"gamma"}.
#' @param phi A coefficient matrix of the VAR(1) process used to generate the
#'   peak magnitudes of the modulating signals.
#' @param sigma The innovation variance-covariance matrix of the VAR(1) process.
#' @param range A numeric vector of length two specifying the minimum and
#'   maximum segment lengths associated with successive peaks of the
#'   modulating signals (in Hz).
#' @param peak.scale A numeric vector of length two specifying the minimum
#'   and maximum magnitudes used to rescale the modulating signals.
#' @param prob.mid A numeric scalar between 0 and 1 specifying the
#'   probability mass assigned to \code{mean(range)} when simulating the
#'   time intervals between successive peaks. Probabilities decrease away
#'   from the center so that the total probability sums to one.
#' @param nsec.burnin A positive numeric scalar specifying the duration
#'   in seconds of the burn-in period to discard.
#'
#' @return A numeric matrix containing the simulated causal band-specific
#'   oscillations. Rows correspond to time points and columns correspond to
#'   the frequency bands specified in \code{fband}.
#' @examples
#' # Simulating causal frequency bands (theta--gamma coupling)
#' set.seed(1)
#' phi_tg <- matrix(c(-0.1,    0,    0,    0,
#'                    2.5,  0.1,    0,    0,
#'                    0,    0, -0.1,    0,
#'                    2.5,    0,  2.5,  0.1), 4, byrow = TRUE)
#' Z <- causal.fband(nsec = 10,samprate = 100,fband = c("theta","theta","gamma","gamma"),
#'                   phi = phi_tg,sigma = diag(rep(1,4)),range = c(20,30))
#' plot.ts(Z[1:200,1])
#' plot.ts(Z[1:200,2])
#' plot.ts(Z[1:200,3])
#' plot.ts(Z[1:200,4])
#' @export
causal.fband <- function(nsec,samprate=100,fband,
                         phi,sigma = diag(rep(1,length(fband))),
                         range = c(45,55),peak.scale = c(1,10),prob.mid = 0.6,
                         nsec.burnin = 5){

  Tlength <- nsec*samprate
  Z_t <- sim.fband(nsec = nsec, samprate = samprate,fband = fband,minlength = mean(range))

  peaks_xy <- tsDyn::VAR.sim(B=phi, n=ceiling(Tlength/range[1])+ceiling(1.1*nsec.burnin*samprate),
                      include="none", varcov=sigma)[-c(1:(nsec.burnin*samprate)),]

  if(all(!is.na(peak.scale))){
    for(f in 1:length(fband)){
      peaks_xy[,f] <- (peak.scale[2]-peak.scale[1])*(peaks_xy[,f] - min(peaks_xy[,f]))/(max(peaks_xy[,f]) - min(peaks_xy[,f])) + peak.scale[1]
    }
  }

  Y_t <- array(NA,dim = c(Tlength,length(fband)))

  if(length(unique(range))==2){
    peaks_lim <- (range[1]):(range[2])
    peaks_side <- (length(peaks_lim) - 1)/2
    mid <- ceiling((sum(1:peaks_side)*2)/(1 - prob.mid)) - sum(1:peaks_side)*2
    probs <- c(1:peaks_side,mid,peaks_side:1)
    probs <- probs/sum(probs)
  }

  for(f in 1:length(fband)){
    if(length(unique(range))==2){
      t_peaks <- sample(x = peaks_lim,prob = probs,size = length(peaks_xy[,f]),replace = TRUE)
    } else {
      t_peaks <- rep(mean(range),length(peaks_xy[,f]))
    }
    Y_t[,f] <- scale(aux_modulating.signal(time = t_peaks,peaks = peaks_xy[,f])[1:Tlength]*Z_t[,f])
  }

  return(Y_t)
}




#' Vine Copula-Based Inference for Transfer Entropy
#'
#' Estimates transfer entropy between two time series, \code{x} and \code{y},
#' in both directions: from \code{x} to \code{y} and from \code{y} to \code{x}.
#' Estimation is based on the vine copula representation of Redondo et al.
#' (2025).
#'
#' @param x A numeric vector containing the first time series.
#' @param y A numeric vector containing the second time series.
#' @param lx A positive integer specifying the number of lags of \code{x}
#'   used to estimate transfer entropy.
#' @param ly A positive integer specifying the number of lags of \code{y}
#'   used to estimate transfer entropy.
#' @param nboot A non-negative integer specifying the number of bootstrap
#'   replicates used to compute p-values. The bootstrap is based on the
#'   null-copula resampling scheme introduced in Redondo et al. (2025).
#'   If \code{nboot = 0}, the bootstrap procedure is skipped and only
#'   point estimates are returned.
#' @param cop.fam A character vector specifying the copula families used
#'   to estimate the vine copula models. Options follow those supported by
#'   the \pkg{rvinecopulib} package.
#' @param margins A character scalar specifying the method used to estimate
#'   the marginal distributions. Possible values are \code{"GEV"}, for a
#'   generalized extreme value distribution, and \code{"ECDF"}, for the
#'   empirical distribution function.
#' @param nMC A positive integer specifying the number of Monte Carlo samples
#'   used to approximate the transfer entropy.
#'
#' @return A list containing transfer entropy estimates in both directions
#'   and, if \code{nboot > 0}, the corresponding bootstrap p-values.
#' @examples
#' # Registering number of cores to be used for bootstrap (only run if nboot>0)
#' # no_cores <- parallel::detectCores() - 1
#' # cl <- parallel::makeCluster(no_cores)
#' # doParallel::registerDoParallel(cl)
#'
#' set.seed(1)
#' phi_s <- matrix(c(-0.1, 0.6,0,0.1), 2, byrow = TRUE)
#' Z <- tsDyn::VAR.sim(B = phi_s,n = 500,lag = 1,include = "none")
#' copTE(x = Z[,1],y = Z[,2],lx = 1,ly = 1,nboot = 0)
#' @export
copTE <- function(x,y,lx = 2,ly = 2,nboot = 0,cop.fam = "gaussian",margins = "ECDF",nMC = 10^4){

  if(margins == "GEV"){
    gev1 <- extRemes::fevd(x,method = "MLE", type = "GEV")
    gev2 <- extRemes::fevd(y,method = "MLE", type = "GEV")

    x <- extRemes::pevd(x,loc = gev1$results$par[1],scale = gev1$results$par[2],shape = gev1$results$par[3])
    y <- extRemes::pevd(y,loc = gev2$results$par[1],scale = gev2$results$par[2],shape = gev2$results$par[3])
  } else if(margins == "ECDF"){
    x <- rvinecopulib::pseudo_obs(x)
    y <- rvinecopulib::pseudo_obs(y)
  }

  Mydata <- cbind(x)

  for(i in 1:lx){
    xlag <- c(rep(NA,i),x[-(length(x) - (1:i - 1))])
    Mydata <- cbind(Mydata,xlag)
  }

  for(i in ly:1){
    ylag <- c(rep(NA,i),y[-(length(y) - (1:i - 1))])
    Mydata <- cbind(Mydata,ylag)
  }

  Mydata <- cbind(Mydata,y)
  u <- Mydata[-(1:(max(lx,ly))),]
  n <- nrow(u)
  d <- ncol(u)

  cop.np  <- rvinecopulib::vinecop(u, family_set=cop.fam,
                     structure = rvinecopulib::dvine_structure(1:d), psi0=0.9,
                     show_trace=FALSE, cores=1,selcrit = "mbicv")
  cop.np.pcs <- rvinecopulib::get_all_pair_copulas(cop.np)

  uMC <- rvinecopulib::rvinecop(nMC,cop.np,cores=1)

  pseudo <- array(NA,dim = c(nMC,2,d-2,d-1))
  for(i in 1:(d-2)){
    for(j in 1:(d-i)){
      cop.sub <- rvinecopulib::get_pair_copula(cop.np,tree = i,edge = j)
      if(i == 1){
        pseudo[,1,i,j] <- rvinecopulib::hbicop(uMC[,c(j,(j+1))],cond_var = 2,cop.sub)
        pseudo[,2,i,j] <- rvinecopulib::hbicop(uMC[,c(j,(j+1))],cond_var = 1,cop.sub)
      } else {
        pseudo[,1,i,j] <- rvinecopulib::hbicop(cbind(pseudo[,1,i-1,j],pseudo[,2,i-1,j+1]),
                                 cond_var = 2,cop.sub)
        pseudo[,2,i,j] <- rvinecopulib::hbicop(cbind(pseudo[,1,i-1,j],pseudo[,2,i-1,j+1]),
                                 cond_var = 1,cop.sub)
      }
    }
  }

  te_comp <- 1
  for(i in 1:ly){
    cop.sub <- rvinecopulib::get_pair_copula(cop.np,tree = lx+i,edge = 1)
    te_comp <- te_comp*rvinecopulib::dbicop(cbind(pseudo[,1,lx+i-1,1],pseudo[,2,lx+i-1,2]),cop.sub)
    cop.np.pcs[[lx+i]][[1]] <- rvinecopulib::bicop_dist(family = "indep")
  }

  mean_te.ytox <- mean(log(te_comp),na.rm = TRUE)

  te_comp <- 1
  for(i in 1:lx){
    cop.sub <- rvinecopulib::get_pair_copula(cop.np,tree = ly+i,edge = d-(ly+i))
    te_comp <- te_comp*rvinecopulib::dbicop(cbind(pseudo[,1,ly+i-1,d-(ly+i)],pseudo[,2,ly+i-1,d-(ly+i)+1]),cop.sub)
    cop.np.pcs[[ly+i]][[d-(ly+i)]] <- rvinecopulib::bicop_dist(family = "indep")
  }

  mean_te.xtoy <- mean(log(te_comp),na.rm = TRUE)

  mean_te <- c(mean_te.ytox,mean_te.xtoy)
  names(mean_te) <- c("Y_to_X","X_to_Y")
  mean_te[mean_te<0] <- 0

  output <- list(estimate = mean_te,pval = rep(NA,2))

  cop.np.null <- rvinecopulib::vinecop_dist(cop.np.pcs,rvinecopulib::dvine_structure(1:d))

  if(nboot > 0){

    mean_te_bs <- foreach::foreach(b=1:nboot,.combine='rbind.data.frame',.errorhandling='remove',
                          .packages=c("rvinecopulib")) %dopar% {

                            u_bs <- rvinecopulib::rvinecop(n,cop.np.null,cores=1)

                            cop.np  <- rvinecopulib::vinecop(u_bs, family_set=cop.fam,
                                               structure = rvinecopulib::dvine_structure(1:d), psi0=0.9,
                                               show_trace=FALSE, cores=1,selcrit = "mbicv")

                            uMC <- rvinecopulib::rvinecop(nMC,cop.np,cores=1)

                            pseudo <- array(NA,dim = c(nMC,2,d-2,d-1))
                            for(i in 1:(d-2)){
                              for(j in 1:(d-i)){
                                cop.sub <- rvinecopulib::get_pair_copula(cop.np,tree = i,edge = j)
                                if(i == 1){
                                  pseudo[,1,i,j] <- rvinecopulib::hbicop(uMC[,c(j,(j+1))],cond_var = 2,cop.sub)
                                  pseudo[,2,i,j] <- rvinecopulib::hbicop(uMC[,c(j,(j+1))],cond_var = 1,cop.sub)
                                } else {
                                  pseudo[,1,i,j] <- rvinecopulib::hbicop(cbind(pseudo[,1,i-1,j],pseudo[,2,i-1,j+1]),
                                                           cond_var = 2,cop.sub)
                                  pseudo[,2,i,j] <- rvinecopulib::hbicop(cbind(pseudo[,1,i-1,j],pseudo[,2,i-1,j+1]),
                                                           cond_var = 1,cop.sub)
                                }
                              }
                            }

                            te_comp <- 1
                            for(i in 1:ly){
                              cop.sub <- rvinecopulib::get_pair_copula(cop.np,tree = lx+i,edge = 1)
                              te_comp <- te_comp*rvinecopulib::dbicop(cbind(pseudo[,1,lx+i-1,1],pseudo[,2,lx+i-1,2]),cop.sub)
                            }

                            mean_te.ytox_bs <- mean(log(te_comp),na.rm = TRUE)

                            te_comp <- 1
                            for(i in 1:lx){
                              cop.sub <- rvinecopulib::get_pair_copula(cop.np,tree = ly+i,edge = d-(ly+i))
                              te_comp <- te_comp*rvinecopulib::dbicop(cbind(pseudo[,1,ly+i-1,d-(ly+i)],pseudo[,2,ly+i-1,d-(ly+i)+1]),cop.sub)
                            }

                            mean_te.xtoy_bs <- mean(log(te_comp),na.rm = TRUE)

                            mean_te_bs_entry <- data.frame(Y_to_X = mean_te.ytox_bs,X_to_Y = mean_te.xtoy_bs)
                          }

    output$pval <- c(sum(mean_te_bs[,1] >= mean_te[1]),sum(mean_te_bs[,2] >= mean_te[2]))/nrow(mean_te_bs)
    names(output$pval) <- c("Y_to_X","X_to_Y")
  }
  return(output)
}

#' Vine Copula-Based Inference for Spectral Transfer Entropy
#'
#' Estimates spectral transfer entropy between the \code{fbandx} oscillation
#' of time series \code{x} and the \code{fbandy} oscillation of time series
#' \code{y}, in both directions: from \code{x} to \code{y} and from
#' \code{y} to \code{x}. Estimation is based on the vine copula representation
#' of Redondo et al. (2025).
#'
#' Available frequency bands are:
#' \itemize{
#'   \item \code{"delta"}: 0.5--4 Hz
#'   \item \code{"theta"}: 4--8 Hz
#'   \item \code{"alpha"}: 8--12 Hz
#'   \item \code{"beta"}: 12--30 Hz
#'   \item \code{"gamma"}: 30--45 Hz
#' }
#'
#' @param x A numeric vector containing the first time series.
#' @param y A numeric vector containing the second time series.
#' @param samprate A positive numeric scalar giving the sampling frequency
#'   in Hz.
#' @param blength A positive integer specifying the block size used to
#'   aggregate oscillation magnitudes into a block maxima series.
#' @param fbandx A character scalar specifying the frequency band to extract
#'   from \code{x}. Possible values are \code{"delta"}, \code{"theta"},
#'   \code{"alpha"}, \code{"beta"}, and \code{"gamma"}.
#' @param fbandy A character scalar specifying the frequency band to extract
#'   from \code{y}. Possible values are \code{"delta"}, \code{"theta"},
#'   \code{"alpha"}, \code{"beta"}, and \code{"gamma"}.
#' @param lx A positive integer specifying the number of lags of \code{x}
#'   used to estimate spectral transfer entropy.
#' @param ly A positive integer specifying the number of lags of \code{y}
#'   used to estimate spectral transfer entropy.
#' @param nboot A non-negative integer specifying the number of bootstrap
#'   replicates used to compute p-values. The bootstrap is based on the
#'   null-copula resampling scheme introduced in Redondo et al. (2025).
#'   If \code{nboot = 0}, the bootstrap procedure is skipped and only
#'   point estimates are returned.
#' @param cop.fam A character vector specifying the copula families used
#'   to estimate the vine copula models. Options follow those supported by
#'   the \pkg{rvinecopulib} package.
#' @param margins A character scalar specifying the method used to estimate
#'   the marginal distributions of the block maxima series. Possible values
#'   are \code{"GEV"}, for a generalized extreme value distribution, and
#'   \code{"ECDF"}, for the empirical distribution function.
#' @param nMC A positive integer specifying the number of Monte Carlo samples
#'   used to approximate the spectral transfer entropy.
#'
#' @return A list containing spectral transfer entropy estimates in both
#'   directions and, if \code{nboot > 0}, the corresponding bootstrap
#'   p-values.
#' @examples
#' # Registering number of cores to be used for bootstrap (only run if nboot>0)
#' # no_cores <- parallel::detectCores() - 1
#' # cl <- parallel::makeCluster(no_cores)
#' # doParallel::registerDoParallel(cl)
#'
#' # Simulating causal frequency bands (theta--gamma coupling)
#' set.seed(1)
#' phi_tg <- matrix(c(-0.1,    0,    0,    0,
#'                    2.5,  0.1,    0,    0,
#'                    0,    0, -0.1,    0,
#'                    2.5,    0,  2.5,  0.1), 4, byrow = TRUE)
#'
#' Z <- causal.fband(nsec = 50,samprate = 100,fband = c("theta","theta","gamma","gamma"),
#'                   phi = phi_tg,sigma = diag(rep(1,4)),range = c(45,55))
#'
#' X <- 0.45*Z[,1] + 0.5*Z[,3] + 0.05*rnorm(10000)
#' Y <- 0.5*Z[,2] + 0.45*Z[,4] + 0.05*rnorm(10000)
#' plot.ts(X)
#' plot.ts(Y)
#'
#' STE(x = X,y = Y,samprate = 100,blength = 50,
#'     fbandx = "theta",fbandy = "theta",lx = 1,ly = 1,nboot = 0)
#' STE(x = X,y = Y,samprate = 100,blength = 50,
#'     fbandx = "theta",fbandy = "gamma",lx = 1,ly = 1,nboot = 0)
#' STE(x = X,y = Y,samprate = 100,blength = 50,
#'     fbandx = "gamma",fbandy = "gamma",lx = 1,ly = 1,nboot = 0)
#' @export
STE <- function(x,y,samprate,blength,fbandx,fbandy,lx = 1,ly = 1,
                nboot = 0,cop.fam = "gaussian",margins = "ECDF",nMC = 10^4){
  Tlength <- nrow(cbind(x,y))
  overlap <- blength
  n_block <- floor((Tlength-blength)/overlap + 1)
  index <- overlap*rep(1:n_block-1,each = blength) + rep(1:blength,times = n_block)

  dat_x <- apply(matrix(abs(extract.fband(x = x,samprate = samprate,fband = fbandx)[[1]][index]),ncol = n_block),MARGIN = 2,FUN = max)
  dat_y <- apply(matrix(abs(extract.fband(x = y,samprate = samprate,fband = fbandy)[[1]][index]),ncol = n_block),MARGIN = 2,FUN = max)

  output <- copTE(x = dat_x,y = dat_y,margins = margins,lx = lx,ly = lx,nboot = nboot,cop.fam = cop.fam,nMC = nMC)

  return(output)
}




