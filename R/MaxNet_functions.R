# This function is based the function of same name from the package 'enmSdm'
# This function is adapted to allow parsing the desired combination of feature classes. 
# The original function from the package forces all combination of feature classes, which sometimes is not what the researchers wants.

#' Calibrate a MaxNet (MaxEnt) model using AICc
#'
#' This function calculates the "best" MaxEnt model using AICc across all possible combinations of a set of master regularization parameters and feature classes. The "best" model has the lowest AICc, with ties broken by number of features (fewer is better), regularization multiplier (higher better), then finally the number of coefficients (fewer better). The function can return the best model (default), a list of models created using all possible combinations of feature classes and regularization multipliers, and/or a data frame with tuning statistics for each model. Models in the list and in the data frame are sorted from best to worst.
#'
#' @param data  Data frame or matrix. Contains a column indicating whether each row is a presence (1) or background (0) site, plus columns for environmental predictors.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh}) or any combination of \code{lpqht}, where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quadratic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical.  If \code{TRUE} (default) then test all possible combinations of classes (note that all tested models will at least have linear features). If \code{FALSE} then use the classes provided (these will not vary between models).
#' @param dropOverparam Logical, if \code{TRUE} (default), drop models if they have more coefficients than training occurrences. It is possible for no models to fulfill this criterion, in which case no models will be returned.
#' @param anyway Logical. Same as \code{dropOverparam} (included for backwards compatibility. If \code{NULL} (default), then the value of \code{dropOverparam} will take precedence. If \code{TRUE} or \code{FALSE} then \code{anyway} will override the value of \code{dropOverparam}.
#' @param out Character or character vector. Indicates type of value returned. Values can be \code{'model'} (default; return model with lowest AICc), \code{'models'} (return a list of all models), and/or \code{'tuning'} (return a data frame with AICc for each model). If more than one value is specified, then the output will be a list with elements named "model", "models", and/or "tuning". If \code{'models'} is specified, they will only be produced if \code{select = TRUE}. The models will appear in the list in same order as they appear in the tuning table (i.e., model with the lowest AICc first, second-lowest next, etc.). If just one value is specified, the output will be either an object of class \code{MaxEnt}, a list with objects of class \code{MaxEnt}, or a data frame.
#' @param forceLinear Logical. If \code{TRUE} (default) then require any tested models to include at least linear features.
#' @param cores Integer >= 1. Number of cores to use. Default is 1.
#' @param verbose Logical. If \code{TRUE} report the AICc table.
#' @param ... Extra arguments. Not used.
#' @return If \code{out = 'model'} this function returns an object of class \code{MaxEnt}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters, log likelihood, and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{MaxEnt} object and the data frame.
#' @details This function is a wrapper for \code{\link[maxnet]{maxnet}}.
#' @seealso \code{\link[maxnet]{maxnet}}, \code{\link[dismo]{maxent}}, \code{\link{trainMaxEnt}}
#' @references
#' Phillips, S.J., Anderson, R.P., Dudík, M. Schapire, R.E., and Blair, M.E.  2017.  Opening the black box: An open-source release of Maxent. \emph{Ecography} 40:887-893.
#' Warren, D.L. and S.N. Siefert. 2011. Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria. \emph{Ecological Applications} 21:335-342.

# spname = sptogo
# data=PresAbsFull
# verbose=TRUE
# resp = names(data)[1]
# preds = names(data)[2:ncol(data)]
# kfolds = NULL
# regMult = seq(0.5, 4, by = 0.5)
# classes = c("l","h","lq","lqh","lqhp")
# out = c("model", "tuning")
# cores = n_jobs
# output_dir = output_dir
# check_if_exist = FALSE
# check_if_exist = TRUE
# dropOverparam = TRUE
# anyway = TRUE

trainMaxNet <- function (spname, data, resp = names(data)[1], preds = names(data)[2:ncol(data)], 
                         kfolds = NULL,
                         regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10), classes = c("l","h","lq","lqh","lqhp"), 
                         dropOverparam = TRUE, anyway = TRUE, 
                         out = "model", cores = 1, verbose = FALSE,
                         output_dir = NULL, check_if_exist = FALSE) 
{
    output_file <- here::here(output_dir,paste0(spname,"_model.RDS"))
    if(check_if_exist){
        test <- file.exists(output_file)
    } else {
        test <- FALSE
    }
    if(test){
        cat("SDMs were fitted for this species. Loading SDMs from file...")
        readRDS(output_file)
    } else {
        cat("Fitting SDMs")
        
        # response and predictors
        if (inherits(data, 'data.table')) data <- as.data.frame(data)
        if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
        if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]
        
        # get response and predictors
        presBg <- data[ , resp, drop=TRUE]
        data <- data[ , preds, drop=FALSE]
        
        ## collate all presences
        allPres <- data[presBg == 1, , drop=FALSE]
        
        ## collate all background sites
        allBg <- data[presBg == 0, , drop=FALSE]
        
        classesToTest <- classes
        if (any("p" %in% classesToTest) & ncol(data) == 1) {
            product <- FALSE
            warning("Data has only one variable so forcing product features to FALSE.")
        }
        tuning <- expand.grid(regMult, classesToTest)
        names(tuning) <- c("regMult","classes")
        if (cores > 1) {
            cat("Fitting", nrow(tuning), "MaxNet models in parallel")
            `%makeWork%` <- foreach::`%dopar%`
            cl <- parallel::makeCluster(cores)
            doParallel::registerDoParallel(cl)
        } else {
            `%makeWork%` <- foreach::`%do%`
        }
        mcOptions <- list(preschedule = TRUE, set.seed = FALSE, silent = FALSE)
        work <- foreach::foreach(
            i = 1:nrow(tuning), .options.multicore = mcOptions, 
            .export = c("trainMaxNetWorker",
                        "predictMaxNet",
                        "crossvalSDM",
                        "evalSDM",
                        "TSS",
                        "sediWeighted",
                        "contBoyce"),
            .inorder = FALSE) %makeWork% {
                trainMaxNetWorker(i = i, kfolds = kfolds, tuning = tuning, presBg = presBg, data = data, out = out)
            }
        
        if (cores > 1) { parallel::stopCluster(cl) }
        if(class(work) == "try-error"){ stop("Someting went wrong with the paralelization of models") }
        tuning <- data.frame()
        for (i in seq_along(work)) {
            tuning <- rbind(tuning, work[[i]]$thisTuning)
        }
        modelOrder <- order(tuning$AICc, tuning$numClasses, tuning$regMult, 
                            tuning$numCoeff, decreasing = c(FALSE, FALSE, TRUE, FALSE))
        tuning <- tuning[modelOrder, ]
        if (nrow(tuning) > 0) {
            tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
            tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
            tuning$aicWeight <- tuning$relLike/sum(tuning$relLike)
            rownames(tuning) <- 1:nrow(tuning)
        }
        if (verbose) {
            omnibus::say("")
            print(tuning, digits = 3)
            omnibus::say("")
        }
        if (length(out) > 1) {
            output <- list()
            if ("models" %in% out) {
                models <- list()
                for (i in seq_along(work)) {
                    models[[i]] <- work[[i]]$model
                }
                models <- models[modelOrder]
                if (dropOverparam) {
                    topModel <- models[[1]]
                    topTuning <- tuning[1, , drop = FALSE]
                    overparamModels <- which(tuning$n < tuning$numCoeff)
                    if (length(overparamModels) > 0) {
                        tuning <- tuning[-overparamModels, ]
                        models <- models[-overparamModels]
                    }
                    if (length(models) == 0 & anyway) {
                        warning("No models had fewer coefficients than predictors. Returning best model anyway.", 
                                immediate. = TRUE)
                        model <- topModel
                        tuning <- topTuning
                    }
                    else if (length(models) == 0 & !anyway) {
                        warning("No models had fewer coefficients than predictors. NA model(s) returned.", 
                                immediate. = TRUE)
                        model <- NULL
                    }
                    else {
                        model <- models[[1]]
                    }
                }
                output$models <- models
            }
            if ("model" %in% out) 
                models <- list()
            for (i in seq_along(work)) {
                models[[i]] <- work[[i]]$model
            }
            model <- models[modelOrder]
            topModel <- models[[1]]
            output$model <- topModel
            if ("tuning" %in% out) 
                output$tuning <- tuning
            if(!is.null(output_dir)){
                if(!dir.exists(output_dir)){
                    dir.create(output_dir,recursive = T)
                }
                saveRDS(output, output_file)
            }
            return(output)
        } else if (out == "models") {
            return(models)
        } else if (out == "model") {
            return(model)
        } else if (out == "tuning") {
            return(tuning)
        }
    }
}

#######################
### worker function ###
#######################

trainMaxNetWorker <- function(
        i,								# iterator
        kfolds,                         # n folds for cross validation
        tuning,							# tuning data frame
        presBg,						# vector with 1 (present) / 0 (background)
        data,							# df with all presence/background environmental data
        out
) {
    
    thisRegMult <- tuning$regMult[i]
    thisClasses <- tuning$classes[i]
    
    allPres <- data[presBg == 1, , drop = FALSE]
    allBg <- data[presBg == 0, , drop = FALSE]
    
    # add dummy column if doing univariate model to avoid error in maxnet.default.regularizationMOD
    if (ncol(data) == 1 & thisClasses == 'l') {
        
        thisData <- data
        thisPres <- allPres
        thisBg <- allBg
        
        thisData$DUMMY <- rep(1, nrow(thisData))
        thisPres$DUMMY <- rep(1, nrow(allPres))
        thisBg$DUMMY <- rep(1, nrow(allBg))
        
    } else {
        
        thisData <- data
        thisPres <- allPres
        thisBg <- allBg
        
    }
    
    # train model
    model <- maxnet::maxnet(
        p=as.vector(presBg),
        data=thisData,
        f=maxnet::maxnet.formula(p=as.vector(presBg), data=thisData, classes=thisClasses),
        regmult=thisRegMult
    )
    
    ## predict to all data
    predFull <- predictMaxNet(
        model=model,
        newdata=thisData)
    
    if(all(is.na(predFull))){
        # impossible to fit the model...
        # do nothing...
        
    } else {
        ## predict to training (and maybe test presences)
        predPres <- predictMaxNet(
            model=model,
            newdata=thisPres)
        
        ## predict to background
        predBg <- predictMaxNet(
            model=model,
            newdata=thisBg)
        
        # AICc
        rawSum <- sum(c(predPres, predBg), na.rm=TRUE)
        ## log likelihood
        ll <- sum(log(predPres / rawSum), na.rm=TRUE)
        ## number of parameters
        numCoeff <- length(model$betas)
        AICc <- -2 * ll + 2 * numCoeff + (2 * numCoeff * (numCoeff + 1)) / (sum(presBg) - numCoeff - 1)
        
        if(!is.null(kfolds)){
            ## predicted cross validation
            predCV <- crossvalSDM(kfolds = 5, traindata = thisData, pa = as.vector(presBg),
                                  classes = thisClasses, regmult = thisRegMult)
            
            ## evaluate SDMs
            evalCV <- evalSDM(observation = as.vector(presBg), predictions = predCV)
            names(evalCV) <- paste(names(evalCV), "CV", sep = "_")
            
            thisTuning <- cbind(thisTuning, evalCV)
        }else{
            thisTuning=data.frame(
                regMult=thisRegMult,
                classes=thisClasses,
                numClasses=nchar(as.character(thisClasses)),
                numCoeff=numCoeff,
                logLik=ll,
                AICc=AICc)
        }
        
        evalFull <- evalSDM(observation = as.vector(presBg), predictions = predFull)
        names(evalFull) <- paste(names(evalFull), "Full", sep = "_")
        
        thisTuning <- cbind(thisTuning, evalFull)
        output <- list()
        if ("model" %in% out) 
            output$model <- model
        if ("tuning" %in% out) 
            output$thisTuning <- thisTuning
        
        return(output)
    }
}


###############################
# modified from 'mecofun' package
#' crossvalSDM 
#'
#' A function for deriving cross-validated predictions. The function partitions the data into k folds, determines the model algorithm, updates the model for the new training data and makes predictions to the hold-out data using this algorithm.
crossvalSDM <- function(kfolds, traindata, pa, classes = NULL, regmult = NULL) {
    
    foldsPres <- which(pa == 1)
    foldsPres <- cut(1:length(foldsPres),breaks=kfolds,labels=FALSE)
    foldsAbs <- which(pa == 0)
    foldsAbs <- cut(1:length(foldsAbs),breaks=kfolds,labels=FALSE)
    ks <- c(foldsPres,foldsAbs)
    ks <- ks[sample(1:length(ks))] # force random positioning of folds
    
    preds <- mclapply(seq_len(kfolds), function(i) {
        data_train <- traindata[ks!=i,]
        data_test <- traindata[ks==i,]
        pa_train <- pa[ks!=i]
        
        # train the model with the training data
        modtmp <- maxnet::maxnet(p = pa_train, data = data_train, 
                                 f=maxnet::maxnet.formula(p=as.vector(pa_train), data=data_train, 
                                                          classes=ifelse(is.null(classes), "default", classes)),
                                 regfun=maxnet::maxnet.default.regularization,
                                 regmult=ifelse(is.null(regmult), 1, regmult))
        
        # make predictions for k-fold using the test data
        predictMaxNet(modtmp, data_test)
    }, mc.cores = kfolds)
    
    # feed
    cross_val_preds = rep(NA,length(ks))
    for(i in seq_len(kfolds)){
        cross_val_preds[ks==i] <- preds[[i]]
    }
    
    cross_val_preds
}


###############################
# This function is based the function of same name from the package 'enmSdm'
#' Predictions from a MaxNet model
#'
#' This function is the same as the \code{predict} function in the \pkg{maxnet} package, except that:
#' \itemize{
#'		\item	If the input is a data frame, the output is a vector as output (not a single-column matrix);
#'		\item	If the input is a \code{SpatRaster} or \code{Raster*}, the output is of the same type;
#'		\item	The default output is on the cloglog scale
#'		\item   The function can be explicitly called (versus doing, say, \code{maxnet:::predict.maxnet}, which does not work even when that would be really useful...).
#' }
#'
#' @param model	Object of class \code{maxnet}.
#' @param newdata	Object of class \code{data.frame}, \code{SpatRaster} (\pkg{terra} package), or any of \code{Raster}, \code{RasterBrick}, or \code{RasterStack} (\pkg{raster} package).
#' @param clamp		If \code{TRUE} (default), predict outside the range of training data by 'clamping' values to the last value.
#' @param type		One of:
#' \itemize{
#'		\item		\code{cloglog} (default): Predictions are on a complementary log-log scale.
#'		\item		\code{logistic}: Predictions are on a logistic scale (and thus technically the same to several decimal places as predictions from MaxEnt <=3.3.3k, except for differences in default features).
#'		\item		\code{link}: Predictions are on the scale of the predictors.
#'		\item		\code{exponential}: Predictions are on an exponential ('raw') scale.
#'	}
#' @param ... Other arguments (unused).
#' @return Numeric vector.
#' @seealso \code{\link[raster]{predict}} from the \pkg{raster} package, \code{\link[terra]{predict}} from the \pkg{terra} package, and \code{\link[maxnet]{maxnet}} (see the \code{predict} function therein)
#' @export

predictMaxNet <- function(model, newdata, clamp=TRUE, type='cloglog', ...) {
    
    if (inherits(newdata, 'Raster')) {
        out <- raster::predict(newdata, model, fun=predictMaxNet, ...)
    } else if (inherits(newdata, 'SpatRaster')) {
        out <- terra::predict(newdata, model, fun=predictMaxNet, ...)
    } else {
        
        if (clamp) {
            for (v in intersect(names(model$varmax), names(newdata))) {
                newdata[ , v] <- pmin(pmax(newdata[ , v], model$varmin[v]), model$varmax[v])
            }
        }
        terms <- sub('hinge\\((.*)\\):(.*):(.*)$', 'maxnet:::hingeval(\\1,\\2,\\3)', names(model$betas))
        terms <- sub('categorical\\((.*)\\):(.*)$', 'maxnet:::categoricalval(\\1,\'\\2\')', terms)
        terms <- sub('thresholds\\((.*)\\):(.*)$', 'maxnet:::thresholdval(\\1,\\2)', terms)
        f <- formula(paste('~', paste(terms, collapse=' + '), '-1'))
        mm <- model.matrix(f, data.frame(newdata))
        if (clamp) mm <- t(pmin(pmax(t(mm), model$featuremins[names(model$betas)]), 
                                model$featuremaxs[names(model$betas)]))
        link <- (mm %*% model$betas) + model$alpha
        if (type=='cloglog') out <- 1 - exp(0 - exp(model$entropy + link))
        if (type=='logistic') out <- 1 / (1 + exp(-model$entropy - link))
        if (type=='exponential') out <- exp(link)
        if (type=='link') out <- link
        out <- out[ , 1]
    }
    out
}

############
# modified from the package 'mecofun'
#' evalSDM
#'
#' Evaluate SDM perfomance
#' @param observation vector containing the observed response 
#' @param predictions vector containing the predictions
#' @param thresh.method a string indicating which method to use for optimising the binarising threshold (see ?PresenceAbsence::optimal.thresholds. Defaults to "MaxSens+Spec" (the maximum of sensitivity+specificity).
#' @param req.sens additional argument to PresenceAbsence::optimal.thresholds()
#' @param req.spec additional argument to PresenceAbsence::optimal.thresholds()
#' @param FPC additional argument to PresenceAbsence::optimal.thresholds()
#' @param FNC additional argument to PresenceAbsence::optimal.thresholds()
#' @return A dataframe with performance statistics.
#' @examples 
#' data(Anguilla_train)
#' m1 <- glm(Angaus ~ poly(SegSumT,2), data=Anguilla_train, family='binomial')
#' preds_cv <- crossvalSDM(m1, kfold=5, traindat=Anguilla_train, colname_species = 'Angaus', colname_pred = 'SegSumT')
#' evalSDM(Anguilla_train$Angaus, preds_cv)

evalSDM <- function(observation, predictions, thresh.method='MaxSens+Spec'){
    thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                             obs = observation,
                             pred = predictions)
    
    thresh <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat)
    cmx.opt <- PresenceAbsence::cmx(DATA= thresh.dat, threshold=thresh[thresh$Method==thresh.method,2])
    sedi <- max(sediWeighted(pres = predictions[which(observation==1)], contrast = predictions[which(observation==0)]))
    boyce <- contBoyce(pres = predictions[which(observation==1)], contrast = predictions[which(observation==0)])
    
    data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
               TSS = TSS(cmx.opt), 
               Kappa = PresenceAbsence::Kappa(cmx.opt, st.dev=F),
               Sens = PresenceAbsence::sensitivity(cmx.opt, st.dev=F),
               Spec = PresenceAbsence::specificity(cmx.opt, st.dev=F),
               PCC = PresenceAbsence::pcc(cmx.opt, st.dev=F),
               SEDI = sedi,
               boyce = boyce,
               thresh = thresh[thresh$Method==thresh.method,2],
               thresh_method=thresh.method)
}
#from the package 'mecofun'
TSS = function(cmx){
    PresenceAbsence::sensitivity(cmx, st.dev=F) + 
        PresenceAbsence::specificity(cmx, st.dev=F) - 1
}

############
# modified from the package 'enmSDM'

#' Symmetric extremal dependence index (SEDI) from {enmSdm}
#'
#' This function calculates the symmetric extremal dependence index (SEDI) using a threshold applied to continuous data to demarcate "presence" from "contrast".
#' @param pres Numeric vector. Predicted values at presence sites.
#' @param contrast Numeric vector. Predicted values at absence/background sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param thresholds Numeric vector, Values at which to threshold predictions for calculation of SEDI.
#' @param delta Positive numeric >0 in the range [0, 1] and usually very small. This value is used only if calculating the SEDI threshold when any true positive rate or false negative rate is 0 or the false negative rate is 1. Since SEDI uses log(x) and log(1 - x), values of 0 and 1 will produce \code{NA}s. To obviate this, a small amount can be added to rates that equal 0 and subtracted from rates that equal 1.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and absence/background predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param bgWeight Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.
#' @param ... Other arguments (unused).
#' @return Numeric value.
#' @references Wunderlich, R.F., Lin, Y-P., Anthony, J., and Petway, J.R.  2019.  Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribution models.  Nature Conservation 35:97-116.
#' @seealso \code{\link[stats]{cor}}, \code{\link{fpb}}, \code{\link{aucWeighted}}, \code{link[enmSdm]{contBoyce}}, \code{link[enmSdm]{contBoyce2x}}, \code{link[enmSdm]{thresholdWeighted}}, \code{link[enmSdm]{thresholdStats}}
#' @examples
#' set.seed(123)
#' pres <- sqrt(runif(100))
#' contrast <- runif(10000)
#' hist(contrast, col='gray', xlim=c(0, 1), breaks=20)
#' hist(pres, col='green', breaks=20, add=TRUE)
#' max(sediWeighted(pres, contrast), na.rm=TRUE)
#' # SEDI is fairly insensitive to weighting (as per Wunderlich et al. 2019)
#' presWeight <- c(rep(1, 50), rep(1000, 50))
#' max(sediWeighted(pres, contrast, presWeight=presWeight), na.rm=TRUE)

sediWeighted <- function(
        pres,
        contrast,
        presWeight = rep(1, length(pres)),
        contrastWeight = rep(1, length(contrast)),
        thresholds = seq(0, 1, by=0.01),
        delta = 0.001,
        na.rm = FALSE,
        bg = NULL,
        bgWeight = NULL,
        ...
) {
    
    if (missing(contrast) & !is.null(bg)) contrast <- bg
    if (missing(contrastWeight) & !is.null(bgWeight)) contrastWeight <- bgWeight
    
    # if all NAs
    if (all(is.na(pres)) | all(is.na(contrast)) | all(is.na(presWeight)) | all(is.na(contrastWeight))) return(NA)
    
    # catch errors
    if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
    if (length(contrastWeight) != length(contrast)) stop('You must have the same number of absence/background predictions and absence/background weights ("contrast" and "contrastWeight").')
    
    # remove NAs
    if (na.rm) {
        
        cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
        pres <- cleanedPres[[1]]
        presWeight <- cleanedPres[[2]]
        
        cleanedContrast <- omnibus::naOmitMulti(contrast, contrastWeight)
        contrast <- cleanedContrast[[1]]
        contrastWeight <- cleanedContrast[[2]]
        
    }
    
    sumPresWeights <- sum(presWeight)
    sumContrastWeights <- sum(contrastWeight)
    
    # SEDI, true positive rate, true negative rate, false negative rate
    sedi <- tpr <- tnr <- fnr <- rep(NA, length(thresholds))
    
    for (i in seq_along(thresholds)) {
        
        thisThresh <- thresholds[i]
        
        # which presences/contrast sites are CORRECTLY predicted at this threshold
        whichCorrectPres <- which(pres >= thisThresh)
        whichCorrectContrast <- which(contrast < thisThresh)
        
        numCorrectPres <- length(whichCorrectPres)
        numCorrectContrast <- length(whichCorrectContrast)
        
        anyCorrectPres <- (numCorrectPres > 0)
        anyCorrectContrast <- (numCorrectContrast > 0)
        
        # which presences/contrast sites are INCORRECTLY predicted at this threshold
        whichIncorrectPres <- which(pres < thisThresh)
        whichIncorrectContrast <- which(contrast >= thisThresh)
        
        numIncorrectPres <- length(whichIncorrectPres)
        numIncorrectContrast <- length(whichIncorrectContrast)
        
        anyIncorrectPres <- (numIncorrectPres > 0)
        anyIncorrectContrast <- (numIncorrectContrast > 0)
        
        # weights of CORRECTLY predicted predictions
        correctPresWeights <- if (anyCorrectPres) {
            sum(presWeight[whichCorrectPres])
        } else {
            0
        }
        
        correctContrastWeights <- if (anyCorrectContrast) {
            sum(contrastWeight[whichCorrectContrast])
        } else {
            0
        }
        
        # weights of INCORRECTLY predicted predictions
        incorrectPresWeights <- if (anyIncorrectPres) {
            sum(presWeight[whichIncorrectPres])
        } else {
            0
        }
        
        incorrectContrastWeights <- if (anyIncorrectContrast) {
            sum(contrastWeight[whichIncorrectContrast])
        } else {
            0
        }
        
        # true positive/negative rates
        tpr[i] <- correctPresWeights / sumPresWeights
        tnr[i] <- correctContrastWeights / sumContrastWeights
        
        # false positive/negative rates
        fnr[i] <- incorrectContrastWeights / sumContrastWeights
        
    }
    
    # SEDI
    if (any(tpr == 0)) tpr[tpr == 0] <- delta
    if (any(fnr == 0)) fnr[fnr == 0] <- delta
    if (any(tpr == 1)) tpr[tpr == 1] <- 1 - delta
    if (any(fnr == 1)) fnr[fnr == 1] <- 1 - delta
    
    numer <- log10(fnr) - log10(tpr) - log10(1 - fnr) + log10(1 - tpr)
    denom <- log10(fnr) + log10(tpr) + log10(1 - fnr) + log10(1 - tpr)
    sedi <- numer / denom
    
    sedi
    
}

############
# modified from the package 'enmSDM'

#' Continuous Boyce Index (CBI) with weighting from {enmSdm}
#'
#' This function calculates the Continuous Boyce Index (CBI), a measure of model accuracy for presence-only test data.
#' @param pres Numeric vector. Predicted values at presence sites.
#' @param contrast Numeric vector. Predicted values at absence/background sites.
#' @param numBins Positive integer. Number of (overlapping) bins into which to divide predictions.
#' @param binWidth Positive numeric. Size of a bin. Each bin will be \code{binWidth * (max - min)} where \code{max} and \code{min}. If \code{autoWindow} is \code{FALSE} (the default) then \code{min} is 0 and \code{max} is 1. If \code{autoWindow} is \code{TRUE} then \code{min} and \code{max} are the maximum and minimum value of all predictions in the background and presence sets (i.e., not necessarily 0 and 1).
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param autoWindow Logical. If \code{FALSE} calculate bin boundaries starting at 0 and ending at 1 + epsilon (where epsilon is a very small number to assure inclusion of cases that equal 1 exactly). If \code{TRUE} (default) then calculate bin boundaries starting at minimum predicted value and ending at maximum predicted value.
#' @param method Character. Type of correlation to calculate. The default is \code{'spearman'}, the Spearman rank correlation coefficient used by Boyce et al. (2002) and Hirzel et al. (2006), which is the traditional CBI. In contrast, \code{'pearson'} or \code{'kendall'} can be used instead.  See [stats::cor()] for more details.
#' @param dropZeros Logical. If TRUE then drop all bins in which the frequency of presences is 0.
#' @param graph Logical. If TRUE then plot P vs E and P/E versus bin.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param bgWeight Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.
#' @param ... Other arguments (not used).
#' @return Numeric value.
#' @details CBI is the Spearman rank correlation coefficient between the proportion of sites in each prediction class and the expected proportion of predictions in each prediction class based on the proportion of the landscape that is in that class. Values >0 indicate the model's output is positively correlated with the true probability of presence.  Values <0 indicate it is negatively correlated with the true probability of presence.
#' @references Boyce, M.S., Vernier, P.R., Nielsen, S.E., and Schmiegelow, F.K.A.  2002.  Evaluating resource selection functions.  \emph{Ecological Modeling} 157:281-300)
#' @references Hirzel, A.H., Le Lay, G., Helfer, V., Randon, C., and Guisan, A.  2006.  Evaluating the ability of habitat suitability models to predict species presences.  \emph{Ecological Modeling} 199:142-152.
#' @seealso \code{\link[stats]{cor}}, \code{\link{fpb}}, \code{\link{aucWeighted}}, \code{link[enmSdm]{contBoyce2x}}
#' @examples
#' set.seed(123)
#' pres <- sqrt(runif(100))
#' contrast <- runif(1000)
#' contBoyce(pres, contrast)
#' contBoyce2x(pres, contrast)
#' presWeight <- c(rep(1, 10), rep(0.5, 90))
#' contBoyce(pres, contrast, presWeight=presWeight)
#' contBoyce2x(pres, contrast, presWeight=presWeight)
#' \dontrun{
#' # compare stability of CBI calculated with ecospat.boyce() in ecospat package
#' library(ecospat)
#' set.seed(123)
#' results <- data.frame()
#' for (perform in c(1, 1.5, 2)) {
#' 	for (i in 1:30) {
#'
#'    pres <- runif(100)^(1 / perform)
#'    contrast <- runif(1000)
#'
#'    cbi_enmSdm <- contBoyce(pres, contrast)
#'    cbi_ecospat <- ecospat.boyce(contrast, pres, PEplot=FALSE)$Spearman.cor
#'
#'    results <- rbind(
#'      results,
#'      data.frame(
#'        performance = rep(perform, 2),
#'        method = c('enmSdm', 'ecospat'),
#'        cbi = c(cbi_enmSdm, cbi_ecospat)
#'      )
#'    )
#'
#' 	}
#'
#' }
#'
#' results$performance[results$performance == 1] <- 'poor'
#' results$performance[results$performance == 1.5] <- 'OK'
#' results$performance[results$performance == 2] <- 'good'
#'
#' results$category <- paste0(results$method, '\n', results$performance)
#'
#' par(mfrow=c(1, 2))
#' boxplot(cbi ~ category,
#' 	data=results,
#' 	ylab='CBI',
#' 	main='CBI of poor, OK, and good models',
#' 	border=c(rep('darkred', 3),
#' 	rep('darkblue', 3))
#' )
#' plot(results$cbi,
#' 	pch=rep(c(21, 22, 23, 24), each=2),
#' 	contrast=ifelse(results$method == 'ecospat', 'darkred', 'cornflowerblue'),
#' 	main='Pairs of CBIs',
#' 	ylab='CBI'
#' )
#' legend('bottomright', fill=c('darkred', 'cornflowerblue'), legend=c('ecospat', 'enmSdm'))
#' }
#' @export

contBoyce <- function(
        pres,
        contrast,
        numBins = 101,
        binWidth = 0.1,
        presWeight = rep(1, length(pres)),
        contrastWeight = rep(1, length(contrast)),
        autoWindow = TRUE,
        method = 'spearman',
        dropZeros = TRUE,
        graph = FALSE,
        na.rm = FALSE,
        bg = NULL,
        bgWeight = NULL,
        ...
) {
    
    if (missing(contrast) & !is.null(bg)) contrast <- bg
    if (missing(contrastWeight) & !is.null(bgWeight)) contrastWeight <- bgWeight
    
    # if all NAs
    if (all(is.na(pres)) | all(is.na(contrast)) | all(is.na(presWeight)) | all(is.na(contrastWeight))) return(NA)
    
    # catch errors
    if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
    if (length(contrastWeight) != length(contrast)) stop('You must have the same number of absence/background predictions and absence/background weights ("contrast" and "contrastWeight").')
    
    # # remove NAs
    # if (na.rm) {
    
    # cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
    # pres <- cleanedPres[[1]]
    # presWeight <- cleanedPres[[2]]
    
    # cleanedContrast <- omnibus::naOmitMulti(contrast, contrastWeight)
    # contrast <- cleanedContrast[[1]]
    # contrastWeight <- cleanedContrast[[2]]
    
    # }
    
    # right hand side of each class (assumes max value is >0)
    lowest <- if (autoWindow) { min(c(pres, contrast), na.rm=na.rm) } else { 0 }
    highest <- if (autoWindow) { max(c(pres, contrast), na.rm=na.rm) + omnibus::eps() } else { 1 + omnibus::eps() }
    
    windowWidth <- binWidth * (highest - lowest)
    
    lows <- seq(lowest, highest - windowWidth, length.out=numBins)
    highs <- seq(lowest + windowWidth + omnibus::eps(), highest, length.out=numBins)
    
    ##########
    ## MAIN ##
    ##########
    
    ## initiate variables to store predicted/expected (P/E) values
    freqPres <- freqContrast <- rep(NA, length(numBins))
    
    ### tally proportion of test presences/background sites in each class
    for (countClass in 1:numBins) {
        
        # number of presence predictions in this class
        presInBin <- pres >= lows[countClass] & pres < highs[countClass]
        presInBin <- presInBin * presWeight
        freqPres[countClass] <- sum(presInBin, na.rm=na.rm)
        
        # number of background predictions in this class
        bgInBin <- contrast >= lows[countClass] & contrast < highs[countClass]
        bgInBin <- bgInBin * contrastWeight
        freqContrast[countClass] <- sum(bgInBin, na.rm=na.rm)
        
    } # next predicted value class
    
    # mean bin prediction
    meanPred <- rowMeans(cbind(lows, highs))
    
    # add small number to each bin that has 0 background frequency but does have a presence frequency > 0
    if (any(freqPres > 0 & freqContrast == 0)) {
        smallValue <- min(0.5 * c(presWeight[presWeight > 0], contrastWeight[contrastWeight > 0]))
        freqContrast[freqPres > 0 & freqContrast == 0] <- smallValue
    }
    
    # remove classes with 0 presence frequency
    if (dropZeros && 0 %in% freqPres) {
        zeros <- which(freqPres == 0)
        meanPred[zeros] <- NA
        freqPres[zeros] <- NA
        freqContrast[zeros] <- NA
    }
    
    # remove classes with 0 background frequency
    if (any(0 %in% freqContrast)) {
        zeros <- which(freqContrast == 0)
        meanPred[zeros] <- NA
        freqPres[zeros] <- NA
        freqContrast[zeros] <- NA
    }
    
    P <- freqPres / sum(presWeight, na.rm=TRUE)
    E <- freqContrast / sum(contrastWeight, na.rm=TRUE)
    PE <- P / E
    
    # plot
    if (graph) {
        graphics::par(mfrow=c(1, 2))
        lims <- c(0, max(P, E, na.rm=TRUE))
        plot(E, P, col='white', xlab='Expected', ylab='Predicted', main='P/E\nNumbered from lowest to highest class', xlim=lims, ylim=lims)
        graphics::text(E, P, labels=1:numBins, col=1:20)
        plot(meanPred, PE, type='l', xlab='Mean Prediction in Bin', ylab='P/E Ratio', main='CBI\nNumbered from lowest to highest class')
        graphics::text(meanPred, PE, labels=1:numBins, col='blue')
    }
    
    # remove NAs
    out <- omnibus::naOmitMulti(meanPred, PE)
    meanPred <- out[[1]]
    PE <- out[[2]]
    
    # calculate continuous Boyce index (cbi)
    cbi <- stats::cor(x=meanPred, y=PE, method=method)
    cbi
    
}