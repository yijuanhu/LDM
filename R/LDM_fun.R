Rarefy = function (otu.tab, depth = min(rowSums(otu.tab)))
{
    otu.tab <- as.matrix(otu.tab)
    ind <- (rowSums(otu.tab) < depth)
    sam.discard <- rownames(otu.tab)[ind]
    otu.tab <- otu.tab[!ind, ]
    rarefy <- function(x, depth) {
        y <- sample(rep(1:length(x), x), depth)
        y.tab <- table(y)
        z <- numeric(length(x))
        z[as.numeric(names(y.tab))] <- y.tab
        z
    }
    otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
    rownames(otu.tab.rff) <- rownames(otu.tab)
    colnames(otu.tab.rff) <- colnames(otu.tab)
    return(list(otu.tab.rff = otu.tab.rff, discard = sam.discard))
} # Rarefy


#' Averaging the squared distance matrices each calculated from a rarefied OTU table
#' 
#' This function computes a distance matrix for each rarefied OTU table, square the distance matrix (in an element-wise manner), 
#' and then average the squared distance matrices.
#' 
#' @param otu.table the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' @param dist.method method for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package. The default is "jaccard". 
#' For more details, see the \code{dist.method} argument in the \code{ldm} function.
#' @param tree the phylogeneic tree. The default is NULL.
#' @param scale.otu.table a logical variable indicating whether to scale the rows of the OTU table. 
#' For count data, this corresponds to dividing by the library size to give relative frequencies. 
#' The default is FALSE.
#' @param n.rarefy number of rarefactions. The default is 100.
#' @param binary the "binary" parameter in \code{vegdist}. The default is TRUE.
#' @param seed a single-value integer seed for the random process of drawing rarefaction replicates. 
#' The seed is user supplied or internally generated. The default is 123.
#' @return a single matrix object
#'   \item{D2.avg}{The average of the squared distance matrices.}

#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gsatten@emory.edu>
#' @export
#' @examples
#' data(throat.otu.tab5)
#' dist.avg.D2 <- avgdist.squared(throat.otu.tab5, dist.method="jaccard", n.rarefy=100)

avgdist.squared = function(otu.table, dist.method="jaccard", tree=NULL, scale.otu.table=FALSE, n.rarefy=100, binary=TRUE, seed=123) {
    
    if (is.null(seed)) {
        seed = sample(1:10^6, 1)
    }
    set.seed(seed)
    
    n.sam = nrow(otu.table)
    D2.avg = matrix(0, n.sam, n.sam)
    
    for (r in 1:n.rarefy) {
        otu.rarefy = Rarefy(otu.table)$otu.tab.rff 
        if (binary) otu.rarefy = (otu.rarefy>0)*1
        
        dist <- calculate.dist(dist.method=dist.method, otu.table=otu.rarefy, tree=tree, scale.otu.table=scale.otu.table, binary=binary)
        D2.avg = D2.avg + dist^2
    }
    D2.avg = D2.avg/n.rarefy
    
    return(D2.avg)
    
} # avgdist2


gower = function (d, square=TRUE, center=TRUE) 
{
    
    a <- as.matrix(d)
    
    #--------------------------------------------------------------
    # squaring the matrix 
    #--------------------------------------------------------------
    
    if (square) a <- a^2
    
    #--------------------------------------------------------------
    # centering rows and columns of the (squared) matrix of d
    #--------------------------------------------------------------
    
    if (center) {

        a =  sweep(a, 1, rowMeans(a) )
        a = -sweep(a, 2, colMeans(a) )/2

    }
    
    return(a)
    
}# gower


fdr.Sandev = function(p.otu) {
    
    m = length(p.otu)
    
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))
    
    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)   
    qval.orig = qval.sort[mat]
    results = qval.orig
    return(results)
    
} # fdr.Sandev


#' @importFrom GUniFrac GUniFrac
#' @importFrom vegan vegdist
calculate.dist <- function(otu.table, tree=NULL, dist.method="bray", 
                           binary=FALSE, rarefy=0, scale.otu.table=TRUE) {
    
    rowsum = rowSums(otu.table)
    if (min(rowsum)==0) {
        warning("There exists sample(s) with zero reads at every OTU!")
    }
    
    rowsum[which(rowsum==0)] = 1
    if (scale.otu.table) freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowsum ) )
    else freq.table <- otu.table
    
    if (grepl(dist.method, "wt-unifrac")) {
        dist <- GUniFrac::GUniFrac(otu.table, tree, alpha=c(1))$unifrac[,,"d_1"]
    } else if (grepl(dist.method, "unwt-unifrac")) {
        dist <- GUniFrac::GUniFrac(otu.table, tree, alpha=c(1))$unifrac[,,"d_UW"]
    } else if (grepl(dist.method, "hellinger")) {
        dist <- 0.5*dist( x=sqrt(freq.table), method='euclidean')
    } else {
        dist <- vegan::vegdist(x=freq.table, method=dist.method, binary=binary)
    }
    
    # if (dist.method==tolower("bray")) {
    #     dist <- 0.5*as.matrix( dist( x=freq.table, method='manhattan') )
    # } else if (dist.method==tolower("Jaccard")) {
    # dist <- vegan::vegdist(x=otu.table, method="jaccard", binary=TRUE)
    # }
    
    dist <- as.matrix(dist)
    
    return(dist)
    
} # calculate.dist


#' Adjusting data (distance matrix and OTU table) by covariates
#' 
#' This function produces adjusted distance matrix and OTU table (if provided) 
#' after removing the effects of covariates (e.g., confounders). 
#' Observations with any missing data are removed.
#' 
#' @param formula a symbolic description of the covariate model in the form \code{ ~ model}, 
#' where \code{model} is specified in the same way as for \code{lm} or \code{glm}. For example, 
#' \code{~ a + b} specifies a model with the main effects of covariates \code{a} and \code{b}, and 
#' \code{~ a*b}, equivalently \code{~ a + b + a:b}, specifies a model with the main effects of 
#' \code{a} and \code{b} as well as their interaction.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the covariates. 
#' If not found in \code{data}, the covariates are taken from environment (formula), 
#' typically the environment from which \code{adjust.data.by.covariates} is called. 
#' The default is .GlobalEnv.
#' @param otu.table the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' If provided, an adjusted (and column-centered) OTU table at the frequency (i.e., relative abundance) scale 
#' and an adjusted (and columnn-centered) OTU table at the arcsin-root-transformed frequency scale are output. If provided, 
#' it is also used for calculating the distance matrix unless the distance matrix is directly 
#' imported through \code{dist}.
#' The default is NULL.
#' @param tree a phylogenetic tree. Only used for calculating a
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of
#'   requested distance does not require a phylogenetic tree, or if the distance
#'   matrix is directly imported through \code{dist}. The default is NULL.
#' @param dist.method method for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package
#'  (i.e., "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", 
#'  "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis")
#'   as well as "hellinger" and "wt-unifrac". 
#'   The default is "bray". 
#'   For more details, see the \code{dist.method} argument in the \code{ldm} function.
#' @param binary the "binary" parameter in \code{vegdist}. The default is FALSE.
#' @param dist a distance matrix. Can be either an object of class "dist" or "matrix".
#'   The elements of the distance matrix will be squared and then the matrix will be centered if the default choices 
#'   \code{square.dist=TRUE} and \code{center.dist=TRUE} are used. If \code{dist=NULL}, the distance matrix is 
#'   calculated from the \code{otu.table}, using the value of \code{dist.method} (and \code{tree} if required). 
#'   The default is NULL.
#' @param square.dist a logical variable indicating whether to square the 
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale the rows of the OTU table 
#'   for the frequency scale.  For count data, this corresponds to dividing by the library size to give 
#'   relative frequencies. The default is TRUE. 
#' @param center.otu.table a logical variable indicating whether to center the 
#'   columns of the OTU table. The OTU table should be centered if the distance 
#'   matrix has been centered. Applied to both OTU tables at frequency and transformed scales. The default is TRUE.
#' @param freq.scale.only a logical variable indicating whether to provide adjusted frequency-scale OTU table only 
#' (not adjusted OTU table at the arcsin-root transformed frequency scale). The default is FALSE.
#' @return a list consisting of 
#'   \item{adj.dist}{the (squared/centered) distance matrix
#'   after adjustment of covariates.}
#'   \item{y.freq}{the (column-centered) frequency-scale OTU table after adjustment of covariates.} 
#'   \item{y.tran}{the (column-centered) arcsin-root-transformed 
#'   OTU table after adjustment of covariates.} 
#' @keywords microbiome PCA ordination distance
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gsatten@emory.edu>
#' @export
#' @examples
#' data(throat.meta)
#' data(throat.otu.tab5)
#' adj.data <- adjust.data.by.covariates(formula= ~ Sex + AntibioticUse, data=throat.meta,
#'                                       otu.table=throat.otu.tab5, dist.method="bray")


adjust.data.by.covariates = function(formula=NULL, data=.GlobalEnv, 
                                     otu.table=NULL, tree=NULL, dist.method="bray", binary=FALSE, dist=NULL, 
                                     square.dist=TRUE, center.dist=TRUE, 
                                     scale.otu.table=TRUE, center.otu.table=TRUE,
                                     freq.scale.only=FALSE) {
    
    #------------------------
    # covariates (e.g., confounders)
    #------------------------
    
    options(na.action=na.pass)
    m1 = model.matrix(object=formula, data=data)
    
    if (!is.null(dist)) {
        dist = as.matrix(dist)
        if (dim(m1)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between covariates and dist' )
    }
    if (!is.null(otu.table)) {
        # remove zero OTUs
        w = which(colSums(abs(otu.table))>0)
        if (length(w) < ncol(otu.table)) {
            warning(paste(ncol(otu.table)-length(w), 'OTU(s) with zero counts in all samples are removed', sep=" "))
            otu.table = otu.table[,w,drop=FALSE]
        }
        
        if (dim(m1)[1] != dim(otu.table)[1]) 
            otu.table <- t(otu.table)
        if (dim(m1)[1] != dim(otu.table)[1]) stop( 'numbers of observations mismatch between covariates and otu.table' )
    }
    
    if (any(is.na(m1))) {
        w = rowAnys(is.na(m1))
        warning(paste(sum(w), 'observation(s) with any missing data are removed', sep=" "))
        
        w = !w
        if (!is.null(dist)) dist = dist[w, w]
        if (!is.null(otu.table)) otu.table = otu.table[w,]
        m1 = m1[w,]
    }
    
    center.m1 = TRUE
    if (center.m1) m1 = scale( m1, center=TRUE, scale=FALSE )
    
    #---------------------------------------
    # checking negative values in otu.table
    #---------------------------------------
    
    if (!is.null(otu.table)) {
        neg.exist = any(otu.table<0)
        if (neg.exist) {
            if (scale.otu.table == TRUE) {
                stop("The OTU table has negative values, so it does not make sense to use 'scale.otu.table=TRUE'")
            }
            if (freq.scale.only == FALSE) {
                stop("The OTU table has negative values, which cannot be arcsin-root transformed by 'freq.scale.only=FALSE'")
            }
        }
    }
    
    #------------------------
    # dist matrix
    #------------------------
    
    if (is.null(dist) & is.null(otu.table)) {
        stop( 'must specify one of dist and otu.table' )
    }
    
    if (!is.null(otu.table) & !is.null(dist)) {
        if (dim(otu.table)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between otu.table and dist' )
    }
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.method=dist.method, otu.table=otu.table, tree=tree, scale.otu.table=scale.otu.table, binary=binary)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)
    
    
    #---------------------
    # calculate d.resid
    #---------------------
    
    tol.d=10^-8
    svd.m1 = svd(m1)
    use = (svd.m1$d>tol.d)
    
    hat.matrix = tcrossprod(svd.m1$u[, use]) # svd.m1$u[, use] %*% t( svd.m1$u[, use] )
    hat.matrix.bar = diag(dim(d.gower)[1]) - hat.matrix
    
    d.resid = hat.matrix.bar %*% d.gower
    d.resid = d.resid %*% hat.matrix.bar
    
    #---------------------
    # calculate adj.otu.table
    #---------------------
    
    y.freq = NULL
    y.tran = NULL
    
    if (!is.null(otu.table)) {
        
        rowsum = rowSums(otu.table)
        if (min(rowsum)==0) {
            warning("There exists sample(s) with zero reads at every OTU!")
        }
        
        # freq
        if (scale.otu.table) {
            rowsum[which(rowsum==0)] = 1
            freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowsum ) )
        } else {
            freq.table <- otu.table
        }
        y.freq <- scale( freq.table, center=center.otu.table, scale=FALSE )
        
        y.tran <- NULL
        if (!freq.scale.only) {
            # arcsin
            theta <- asin(sqrt(freq.table))
            y.tran <- scale( theta, center=center.otu.table, scale=FALSE)
        }
        
        # # check discordant centering (turn off)
        # max.center.gower <- max(abs(rowMeans(d.gower)))
        # max.center.yfreq  <- max(abs(colMeans(y.freq)))
        # if ( (max.center.gower - 10^-6) * (max.center.yfreq - 10^-6) < 0) {
        #     stop( 'discordant centering of the OTU table and distance matrix' )
        # }
        
        # x.model
        x.model = svd.m1$u[, use]
        
        # adj.otu.table
        x1.tilda.freq = crossprod(x.model, y.freq) # t(x.model) %*% y.freq
        x1.tilda.freq = x.model %*% x1.tilda.freq
        y.freq = y.freq - x1.tilda.freq
        
        if (!freq.scale.only) {
            x1.tilda.tran = crossprod(x.model, y.tran) # t(x.model) %*% y.tran
            x1.tilda.tran = x.model %*% x1.tilda.tran
            y.tran = y.tran - x1.tilda.tran
        }
    }
    
    res <- list(y.freq=y.freq,
                y.tran=y.tran,
                adj.dist=d.resid)
    
    return(res)
    
} # adjust.data.by.covariates


#' Testing hypotheses about the microbiome using a linear decomposition model (LDM)
#' 
#' This function allows you to 
#' 1. simultaneously test the global association with the overall  
#' microbiome composition and individual OTU associations to give coherent 
#' results;
#' 2. test hypotheses based on data at both the frequency (i.e., relative abundance) and arcsine-root-transformed frequency scales, 
#' and perform an ``omnibus" test that combines results from analyses conducted on the two scales;
#' 3. test presence-absence associations based on infinite number of rarefaction replicates;
#' 4. handle complex design features such as 
#' confounders, interactions, and clustered data (with between- and within-cluster covariates);
#' 5. test associations with a survival outcome (i.e., censored survival times);
#' 6. perform mediation analysis of the microbiome;
#' 7. perform the omnibus test LDM-omni3 that combines results from analyses conducted on the frequency, arcsine-root-transformed, and presence-absence scales.
#' 
#' The formula has the form 
#' 
#' \code{otu.table ~ (first set of covariates) + (second set of covariates)
#' ... + (last set of covariates)} 
#' 
#' or 
#' 
#' \code{otu.table | confounders ~ (first set of covariates) + (second set of covariates)
#' ... + (last set of covariates)} 
#' 
#' where \code{otu.table} is
#' the OTU table with rows for samples and columns for OTUs and each set of 
#' covariates are enclosed in parentheses. The covariates in each submodel (set of covariates) are tested jointly,
#' after projecting off terms in submodels that appear earlier in the model.
#' 
#' For example, given OTU table \code{y} and a data frame \code{metadata} that contains 4 covariates, 
#' \code{a}, \code{b}, \code{c} and \code{d},  
#' some valid formulas would be:
#' 
#' \code{y ~ a + b + c + d} ### no confounders, 4 submodels (i.e., sets of covariates)
#' 
#' \code{y ~ (a+b) + (c+d)} ### no confounders, 2 submodels each having 
#' 2 covariates
#' 
#' \code{y | b ~ (a+c) + d} ### \code{b} is a confounder, submodel 1 is 
#' \code{(a+c)}, and submodel 2 is \code{d}
#' 
#' \code{y | b+c ~ a*d}     ### there are 2 confounders \code{b} 
#' and \code{c}; there is 1 submodel consisting of the three terms \code{a}, \code{d}, and \code{a:d} (interaction). 
#' This example is equivalent to \code{y | b+c ~ (a+d+a:d)}
#' 
#' \code{y | as.factor(b) ~ (a+d) + a:d}  ### the confounder 
#' \code{b} will be treated as a factor variable, submodel 1 will have the main 
#' effects \code{a} and \code{d}, and submodel 2 will have only the interaction 
#' between \code{a} and \code{d}
#' 
#' \code{y | as.factor(b) ~ (a) + (d) + (a:d)} ### there are 3 submodels \code{a}, \code{d}, and \code{a:d}.
#' Putting paratheses around a single variable is allowed but not necessary.
#'
#' Submodels that combine character and numeric values are allowed; character-valued variables are coerced into factor 
#' variables.  Confounders are distinguished from other covariates as test statistics are not calculated for confounders
#' (which are included for scientific reasons, not by virtue of significance test results); 
#' consequently they also do not contribute to stopping criteria.  If tests of confounders are desired, confounders should
#' put on the right hand side of the formula as the first submodel.
#' 
#' For testing mediation effects of the microbiome that mediate the effect of the exposure(s) on the outcome(s), 
#' the formula takes the specific form:
#' 
#' \code{otu.table ~ exposure + outcome}
#' 
#' or most generally
#' 
#' \code{otu.table | (set of confounders) ~ (set of exposures) + (set of outcomes)}
#' 
#' in which there should be exactly two terms on the right hand side of the regression, 
#' corresponding to the exposure(s) and the outcome(s), the outcome(s) must appear after the exposure(s), 
#' and the covariates or confounders must appear after \code{|}.
#' 
#' LDM uses two sequential stopping criteria. For the global test, LDM uses the 
#' stopping rule of Besag and Clifford (1991), which stops permutation when a 
#' pre-specified minimum number (default=100) of rejections (i.e., the permutation 
#' statistic exceeded the observed test statistic) has been reached. For the 
#' OTU-specific tests, LDM uses the stopping rule of Sandve et al. (2011), 
#' which stops permutation when every OTU test has either reached the pre-specified 
#' number (default=100) of rejections or yielded a q-value that is below the 
#' nominal FDR level (default=0.1). As a convention, we call a test "stopped"
#' if the corresponding stopping criterion has been satisfied. Although all tests 
#' are always terminated if a pre-specified maximum number (see description of \code{n.perm.max} in Arguments list) of 
#' permutations have been generated, some tests may not have "stopped".  This typically occurs when
#' the relevant p-value is small or near the cutoff for inclusion in a list of significant findings; 
#' for global tests meeting the stopping criterion is not critical, but 
#' caution is advised when interpreting OTU-level tests that have not stopped as additional OTUs may be found 
#' with a larger number of permutations.
#' 
#' 
#' @param formula a symbolic description of the 
#'   model to be fitted. The details of model specification are given under 
#'   "Details".
#' @param other.surv.resid a vector of data, usually the Martingale or deviance residuals from fitting the Cox model to the survival outcome (if it is the outcome of interest) and other covariates.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the covariates of interest and 
#' confounding covariates. 
#' If not found in \code{data}, the covariates are taken from environment(formula), 
#' typically the environment from which \code{ldm} is called. The default is .GlobalEnv.
#' @param tree a phylogenetic tree. Only used for calculating a 
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of 
#'   the requested distance does not involve a phylogenetic tree, or if the 
#'   distance matrix is directly imported through \code{dist}.
#' @param dist.method method for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package
#'  (i.e., "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", 
#'  "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
#'  "chao", "cao", "mahalanobis") as well as "hellinger" and "wt-unifrac". 
#'  The Hellinger distance measure (\code{dist.method="hellinger"}) takes the form
#'  \code{0.5*E}, where E is the Euclidean distance between the square-root-transformed 
#'  frequency data. The weighted UniFrac distance (\code{dist.method="wt-unifrac"}) 
#'  is calculated by interally calling \code{GUniFrac} in the \code{GUniFrac} package.
#'   Not used when anything other than \code{dist=NULL} is specified for \code{dist}.
#'   The default is "bray".
#' @param dist a distance matrix. Can be an object of class either "dist" or "matrix".
#'   The elements of the distance matrix will be squared and then the matrix will be centered if the default choices 
#'   \code{square.dist=TRUE} and \code{center.otu.table=TRUE} are used. If \code{dist=NULL}, the distance matrix is 
#'   calculated from the \code{otu.table}, using the value of \code{dist.method} (and \code{tree} if required). 
#'   The default is NULL.
#' @param cluster.id character or factor variable that identifies clusters. The default value
#'   cluster.id=NULL if the observations are not clustered (i.e., are independent).
#' @param strata a character or factor variable that defines strata (groups), within which to constrain permutations. 
#'   The default is NULL.
#' @param how a permutation control list, for users who want to specify their own call to the \code{how} function from the \code{permute} package.  
#'   The default is NULL.
#' @param perm.within.type a character string that takes values "free", "none", "series", or "grid".  
#'   The default is "free" (for random permutations).
#' @param perm.between.type a character string that takes values "free", "none", or "series".  
#'   The default is "none".
#' @param perm.within.nrow a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See documentation for permute package for additional details.
#' @param perm.within.ncol a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See documentation for permute package for additional details.
#' @param n.perm.max the maximum number of permutations. The default is NULL, in which case a maximum of
#'   5000 permutations are used for the global test and a maximum of \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal}) 
#'   are used for the OTU test, where \code{n.otu} is the number of OTUs.  If a numeric value for \code{n.perm.max} is specified, 
#'   this value is used for both global and OTU-level tests.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. The 
#'   default is 100.
#' @param seed a user-supplied integer seed for the random number generator in the 
#'   permutation procedure. The default is NULL; with the default value, an integer seed will be 
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case 
#'   the user wants to reproduce the permutation replicates.
#' @param fdr.nominal the nominal FDR value. The default is 0.1.
#' @param square.dist a logical variable indicating whether to square the 
#'   distance matrix. The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale the rows of the OTU table.  For count 
#'   data, this corresponds to dividing by the library size to give frequencies (i.e., relative abundances).  Does not affect the tran scale.  
#'   The default is TRUE. 
#' @param center.otu.table a logical variable indicating whether to center the 
#'   columns of the OTU table. The OTU table should be centered if the distance 
#'   matrix has been centered. Applied to both the frequency and transformed scales.  The default is TRUE.
#' @param freq.scale.only a logical variable indicating whether to perform analysis of the frequency-scale data only 
#' (not the arcsin-root transformed frequency data and the omnibus test). The default is FALSE.
#' @param binary a logical value indicating whether to perform presence-absence
#'   analysis. The default is FALSE (analyzing relative abundance data).
#' @param n.rarefy an integer-valued number of rarefactions. The value "all" is also allowed, 
#' and requests the LDM-A method that essentially aggregate information from all rarefactions.
#'  The default is 0 (no rarefaction).
#' @param n.cores The number of cores to use in parallel computing, i.e., at most how many child processes will be run simultaneously. 
#' The default is 4.
#' @param test.mediation a logical value indicating whether to perform the mediation analysis. The default is FALSE. 
#' If TRUE, the formula takes the specific form \code{otu.table ~ exposure + outcome} or most generally
#' \code{otu.table | (set of confounders) ~ (set of exposures) + (set of outcomes)}.
#' @param test.omni3 a logical value indicating whether to perform the new omnibus test (LDM-omni3). The default is FALSE. 
#' @param comp.anal a logical value indicating whether the centered-log-ratio taxa count data are used (LDM-clr). The default is FALSE. 
#' @param comp.anal.adjust a character string that takes value "median" or "mode" to choose the estimator for the beta mean (Hu and Satten, 2023). The default is "median".
#' @return a list consisting of 
#'   \item{x}{the (orthonormal) design matrix X as defined in Hu and Satten (2020)} 
#'   \item{dist}{the (squared/centered) distance matrix} 
#'   \item{mean.freq}{the mean relative abundance of OTUs (the column means of the frequency-scale OTU table)} 
#'   \item{y.freq}{the frequency-scale OTU table, scaled and centered if so specified} 
#'   \item{d.freq}{a vector of the non-negative diagonal elements of \code{D} that satisfies
#'   \code{x^T y.freq = D v^T}}
#'   \item{v.freq}{the v matrix with unit columns that satisfies
#'   \code{x^T y.freq = D v^T}}
#'   \item{y.tran}{the (column-centered) arcsin-root-transformed 
#'   OTU table} 
#'   \item{d.tran}{a vector of the non-negative diagonal elements of \code{D} that satisfies
#'   \code{x^T y.tran = D v^T}}
#'   \item{v.tran}{the v matrix with unit columns that satisfies
#'   \code{x^T y.tran = D v^T}}
#'   \item{low}{a vector of lower indices for confounders (if there is any) and submodels}
#'   \item{up}{a vector of upper indices for confounders (if there is any) and submodels}
#'   \item{beta}{a matrix of effect sizes of every trait on every OTU}
#'   \item{phi}{a matrix of probabilities that the rarefied count of an OTU in a sample is non-zero}
#'   \item{VE.global.freq.confounders}{Variance explained (VE) by confounders, based on the frequency-scale data}
#'   \item{VE.global.freq.submodels}{VE by each submodel, based on the frequency-scale data}
#'   \item{VE.global.freq.residuals}{VE by each component in the residual distance, based on the frequency-scale data}
#'   \item{VE.otu.freq.confounders}{Contribution of each OTU to VE by confounders, based on the frequency-scale data}
#'   \item{VE.otu.freq.submodel}{Contribution of each OTU to VE by each submodel, based on the frequency-scale data}
#'   \item{VE.global.tran.confounders}{VE by confounders, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.global.tran.submodels}{VE by each submodel, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.global.tran.residuals}{VE by each component in the residual distance, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.otu.tran.confounders}{Contribution of each OTU to VE by confounders, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.otu.tran.submodels}{Contribution of each OTU to VE by each submodel, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.df.confounders}{Degree of freedom (i.e., number of components) associated with the VE for confounders}
#'   \item{VE.df.submodels}{Degree of freedom (i.e., number of components) associated with the VE for each submodel}
#'   \item{F.global.freq}{F statistics for testing each submodel, based on
#'   the frequency-scale data} 
#'   \item{F.global.tran}{F statistics for testing each submodel, based on 
#'   the arcsin-root-transformed frequency data} 
#'   \item{F.otu.freq}{F statistics for testing each OTU for each submodel, based on the frequency-scale data} 
#'   \item{F.otu.tran}{F statistics for testing each OTU for each submodel, based on the arcsin-root-transformed data} 
#'   \item{p.global.freq}{p-values for the global test of each set of covariates
#'   based on the frequency-scale data} 
#'   \item{p.global.tran}{p-values for the global test of each set of covariates
#'   based on the arcsin-root-transformed frequency data} 
#'   \item{p.global.pa}{p-values for the global test of each set of covariates
#'   based on the presence-absence data} 
#'   \item{p.global.omni}{p-values for the global test of each set of covariates 
#'   based on the omnibus statistics in LDM-omni, which are the minima of the p-values obtained 
#'   from the frequency scale and the arcsin-root-transformed frequency data 
#'   as the final test statistics, and use the corresponding minima from the 
#'   permuted data to simulate the null distributions} 
#'   \item{p.global.harmonic}{p-values for the global test of each set of covariates
#'   based on the Harmonic-mean p-value combination method applied to the OTU-level omnibus p-values} 
#'   \item{p.global.fisher}{p-values for the global test of each set of covariates
#'   based on the Fisher p-value combination method applied to the OTU-level omnibus p-values} 
#'   \item{p.global.omni3}{p-values for the global test of each set of covariates 
#'   based on the omnibus test LDM-omni3} 
#'   \item{p.global.freq.OR, p.global.tran.OR, p.global.pa.OR, p.global.omni.OR, p.global.harmonic.OR, p.global.fisher.OR, p.global.omni3.OR}{global p-values for testing \code{other.surv.resid}}
#'   \item{p.global.freq.com, p.global.tran.com, p.global.pa.com, p.global.omni.com, p.global.harmonic.com, p.global.fisher.com, p.global.omni3.com}{global p-values from the combination test that combines the results from analyzing both the Martingale and deviance residuals from a Cox model (one of them is supplied by \code{other.surv.resid})}
#'   \item{p.otu.freq}{p-values for the OTU-specific tests based on the 
#'   frequency scale data} 
#'   \item{p.otu.tran}{p-values for the OTU-specific tests based on the 
#'   arcsin-root-transformed frequency data}
#'   \item{p.otu.pa}{p-values for the OTU-specific tests based on the 
#'   presence-absence data} 
#'   \item{p.otu.omni}{p-values for the OTU-specific tests based on the 
#'   omnibus test LDM-omni} 
#'   \item{p.otu.omni3}{p-values for the OTU-specific tests based on the 
#'   omnibus test LDM-omni3} 
#'   \item{q.otu.freq}{q-values (i.e., FDR-adjusted p-values) 
#'   for the OTU-specific tests based on the frequency scale data} 
#'   \item{q.otu.tran}{q-values for the OTU-specific tests based on 
#'   the arcsin-root-transformed frequency data} 
#'   \item{q.otu.pa}{q-values (i.e., FDR-adjusted p-values) 
#'   for the OTU-specific tests based on the presence-absence data} 
#'   \item{q.otu.omni}{q-values for the OTU-specific tests based on the 
#'   omnibus test LDM-omni}
#'   \item{q.otu.omni3}{q-values for the OTU-specific tests based on the 
#'   omnibus test LDM-omni3} 
#'   \item{p.otu.freq.OR, p.otu.tran.OR, p.otu.pa.OR, p.otu.omni.OR, p.otu.omni3.OR, q.otu.freq.OR, q.otu.tran.OR, q.otu.pa.OR, q.otu.omni.OR, q.otu.omni3.OR}{OTU-level p-values and q-values for testing \code{other.surv.resid}}
#'   \item{p.otu.freq.com, p.otu.tran.com, p.otu.pa.com, p.otu.omni.com, p.otu.omni3.com, q.otu.freq.com, q.otu.tran.com, q.otu.pa.com, q.otu.omni.com, q.otu.omni3.com}{OTU-level p-values and q-values from the combination tests that combine the results from analyzing both the Martingale and deviance residuals from a Cox model (one of them is supplied by \code{other.surv.resid})}
#'   \item{detected.otu.freq}{detected OTUs (whose names are found in the column names of the OTU table) at the nominal FDR, based on the frequency scale data} 
#'   \item{detected.otu.tran}{detected OTUs based on the arcsin-root-transformed frequency data} 
#'   \item{detected.otu.pa}{detected OTUs based on the presence-absence data} 
#'   \item{detected.otu.omni}{detected OTU based on the 
#'   omnibus test LDM-omni} 
#'   \item{detected.otu.omni3}{detected OTU based on the 
#'   omnibus test LDM-omni3} 
#'   \item{detected.otu.freq.OR, detected.otu.tran.OR, detected.otu.pa.OR, detected.otu.omni.OR, detected.otu.omni3.OR}{detected OTUs for \code{other.surv.resid}}
#'   \item{detected.otu.freq.com, detected.otu.tran.com, detected.otu.pa.com, detected.otu.omni.com, detected.otu.omni3.com}{detected OTUs by the combination tests that combines the Martingale and deviance residuals from a Cox model (one of them is supplied by \code{other.surv.resid})}
#'   \item{med.p.global.freq, med.p.global.tran, med.p.global.omni, med.p.global.pa, med.p.global.harmonic, med.p.global.fisher, med.p.global.omni3}{p-values for the global tests of the overall mediation effect by the microbiome}
#'   \item{med.p.global.freq.OR, med.p.global.tran.OR, med.p.global.omni.OR, med.p.global.pa.OR, med.p.global.harmonic.OR, med.p.global.fisher.OR, med.p.global.omni3.OR}{p-values for the global tests of the overall mediation effect by the microbiome, when the outcome is \code{other.surv.resid}}
#'   \item{med.p.global.freq.com, med.p.global.tran.com, med.p.global.omni.com, med.p.global.pa.com, med.p.global.harmonic.com, med.p.global.fisher.com, med.p.global.omni3.com}{p-values for the global tests of the overall mediation effect by the microbiome, combining the results from analyzing both the Martingale and deviance residuals as outcomes}
#'   \item{med.detected.otu.freq, med.detected.otu.tran, med.detected.otu.omni, med.detected.otu.pa, med.detected.otu.omni3}{detected mediating OTUs}
#'   \item{med.detected.otu.freq.OR, med.detected.otu.tran.OR, med.detected.otu.omni.OR, med.detected.otu.pa.OR, med.detected.otu.omni3.OR}{detected mediating OTUs for the outcome \code{other.surv.resid}}
#'   \item{med.detected.otu.freq.com, med.detected.otu.tran.com, med.detected.otu.omni.com, med.detected.otu.pa.com, med.detected.otu.omni3.com}{detected mediating OTUs, combining the results from analyzing both the Martingale and deviance residuals as outcomes}
#'   \item{n.perm.completed}{number of permutations completed} 
#'   \item{global.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all global tests} 
#'   \item{otu.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all OTU-specific tests}
#'   \item{seed}{the seed that is user supplied or internally generated, stored in case 
#'   the user wants to reproduce the permutation replicates}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gsatten@emory.edu>
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom parallel mclapply
#' @importFrom permute how shuffleSet Plots Within
#' @importFrom utils tail
#' @importFrom stats as.formula cor model.frame model.matrix na.omit na.pass p.adjust pf
#' @import matrixStats
#' @export
#' @references Hu YJ, Satten GA (2020). Testing hypotheses about the microbiome using the linear decomposition model (LDM) 
#'   Bioinformatics, 36(14), 4106-4115.
#' @references Hu YJ, Lane A, and Satten GA (2021). A rarefaction-based extension of the LDM for testing presence-absence associations in the microbiome. Bioinformatics, 37(12):1652-1657.
#' @references Zhu Z, Satten GA, Caroline M, and Hu YJ (2020). Analyzing matched sets of microbiome data using the LDM and PERMANOVA. Microbiome, 9(133), https://doi.org/10.1186/s40168-021-01034-9.
#' @references Zhu Z, Satten GA, and Hu YJ (2022). Integrative analysis of relative abundance data and presence-absence data of the microbiome using the LDM. Bioinformatics, doi.org/10.1093/bioinformatics/btac181.
#' @references Yue Y and Hu YJ (2021) A new approach to testing mediation of the microbiome using the LDM. bioRxiv, https://doi.org/10.1101/2021.11.12.468449.
#' @references Hu Y, Li Y, Satten GA, and Hu YJ (2022) Testing microbiome associations with censored survival outcomes at both the community and individual taxon levels. bioRxiv, doi.org/10.1101/2022.03.11.483858.
#' @examples
#'data(throat.otu.tab5)
#'data(throat.meta)
#'res.ldm <- ldm(formula=throat.otu.tab5 | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
#'               data=throat.meta, seed=67817, fdr.nominal=0.1) 

ldm = function( formula, other.surv.resid=NULL, data=.GlobalEnv, tree=NULL, dist.method="bray", dist=NULL, 
                     cluster.id=NULL, strata=NULL, how=NULL,
                     perm.within.type="free", perm.between.type="none",
                     perm.within.ncol=0, perm.within.nrow=0,
                     n.perm.max=NULL, n.rej.stop=100, seed=NULL, 
                     fdr.nominal=0.1,
                     square.dist=TRUE, 
                     scale.otu.table=TRUE, center.otu.table=TRUE,
                     freq.scale.only=FALSE, binary=FALSE, n.rarefy=0,
                     test.mediation=FALSE,
                     test.omni3=FALSE,
                     comp.anal=FALSE, comp.anal.adjust="median",
                     n.cores=4) {  
    
    #------------------------
    # form.call
    #------------------------
    options(na.action=na.omit) # fixed a bug here
    object=formula
    #
    #   extract cluster.id from dataframe
    #
    cl=match.call()
    mf=match.call(expand.dots=FALSE)
    m=match( x='cluster.id', table=names(mf) )
    mf.string=as.character( mf[c(1L,m)] )
    cluster.name=mf.string[2]
    if (cluster.name=='NULL') {
        cluster.id=NULL
    } else {   
        loc.dollar=utils::tail( gregexpr('\\$', cluster.name)[[1]] , n=1 )
        if (loc.dollar<0)  {
            cluster.id=getElement(data,cluster.name)
            if( is.null(cluster.id) ) cluster.id=get(cluster.name)
        } else {   
            df.name=get( substr(cluster.name, start=1, stop=loc.dollar-1) )
            var.name=substr(cluster.name, start=loc.dollar+1, stop=nchar(cluster.name))            
            cluster.id= getElement(df.name,var.name) 
        }
    }
    #        
    #   extract model from formula    
    #    
    obj=toString(object)
    obj=gsub('\\s','',obj)
    prefix=' ~ + 0 + '
    loc.comma=gregexpr(',',obj)[[1]]
    start.terms=loc.comma[2]
    terms=substr(obj,start=start.terms+1, stop=nchar(obj))
    #
    #   find n.obs and full set of rownames
    #   
    if (class(data)=='data.frame') {
        row.names=rownames(data)
        n.obs=length(row.names)
    } else {   
        df=model.frame( as.formula(paste('~',terms)) , na.action=na.pass )
        row.names=rownames(df)
        n.obs=length(row.names)
    }
    #
    #   check for missing values in cluster.id
    #        
    
    if (is.null(cluster.id)) {
        use.rows=row.names
    } else {   
        use=!is.na(cluster.id)
        use.rows=row.names[use]
    }
    #
    #   check for and extract confounders
    #
    model=list()
    j=1
    loc.bar=regexpr('\\|',obj)[1]
    loc.minus=regexpr('-',obj)[1]
    loc.delim=max( loc.bar, loc.minus)
    if (loc.delim>0) {
        end.confound=loc.comma[2]
        c=substr(obj,start=loc.delim+1, stop=end.confound-1)
        conf=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        model[[j]]=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        #       use.rows=intersect( use.rows, rownames(conf) )
        use.rows=rownames(model[[1]]) 
        j=j+1
    } else {
        conf=NULL
    }     
    #
    #   extract model terms
    #
    #   j=1
    continue=TRUE
    while (continue) {
        if (substr(terms,1,1)=='(') {
            stop=regexpr(')\\+',terms)[1]
        } else {
            stop=regexpr('\\+',terms)[1] - 1
        }          
        
        if (stop<=0) stop=nchar(terms) 
        m=substr(terms, start=1, stop=stop)
        model[[j]]=model.matrix( as.formula( paste(prefix,m) ) , data=data)
        use.rows=intersect( use.rows, rownames(model[[j]]) )
        #        if (j==1) {
        #            use.rows=rownames(model[[1]])
        #            }
        #        else {
        #            use.rows=intersect( use.rows, rownames(model[[j]]) )
        #            }         
        if (stop+2<=nchar(terms)) {
            terms=substr(terms, start=stop+2, stop=nchar(terms))
            j=j+1
        } else {
            continue=FALSE
        }             
    }   
    n.model=j    
    #
    #  extract OTU table
    #      
    if (is.null(conf)) loc.delim=loc.comma[2]
    otu.name=substr(obj, start=loc.comma[1]+1, stop=loc.delim-1)
    #   loc.dollar=regexpr('\\$', otu.name)[1]
    loc.dollar=utils::tail( gregexpr('\\$', otu.name)[[1]] , n=1 )
    if (loc.dollar<0)  {
        if (class(data)=='data.frame') {
            otu.table=getElement(data, otu.name)
            if (is.null(otu.table)) otu.table= get(otu.name) 
            otu.table=as.matrix(otu.table)
        } else {
            otu.table=as.matrix( get(otu.name) )
        }
    } else {
        df.name=get( substr(otu.name, start=1, stop=loc.dollar-1) )
        var.name=substr(otu.name, start=loc.dollar+1, stop=nchar(otu.name))
        otu.table=as.matrix( getElement(df.name,var.name) )
    }        
    
    #---------------------------------------
    # checking negative values in otu.table
    #---------------------------------------
    
    neg.exist = any(otu.table<0)
    if (neg.exist) {
        if (scale.otu.table == TRUE) {
            stop("The OTU table has negative values, so it does not make sense to use 'scale.otu.table=TRUE'")
        }
        if (freq.scale.only == FALSE) {
            stop("The OTU table has negative values, which cannot be arcsin-root transformed by 'freq.scale.only=FALSE'")
        }
    }
    
    #    if (is.null(otu.table)) otu.table=as.matrix( getElement(.GlobalEnv,otu.name) )
    if ( nrow(otu.table) != n.obs ) {
        if (ncol(otu.table)==n.obs ) {
            otu.table=t(otu.table)
        } else {   
            stop('OTU table and covariates have different number of observations')
        }
    }   
    
    if (!is.null(dist)) {
        dist <- as.matrix(dist)
        if (dim(otu.table)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between the OTU table and dist' )
    }
    
    #
    #   remove rows having NA 
    #    
    for (j in 1:n.model) {
        keep =  rownames( model[[j]] ) %in% use.rows
        model[[j]]=model[[j]][keep,,drop=FALSE]
    }
    if (!is.null(conf)) {
        keep =  rownames(conf) %in% use.rows 
        conf=conf[keep,,drop=FALSE]
    }
    keep=row.names %in% use.rows    
    otu.table=otu.table[keep,,drop=FALSE]    
    if (!is.null(dist)) dist=dist[keep,keep]
    if (!is.null(cluster.id)) cluster.id=cluster.id[keep]
    
    # transpose
    if (dim(model[[1]])[1] != dim(otu.table)[1]) 
        otu.table <- t(otu.table)
    if (dim(model[[1]])[1] != dim(otu.table)[1]) stop( 'numbers of observations mismatch between covariates and the OTU table' )
    
    # OTU names
    if (is.null(colnames(otu.table))) { 
        colnames(otu.table) = 1:ncol(otu.table)
    }
    
    # remove zero OTUs
    w = which(colSums(abs(otu.table))>0)
    if (length(w) < ncol(otu.table)) {
        warning(paste(ncol(otu.table)-length(w), 'OTU(s) with zero counts in all samples are removed', sep=" "))
        otu.table = otu.table[,w,drop=FALSE]
    }
    
    
    rowsum = rowSums(otu.table)
    if (min(rowsum)==0) {
        warning("There exists sample(s) with zero reads at every OTU!")
    }
    
    
    #------------------------
    # setup permutation
    #------------------------
    
    if (class(how)=='how') {
        CTRL=how                   # user-provided how list
    } else {
        if (is.null(cluster.id)) {
            if (is.null(perm.within.type) & is.null(perm.between.type)) {
                # default when no unclustered data has no type specified is 'free'
                perm.within.type='free'    
            }
            if (is.null(strata)) {
                # setup for unclustered permutation
                CTRL = permute::how( within=permute::Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))  
            } else {
                # setup for unclustered, stratified permutation
                strata=as.factor(strata)
                CTRL = permute::how( blocks=strata, within=permute::Within(type=perm.within.type, 
                                                         nrow=perm.within.nrow, 
                                                         ncol=perm.within.ncol))  
            }    
        } else {        
            cluster.id=as.factor(cluster.id)
            if (is.null(strata)) {            
                #  clustered but unstratified data
                CTRL = permute::how( plots=permute::Plots(cluster.id, type=perm.between.type ), 
                            within=permute::Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            } else {
                #   clustered and stratified data
                strata=as.factor(strata)             
                CTRL = permute::how( blocks=strata, 
                            plots=permute::Plots(cluster.id, type=perm.between.type ), 
                            within=permute::Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            }
        }
    }    
    
    
    #------------------------
    # deciding methods
    #------------------------
    
    all.rarefy = (tolower(n.rarefy)=="all")
    no.rarefy = (n.rarefy==0)
    if (no.rarefy | all.rarefy) n.rarefy=1
    
    if (all.rarefy) {
        scale.otu.table = FALSE
        center.otu.table = FALSE
        freq.scale.only=TRUE
        binary=TRUE
    } 
    
    if (binary) {
        scale.otu.table=FALSE
        freq.scale.only=TRUE
        
        if(all(otu.table>0)) stop("All read counts are greater than zero, so the presence-absence analysis requested by 'binary=TRUE' cannot be performed!")
    }
    
    #------------------------
    # setup model
    #------------------------
    
    OR = other.surv.resid # e.g., Deviance residual
    
    n.obs = dim(model[[1]])[1]
    n.otu = ncol(otu.table)
    
    adjust.for.confounders = !is.null(conf)
    n.var = length(model)
    
    if (!all.rarefy | test.omni3) {
        
        n.var1 = n.var - as.numeric(adjust.for.confounders)
        
        center.vars=center.otu.table
        
        index = rep(0, n.var)
        
        for (i in 1:n.var) {
            m.i = model[[i]]
            if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
            
            if (i==1) {
                m = m.i
                index[i] = dim(m.i)[2] 
            } else {
                m = cbind(m, m.i)   
                index[i] = index[i-1] + dim(m.i)[2]    
            }
        }
        
        if (!is.null(OR)) {
            m.OR = m
            if (center.vars) OR = scale( OR, center=TRUE, scale=FALSE )
            m.OR[,ncol(m.OR)] = OR
        }
    }
    
    if (all.rarefy | test.omni3) {
        
        model.pa = model
        adjust.for.confounders.pa = adjust.for.confounders
        
        if (adjust.for.confounders.pa) {
            model.pa[[1]] = cbind(rep(1, n.obs), model.pa[[1]])
        } else {
            adjust.for.confounders.pa = TRUE
            for (i in n.var:1) {
                model.pa[[i+1]] = model.pa[[i]]
            }
            model.pa[[1]] = matrix(1, nrow=n.obs, ncol=1)
        }
        n.var.pa = length(model.pa)
        n.var1 = n.var.pa - as.numeric(adjust.for.confounders.pa)
        
        index.pa = rep(0, n.var.pa)
        
        for (i in 1:n.var.pa) {
            m.i.pa = model.pa[[i]]
            
            if (i==1) {
                m.pa = m.i.pa
                index.pa[i] = dim(m.i.pa)[2] 
            } else {
                m.pa = cbind(m.pa, m.i.pa)   
                index.pa[i] = index.pa[i-1] + dim(m.i.pa)[2]    
            }
        }
        if (!is.null(OR)) {
            m.pa.OR = m.pa
            if (center.vars) OR = scale( OR, center=TRUE, scale=FALSE )
            m.pa.OR[,ncol(m.pa.OR)] = OR
        }
    }    
    
    #---------------------
    # rarefaction or not?
    #---------------------
    
    if (is.null(seed)) {
        seed = sample(1:10^6, 1)
    }
    set.seed(seed)
    
    
    ldm.obs.freq = NULL
    ldm.obs.tran = NULL
    ldm.obs.pa = NULL
    ldm.obs.freq.OR = NULL
    ldm.obs.tran.OR = NULL
    ldm.obs.pa.OR = NULL
    
    if (!all.rarefy | test.omni3) {
        ss.tot.freq = array(NA, dim=c(n.var1, n.otu, n.rarefy))
        resid.freq = array(NA, dim=c(n.obs, n.otu, n.var1, n.rarefy))
        ss.tot.tran = NULL
        resid.tran = NULL
        if (!freq.scale.only) {
            ss.tot.tran = array(NA, dim=c(n.var1, n.otu, n.rarefy))
            resid.tran = array(NA, dim=c(n.obs, n.otu, n.var1, n.rarefy))
        }
        if (!is.null(OR)) {
            ss.tot.freq.OR = array(NA, dim=c(n.var1, n.otu, n.rarefy))
            resid.freq.OR = array(NA, dim=c(n.obs, n.otu, n.var1, n.rarefy))
            ss.tot.tran.OR = NULL
            resid.tran.OR = NULL
            if (!freq.scale.only) {
                ss.tot.tran.OR = array(NA, dim=c(n.var1, n.otu, n.rarefy))
                resid.tran.OR = array(NA, dim=c(n.obs, n.otu, n.var1, n.rarefy))
            }
        }
        
        for (r in 1:n.rarefy) {
            
            if (!no.rarefy) {
                otu.rarefy= Rarefy(otu.table)$otu.tab.rff
            } else {
                otu.rarefy = otu.table
            }
            if (binary) {
                otu.rarefy = (otu.rarefy>0)*1
            } else {
                const <- max(abs(rowSums(otu.rarefy)))
                if (const > 1e+08) otu.rarefy = otu.rarefy/const
            }
            
            #------------------------
            # dist matrix
            #------------------------
            
            if (is.null(dist)) {
                if (ncol(otu.rarefy)==1) {
                    dist <- diag(nrow(otu.rarefy))
                } else {
                    dist <- calculate.dist(dist.method=dist.method, otu.table=otu.rarefy, tree=tree, scale.otu.table=scale.otu.table, binary=binary)
                }
            }
            d.gower <- gower(d=dist, square=square.dist, center=center.otu.table)
            
            #------------------------
            # data matrix y.freq, y.tran
            #------------------------
            
            if (scale.otu.table) {
                rowsum = rowSums(otu.rarefy)
                rowsum[which(rowsum==0)] = 1
                freq.table <- t( scale( t(otu.rarefy), center=FALSE, scale=rowsum ) )
            } else {
                freq.table <- otu.rarefy
            }
            mean.freq <- colMeans(freq.table)
            y.freq <- scale( freq.table, center=center.otu.table, scale=FALSE )
            
            y.tran <- NULL
            if (!freq.scale.only) {
                theta <- asin(sqrt(freq.table))
                y.tran <- scale( theta, center=center.otu.table, scale=FALSE)
            }
            
            if (r==1) {
                mean.gower <- mean(abs(d.gower))
                mean.yfreq <- mean(abs(y.freq))
                max.center.gower <- max(abs(rowMeans(d.gower)))
                max.center.yfreq  <- max(abs(colMeans(y.freq)))
                if ( (max.center.gower - 10^-8*mean.gower) * (max.center.yfreq - 10^-8*mean.yfreq) < 0) {
                    stop( 'discordant centering of the OTU table and distance matrix' )
                }
            }
            
            #---------------------
            # model fitting
            #---------------------
            
            fit.ldm = calculate.x.and.resid( d.gower=d.gower, y.freq=y.freq, y.tran=y.tran, 
                                             index=index, m=m, adjust.for.confounders=adjust.for.confounders)  
            
            ss.tot.freq[,,r] = fit.ldm$ss.tot.freq
            resid.freq[,,,r] = fit.ldm$resid.freq
            if (!freq.scale.only) {
                ss.tot.tran[,,r] = fit.ldm$ss.tot.tran
                resid.tran[,,,r] = fit.ldm$resid.tran
            }
            if (!is.null(OR)) {
                fit.ldm.OR = calculate.x.and.resid( d.gower=d.gower, y.freq=y.freq, y.tran=y.tran, 
                                                    index=index, m=m.OR, adjust.for.confounders=adjust.for.confounders)  
                ss.tot.freq.OR[,,r] = fit.ldm.OR$ss.tot.freq
                resid.freq.OR[,,,r] = fit.ldm.OR$resid.freq
                if (!freq.scale.only) {
                    ss.tot.tran.OR[,,r] = fit.ldm.OR$ss.tot.tran
                    resid.tran.OR[,,,r] = fit.ldm.OR$resid.tran
                }
            }
            
            if (r==1) {
                x.design = fit.ldm$x
                if (!is.null(OR)) x.design.OR = fit.ldm.OR$x
                low = fit.ldm$low
                up = fit.ldm$up
                
                ndf = fit.ldm$ndf
            }
            
        }# rarefaction
        
        
        #---------------------
        # observed statistic
        #---------------------
        
        ldm.obs.freq = ldm.stat(x=x.design, low=low, up=up, resid=resid.freq, ss.tot=ss.tot.freq, adjust.for.confounders=adjust.for.confounders, comp.anal=comp.anal, comp.anal.adjust=comp.anal.adjust)
        if (!freq.scale.only) ldm.obs.tran = ldm.stat(x=x.design, low=low, up=up, resid=resid.tran, ss.tot=ss.tot.tran, adjust.for.confounders=adjust.for.confounders)
        if (!is.null(OR)) {
            ldm.obs.freq.OR = ldm.stat(x=x.design.OR, low=low, up=up, resid=resid.freq.OR, ss.tot=ss.tot.freq.OR, adjust.for.confounders=adjust.for.confounders, comp.anal=comp.anal, comp.anal.adjust=comp.anal.adjust)
            if (!freq.scale.only) ldm.obs.tran.OR = ldm.stat(x=x.design.OR, low=low, up=up, resid=resid.tran.OR, ss.tot=ss.tot.tran.OR, adjust.for.confounders=adjust.for.confounders)
        }
        
    } # if (!all.rarefy | test.omni3)
    
    if (all.rarefy | test.omni3) { 
        
        d.gower = NULL
        y.freq = NULL
        y.tran = NULL
        mean.freq = NULL
        
        #---------------------
        # model fitting
        #---------------------
        
        fit.ldm.pa = calculate.x.and.resid.allrarefy( y=otu.table, index=index.pa, m=m.pa, adjust.for.confounders=adjust.for.confounders.pa)  
        if (!is.null(OR)) fit.ldm.pa.OR = calculate.x.and.resid.allrarefy( y=otu.table, index=index.pa, m=m.pa.OR, adjust.for.confounders=adjust.for.confounders.pa)  
        
        #---------------------
        # observed statistic
        #---------------------
        
        ldm.obs.pa = ldm.stat.allrarefy(x=fit.ldm.pa$x, low=fit.ldm.pa$low, up=fit.ldm.pa$up, 
                                        resid=fit.ldm.pa$resid, ss.tot=fit.ldm.pa$ss.tot, 
                                        P.resid=fit.ldm.pa$P.resid, ss.tot.1=fit.ldm.pa$ss.tot.1, 
                                        phi_1phi=fit.ldm.pa$phi_1phi,
                                        adjust.for.confounders=adjust.for.confounders.pa)
        if (!is.null(OR)) {
            ldm.obs.pa.OR = ldm.stat.allrarefy(x=fit.ldm.pa.OR$x, low=fit.ldm.pa.OR$low, up=fit.ldm.pa.OR$up, 
                                               resid=fit.ldm.pa.OR$resid, ss.tot=fit.ldm.pa.OR$ss.tot, 
                                               P.resid=fit.ldm.pa.OR$P.resid, ss.tot.1=fit.ldm.pa.OR$ss.tot.1, 
                                               phi_1phi=fit.ldm.pa.OR$phi_1phi,
                                               adjust.for.confounders=adjust.for.confounders.pa)
        }
    }  # if (all.rarefy | test.omni3)
    
    p.otu.freq = NULL
    p.otu.tran = NULL
    p.otu.omni = NULL
    p.otu.pa = NULL
    p.otu.omni3 = NULL
    q.otu.freq = NULL
    q.otu.tran = NULL
    q.otu.omni = NULL
    q.otu.pa = NULL
    q.otu.omni3 = NULL
    
    p.otu.freq.OR = NULL
    p.otu.tran.OR = NULL
    p.otu.omni.OR = NULL
    p.otu.pa.OR = NULL
    p.otu.omni3.OR = NULL
    q.otu.freq.OR = NULL
    q.otu.tran.OR = NULL
    q.otu.omni.OR = NULL
    q.otu.pa.OR = NULL
    q.otu.omni3.OR = NULL
    
    p.otu.freq.com = NULL
    p.otu.tran.com = NULL
    p.otu.omni.com = NULL
    p.otu.pa.com = NULL
    p.otu.omni3.com = NULL
    q.otu.freq.com = NULL
    q.otu.tran.com = NULL
    q.otu.omni.com = NULL
    q.otu.pa.com = NULL
    q.otu.omni3.com = NULL
    
    
    p.global.freq = NULL
    p.global.tran = NULL
    p.global.pa = NULL
    p.global.harmonic = NULL
    p.global.fisher = NULL
    p.global.omni = NULL
    p.global.omni3 = NULL
    
    p.global.freq.OR = NULL
    p.global.tran.OR = NULL
    p.global.pa.OR = NULL
    p.global.harmonic.OR = NULL
    p.global.fisher.OR = NULL
    p.global.omni.OR = NULL
    p.global.omni3.OR = NULL
    
    p.global.freq.com = NULL
    p.global.tran.com = NULL
    p.global.pa.com = NULL
    p.global.harmonic.com = NULL
    p.global.fisher.com = NULL
    p.global.omni.com = NULL
    p.global.omni3.com = NULL
    
    
    med.p.global.freq = NULL
    med.p.global.tran = NULL
    med.p.global.omni = NULL
    med.p.global.pa = NULL
    med.p.global.harmonic = NULL
    med.p.global.fisher = NULL
    med.p.global.omni3 = NULL
    
    med.p.global.freq.OR = NULL
    med.p.global.tran.OR = NULL
    med.p.global.omni.OR = NULL
    med.p.global.pa.OR = NULL
    med.p.global.harmonic.OR = NULL
    med.p.global.fisher.OR = NULL
    med.p.global.omni3.OR = NULL
    
    med.p.global.freq.com = NULL
    med.p.global.tran.com = NULL
    med.p.global.omni.com = NULL
    med.p.global.pa.com = NULL
    med.p.global.harmonic.com = NULL
    med.p.global.fisher.com = NULL
    med.p.global.omni3.com = NULL
    
    
    n.perm.completed = NULL
    n.global.perm.completed = NULL
    n.otu.perm.completed = NULL
    
    global.tests.stopped = NULL
    otu.tests.stopped    = NULL
    
    
    if (ifelse(is.null(n.perm.max), 1, (n.perm.max>0))) {
        
        
        ##############################################################################
        
        parallel.perm <- function(i) {
            
            ldm.perm.freq = NULL
            ldm.perm.tran = NULL
            ldm.perm.pa = NULL
            
            ldm.perm.freq.OR = NULL
            ldm.perm.tran.OR = NULL
            ldm.perm.pa.OR = NULL
            
            if (!all.rarefy | test.omni3) {
                ldm.perm.freq = ldm.stat(x=x.design[perm[i,], ], low=low, up=up, resid=resid.freq[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.freq[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders, comp.anal=comp.anal, comp.anal.adjust=comp.anal.adjust, comp.effect=comp.effect.iperm[,i,drop=FALSE])
                if (!freq.scale.only) ldm.perm.tran = ldm.stat(x=x.design[perm[i,], ], low=low, up=up, resid=resid.tran[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.tran[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders)
                if (!is.null(OR)) {
                    ldm.perm.freq.OR = ldm.stat(x=x.design.OR[perm[i,], ], low=low, up=up, resid=resid.freq.OR[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.freq.OR[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders, comp.anal=comp.anal, comp.anal.adjust=comp.anal.adjust, comp.effect=comp.effect.iperm.OR[,i,drop=FALSE])
                    if (!freq.scale.only) ldm.perm.tran.OR = ldm.stat(x=x.design.OR[perm[i,], ], low=low, up=up, resid=resid.tran.OR[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.tran.OR[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders)
                }
            } 
            
            if (all.rarefy | test.omni3) {
                ldm.perm.pa = ldm.stat.allrarefy(x=fit.ldm.pa$x[perm[i,], ], low=fit.ldm.pa$low, up=fit.ldm.pa$up, 
                                                 resid=fit.ldm.pa$resid[,otu.smallp,,drop=FALSE], ss.tot=fit.ldm.pa$ss.tot[,otu.smallp,drop=FALSE], 
                                                 P.resid=fit.ldm.pa$P.resid, ss.tot.1=fit.ldm.pa$ss.tot.1[,otu.smallp,drop=FALSE], 
                                                 phi_1phi=fit.ldm.pa$phi_1phi[,otu.smallp,drop=FALSE],
                                                 adjust.for.confounders=adjust.for.confounders.pa)
                if (!is.null(OR)) {
                    ldm.perm.pa.OR = ldm.stat.allrarefy(x=fit.ldm.pa.OR$x[perm[i,], ], low=fit.ldm.pa.OR$low, up=fit.ldm.pa.OR$up, 
                                                        resid=fit.ldm.pa.OR$resid[,otu.smallp,,drop=FALSE], ss.tot=fit.ldm.pa.OR$ss.tot[,otu.smallp,drop=FALSE], 
                                                        P.resid=fit.ldm.pa.OR$P.resid, ss.tot.1=fit.ldm.pa.OR$ss.tot.1[,otu.smallp,drop=FALSE], 
                                                        phi_1phi=fit.ldm.pa.OR$phi_1phi[,otu.smallp,drop=FALSE],
                                                        adjust.for.confounders=adjust.for.confounders.pa)
                }
            }
            
            list(ldm.perm.freq$ve.global, ldm.perm.freq$ve.otu, 
                 ldm.perm.tran$ve.global, ldm.perm.tran$ve.otu, 
                 ldm.perm.pa$ve.global, ldm.perm.pa$ve.otu,
                 ldm.perm.freq.OR$ve.global, ldm.perm.freq.OR$ve.otu, 
                 ldm.perm.tran.OR$ve.global, ldm.perm.tran.OR$ve.otu, 
                 ldm.perm.pa.OR$ve.global, ldm.perm.pa.OR$ve.otu,
                 ldm.perm.freq$comp.effect, ldm.perm.freq.OR$comp.effect)
            
        } # parallel.perm
        
        ##############################################################################
        
        
        tol.eq = 10^-8
        
        n.global.perm.max = 5000
        n.global.perm.min = 1000
        if (test.mediation) {
            n.global.perm.max = 10000
            n.global.perm.min = 10000
        }
        n.perm.step = 1000
        
        
        if (is.null(n.perm.max)) {
            n.otu.perm.max = n.otu * n.rej.stop * (1/fdr.nominal)
            n.perm.max = max(n.global.perm.max, n.otu.perm.max, na.rm=TRUE)
        } else {
            n.global.perm.max = n.perm.max
            n.otu.perm.max = n.perm.max
        }
        
        n.perm.completed = 0
        n.global.perm.completed = n.global.perm.max # read max if not stop early
        n.otu.perm.completed = n.otu.perm.max
        
        
        otu.smallp = 1:n.otu
        n.otu.smallp = n.otu
        
        
        # saving global stat
        global.tests.stopped = FALSE
        global.tests.done = FALSE # stopped or reach max
        
        if (!all.rarefy | test.omni3) {
            global.freq = ldm.obs.freq$ve.global
            global.freq.perm = array(NA, dim=c(n.var1, n.rarefy, n.global.perm.max))
            n.global.freq = 0
            if (!freq.scale.only) {
                global.tran = ldm.obs.tran$ve.global
                global.tran.perm = array(NA, dim=c(n.var1, n.rarefy, n.global.perm.max))
                n.global.tran = 0
            }
            if (!is.null(OR)) {
                global.freq.OR = ldm.obs.freq.OR$ve.global
                global.freq.perm.OR = array(NA, dim=c(n.var1, n.rarefy, n.global.perm.max))
                n.global.freq.OR = 0
                if (!freq.scale.only) {
                    global.tran.OR = ldm.obs.tran.OR$ve.global
                    global.tran.perm.OR = array(NA, dim=c(n.var1, n.rarefy, n.global.perm.max))
                    n.global.tran.OR = 0
                }
            }
        }
        
        if (all.rarefy | test.omni3) {
            global.pa = ldm.obs.pa$ve.global
            global.pa.perm = array(NA, dim=c(n.var1, n.rarefy, n.global.perm.max))
            n.global.pa = 0
            if (!is.null(OR)) {
                global.pa.OR = ldm.obs.pa.OR$ve.global
                global.pa.perm.OR = array(NA, dim=c(n.var1, n.rarefy, n.global.perm.max))
                n.global.pa.OR = 0
            }
        }
        
        p.global.freq.tmp = NULL
        p.global.tran.tmp = NULL
        p.global.freq.null = NULL
        p.global.tran.null = NULL
        p.global.freq.tmp.OR = NULL
        p.global.tran.tmp.OR = NULL
        p.global.freq.null.OR = NULL
        p.global.tran.null.OR = NULL
        p.global.pa.tmp = NULL
        p.global.pa.null = NULL
        p.global.pa.tmp.OR = NULL
        p.global.pa.null.OR = NULL
        
        if (test.mediation) {
            if (!all.rarefy | test.omni3) {
                med.n.global.freq = 0
                if (!freq.scale.only) med.n.global.tran = 0 
                if (!is.null(OR)) {
                    med.n.global.freq.OR = 0
                    if (!freq.scale.only) med.n.global.tran.OR = 0 
                }
            }
            if (all.rarefy | test.omni3) {
                med.n.global.pa = 0
                if (!is.null(OR)) med.n.global.pa.OR = 0
            }
            
            med.p.global.freq.tmp = NULL
            med.p.global.tran.tmp = NULL
            med.p.global.freq.null = NULL
            med.p.global.tran.null = NULL
            med.p.global.freq.tmp.OR = NULL
            med.p.global.tran.tmp.OR = NULL
            med.p.global.freq.null.OR = NULL
            med.p.global.tran.null.OR = NULL
            med.p.global.pa.tmp = NULL
            med.p.global.pa.null = NULL
            med.p.global.pa.tmp.OR = NULL
            med.p.global.pa.null.OR = NULL
        }
        
        # saving OTU stat
        otu.tests.stopped = FALSE
        
        if (!all.rarefy | test.omni3) {
            otu.freq = ldm.obs.freq$ve.otu
            otu.freq.perm = array(NA, c(n.var1, n.otu, n.global.perm.max))
            
            comp.effect.iperm = NULL
            comp.effect.iperm.OR = NULL
            if (comp.anal) {
                comp.effect.perm = matrix(NA, up[n.var], n.global.perm.max)
                if (!is.null(OR)) comp.effect.perm.OR = matrix(NA, up[n.var], n.global.perm.max)
                last = 1
            }
            
            n.otu.freq = 0
            Aset.freq = matrix(TRUE, n.var1, n.otu)
            p.otu.freq = matrix(NA, n.var1, n.otu)
            
            if (!freq.scale.only) {
                otu.tran = ldm.obs.tran$ve.otu
                otu.tran.perm = array(NA, c(n.var1, n.otu, n.global.perm.max))
                n.otu.tran = 0
                Aset.tran = matrix(TRUE, n.var1, n.otu)
                p.otu.tran = matrix(NA, n.var1, n.otu)
                
                Aset.omni = matrix(TRUE, n.var1, n.otu)
                p.otu.omni = matrix(NA, n.var1, n.otu)
            }
            
            if (!is.null(OR)) {
                otu.freq.OR = ldm.obs.freq.OR$ve.otu
                otu.freq.perm.OR = array(NA, c(n.var1, n.otu, n.global.perm.max))
                n.otu.freq.OR = 0
                Aset.freq.OR = matrix(TRUE, n.var1, n.otu)
                p.otu.freq.OR = matrix(NA, n.var1, n.otu)
                
                Aset.freq.com = matrix(TRUE, n.var1, n.otu)
                p.otu.freq.com = matrix(NA, n.var1, n.otu)
                
                if (!freq.scale.only) {
                    otu.tran.OR = ldm.obs.tran.OR$ve.otu
                    otu.tran.perm.OR = array(NA, c(n.var1, n.otu, n.global.perm.max))
                    n.otu.tran.OR = 0
                    Aset.tran.OR = matrix(TRUE, n.var1, n.otu)
                    p.otu.tran.OR = matrix(NA, n.var1, n.otu)
                    
                    Aset.tran.com = matrix(TRUE, n.var1, n.otu)
                    p.otu.tran.com = matrix(NA, n.var1, n.otu)
                    
                    Aset.omni.OR = matrix(TRUE, n.var1, n.otu)
                    p.otu.omni.OR = matrix(NA, n.var1, n.otu)
                    
                    Aset.omni.com = matrix(TRUE, n.var1, n.otu)
                    p.otu.omni.com = matrix(NA, n.var1, n.otu)
                }
            }
        }
        
        if (all.rarefy | test.omni3) {
            otu.pa = ldm.obs.pa$ve.otu
            otu.pa.perm = array(NA, c(n.var1, n.otu, n.global.perm.max))
            n.otu.pa = 0
            Aset.pa = matrix(TRUE, n.var1, n.otu)
            p.otu.pa = matrix(NA, n.var1, n.otu)
            
            if (!is.null(OR)) {
                otu.pa.OR = ldm.obs.pa.OR$ve.otu
                otu.pa.perm.OR = array(NA, c(n.var1, n.otu, n.global.perm.max))
                n.otu.pa.OR = 0
                Aset.pa.OR = matrix(TRUE, n.var1, n.otu)
                p.otu.pa.OR = matrix(NA, n.var1, n.otu)
                
                Aset.pa.com = matrix(TRUE, n.var1, n.otu)
                p.otu.pa.com = matrix(NA, n.var1, n.otu)
            }
        }
        
        if (test.omni3) {
            Aset.omni3 = matrix(TRUE, n.var1, n.otu)
            p.otu.omni3 = matrix(NA, n.var1, n.otu)
            if (!is.null(OR)) {
                Aset.omni3.OR = matrix(TRUE, n.var1, n.otu)
                p.otu.omni3.OR = matrix(NA, n.var1, n.otu)   
                Aset.omni3.com = matrix(TRUE, n.var1, n.otu)
                p.otu.omni3.com = matrix(NA, n.var1, n.otu)   
            }
        }
        
        # checking if the process is stabilized
        
        n.stable.max = 10
        if (!all.rarefy | test.omni3) {
            n.stable.freq = 0
            if (!is.null(OR)) {
                n.stable.freq.OR = 0
                n.stable.freq.com = 0
            }
            if (!freq.scale.only) {
                n.stable.tran = 0
                n.stable.omni = 0
                if (!is.null(OR)) {
                    n.stable.tran.OR = 0
                    n.stable.omni.OR = 0
                    n.stable.tran.com = 0
                    n.stable.omni.com = 0
                }
            }
        }
        if (all.rarefy | test.omni3) {
            n.stable.pa = 0
            if (!is.null(OR)) {
                n.stable.pa.OR = 0
                n.stable.pa.com = 0
            }
        }
        if (test.omni3) {
            n.stable.omni3 = 0
            if (!is.null(OR)) {
                n.stable.omni3.OR = 0
                n.stable.omni3.com = 0   
            }
        }
        
        n.perm.block = 1000
        nblock = ceiling(n.perm.max/n.perm.block)
        i.block = 0
        
        i.sim.low = 0
        i.sim.up = 0
        
        set.seed(seed) 
        
        while (n.perm.completed < n.perm.max) {
            
            i.block = i.block + 1
            
            if (n.perm.completed == 100000) {
                n.perm.block = n.perm.block*10
            }
            
            #######################   reduce OTUs to speed-up computation   ########################
            
            if (global.tests.done & (n.perm.completed %% n.perm.step == 0)) {
                
                if (test.omni3) {
                    column.sums = colSums(Aset.freq + Aset.tran + Aset.omni + Aset.pa + Aset.omni3)
                    if (!is.null(OR)) column.sums = column.sums + colSums(Aset.freq.OR + Aset.tran.OR + Aset.omni.OR + Aset.pa.OR + Aset.omni3.OR)
                } else if (all.rarefy) {
                    column.sums = colSums(Aset.pa)
                    if (!is.null(OR)) column.sums = column.sums + colSums(Aset.pa.OR)
                } else if (freq.scale.only) {
                    column.sums = colSums(Aset.freq) 
                    if (!is.null(OR)) column.sums = column.sums + colSums(Aset.freq.OR) 
                } else {
                    column.sums = colSums(Aset.freq + Aset.tran + Aset.omni)
                    if (!is.null(OR)) column.sums = column.sums + colSums(Aset.freq.OR + Aset.tran.OR + Aset.omni.OR)
                }
                w.otu.smallp = which(column.sums > 0)
                
                cat("number of OTUs do not meet early stopping criterion:", length(w.otu.smallp), "\n")
                
                if (length(w.otu.smallp) != n.otu.smallp) {
                    
                    otu.smallp = otu.smallp[w.otu.smallp]
                    n.otu.smallp = length(otu.smallp)
                    
                    if (!all.rarefy | test.omni3) {
                        n.otu.freq = n.otu.freq[,w.otu.smallp,drop=FALSE]
                        Aset.freq = Aset.freq[,w.otu.smallp,drop=FALSE]
                        if (!is.null(OR)) {
                            n.otu.freq.OR = n.otu.freq.OR[,w.otu.smallp,drop=FALSE]
                            Aset.freq.OR = Aset.freq.OR[,w.otu.smallp,drop=FALSE]
                            Aset.freq.com = Aset.freq.com[,w.otu.smallp,drop=FALSE]
                        }
                        if (!freq.scale.only) {
                            n.otu.tran = n.otu.tran[,w.otu.smallp,drop=FALSE]
                            Aset.tran = Aset.tran[,w.otu.smallp,drop=FALSE]
                            Aset.omni = Aset.omni[,w.otu.smallp,drop=FALSE]
                            if (!is.null(OR)) {
                                n.otu.tran.OR = n.otu.tran.OR[,w.otu.smallp,drop=FALSE]
                                Aset.tran.OR = Aset.tran.OR[,w.otu.smallp,drop=FALSE]
                                Aset.omni.OR = Aset.omni.OR[,w.otu.smallp,drop=FALSE]      
                                Aset.tran.com = Aset.tran.com[,w.otu.smallp,drop=FALSE]
                                Aset.omni.com = Aset.omni.com[,w.otu.smallp,drop=FALSE]   
                            }
                        }
                    }
                    if (all.rarefy | test.omni3) {
                        n.otu.pa = n.otu.pa[,w.otu.smallp,drop=FALSE]
                        Aset.pa = Aset.pa[,w.otu.smallp,drop=FALSE]
                        if (!is.null(OR)) {
                            n.otu.pa.OR = n.otu.pa.OR[,w.otu.smallp,drop=FALSE]
                            Aset.pa.OR = Aset.pa.OR[,w.otu.smallp,drop=FALSE]
                            Aset.pa.com = Aset.pa.com[,w.otu.smallp,drop=FALSE]
                        }
                    }
                    if (test.omni3) {
                        Aset.omni3 = Aset.omni3[,w.otu.smallp,drop=FALSE]
                        if (!is.null(OR)) {
                            Aset.omni3.OR = Aset.omni3.OR[,w.otu.smallp,drop=FALSE]
                            Aset.omni3.com = Aset.omni3.com[,w.otu.smallp,drop=FALSE]
                        }
                    }
                } # if (length(w.otu.smallp) != n.otu.smallp)
                
                if (!all.rarefy | test.omni3) {
                    tmp = otu.freq.perm[,w.otu.smallp,1:n.perm.completed,drop=FALSE]
                    otu.freq.perm = array(NA, c(n.var1, n.otu.smallp, n.perm.completed+n.perm.block))
                    otu.freq.perm[,,1:n.perm.completed] = tmp
                    if (!is.null(OR)) {
                        tmp = otu.freq.perm.OR[,w.otu.smallp,1:n.perm.completed,drop=FALSE]
                        otu.freq.perm.OR = array(NA, c(n.var1, n.otu.smallp, n.perm.completed+n.perm.block))
                        otu.freq.perm.OR[,,1:n.perm.completed] = tmp
                    }
                    if (!freq.scale.only) {
                        tmp = otu.tran.perm[,w.otu.smallp,1:n.perm.completed,drop=FALSE]
                        otu.tran.perm = array(NA, c(n.var1, n.otu.smallp, n.perm.completed+n.perm.block))
                        otu.tran.perm[,,1:n.perm.completed] = tmp
                        if (!is.null(OR)) {
                            tmp = otu.tran.perm.OR[,w.otu.smallp,1:n.perm.completed,drop=FALSE]
                            otu.tran.perm.OR = array(NA, c(n.var1, n.otu.smallp, n.perm.completed+n.perm.block))
                            otu.tran.perm.OR[,,1:n.perm.completed] = tmp
                        }
                    }
                }
                if (all.rarefy | test.omni3) {
                    tmp = otu.pa.perm[,w.otu.smallp,1:n.perm.completed,drop=FALSE]
                    otu.pa.perm = array(NA, c(n.var1, n.otu.smallp, n.perm.completed+n.perm.block))
                    otu.pa.perm[,,1:n.perm.completed] = tmp
                    if (!is.null(OR)) {
                        tmp = otu.pa.perm.OR[,w.otu.smallp,1:n.perm.completed,drop=FALSE]
                        otu.pa.perm.OR = array(NA, c(n.var1, n.otu.smallp, n.perm.completed+n.perm.block))
                        otu.pa.perm.OR[,,1:n.perm.completed] = tmp
                    }
                }
            } # reduce OTUs
            
            #####################   calculate test statistics   ######################
            ######################   (parallel permutations)   #######################
            
            perm = permute::shuffleSet(n.obs, n.perm.block, CTRL)
            
            if (n.otu.smallp > 20 & n.cores > 1) { 
                
                i.sim.low = i.sim.up + 1
                i.sim.up = i.sim.up + n.perm.block
                i.sim = i.sim.low:i.sim.up
                
                if (comp.anal & n.otu.smallp != n.otu) {
                    i.comp.effect <- sample(1:last, size=n.perm.block)
                    comp.effect.iperm <- comp.effect.perm[,i.comp.effect,drop=FALSE]
                    if (!is.null(OR)) comp.effect.iperm.OR <- comp.effect.perm.OR[,i.comp.effect,drop=FALSE]
                }
                
                if (Sys.info()[['sysname']] == 'Windows') {
                    parallel.stat = BiocParallel::bplapply(1:n.perm.block, parallel.perm, BPPARAM = BiocParallel::MulticoreParam(workers=n.cores))
                } else {
                    parallel.stat = parallel::mclapply(1:n.perm.block, parallel.perm, mc.cores = n.cores)
                }
                
                if (comp.anal & n.otu.smallp == n.otu) {
                    comp.effect.perm[,i.sim] = sapply(parallel.stat, "[[", 13)
                    if (!is.null(OR)) comp.effect.perm.OR[,i.sim] = sapply(parallel.stat, "[[", 14)
                    last = i.sim.up
                } 
                
                if (!global.tests.done) {
                    
                    if (!all.rarefy | test.omni3) {
                        global.freq.block = sapply(parallel.stat, "[[", 1)
                        global.freq.perm[,,i.sim] = array(global.freq.block, dim=c(n.var1, n.rarefy, n.perm.block))
                        n.global.freq = n.global.freq + matrix(rowSums((matrix(global.freq.block, ncol=n.perm.block) > as.vector(global.freq) + tol.eq) + (matrix(global.freq.block, ncol=n.perm.block) > as.vector(global.freq) - tol.eq)), nrow=n.var1)
                        if (!is.null(OR)) {
                            global.freq.block.OR = sapply(parallel.stat, "[[", 7)
                            global.freq.perm.OR[,,i.sim] = array(global.freq.block.OR, dim=c(n.var1, n.rarefy, n.perm.block))
                            n.global.freq.OR = n.global.freq.OR + matrix(rowSums((matrix(global.freq.block.OR, ncol=n.perm.block) > as.vector(global.freq.OR) + tol.eq) + (matrix(global.freq.block.OR, ncol=n.perm.block) > as.vector(global.freq.OR) - tol.eq)), nrow=n.var1)
                        }
                        if (!freq.scale.only) {
                            global.tran.block = sapply(parallel.stat, "[[", 3)
                            global.tran.perm[,,i.sim] = global.tran.block
                            n.global.tran = n.global.tran + matrix(rowSums((matrix(global.tran.block, ncol=n.perm.block) > as.vector(global.tran) + tol.eq) + (matrix(global.tran.block, ncol=n.perm.block) > as.vector(global.tran) - tol.eq)), nrow=n.var1)
                            if (!is.null(OR)) {
                                global.tran.block.OR = sapply(parallel.stat, "[[", 9)
                                global.tran.perm.OR[,,i.sim] = global.tran.block.OR
                                n.global.tran.OR = n.global.tran.OR + matrix(rowSums((matrix(global.tran.block.OR, ncol=n.perm.block) > as.vector(global.tran.OR) + tol.eq) + (matrix(global.tran.block.OR, ncol=n.perm.block) > as.vector(global.tran.OR) - tol.eq)), nrow=n.var1)
                            }
                        }
                    }
                    if (all.rarefy | test.omni3) {
                        global.pa.block = sapply(parallel.stat, "[[", 5)
                        global.pa.perm[,,i.sim] = array(global.pa.block, dim=c(n.var1, n.rarefy, n.perm.block))
                        n.global.pa = n.global.pa + matrix(rowSums((matrix(global.pa.block, ncol=n.perm.block) > as.vector(global.pa) + tol.eq) + (matrix(global.pa.block, ncol=n.perm.block) > as.vector(global.pa) - tol.eq)), nrow=n.var1)
                        if (!is.null(OR)) {
                            global.pa.block.OR = sapply(parallel.stat, "[[", 11)
                            global.pa.perm.OR[,,i.sim] = array(global.pa.block.OR, dim=c(n.var1, n.rarefy, n.perm.block))
                            n.global.pa.OR = n.global.pa.OR + matrix(rowSums((matrix(global.pa.block.OR, ncol=n.perm.block) > as.vector(global.pa.OR) + tol.eq) + (matrix(global.pa.block.OR, ncol=n.perm.block) > as.vector(global.pa.OR) - tol.eq)), nrow=n.var1)
                        }
                    }
                }
                
                if (!otu.tests.stopped | (test.mediation & n.perm.completed < n.global.perm.min)) {
                    
                    if (!all.rarefy | test.omni3) {
                        otu.freq.block = sapply(parallel.stat, "[[", 2)
                        otu.freq.perm[,,i.sim] = array(otu.freq.block, dim=c(n.var1, nrow(otu.freq.block)/n.var1, n.perm.block))
                        n.otu.freq = n.otu.freq + matrix(rowSums((otu.freq.block>as.vector(otu.freq[,otu.smallp,drop=FALSE])+tol.eq) + (otu.freq.block>as.vector(otu.freq[,otu.smallp,drop=FALSE])-tol.eq)), nrow=n.var1)

                        if (!is.null(OR)) {
                            otu.freq.block.OR = sapply(parallel.stat, "[[", 8)
                            otu.freq.perm.OR[,,i.sim] = array(otu.freq.block.OR, dim=c(n.var1, nrow(otu.freq.block.OR)/n.var1, n.perm.block))
                            n.otu.freq.OR = n.otu.freq.OR + matrix(rowSums((otu.freq.block.OR>as.vector(otu.freq.OR[,otu.smallp,drop=FALSE])+tol.eq) + (otu.freq.block.OR>as.vector(otu.freq.OR[,otu.smallp,drop=FALSE])-tol.eq)), nrow=n.var1)
                        }
                        if (!freq.scale.only) {
                            otu.tran.block = sapply(parallel.stat, "[[", 4)
                            otu.tran.perm[,,i.sim] = array(otu.tran.block, dim=c(n.var1, nrow(otu.tran.block)/n.var1, n.perm.block))
                            n.otu.tran = n.otu.tran + matrix(rowSums((otu.tran.block>as.vector(otu.tran[,otu.smallp,drop=FALSE])+tol.eq) + (otu.tran.block>as.vector(otu.tran[,otu.smallp,drop=FALSE])-tol.eq)), nrow=n.var1)
                            if (!is.null(OR)) {
                                otu.tran.block.OR = sapply(parallel.stat, "[[", 10)
                                otu.tran.perm.OR[,,i.sim] = array(otu.tran.block.OR, dim=c(n.var1, nrow(otu.tran.block.OR)/n.var1, n.perm.block))
                                n.otu.tran.OR = n.otu.tran.OR + matrix(rowSums((otu.tran.block.OR>as.vector(otu.tran.OR[,otu.smallp,drop=FALSE])+tol.eq) + (otu.tran.block.OR>as.vector(otu.tran.OR[,otu.smallp,drop=FALSE])-tol.eq)), nrow=n.var1)
                            }
                        }
                    }
                    if (all.rarefy | test.omni3) {
                        otu.pa.block = sapply(parallel.stat, "[[", 6)
                        otu.pa.perm[,,i.sim] = array(otu.pa.block, dim=c(n.var1, nrow(otu.pa.block)/n.var1, n.perm.block))
                        n.otu.pa = n.otu.pa + matrix(rowSums((otu.pa.block>as.vector(otu.pa[,otu.smallp,drop=FALSE])+tol.eq) + (otu.pa.block>as.vector(otu.pa[,otu.smallp,drop=FALSE])-tol.eq)), nrow=n.var1)
                        if (!is.null(OR)) {
                            otu.pa.block.OR = sapply(parallel.stat, "[[", 12)
                            otu.pa.perm.OR[,,i.sim] = array(otu.pa.block.OR, dim=c(n.var1, nrow(otu.pa.block.OR)/n.var1, n.perm.block))
                            n.otu.pa.OR = n.otu.pa.OR + matrix(rowSums((otu.pa.block.OR>as.vector(otu.pa.OR[,otu.smallp,drop=FALSE])+tol.eq) + (otu.pa.block.OR>as.vector(otu.pa.OR[,otu.smallp,drop=FALSE])-tol.eq)), nrow=n.var1)
                        }
                    }
                }
                
            } else { # no parallel computing
                
                for (i in 1:n.perm.block) {
                    
                    i.sim = n.perm.completed + i
                    
                    ldm.perm.tran = NULL
                    
                    if (!all.rarefy | test.omni3) {
                        
                        if (comp.anal & n.otu.smallp != n.otu) {
                            i.comp.effect <- sample(1:last, size=1)
                            comp.effect.iperm <- comp.effect.perm[,i.comp.effect,drop=FALSE]
                            if (!is.null(OR)) comp.effect.iperm.OR <- comp.effect.perm.OR[,i.comp.effect,drop=FALSE]
                        }
                        
                        ldm.perm.freq = ldm.stat(x=x.design[perm[i,], ], low=low, up=up, resid=resid.freq[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.freq[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders, comp.anal=comp.anal, comp.anal.adjust=comp.anal.adjust, comp.effect=comp.effect.iperm)
                        if (!is.null(OR)) ldm.perm.freq.OR = ldm.stat(x=x.design.OR[perm[i,], ], low=low, up=up, resid=resid.freq.OR[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.freq.OR[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders, comp.anal=comp.anal, comp.anal.adjust=comp.anal.adjust, comp.effect=comp.effect.iperm.OR)
                        
                        if (comp.anal & n.otu.smallp == n.otu) {
                            comp.effect.perm[,i.sim] = ldm.perm.freq$comp.effect
                            if (!is.null(OR)) comp.effect.perm.OR[,i.sim] = ldm.perm.freq.OR$comp.effect
                            last = i.sim
                        } 
                        
                        if (!freq.scale.only) {
                            ldm.perm.tran = ldm.stat(x=x.design[perm[i,], ], low=low, up=up, resid=resid.tran[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.tran[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders)
                            if (!is.null(OR)) ldm.perm.tran.OR = ldm.stat(x=x.design.OR[perm[i,], ], low=low, up=up, resid=resid.tran.OR[,otu.smallp,,,drop=FALSE], ss.tot=ss.tot.tran.OR[,otu.smallp,,drop=FALSE], adjust.for.confounders=adjust.for.confounders)
                        } 
                    }
                    if (all.rarefy | test.omni3) {
                        ldm.perm.pa = ldm.stat.allrarefy(x=fit.ldm.pa$x[perm[i,], ]  , low=fit.ldm.pa$low, up=fit.ldm.pa$up, 
                                                         resid=fit.ldm.pa$resid[,otu.smallp,,drop=FALSE], ss.tot=fit.ldm.pa$ss.tot[,otu.smallp,drop=FALSE], 
                                                         P.resid=fit.ldm.pa$P.resid, ss.tot.1=fit.ldm.pa$ss.tot.1[,otu.smallp,drop=FALSE], 
                                                         phi_1phi=fit.ldm.pa$phi_1phi[,otu.smallp,drop=FALSE],
                                                         adjust.for.confounders=adjust.for.confounders.pa)
                        if (!is.null(OR)) ldm.perm.pa.OR = ldm.stat.allrarefy(x=fit.ldm.pa.OR$x[perm[i,], ]  , low=fit.ldm.pa.OR$low, up=fit.ldm.pa.OR$up, 
                                                                              resid=fit.ldm.pa.OR$resid[,otu.smallp,,drop=FALSE], ss.tot=fit.ldm.pa.OR$ss.tot[,otu.smallp,drop=FALSE], 
                                                                              P.resid=fit.ldm.pa.OR$P.resid, ss.tot.1=fit.ldm.pa.OR$ss.tot.1[,otu.smallp,drop=FALSE], 
                                                                              phi_1phi=fit.ldm.pa.OR$phi_1phi[,otu.smallp,drop=FALSE],
                                                                              adjust.for.confounders=adjust.for.confounders.pa)
                    }
                    
                    if (!global.tests.done) {
                        if (!all.rarefy | test.omni3) {
                            global.freq.perm[,,i.sim] = ldm.perm.freq$ve.global
                            n.global.freq = n.global.freq + (ldm.perm.freq$ve.global > global.freq + tol.eq) + (ldm.perm.freq$ve.global > global.freq - tol.eq)
                            if (!is.null(OR)) {
                                global.freq.perm.OR[,,i.sim] = ldm.perm.freq.OR$ve.global
                                n.global.freq.OR = n.global.freq.OR + (ldm.perm.freq.OR$ve.global > global.freq.OR + tol.eq) + (ldm.perm.freq.OR$ve.global > global.freq.OR - tol.eq)
                            }
                            if (!freq.scale.only) {
                                global.tran.perm[,,i.sim] = ldm.perm.tran$ve.global
                                n.global.tran = n.global.tran + (ldm.perm.tran$ve.global > global.tran + tol.eq) + (ldm.perm.tran$ve.global > global.tran - tol.eq)
                                if (!is.null(OR)) {
                                    global.tran.perm.OR[,,i.sim] = ldm.perm.tran.OR$ve.global
                                    n.global.tran.OR = n.global.tran.OR + (ldm.perm.tran.OR$ve.global > global.tran.OR + tol.eq) + (ldm.perm.tran.OR$ve.global > global.tran.OR - tol.eq)
                                }
                            }
                        }
                        if (all.rarefy | test.omni3) {
                            global.pa.perm[,,i.sim] = ldm.perm.pa$ve.global
                            n.global.pa = n.global.pa + (ldm.perm.pa$ve.global > global.pa + tol.eq) + (ldm.perm.pa$ve.global > global.pa - tol.eq)
                            if (!is.null(OR)) {
                                global.pa.perm.OR[,,i.sim] = ldm.perm.pa.OR$ve.global
                                n.global.pa.OR = n.global.pa.OR + (ldm.perm.pa.OR$ve.global > global.pa.OR + tol.eq) + (ldm.perm.pa.OR$ve.global > global.pa.OR - tol.eq)
                            }
                        }
                    }
                    
                    if (!otu.tests.stopped | (test.mediation & i.sim <= n.global.perm.min)) {
                        if (!all.rarefy | test.omni3) {
                            otu.freq.perm[,,i.sim] = ldm.perm.freq$ve.otu
                            n.otu.freq = n.otu.freq + (ldm.perm.freq$ve.otu>otu.freq[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.freq$ve.otu>otu.freq[,otu.smallp,drop=FALSE]-tol.eq)
                            if (!is.null(OR)) {
                                otu.freq.perm.OR[,,i.sim] = ldm.perm.freq.OR$ve.otu
                                n.otu.freq.OR = n.otu.freq.OR + (ldm.perm.freq.OR$ve.otu>otu.freq.OR[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.freq.OR$ve.otu>otu.freq.OR[,otu.smallp,drop=FALSE]-tol.eq)
                            }
                            if (!freq.scale.only) {
                                otu.tran.perm[,,i.sim] = ldm.perm.tran$ve.otu
                                n.otu.tran = n.otu.tran + (ldm.perm.tran$ve.otu>otu.tran[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.tran$ve.otu>otu.tran[,otu.smallp,drop=FALSE]-tol.eq)
                                if (!is.null(OR)) {                                
                                    otu.tran.perm.OR[,,i.sim] = ldm.perm.tran.OR$ve.otu
                                    n.otu.tran.OR = n.otu.tran.OR + (ldm.perm.tran.OR$ve.otu>otu.tran.OR[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.tran.OR$ve.otu>otu.tran.OR[,otu.smallp,drop=FALSE]-tol.eq)
                                }
                            }
                        }
                        if (all.rarefy | test.omni3) {
                            otu.pa.perm[,,i.sim] = ldm.perm.pa$ve.otu
                            n.otu.pa = n.otu.pa + (ldm.perm.pa$ve.otu>otu.pa[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.pa$ve.otu>otu.pa[,otu.smallp,drop=FALSE]-tol.eq)
                            if (!is.null(OR)) {  
                                otu.pa.perm.OR[,,i.sim] = ldm.perm.pa.OR$ve.otu
                                n.otu.pa.OR = n.otu.pa.OR + (ldm.perm.pa.OR$ve.otu>otu.pa.OR[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.pa.OR$ve.otu>otu.pa.OR[,otu.smallp,drop=FALSE]-tol.eq)
                            }
                        }
                    }
                    
                } # for (i in 1:n.perm.block)
            } # no parallel computing
            
            n.perm.completed = n.perm.completed + n.perm.block
            inv.n.perm.completed = 1/n.perm.completed
            inv.n.perm.completed.1 = 1/(n.perm.completed+1)
            
            cat("permutations:", n.perm.completed, "\n")
            
            if (n.perm.completed < n.global.perm.min) next
            if (n.perm.completed >= n.global.perm.min && n.perm.completed %% n.perm.step != 0) next
            
            
            #####################   calculate p value   ######################
            #####################   check if stop early  ######################
            
            # intermediate statistics
            
            p.otu.freq.null = NULL
            p.otu.freq.null.OR = NULL
            p.otu.tran.null = NULL
            p.otu.tran.null.OR = NULL
            p.otu.pa.null = NULL
            p.otu.pa.null.OR = NULL
            
            p.global.freq.tmp = NULL
            p.global.tran.tmp = NULL
            p.global.freq.null = NULL
            p.global.tran.null = NULL
            p.global.freq.tmp.OR = NULL
            p.global.tran.tmp.OR = NULL
            p.global.freq.null.OR = NULL
            p.global.tran.null.OR = NULL
            p.global.pa.tmp = NULL
            p.global.pa.null = NULL
            p.global.pa.tmp.OR = NULL
            p.global.pa.null.OR = NULL
            
            if (test.mediation) {
                med.p.global.freq.tmp = NULL
                med.p.global.tran.tmp = NULL
                med.p.global.freq.null = NULL
                med.p.global.tran.null = NULL
                med.p.global.freq.tmp.OR = NULL
                med.p.global.tran.tmp.OR = NULL
                med.p.global.freq.null.OR = NULL
                med.p.global.tran.null.OR = NULL
                med.p.global.pa.tmp = NULL
                med.p.global.pa.null = NULL
                med.p.global.pa.tmp.OR = NULL
                med.p.global.pa.null.OR = NULL
            }
            
            ################
            # test otu
            ################
            
            if (!all.rarefy | test.omni3) {
                
                if (any(Aset.freq)) {
                    AtoB.freq <- Aset.freq & (n.otu.freq >= n.rej.stop*2)
                    Aset.freq <- Aset.freq & !AtoB.freq
                    p.otu.freq[,otu.smallp][AtoB.freq] <- 0.5*n.otu.freq[AtoB.freq]*inv.n.perm.completed
                    p.otu.freq[,otu.smallp][Aset.freq] <- (0.5*n.otu.freq[Aset.freq]+1)*inv.n.perm.completed.1
                    
                    q.otu.freq <- t(apply(p.otu.freq, 1, fdr.Sandev))
                    if (n.otu == 1) q.otu.freq = matrix(q.otu.freq, ncol=1)
                    
                    Aset.freq.meet.criteria <- rowAlls( ((q.otu.freq[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.freq) | (!Aset.freq) ) 
                    n.stable.freq <- ifelse(Aset.freq.meet.criteria, n.stable.freq + Aset.freq.meet.criteria, 0)
                    Aset.freq.rm.row <- (n.stable.freq >= n.stable.max)                    
                    Aset.freq[Aset.freq.rm.row,] = FALSE
                }
                if (!is.null(OR)) {
                    if (any(Aset.freq.OR)) {
                        AtoB.freq.OR <- Aset.freq.OR & (n.otu.freq.OR >= n.rej.stop*2)
                        Aset.freq.OR <- Aset.freq.OR & !AtoB.freq.OR
                        p.otu.freq.OR[,otu.smallp][AtoB.freq.OR] <- 0.5*n.otu.freq.OR[AtoB.freq.OR]*inv.n.perm.completed
                        p.otu.freq.OR[,otu.smallp][Aset.freq.OR] <- (0.5*n.otu.freq.OR[Aset.freq.OR]+1)*inv.n.perm.completed.1
                        
                        q.otu.freq.OR <- t(apply(p.otu.freq.OR, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.freq.OR = matrix(q.otu.freq.OR, ncol=1)
                        
                        Aset.freq.meet.criteria.OR <- rowAlls( ((q.otu.freq.OR[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.freq.OR) | (!Aset.freq.OR) )
                        n.stable.freq.OR <- ifelse(Aset.freq.meet.criteria.OR, n.stable.freq.OR + Aset.freq.meet.criteria.OR, 0)
                        Aset.freq.rm.row.OR <- (n.stable.freq.OR >= n.stable.max)                    
                        Aset.freq.OR[Aset.freq.rm.row.OR,] = FALSE
                    }
                    
                    # combination test
                    if (any(Aset.freq.com)) {
                        
                        p.otu.freq.tmp <- 0.5*n.otu.freq 
                        p.otu.freq.OR.tmp <- 0.5*n.otu.freq.OR
                        pmin.otu.freq.com <- pmin(p.otu.freq.tmp, p.otu.freq.OR.tmp)
                        
                        mat = matrix(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                        p.otu.freq.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        
                        mat = matrix(otu.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                        p.otu.freq.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        
                        pmin.otu.freq.null.com <- pmin(p.otu.freq.null, p.otu.freq.null.OR)
                        if (length(dim(pmin.otu.freq.null.com))==2) pmin.otu.freq.null.com <- array(pmin.otu.freq.null.com, c(dim(pmin.otu.freq.null.com), 1))
                        
                        n.otu.freq.com <- rowSums( (pmin.otu.freq.null.com < c(pmin.otu.freq.com) - tol.eq) + 0.5 * (abs(pmin.otu.freq.null.com - c(pmin.otu.freq.com)) < tol.eq), dims=2) 
                        
                        AtoB.freq.com <- Aset.freq.com & (n.otu.freq.com >= n.rej.stop)
                        Aset.freq.com <- Aset.freq.com & !AtoB.freq.com
                        
                        p.otu.freq.com[,otu.smallp][AtoB.freq.com] <- n.otu.freq.com[AtoB.freq.com]*inv.n.perm.completed
                        p.otu.freq.com[,otu.smallp][Aset.freq.com] <- (n.otu.freq.com[Aset.freq.com]+1)*inv.n.perm.completed.1
                        
                        q.otu.freq.com <- t(apply(p.otu.freq.com, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.freq.com = matrix(q.otu.freq.com, ncol=1)
                        
                        Aset.freq.meet.criteria.com <- rowAlls( ((q.otu.freq.com[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.freq.com) | (!Aset.freq.com) ) 
                        n.stable.freq.com <- ifelse(Aset.freq.meet.criteria.com, n.stable.freq.com + Aset.freq.meet.criteria.com, 0)
                        Aset.freq.rm.row.com <- (n.stable.freq.com >= n.stable.max)                    
                        Aset.freq.com[Aset.freq.rm.row.com,] = FALSE
                    }
                }
                
                if (!freq.scale.only) {
                    if (any(Aset.tran)) {
                        AtoB.tran <- Aset.tran & (n.otu.tran >= n.rej.stop*2)
                        Aset.tran <- Aset.tran & !AtoB.tran
                        p.otu.tran[,otu.smallp][AtoB.tran] <- 0.5*n.otu.tran[AtoB.tran]*inv.n.perm.completed
                        p.otu.tran[,otu.smallp][Aset.tran] <- (0.5*n.otu.tran[Aset.tran]+1)*inv.n.perm.completed.1
                        
                        q.otu.tran <- t(apply(p.otu.tran, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.tran = matrix(q.otu.tran, ncol=1)
                        
                        Aset.tran.meet.criteria <- rowAlls( ((q.otu.tran[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.tran) | (!Aset.tran) ) 
                        n.stable.tran <- ifelse(Aset.tran.meet.criteria, n.stable.tran + Aset.tran.meet.criteria, 0)
                        Aset.tran.rm.row <- (n.stable.tran >= n.stable.max)                    
                        Aset.tran[Aset.tran.rm.row,] = FALSE
                    }
                    if (!is.null(OR)) {
                        if (any(Aset.tran.OR)) {
                            AtoB.tran.OR <- Aset.tran.OR & (n.otu.tran.OR >= n.rej.stop*2)
                            Aset.tran.OR <- Aset.tran.OR & !AtoB.tran.OR
                            p.otu.tran.OR[,otu.smallp][AtoB.tran.OR] <- 0.5*n.otu.tran.OR[AtoB.tran.OR]*inv.n.perm.completed
                            p.otu.tran.OR[,otu.smallp][Aset.tran.OR] <- (0.5*n.otu.tran.OR[Aset.tran.OR]+1)*inv.n.perm.completed.1
                            
                            q.otu.tran.OR <- t(apply(p.otu.tran.OR, 1, fdr.Sandev))
                            if (n.otu == 1) q.otu.tran.OR = matrix(q.otu.tran.OR, ncol=1)
                            
                            Aset.tran.meet.criteria.OR <- rowAlls( ((q.otu.tran.OR[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.tran.OR) | (!Aset.tran.OR) ) 
                            n.stable.tran.OR <- ifelse(Aset.tran.meet.criteria.OR, n.stable.tran.OR + Aset.tran.meet.criteria.OR, 0)
                            Aset.tran.rm.row.OR <- (n.stable.tran.OR >= n.stable.max)                    
                            Aset.tran.OR[Aset.tran.rm.row.OR,] = FALSE
                        }
                        
                        # combination test
                        if (any(Aset.tran.com)) {
                            
                            p.otu.tran.tmp <- 0.5*n.otu.tran 
                            p.otu.tran.OR.tmp <- 0.5*n.otu.tran.OR
                            pmin.otu.tran.com <- pmin(p.otu.tran.tmp, p.otu.tran.OR.tmp)
                            
                            mat = matrix(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                            p.otu.tran.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                            
                            mat = matrix(otu.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                            p.otu.tran.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                            
                            pmin.otu.tran.null.com <- pmin(p.otu.tran.null, p.otu.tran.null.OR)
                            if (length(dim(pmin.otu.tran.null.com))==2) pmin.otu.tran.null.com <- array(pmin.otu.tran.null.com, c(dim(pmin.otu.tran.null.com), 1))
                            
                            n.otu.tran.com <- rowSums( (pmin.otu.tran.null.com < c(pmin.otu.tran.com) - tol.eq) + 0.5 * (abs(pmin.otu.tran.null.com - c(pmin.otu.tran.com)) < tol.eq), dims=2) 
                            
                            AtoB.tran.com <- Aset.tran.com & (n.otu.tran.com >= n.rej.stop)
                            Aset.tran.com <- Aset.tran.com & !AtoB.tran.com
                            
                            p.otu.tran.com[,otu.smallp][AtoB.tran.com] <- n.otu.tran.com[AtoB.tran.com]*inv.n.perm.completed
                            p.otu.tran.com[,otu.smallp][Aset.tran.com] <- (n.otu.tran.com[Aset.tran.com]+1)*inv.n.perm.completed.1
                            
                            q.otu.tran.com <- t(apply(p.otu.tran.com, 1, fdr.Sandev))
                            if (n.otu == 1) q.otu.tran.com = matrix(q.otu.tran.com, ncol=1)
                            
                            Aset.tran.meet.criteria.com <- rowAlls( ((q.otu.tran.com[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.tran.com) | (!Aset.tran.com) ) 
                            n.stable.tran.com <- ifelse(Aset.tran.meet.criteria.com, n.stable.tran.com + Aset.tran.meet.criteria.com, 0)
                            Aset.tran.rm.row.com <- (n.stable.tran.com >= n.stable.max)                    
                            Aset.tran.com[Aset.tran.rm.row.com,] = FALSE
                        }
                    }
                }
            } # if (!all.rarefy | test.omni3)
            
            if (all.rarefy | test.omni3) {
                if (any(Aset.pa)) {
                    AtoB.pa <- Aset.pa & (n.otu.pa >= n.rej.stop*2)
                    Aset.pa <- Aset.pa & !AtoB.pa
                    p.otu.pa[,otu.smallp][AtoB.pa] <- 0.5*n.otu.pa[AtoB.pa]*inv.n.perm.completed
                    p.otu.pa[,otu.smallp][Aset.pa] <- (0.5*n.otu.pa[Aset.pa]+1)*inv.n.perm.completed.1
                    
                    q.otu.pa <- t(apply(p.otu.pa, 1, fdr.Sandev))
                    if (n.otu == 1) q.otu.pa = matrix(q.otu.pa, ncol=1)
                    
                    Aset.pa.meet.criteria <- rowAlls( ((q.otu.pa[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.pa) | (!Aset.pa) ) 
                    n.stable.pa <- ifelse(Aset.pa.meet.criteria, n.stable.pa + Aset.pa.meet.criteria, 0)
                    Aset.pa.rm.row <- (n.stable.pa >= n.stable.max)    
                    Aset.pa[Aset.pa.rm.row,] = FALSE
                }
                if (!is.null(OR)) {
                    if (any(Aset.pa.OR)) {
                        AtoB.pa.OR <- Aset.pa.OR & (n.otu.pa.OR >= n.rej.stop*2)
                        Aset.pa.OR <- Aset.pa.OR & !AtoB.pa.OR
                        p.otu.pa.OR[,otu.smallp][AtoB.pa.OR] <- 0.5*n.otu.pa.OR[AtoB.pa.OR]*inv.n.perm.completed
                        p.otu.pa.OR[,otu.smallp][Aset.pa.OR] <- (0.5*n.otu.pa.OR[Aset.pa.OR]+1)*inv.n.perm.completed.1
                        
                        q.otu.pa.OR <- t(apply(p.otu.pa.OR, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.pa.OR = matrix(q.otu.pa.OR, ncol=1)
                        
                        Aset.pa.meet.criteria.OR <- rowAlls( ((q.otu.pa.OR[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.pa.OR) | (!Aset.pa.OR) ) 
                        n.stable.pa.OR <- ifelse(Aset.pa.meet.criteria.OR, n.stable.pa.OR + Aset.pa.meet.criteria.OR, 0)
                        Aset.pa.rm.row.OR <- (n.stable.pa.OR >= n.stable.max)    
                        Aset.pa.OR[Aset.pa.rm.row.OR,] = FALSE
                    }
                    
                    # combination test
                    if (any(Aset.pa.com)) {
                        
                        p.otu.pa.tmp <- 0.5*n.otu.pa 
                        p.otu.pa.OR.tmp <- 0.5*n.otu.pa.OR
                        pmin.otu.pa.com <- pmin(p.otu.pa.tmp, p.otu.pa.OR.tmp)
                        
                        mat = matrix(otu.pa.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                        p.otu.pa.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        
                        mat = matrix(otu.pa.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                        p.otu.pa.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        
                        pmin.otu.pa.null.com <- pmin(p.otu.pa.null, p.otu.pa.null.OR)
                        if (length(dim(pmin.otu.pa.null.com))==2) pmin.otu.pa.null.com <- array(pmin.otu.pa.null.com, c(dim(pmin.otu.pa.null.com), 1))
                        
                        n.otu.pa.com <- rowSums( (pmin.otu.pa.null.com < c(pmin.otu.pa.com) - tol.eq) + 0.5 * (abs(pmin.otu.pa.null.com - c(pmin.otu.pa.com)) < tol.eq), dims=2) 
                        
                        AtoB.pa.com <- Aset.pa.com & (n.otu.pa.com >= n.rej.stop)
                        Aset.pa.com <- Aset.pa.com & !AtoB.pa.com
                        
                        p.otu.pa.com[,otu.smallp][AtoB.pa.com] <- n.otu.pa.com[AtoB.pa.com]*inv.n.perm.completed
                        p.otu.pa.com[,otu.smallp][Aset.pa.com] <- (n.otu.pa.com[Aset.pa.com]+1)*inv.n.perm.completed.1
                        
                        q.otu.pa.com <- t(apply(p.otu.pa.com, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.pa.com = matrix(q.otu.pa.com, ncol=1)
                        
                        Aset.pa.meet.criteria.com <- rowAlls( ((q.otu.pa.com[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.pa.com) | (!Aset.pa.com) ) 
                        n.stable.pa.com <- ifelse(Aset.pa.meet.criteria.com, n.stable.pa.com + Aset.pa.meet.criteria.com, 0)
                        Aset.pa.rm.row.com <- (n.stable.pa.com >= n.stable.max)                    
                        Aset.pa.com[Aset.pa.rm.row.com,] = FALSE
                    }
                }
            } # if (all.rarefy | test.omni3)
            
            if (test.omni3 | (!all.rarefy & !freq.scale.only)) {
                
                #################################
                # Dependent sequential stop rule
                #################################
                
                any.Aset.omni.com = ifelse(!is.null(OR), any(Aset.omni, Aset.omni.OR, Aset.omni.com), any(Aset.omni))
                any.Aset.omni3 = ifelse(test.omni3, any(Aset.omni, Aset.omni3), any(Aset.omni))
                any.Aset.omni3.com = ifelse(test.omni3 & !is.null(OR), any(Aset.omni, Aset.omni3, Aset.omni.OR, Aset.omni3.OR, Aset.omni3.com), any(Aset.omni))
                
                if (any.Aset.omni.com | any.Aset.omni3 | any.Aset.omni3.com) {
                    
                    p.otu.freq.tmp <- 0.5*n.otu.freq 
                    p.otu.tran.tmp <- 0.5*n.otu.tran
                    pmin.otu.omni <- pmin(p.otu.freq.tmp, p.otu.tran.tmp)
                    
                    if (is.null(p.otu.freq.null)) {
                        mat = matrix(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                        p.otu.freq.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                    }
                    if (is.null(p.otu.tran.null)) {
                        mat = matrix(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                        p.otu.tran.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                    }
                    pmin.otu.omni.null <- pmin(p.otu.freq.null, p.otu.tran.null)
                    if (length(dim(pmin.otu.omni.null))==2) pmin.otu.omni.null <- array(pmin.otu.omni.null, c(dim(pmin.otu.omni.null), 1))
                    
                    n.otu.omni <- rowSums( (pmin.otu.omni.null < c(pmin.otu.omni) - tol.eq) + 0.5 * (abs(pmin.otu.omni.null - c(pmin.otu.omni)) < tol.eq), dims=2) 
                    
                    AtoB.omni <- Aset.omni & (n.otu.omni >= n.rej.stop)
                    Aset.omni <- Aset.omni & !AtoB.omni
                    
                    p.otu.omni[,otu.smallp][AtoB.omni] <- n.otu.omni[AtoB.omni]*inv.n.perm.completed
                    p.otu.omni[,otu.smallp][Aset.omni] <- (n.otu.omni[Aset.omni]+1)*inv.n.perm.completed.1
                    
                    q.otu.omni <- t(apply(p.otu.omni, 1, fdr.Sandev))
                    if (n.otu == 1) q.otu.omni = matrix(q.otu.omni, ncol=1)
                    
                    Aset.omni.meet.criteria <- rowAlls( ((q.otu.omni[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni) | (!Aset.omni) ) 
                    n.stable.omni <- ifelse(Aset.omni.meet.criteria, n.stable.omni + Aset.omni.meet.criteria, 0)
                    Aset.omni.rm.row <- (n.stable.omni >= n.stable.max)    
                    Aset.omni[Aset.omni.rm.row,] = FALSE
                }
                
                if (!is.null(OR)) {
                    if (any.Aset.omni.com | any.Aset.omni3.com) {
                        
                        p.otu.freq.tmp.OR <- 0.5*n.otu.freq.OR 
                        p.otu.tran.tmp.OR <- 0.5*n.otu.tran.OR
                        pmin.otu.omni.OR <- pmin(p.otu.freq.tmp.OR, p.otu.tran.tmp.OR)
                        
                        if (is.null(p.otu.freq.null.OR)) {
                            mat = matrix(otu.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                            p.otu.freq.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        }
                        if (is.null(p.otu.tran.null.OR)) {
                            mat = matrix(otu.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                            p.otu.tran.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        }
                        pmin.otu.omni.null.OR <- pmin(p.otu.freq.null.OR, p.otu.tran.null.OR)
                        if (length(dim(pmin.otu.omni.null.OR))==2) pmin.otu.omni.null.OR <- array(pmin.otu.omni.null.OR, c(dim(pmin.otu.omni.null.OR), 1))
                        
                        n.otu.omni.OR <- rowSums( (pmin.otu.omni.null.OR < c(pmin.otu.omni.OR) - tol.eq) + 0.5 * (abs(pmin.otu.omni.null.OR - c(pmin.otu.omni.OR)) < tol.eq), dims=2) 
                        AtoB.omni.OR <- Aset.omni.OR & (n.otu.omni.OR >= n.rej.stop)
                        Aset.omni.OR <- Aset.omni.OR & !AtoB.omni.OR
                        
                        p.otu.omni.OR[,otu.smallp][AtoB.omni.OR] <- n.otu.omni.OR[AtoB.omni.OR]*inv.n.perm.completed
                        p.otu.omni.OR[,otu.smallp][Aset.omni.OR] <- (n.otu.omni.OR[Aset.omni.OR]+1)*inv.n.perm.completed.1
                        
                        q.otu.omni.OR <- t(apply(p.otu.omni.OR, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.omni.OR = matrix(q.otu.omni.OR, ncol=1)
                        
                        Aset.omni.meet.criteria.OR <- rowAlls( ((q.otu.omni.OR[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni.OR) | (!Aset.omni.OR) ) 
                        n.stable.omni.OR <- ifelse(Aset.omni.meet.criteria.OR, n.stable.omni.OR + Aset.omni.meet.criteria.OR, 0)
                        Aset.omni.rm.row.OR <- (n.stable.omni.OR >= n.stable.max)    
                        Aset.omni.OR[Aset.omni.rm.row.OR,] = FALSE
                    }
                    
                    # combination test
                    if (any(Aset.omni.com)) {
                        
                        pmin.otu.omni.com <- pmin(pmin.otu.omni, pmin.otu.omni.OR)
                        pmin.otu.omni.null.com <- pmin(pmin.otu.omni.null, pmin.otu.omni.null.OR)
                        if (length(dim(pmin.otu.omni.null.com))==2) pmin.otu.omni.null.com <- array(pmin.otu.omni.null.com, c(dim(pmin.otu.omni.null.com), 1))
                        
                        n.otu.omni.com <- rowSums( (pmin.otu.omni.null.com < c(pmin.otu.omni.com) - tol.eq) + 0.5 * (abs(pmin.otu.omni.null.com - c(pmin.otu.omni.com)) < tol.eq), dims=2) 
                        
                        AtoB.omni.com <- Aset.omni.com & (n.otu.omni.com >= n.rej.stop)
                        Aset.omni.com <- Aset.omni.com & !AtoB.omni.com
                        
                        p.otu.omni.com[,otu.smallp][AtoB.omni.com] <- n.otu.omni.com[AtoB.omni.com]*inv.n.perm.completed
                        p.otu.omni.com[,otu.smallp][Aset.omni.com] <- (n.otu.omni.com[Aset.omni.com]+1)*inv.n.perm.completed.1
                        
                        q.otu.omni.com <- t(apply(p.otu.omni.com, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.omni.com = matrix(q.otu.omni.com, ncol=1)
                        
                        Aset.omni.meet.criteria.com <- rowAlls( ((q.otu.omni.com[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni.com) | (!Aset.omni.com) ) 
                        n.stable.omni.com <- ifelse(Aset.omni.meet.criteria.com, n.stable.omni.com + Aset.omni.meet.criteria.com, 0)
                        Aset.omni.rm.row.com <- (n.stable.omni.com >= n.stable.max)                    
                        Aset.omni.com[Aset.omni.rm.row.com,] = FALSE
                    }
                }
                
                
                if (test.omni3) {
                    
                    if (any.Aset.omni3 | any.Aset.omni3.com) {
                        p.otu.pa.tmp <- 0.5*n.otu.pa
                        pmin.otu.omni3 <- pmin(p.otu.freq.tmp, p.otu.tran.tmp, p.otu.pa.tmp)
                        if (is.null(p.otu.pa.null)) {
                            mat = matrix(otu.pa.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                            p.otu.pa.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                        }
                        pmin.otu.omni3.null <- pmin(p.otu.freq.null, p.otu.tran.null, p.otu.pa.null) # n.var1 x n.otu x n.perm.completed
                        
                        if (length(dim(pmin.otu.omni3.null))==2) pmin.otu.omni3.null <- array(pmin.otu.omni3.null, c(dim(pmin.otu.omni3.null), 1))
                        
                        n.otu.omni3 <- rowSums( (pmin.otu.omni3.null < c(pmin.otu.omni3) - tol.eq) + 0.5 * (abs(pmin.otu.omni3.null - c(pmin.otu.omni3)) < tol.eq), dims=2) 
                        
                        AtoB.omni3 <- Aset.omni3 & (n.otu.omni3 >= n.rej.stop)
                        Aset.omni3 <- Aset.omni3 & !AtoB.omni3
                        
                        p.otu.omni3[,otu.smallp][AtoB.omni3] <- n.otu.omni3[AtoB.omni3]*inv.n.perm.completed
                        p.otu.omni3[,otu.smallp][Aset.omni3] <- (n.otu.omni3[Aset.omni3]+1)*inv.n.perm.completed.1
                        
                        q.otu.omni3 <- t(apply(p.otu.omni3, 1, fdr.Sandev))
                        if (n.otu == 1) q.otu.omni3 = matrix(q.otu.omni3, ncol=1)
                        
                        Aset.omni3.meet.criteria <- rowAlls( ((q.otu.omni3[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni3) | (!Aset.omni3) ) 
                        n.stable.omni3 <- ifelse(Aset.omni3.meet.criteria, n.stable.omni3 + Aset.omni3.meet.criteria, 0)
                        Aset.omni3.rm.row <- (n.stable.omni3 >= n.stable.max)    
                        Aset.omni3[Aset.omni3.rm.row,] = FALSE
                    }
                    
                    if (!is.null(OR)) {
                        if (any.Aset.omni3.com) {
                            p.otu.pa.tmp.OR <- 0.5*n.otu.pa.OR
                            pmin.otu.omni3.OR <- pmin(p.otu.freq.tmp.OR, p.otu.tran.tmp.OR, p.otu.pa.tmp.OR)
                            if (is.null(p.otu.pa.null.OR)) {
                                mat = matrix(otu.pa.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                p.otu.pa.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) 
                            }
                            pmin.otu.omni3.null.OR <- pmin(p.otu.freq.null.OR, p.otu.tran.null.OR, p.otu.pa.null.OR) # n.var1 x n.otu x n.perm.completed
                            if (length(dim(pmin.otu.omni3.null.OR))==2) pmin.otu.omni3.null.OR <- array(pmin.otu.omni3.null.OR, c(dim(pmin.otu.omni3.null.OR), 1))
                            
                            n.otu.omni3.OR <- rowSums( (pmin.otu.omni3.null.OR < c(pmin.otu.omni3.OR) - tol.eq) + 0.5 * (abs(pmin.otu.omni3.null.OR - c(pmin.otu.omni3.OR)) < tol.eq), dims=2) 
                            AtoB.omni3.OR <- Aset.omni3.OR & (n.otu.omni3.OR >= n.rej.stop)
                            Aset.omni3.OR <- Aset.omni3.OR & !AtoB.omni3.OR
                            
                            p.otu.omni3.OR[,otu.smallp][AtoB.omni3.OR] <- n.otu.omni3.OR[AtoB.omni3.OR]*inv.n.perm.completed
                            p.otu.omni3.OR[,otu.smallp][Aset.omni3.OR] <- (n.otu.omni3.OR[Aset.omni3.OR]+1)*inv.n.perm.completed.1
                            
                            q.otu.omni3.OR <- t(apply(p.otu.omni3.OR, 1, fdr.Sandev))
                            if (n.otu == 1) q.otu.omni3.OR = matrix(q.otu.omni3.OR, ncol=1)
                            
                            Aset.omni3.meet.criteria.OR <- rowAlls( ((q.otu.omni3.OR[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni3.OR) | (!Aset.omni3.OR) ) 
                            n.stable.omni3.OR <- ifelse(Aset.omni3.meet.criteria.OR, n.stable.omni3.OR + Aset.omni3.meet.criteria.OR, 0)
                            Aset.omni3.rm.row.OR <- (n.stable.omni3.OR >= n.stable.max)    
                            Aset.omni3.OR[Aset.omni3.rm.row.OR,] = FALSE
                        }
                        
                        # combination test
                        if (any(Aset.omni3.com)) {
                            
                            pmin.otu.omni3.com <- pmin(pmin.otu.omni3, pmin.otu.omni3.OR)
                            pmin.otu.omni3.null.com <- pmin(pmin.otu.omni3.null, pmin.otu.omni3.null.OR)
                            if (length(dim(pmin.otu.omni3.null.com))==2) pmin.otu.omni3.null.com <- array(pmin.otu.omni3.null.com, c(dim(pmin.otu.omni3.null.com), 1))
                            
                            n.otu.omni3.com <- rowSums( (pmin.otu.omni3.null.com < c(pmin.otu.omni3.com) - tol.eq) + 0.5 * (abs(pmin.otu.omni3.null.com - c(pmin.otu.omni3.com)) < tol.eq), dims=2) 
                            
                            AtoB.omni3.com <- Aset.omni3.com & (n.otu.omni3.com >= n.rej.stop)
                            Aset.omni3.com <- Aset.omni3.com & !AtoB.omni3.com
                            
                            p.otu.omni3.com[,otu.smallp][AtoB.omni3.com] <- n.otu.omni3.com[AtoB.omni3.com]*inv.n.perm.completed
                            p.otu.omni3.com[,otu.smallp][Aset.omni3.com] <- (n.otu.omni3.com[Aset.omni3.com]+1)*inv.n.perm.completed.1
                            
                            q.otu.omni3.com <- t(apply(p.otu.omni3.com, 1, fdr.Sandev))
                            if (n.otu == 1) q.otu.omni3.com = matrix(q.otu.omni3.com, ncol=1)
                            
                            Aset.omni3.meet.criteria.com <- rowAlls( ((q.otu.omni3.com[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni3.com) | (!Aset.omni3.com) ) 
                            n.stable.omni3.com <- ifelse(Aset.omni3.meet.criteria.com, n.stable.omni3.com + Aset.omni3.meet.criteria.com, 0)
                            Aset.omni3.rm.row.com <- (n.stable.omni3.com >= n.stable.max)                    
                            Aset.omni3.com[Aset.omni3.rm.row.com,] = FALSE
                        }
                    }
                } # if (test.omni3)
            } # if (test.omni3 | (!all.rarefy & !freq.scale.only)) 
            
            if (test.omni3) {
                meet.sandve.stop = !any(Aset.freq) & !any(Aset.tran) & !any(Aset.pa) & !any(Aset.omni3)
                if (!is.null(OR)) meet.sandve.stop = meet.sandve.stop & !any(Aset.freq.OR) & !any(Aset.tran.OR) & !any(Aset.pa.OR) & !any(Aset.omni3.OR)
            } else if (all.rarefy) {
                meet.sandve.stop = !any(Aset.pa)
                if (!is.null(OR)) meet.sandve.stop = meet.sandve.stop & !any(Aset.pa.OR)
            } else if (freq.scale.only) {
                meet.sandve.stop = !any(Aset.freq)
                if (!is.null(OR)) meet.sandve.stop = meet.sandve.stop & !any(Aset.freq.OR)
            } else {
                meet.sandve.stop = !any(Aset.freq) & !any(Aset.tran) & !any(Aset.omni)
                if (!is.null(OR)) meet.sandve.stop = meet.sandve.stop & !any(Aset.freq.OR) & !any(Aset.tran.OR) & !any(Aset.omni.OR)
            }
            
            if (meet.sandve.stop) {
                otu.tests.stopped = TRUE 
                cat("otu test stopped at permutation", n.perm.completed, "\n")
                n.otu.perm.completed = n.perm.completed
            }
            
            ################
            # test global
            ################
            
            if (!global.tests.done) {
                
                if (test.omni3) {
                    meet.rej.stop = all(n.global.freq >= n.rej.stop*2) & all(n.global.tran >= n.rej.stop*2) & all(n.global.pa >= n.rej.stop*2)
                    if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(n.global.freq.OR >= n.rej.stop*2) & all(n.global.tran.OR >= n.rej.stop*2) & all(n.global.pa.OR >= n.rej.stop*2)
                } else if (all.rarefy) {
                    meet.rej.stop = all(n.global.pa >= n.rej.stop*2)
                    if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(n.global.pa.OR >= n.rej.stop*2)
                } else if (freq.scale.only) {
                    meet.rej.stop = all(n.global.freq >= n.rej.stop*2)
                    if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(n.global.freq.OR >= n.rej.stop*2)
                } else {
                    meet.rej.stop = all(n.global.freq >= n.rej.stop*2) & all(n.global.tran >= n.rej.stop*2)
                    if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(n.global.freq.OR >= n.rej.stop*2) & all(n.global.tran.OR >= n.rej.stop*2)
                }
                
                if (test.mediation) {
                    if (test.omni3) {
                        meet.rej.stop = meet.rej.stop & all(med.n.global.freq >= n.rej.stop*2) & all(med.n.global.tran >= n.rej.stop*2) & all(med.n.global.pa >= n.rej.stop*2)
                        if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(med.n.global.freq.OR >= n.rej.stop*2) & all(med.n.global.tran.OR >= n.rej.stop*2) & all(med.n.global.pa.OR >= n.rej.stop*2)
                    } else if (all.rarefy) {
                        meet.rej.stop = meet.rej.stop & all(med.n.global.pa >= n.rej.stop*2)
                        if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(med.n.global.pa.OR >= n.rej.stop*2)
                    } else if (freq.scale.only) {
                        meet.rej.stop = meet.rej.stop & all(med.n.global.freq >= n.rej.stop*2)
                        if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(med.n.global.freq.OR >= n.rej.stop*2)
                    } else {
                        meet.rej.stop = meet.rej.stop & all(med.n.global.freq >= n.rej.stop*2) & all(med.n.global.tran >= n.rej.stop*2)
                        if (!is.null(OR)) meet.rej.stop = meet.rej.stop & all(med.n.global.freq.OR >= n.rej.stop*2) & all(med.n.global.tran.OR >= n.rej.stop*2)
                    }
                }
                
                if (meet.rej.stop | (n.perm.completed >= n.global.perm.max)) {
                    
                    if (!all.rarefy | test.omni3) {
                        
                        p.global.freq = ifelse((n.global.freq >= n.rej.stop*2), 0.5*n.global.freq*inv.n.perm.completed, (0.5*n.global.freq+1)*inv.n.perm.completed.1)
                        if (!is.null(OR)) {
                            p.global.freq.OR = ifelse((n.global.freq.OR >= n.rej.stop*2), 0.5*n.global.freq.OR*inv.n.perm.completed, (0.5*n.global.freq.OR+1)*inv.n.perm.completed.1)
                            
                            # combination test
                            p.global.freq.tmp <- 0.5*n.global.freq
                            p.global.freq.tmp.OR <- 0.5*n.global.freq.OR
                            p.global.freq.null <- n.perm.completed + 0.5 - array(rowRanks(global.freq.perm[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                            p.global.freq.null.OR <- n.perm.completed + 0.5 - array(rowRanks(global.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                            
                            pmin.global.freq.com <- pmin(p.global.freq.tmp, p.global.freq.tmp.OR)
                            pmin.global.freq.null.com <- pmin(p.global.freq.null, p.global.freq.null.OR)
                            n.global.freq.com <- rowSums( (pmin.global.freq.null.com < c(pmin.global.freq.com) - tol.eq) + 0.5 * (abs(pmin.global.freq.null.com - c(pmin.global.freq.com)) < tol.eq))
                            p.global.freq.com = ifelse((n.global.freq.com >= n.rej.stop), n.global.freq.com*inv.n.perm.completed, (n.global.freq.com+1)*inv.n.perm.completed.1)
                        }
                        
                        
                        if (!freq.scale.only) {
                            
                            p.global.tran = ifelse((n.global.tran >= n.rej.stop*2), 0.5*n.global.tran*inv.n.perm.completed, (0.5*n.global.tran+1)*inv.n.perm.completed.1)
                            if (!is.null(OR)) {
                                p.global.tran.OR = ifelse((n.global.tran.OR >= n.rej.stop*2), 0.5*n.global.tran.OR*inv.n.perm.completed, (0.5*n.global.tran.OR+1)*inv.n.perm.completed.1)
                                
                                # combination test
                                p.global.tran.tmp <- 0.5*n.global.tran
                                p.global.tran.tmp.OR <- 0.5*n.global.tran.OR
                                p.global.tran.null <- n.perm.completed + 0.5 - array(rowRanks(global.tran.perm[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                                p.global.tran.null.OR <- n.perm.completed + 0.5 - array(rowRanks(global.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                                
                                pmin.global.tran.com <- pmin(p.global.tran.tmp, p.global.tran.tmp.OR)
                                pmin.global.tran.null.com <- pmin(p.global.tran.null, p.global.tran.null.OR)
                                n.global.tran.com <- rowSums( (pmin.global.tran.null.com < c(pmin.global.tran.com) - tol.eq) + 0.5 * (abs(pmin.global.tran.null.com - c(pmin.global.tran.com)) < tol.eq))
                                p.global.tran.com = ifelse((n.global.tran.com >= n.rej.stop), n.global.tran.com*inv.n.perm.completed, (n.global.tran.com+1)*inv.n.perm.completed.1)
                            }
                        }
                    } # if (!all.rarefy | test.omni3)
                    
                    if (all.rarefy | test.omni3) {
                        
                        p.global.pa = ifelse((n.global.pa >= n.rej.stop*2), 0.5*n.global.pa*inv.n.perm.completed, (0.5*n.global.pa+1)*inv.n.perm.completed.1)
                        if (!is.null(OR)) {
                            p.global.pa.OR = ifelse((n.global.pa.OR >= n.rej.stop*2), 0.5*n.global.pa.OR*inv.n.perm.completed, (0.5*n.global.pa.OR+1)*inv.n.perm.completed.1)
                            
                            # combination test
                            p.global.pa.tmp <- 0.5*n.global.pa
                            p.global.pa.tmp.OR <- 0.5*n.global.pa.OR
                            p.global.pa.null <- n.perm.completed + 0.5 - rowRanks(matrix(global.pa.perm[,1,1:n.perm.completed],ncol=n.perm.completed))
                            p.global.pa.null.OR <- n.perm.completed + 0.5 - rowRanks(matrix(global.pa.perm.OR[,1,1:n.perm.completed],ncol=n.perm.completed))
                            
                            pmin.global.pa.com <- pmin(p.global.pa.tmp, p.global.pa.tmp.OR)
                            pmin.global.pa.null.com <- pmin(p.global.pa.null, p.global.pa.null.OR)
                            n.global.pa.com <- rowSums( (pmin.global.pa.null.com < c(pmin.global.pa.com) - tol.eq) + 0.5 * (abs(pmin.global.pa.null.com - c(pmin.global.pa.com)) < tol.eq))
                            p.global.pa.com = ifelse((n.global.pa.com >= n.rej.stop), n.global.pa.com*inv.n.perm.completed, (n.global.pa.com+1)*inv.n.perm.completed.1)
                        }
                    }
                    
                    if (test.omni3 | (!all.rarefy & !freq.scale.only)) {
                        
                        if (is.null(p.global.freq.tmp)) p.global.freq.tmp <- 0.5*n.global.freq
                        if (is.null(p.global.tran.tmp)) p.global.tran.tmp <- 0.5*n.global.tran
                        
                        pmin.global.omni <- pmin(p.global.freq.tmp, p.global.tran.tmp)
                        
                        if (is.null(p.global.freq.null)) {
                            p.global.freq.null <- n.perm.completed + 0.5 - array(rowRanks(global.freq.perm[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                        }
                        if (is.null(p.global.tran.null)) {
                            p.global.tran.null <- n.perm.completed + 0.5 - array(rowRanks(global.tran.perm[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                        }
                        
                        pmin.global.omni.null <- pmin(p.global.freq.null, p.global.tran.null)
                        
                        n.global.omni <- rowSums( (pmin.global.omni.null < c(pmin.global.omni) - tol.eq) + 0.5 * (abs(pmin.global.omni.null - c(pmin.global.omni)) < tol.eq), dims=2) 
                        p.global.omni = ifelse((n.global.omni >= n.rej.stop), n.global.omni*inv.n.perm.completed, (n.global.omni+1)*inv.n.perm.completed.1)
                        
                        if (!is.null(OR)) {
                            
                            if (is.null(p.global.freq.tmp.OR)) p.global.freq.tmp.OR <- 0.5*n.global.freq.OR 
                            if (is.null(p.global.tran.tmp.OR)) p.global.tran.tmp.OR <- 0.5*n.global.tran.OR
                            
                            pmin.global.omni.OR <- pmin(p.global.freq.tmp.OR, p.global.tran.tmp.OR)
                            
                            if (is.null(p.global.freq.null.OR)) {
                                p.global.freq.null.OR <- n.perm.completed + 0.5 - array(rowRanks(global.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                            }
                            if (is.null(p.global.tran.null.OR)) {
                                p.global.tran.null.OR <- n.perm.completed + 0.5 - array(rowRanks(global.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], ties.method = "average", dim.=c(n.var1*n.rarefy, n.global.perm.max)), dim = c(n.var1, n.rarefy, n.global.perm.max)) 
                            }
                            
                            pmin.global.omni.null.OR <- pmin(p.global.freq.null.OR, p.global.tran.null.OR)
                            
                            n.global.omni.OR <- rowSums( (pmin.global.omni.null.OR < c(pmin.global.omni.OR) - tol.eq) + 0.5 * (abs(pmin.global.omni.null.OR - c(pmin.global.omni.OR)) < tol.eq), dims=2) 
                            p.global.omni.OR = ifelse((n.global.omni.OR >= n.rej.stop), n.global.omni.OR*inv.n.perm.completed, (n.global.omni.OR+1)*inv.n.perm.completed.1)
                            
                            # combination test
                            pmin.global.omni.com <- pmin(pmin.global.omni, pmin.global.omni.OR)
                            pmin.global.omni.null.com <- pmin(pmin.global.omni.null, pmin.global.omni.null.OR)
                            n.global.omni.com <- rowSums( (pmin.global.omni.null.com < c(pmin.global.omni.com) - tol.eq) + 0.5 * (abs(pmin.global.omni.null.com - c(pmin.global.omni.com)) < tol.eq))
                            p.global.omni.com = ifelse((n.global.omni.com >= n.rej.stop), n.global.omni.com*inv.n.perm.completed, (n.global.omni.com+1)*inv.n.perm.completed.1)
                        }
                        
                        if (test.omni3) {
                            
                            # omni other methods
                            df1 <- rep(NA, n.var1)
                            for (k in 1:n.var1) {
                                k1 = k + as.numeric(adjust.for.confounders)
                                df1[k] <- fit.ldm$up[k1] - fit.ldm$low[k1] + 1
                            }
                            df2 <- n.obs - 1 - fit.ldm$up[n.var1]
                            F.sd <- sqrt(2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4)))
                            thresh <- F.sd/df2/10
                            
                            otu.freq.sd <- array(rowSds(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], dim.=c(n.var1*n.otu, n.perm.completed)), dim = c(n.var1, n.otu))
                            stats <- otu.freq.sd/F.sd
                            panal.otu.freq <- pf(sweep(otu.freq, MARGIN=c(1,2), FUN="/", STATS=stats), df1=df1, df2=df2, lower.tail=FALSE)
                            panal.otu.freq.perm <- pf(sweep(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], MARGIN=c(1,2), FUN="/", STATS=stats), df1=df1, df2=df2, lower.tail=FALSE)
                            panal.otu.freq[otu.freq.sd<thresh] <- 1
                            panal.otu.freq.perm[c(otu.freq.sd<thresh)] <- 1
                            
                            otu.tran.sd <- array(rowSds(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], dim.=c(n.var1*n.otu, n.perm.completed)), dim = c(n.var1, n.otu))
                            stats <- otu.tran.sd/F.sd
                            panal.otu.tran <- pf(sweep(otu.tran, MARGIN=c(1,2), FUN="/", STATS=stats), df1=df1, df2=df2, lower.tail=FALSE)
                            panal.otu.tran.perm <- pf(sweep(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], MARGIN=c(1,2), FUN="/", STATS=stats), df1=df1, df2=df2, lower.tail=FALSE)
                            panal.otu.tran[otu.tran.sd<thresh] <- 1
                            panal.otu.tran.perm[c(otu.tran.sd<thresh)] <- 1
                            
                            otu.pa.sd <- array(rowSds(otu.pa.perm[,,1:n.perm.completed,drop=FALSE], dim.=c(n.var1*n.otu, n.perm.completed)), dim = c(n.var1, n.otu))
                            stats <- otu.pa.sd/F.sd
                            panal.otu.pa <- pf(sweep(otu.pa, MARGIN=c(1,2), FUN="/", STATS=stats), df1=df1, df2=df2, lower.tail=FALSE)
                            panal.otu.pa.perm <- pf(sweep(otu.pa.perm[,,1:n.perm.completed,drop=FALSE], MARGIN=c(1,2), FUN="/", STATS=stats), df1=df1, df2=df2, lower.tail=FALSE)
                            panal.otu.pa[otu.pa.sd<thresh] <- 1
                            panal.otu.pa.perm[c(otu.pa.sd<thresh)] <- 1
                            
                            panalmin.otu.omni3 <- pmin(panal.otu.freq, panal.otu.tran, panal.otu.pa)
                            panalmin.otu.omni3.perm <- pmin(panal.otu.freq.perm, panal.otu.tran.perm, panal.otu.pa.perm)
                            
                            # harmonic mean
                            stat.global.harmonic.obs <- rowSums(1/panalmin.otu.omni3)
                            stat.global.harmonic.perm <- apply(1/panalmin.otu.omni3.perm, c(1,3), sum)
                            n.global.harmonic <- rowSums((stat.global.harmonic.perm > stat.global.harmonic.obs + tol.eq) + (stat.global.harmonic.perm > stat.global.harmonic.obs - tol.eq))
                            p.global.harmonic = ifelse((n.global.harmonic >= n.rej.stop*2), 0.5*n.global.harmonic*inv.n.perm.completed, (0.5*n.global.harmonic+1)*inv.n.perm.completed.1)
                            
                            # fisher
                            stat.global.fisher.obs <- rowSums((1/panalmin.otu.omni3))
                            stat.global.fisher.perm <- apply((1/panalmin.otu.omni3.perm), c(1,3), sum)
                            n.global.fisher <- rowSums((stat.global.fisher.perm > stat.global.fisher.obs + tol.eq) + (stat.global.fisher.perm > stat.global.fisher.obs - tol.eq))
                            p.global.fisher = ifelse((n.global.fisher >= n.rej.stop*2), 0.5*n.global.fisher*inv.n.perm.completed, (0.5*n.global.fisher+1)*inv.n.perm.completed.1)
                            
                            # omni3
                            if (is.null(p.global.pa.tmp)) p.global.pa.tmp <- 0.5*n.global.pa
                            p.global.harmonic.tmp <- 0.5*n.global.harmonic
                            p.global.fisher.tmp <- 0.5*n.global.fisher
                            
                            if (is.null(p.global.pa.null)) p.global.pa.null <- n.perm.completed + 0.5 - rowRanks(matrix(global.pa.perm[,1,1:n.perm.completed],ncol=n.perm.completed))
                            p.global.harmonic.null <- n.perm.completed + 0.5 - rowRanks(stat.global.harmonic.perm)
                            p.global.fisher.null <- n.perm.completed + 0.5 - rowRanks(stat.global.fisher.perm)
                            p.global.freq.null <- array(p.global.freq.null, dim=dim(p.global.fisher.null)) # only rarefy1 is used for omni3
                            p.global.tran.null <- array(p.global.tran.null, dim=dim(p.global.fisher.null))
                            
                            pmin.global.omni3 <- pmin(p.global.freq.tmp, p.global.tran.tmp, p.global.pa.tmp,
                                                      p.global.harmonic.tmp, p.global.fisher.tmp)
                            pmin.global.omni3.null <- pmin(p.global.freq.null, p.global.tran.null, p.global.pa.null,
                                                           p.global.harmonic.null, p.global.fisher.null)
                            
                            n.global.omni3 <- rowSums( (pmin.global.omni3.null < c(pmin.global.omni3) - tol.eq) + 0.5 * (abs(pmin.global.omni3.null - c(pmin.global.omni3)) < tol.eq))
                            p.global.omni3 = ifelse((n.global.omni3 >= n.rej.stop), n.global.omni3*inv.n.perm.completed, (n.global.omni3+1)*inv.n.perm.completed.1)
                            
                            if (!is.null(OR)) {
                                otu.freq.sd.OR <- array(rowSds(otu.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], dim.=c(n.var1*n.otu, n.perm.completed)), dim = c(n.var1, n.otu))
                                stats.OR <- otu.freq.sd.OR/F.sd
                                panal.otu.freq.OR <- pf(sweep(otu.freq.OR, MARGIN=c(1,2), FUN="/", STATS=stats.OR), df1=df1, df2=df2, lower.tail=FALSE)
                                panal.otu.freq.perm.OR <- pf(sweep(otu.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], MARGIN=c(1,2), FUN="/", STATS=stats.OR), df1=df1, df2=df2, lower.tail=FALSE)
                                panal.otu.freq.OR[otu.freq.sd.OR<thresh] <- 1
                                panal.otu.freq.perm.OR[c(otu.freq.sd.OR<thresh)] <- 1
                                
                                otu.tran.sd.OR <- array(rowSds(otu.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], dim.=c(n.var1*n.otu, n.perm.completed)), dim = c(n.var1, n.otu))
                                stats.OR <- otu.tran.sd.OR/F.sd
                                panal.otu.tran.OR <- pf(sweep(otu.tran.OR, MARGIN=c(1,2), FUN="/", STATS=stats.OR), df1=df1, df2=df2, lower.tail=FALSE)
                                panal.otu.tran.perm.OR <- pf(sweep(otu.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], MARGIN=c(1,2), FUN="/", STATS=stats.OR), df1=df1, df2=df2, lower.tail=FALSE)
                                panal.otu.tran.OR[otu.tran.sd.OR<thresh] <- 1
                                panal.otu.tran.perm.OR[c(otu.tran.sd.OR<thresh)] <- 1
                                
                                otu.pa.sd.OR <- array(rowSds(otu.pa.perm.OR[,,1:n.perm.completed,drop=FALSE], dim.=c(n.var1*n.otu, n.perm.completed)), dim = c(n.var1, n.otu))
                                stats.OR <- otu.pa.sd.OR/F.sd
                                panal.otu.pa.OR <- pf(sweep(otu.pa.OR, MARGIN=c(1,2), FUN="/", STATS=stats.OR), df1=df1, df2=df2, lower.tail=FALSE)
                                panal.otu.pa.perm.OR <- pf(sweep(otu.pa.perm.OR[,,1:n.perm.completed,drop=FALSE], MARGIN=c(1,2), FUN="/", STATS=stats.OR), df1=df1, df2=df2, lower.tail=FALSE)
                                panal.otu.pa.OR[otu.pa.sd.OR<thresh] <- 1
                                panal.otu.pa.perm.OR[c(otu.pa.sd.OR<thresh)] <- 1
                                
                                panalmin.otu.omni3.OR <- pmin(panal.otu.freq.OR, panal.otu.tran.OR, panal.otu.pa.OR)
                                panalmin.otu.omni3.perm.OR <- pmin(panal.otu.freq.perm.OR, panal.otu.tran.perm.OR, panal.otu.pa.perm.OR)
                                
                                # harmonic mean
                                stat.global.harmonic.obs.OR <- rowSums(1/panalmin.otu.omni3.OR)
                                stat.global.harmonic.perm.OR <- apply(1/panalmin.otu.omni3.perm.OR, c(1,3), sum)
                                n.global.harmonic.OR <- rowSums((stat.global.harmonic.perm.OR > stat.global.harmonic.obs.OR + tol.eq) + (stat.global.harmonic.perm.OR > stat.global.harmonic.obs.OR - tol.eq))
                                p.global.harmonic.OR = ifelse((n.global.harmonic.OR >= n.rej.stop*2), 0.5*n.global.harmonic.OR*inv.n.perm.completed, (0.5*n.global.harmonic.OR+1)*inv.n.perm.completed.1)
                                
                                # fisher
                                stat.global.fisher.obs.OR <- rowSums(-log(panalmin.otu.omni3.OR))
                                stat.global.fisher.perm.OR <- apply(-log(panalmin.otu.omni3.perm.OR), c(1,3), sum)
                                n.global.fisher.OR <- rowSums((stat.global.fisher.perm.OR > stat.global.fisher.obs.OR + tol.eq) + (stat.global.fisher.perm.OR > stat.global.fisher.obs.OR - tol.eq))
                                p.global.fisher.OR = ifelse((n.global.fisher.OR >= n.rej.stop*2), 0.5*n.global.fisher.OR*inv.n.perm.completed, (0.5*n.global.fisher.OR+1)*inv.n.perm.completed.1)
                                
                                # omni3
                                if (is.null(p.global.pa.tmp.OR)) p.global.pa.tmp.OR <- 0.5*n.global.pa.OR
                                p.global.harmonic.tmp.OR <- 0.5*n.global.harmonic.OR
                                p.global.fisher.tmp.OR <- 0.5*n.global.fisher.OR
                                
                                if (is.null(p.global.pa.null.OR)) p.global.pa.null.OR <- n.perm.completed + 0.5 - rowRanks(matrix(global.pa.perm.OR[,1,1:n.perm.completed],ncol=n.perm.completed))
                                p.global.harmonic.null.OR <- n.perm.completed + 0.5 - rowRanks(stat.global.harmonic.perm.OR)
                                p.global.fisher.null.OR <- n.perm.completed + 0.5 - rowRanks(stat.global.fisher.perm.OR)
                                p.global.freq.null.OR <- array(p.global.freq.null.OR, dim=dim(p.global.fisher.null.OR)) # only rarefy1 is used for omni3
                                p.global.tran.null.OR <- array(p.global.tran.null.OR, dim=dim(p.global.fisher.null.OR))
                                
                                pmin.global.omni3.OR <- pmin(p.global.freq.tmp.OR, p.global.tran.tmp.OR, p.global.pa.tmp.OR,
                                                             p.global.harmonic.tmp.OR, p.global.fisher.tmp.OR)
                                pmin.global.omni3.null.OR <- pmin(p.global.freq.null.OR, p.global.tran.null.OR, p.global.pa.null.OR,
                                                                  p.global.harmonic.null.OR, p.global.fisher.null.OR)
                                n.global.omni3.OR <- rowSums( (pmin.global.omni3.null.OR < c(pmin.global.omni3.OR) - tol.eq) + 0.5 * (abs(pmin.global.omni3.null.OR - c(pmin.global.omni3.OR)) < tol.eq))
                                p.global.omni3.OR = ifelse((n.global.omni3.OR >= n.rej.stop), n.global.omni3.OR*inv.n.perm.completed, (n.global.omni3.OR+1)*inv.n.perm.completed.1)
                                
                                # combination test
                                pmin.global.harmonic.com <- pmin(p.global.harmonic, p.global.harmonic.OR)
                                pmin.global.harmonic.null.com <- pmin(p.global.harmonic.null, p.global.harmonic.null.OR)
                                n.global.harmonic.com <- rowSums( (pmin.global.harmonic.null.com < c(pmin.global.harmonic.com) - tol.eq) + 0.5 * (abs(pmin.global.harmonic.null.com - c(pmin.global.harmonic.com)) < tol.eq))
                                p.global.harmonic.com = ifelse((n.global.harmonic.com >= n.rej.stop), n.global.harmonic.com*inv.n.perm.completed, (n.global.harmonic.com+1)*inv.n.perm.completed.1)
                                
                                pmin.global.fisher.com <- pmin(p.global.fisher, p.global.fisher.OR)
                                pmin.global.fisher.null.com <- pmin(p.global.fisher.null, p.global.fisher.null.OR)
                                n.global.fisher.com <- rowSums( (pmin.global.fisher.null.com < c(pmin.global.fisher.com) - tol.eq) + 0.5 * (abs(pmin.global.fisher.null.com - c(pmin.global.fisher.com)) < tol.eq))
                                p.global.fisher.com = ifelse((n.global.fisher.com >= n.rej.stop), n.global.fisher.com*inv.n.perm.completed, (n.global.fisher.com+1)*inv.n.perm.completed.1)
                                
                                pmin.global.omni3.com <- pmin(pmin.global.omni3, pmin.global.omni3.OR)
                                pmin.global.omni3.null.com <- pmin(pmin.global.omni3.null, pmin.global.omni3.null.OR)
                                n.global.omni3.com <- rowSums( (pmin.global.omni3.null.com < c(pmin.global.omni3.com) - tol.eq) + 0.5 * (abs(pmin.global.omni3.null.com - c(pmin.global.omni3.com)) < tol.eq))
                                p.global.omni3.com = ifelse((n.global.omni3.com >= n.rej.stop), n.global.omni3.com*inv.n.perm.completed, (n.global.omni3.com+1)*inv.n.perm.completed.1)
                                
                            }
                        } # if (test.omni3)
                    } # if (test.omni3 | (!all.rarefy & !freq.scale.only))
                    
                    
                    if (test.omni3) {
                        meet.all.rej.stop = all(n.global.freq >= n.rej.stop*2) & all(n.global.tran >= n.rej.stop*2) & all(n.global.pa >= n.rej.stop*2) & all(n.global.omni3 >= n.rej.stop)
                        if (!is.null(OR)) meet.all.rej.stop = meet.all.rej.stop & all(n.global.freq.OR >= n.rej.stop*2) & all(n.global.tran.OR >= n.rej.stop*2) & all(n.global.pa.OR >= n.rej.stop*2) & all(n.global.omni3.OR >= n.rej.stop)
                    } else if (all.rarefy) {
                        meet.all.rej.stop = all(n.global.pa >= n.rej.stop*2)
                        if (!is.null(OR)) meet.all.rej.stop = meet.all.rej.stop & all(n.global.pa.OR >= n.rej.stop*2)
                    } else if (freq.scale.only) {
                        meet.all.rej.stop = all(n.global.freq >= n.rej.stop*2)
                        if (!is.null(OR)) meet.all.rej.stop = meet.all.rej.stop & all(n.global.freq.OR >= n.rej.stop*2)
                    } else {
                        meet.all.rej.stop = all(n.global.freq >= n.rej.stop*2) & all(n.global.tran >= n.rej.stop*2) & all(n.global.omni >= n.rej.stop)
                        if (!is.null(OR)) meet.all.rej.stop = meet.all.rej.stop & all(n.global.freq.OR >= n.rej.stop*2) & all(n.global.tran.OR >= n.rej.stop*2) & all(n.global.omni.OR >= n.rej.stop)
                    }
                    
                    ######################
                    ##### mediation ######
                    ######################
                    
                    med.meet.all.rej.stop = TRUE
                    
                    if (test.mediation) {
                        
                        if (!all.rarefy | test.omni3) {
                            
                            # freq
                            
                            if (is.null(p.otu.freq.null)) {
                                mat = matrix(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                p.otu.freq.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) # pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                            }
                            p.otu.freq.null <- p.otu.freq.null*inv.n.perm.completed
                            
                            med.T.otu.freq = colMaxs(p.otu.freq)
                            
                            med.T.otu.freq.null1 = pmax(p.otu.freq.null[1,,], p.otu.freq[2,]) # dim: J*B
                            med.T.otu.freq.null2 = pmax(p.otu.freq.null[2,,], p.otu.freq[1,]) # dim: J*B
                            med.T.otu.freq.null3 = pmax(p.otu.freq.null[1,,], p.otu.freq.null[2,,]) # dim: J*B
                            med.T.otu.freq.null = pmin(med.T.otu.freq.null1, med.T.otu.freq.null2, med.T.otu.freq.null3) # dim: J*B
                            
                            med.global.freq = sum(1/med.T.otu.freq)
                            med.global.freq.perm = colSums(1/med.T.otu.freq.null) # length: B
                            
                            med.n.global.freq = sum(med.global.freq.perm > med.global.freq - tol.eq) + sum(med.global.freq.perm > med.global.freq + tol.eq)
                            med.p.global.freq = ifelse((med.n.global.freq >= n.rej.stop*2), 0.5*med.n.global.freq*inv.n.perm.completed, (0.5*med.n.global.freq+1)*inv.n.perm.completed.1)
                            
                            if (!is.null(OR)) {
                                if (is.null(p.otu.freq.null.OR)) {
                                    mat = matrix(otu.freq.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                    p.otu.freq.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) # pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                                }
                                p.otu.freq.null.OR <- p.otu.freq.null.OR*inv.n.perm.completed
                                
                                med.T.otu.freq.OR = colMaxs(p.otu.freq.OR)
                                
                                med.T.otu.freq.null1.OR = pmax(p.otu.freq.null.OR[1,,], p.otu.freq.OR[2,]) # dim: J*B
                                med.T.otu.freq.null2.OR = pmax(p.otu.freq.null.OR[2,,], p.otu.freq.OR[1,]) # dim: J*B
                                med.T.otu.freq.null3.OR = pmax(p.otu.freq.null.OR[1,,], p.otu.freq.null.OR[2,,]) # dim: J*B
                                med.T.otu.freq.null.OR = pmin(med.T.otu.freq.null1.OR, med.T.otu.freq.null2.OR, med.T.otu.freq.null3.OR) # dim: J*B
                                
                                med.global.freq.OR = sum(1/med.T.otu.freq.OR)
                                med.global.freq.perm.OR = colSums(1/med.T.otu.freq.null.OR) # length: B
                                
                                med.n.global.freq.OR = sum(med.global.freq.perm.OR > med.global.freq.OR - tol.eq) + sum(med.global.freq.perm.OR > med.global.freq.OR + tol.eq)
                                med.p.global.freq.OR = ifelse((med.n.global.freq.OR >= n.rej.stop*2), 0.5*med.n.global.freq.OR*inv.n.perm.completed, (0.5*med.n.global.freq.OR+1)*inv.n.perm.completed.1)
                                
                                # combination test
                                med.p.global.freq.tmp <- 0.5*med.n.global.freq
                                med.p.global.freq.tmp.OR <- 0.5*med.n.global.freq.OR
                                med.p.global.freq.null <- n.perm.completed + 0.5 - rank(med.global.freq.perm[1:n.perm.completed])
                                med.p.global.freq.null.OR <- n.perm.completed + 0.5 - rank(med.global.freq.perm.OR[1:n.perm.completed])
                                
                                med.pmin.global.freq.com <- pmin(med.p.global.freq.tmp, med.p.global.freq.tmp.OR)
                                med.pmin.global.freq.null.com <- pmin(med.p.global.freq.null, med.p.global.freq.null.OR)
                                med.n.global.freq.com <- sum( (med.pmin.global.freq.null.com < c(med.pmin.global.freq.com) - tol.eq) + 0.5 * (abs(med.pmin.global.freq.null.com - c(med.pmin.global.freq.com)) < tol.eq))
                                med.p.global.freq.com = ifelse((med.n.global.freq.com >= n.rej.stop), med.n.global.freq.com*inv.n.perm.completed, (med.n.global.freq.com+1)*inv.n.perm.completed.1)
                            }
                            
                            if (!freq.scale.only) {
                                
                                # tran
                                
                                if (is.null(p.otu.tran.null)) {
                                    mat = matrix(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                    p.otu.tran.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) # pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                                }
                                p.otu.tran.null <- p.otu.tran.null*inv.n.perm.completed
                                
                                med.T.otu.tran = colMaxs(p.otu.tran) # I_j
                                
                                med.T.otu.tran.null1 = pmax(p.otu.tran.null[1,,], p.otu.tran[2,]) # dim: J*B
                                med.T.otu.tran.null2 = pmax(p.otu.tran.null[2,,], p.otu.tran[1,]) # dim: J*B
                                med.T.otu.tran.null3 = pmax(p.otu.tran.null[1,,], p.otu.tran.null[2,,]) # dim: J*B
                                med.T.otu.tran.null = pmin(med.T.otu.tran.null1, med.T.otu.tran.null2, med.T.otu.tran.null3) # I_j^(b)
                                
                                med.global.tran = sum(1/med.T.otu.tran) # I
                                med.global.tran.perm = colSums(1/med.T.otu.tran.null) # I^(b)
                                
                                med.n.global.tran = sum(med.global.tran.perm > med.global.tran - tol.eq) + sum(med.global.tran.perm > med.global.tran + tol.eq)
                                med.p.global.tran = ifelse((med.n.global.tran >= n.rej.stop*2), 0.5*med.n.global.tran*inv.n.perm.completed, (0.5*med.n.global.tran+1)*inv.n.perm.completed.1)
                                
                                if (!is.null(OR)) {
                                    if (is.null(p.otu.tran.null.OR)) {
                                        mat = matrix(otu.tran.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                        p.otu.tran.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) # pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                                    }
                                    p.otu.tran.null.OR <- p.otu.tran.null.OR*inv.n.perm.completed
                                    
                                    med.T.otu.tran.OR = colMaxs(p.otu.tran.OR) # I_j
                                    
                                    med.T.otu.tran.null1.OR = pmax(p.otu.tran.null.OR[1,,], p.otu.tran.OR[2,]) # dim: J*B
                                    med.T.otu.tran.null2.OR = pmax(p.otu.tran.null.OR[2,,], p.otu.tran.OR[1,]) # dim: J*B
                                    med.T.otu.tran.null3.OR = pmax(p.otu.tran.null.OR[1,,], p.otu.tran.null.OR[2,,]) # dim: J*B
                                    med.T.otu.tran.null.OR = pmin(med.T.otu.tran.null1.OR, med.T.otu.tran.null2.OR, med.T.otu.tran.null3.OR) # I_j^(b)
                                    
                                    med.global.tran.OR = sum(1/med.T.otu.tran.OR) # I
                                    med.global.tran.perm.OR = colSums(1/med.T.otu.tran.null.OR) # I^(b)
                                    
                                    med.n.global.tran.OR = sum(med.global.tran.perm.OR > med.global.tran.OR - tol.eq) + sum(med.global.tran.perm.OR > med.global.tran.OR + tol.eq)
                                    med.p.global.tran.OR = ifelse((med.n.global.tran.OR >= n.rej.stop*2), 0.5*med.n.global.tran.OR*inv.n.perm.completed, (0.5*med.n.global.tran.OR+1)*inv.n.perm.completed.1)
                                    
                                    # combination test
                                    med.p.global.tran.tmp <- 0.5*med.n.global.tran
                                    med.p.global.tran.tmp.OR <- 0.5*med.n.global.tran.OR
                                    med.p.global.tran.null <- n.perm.completed + 0.5 - rank(med.global.tran.perm[1:n.perm.completed])
                                    med.p.global.tran.null.OR <- n.perm.completed + 0.5 - rank(med.global.tran.perm.OR[1:n.perm.completed])
                                    
                                    med.pmin.global.tran.com <- pmin(med.p.global.tran.tmp, med.p.global.tran.tmp.OR)
                                    med.pmin.global.tran.null.com <- pmin(med.p.global.tran.null, med.p.global.tran.null.OR)
                                    med.n.global.tran.com <- sum( (med.pmin.global.tran.null.com < c(med.pmin.global.tran.com) - tol.eq) + 0.5 * (abs(med.pmin.global.tran.null.com - c(med.pmin.global.tran.com)) < tol.eq))
                                    med.p.global.tran.com = ifelse((med.n.global.tran.com >= n.rej.stop), med.n.global.tran.com*inv.n.perm.completed, (med.n.global.tran.com+1)*inv.n.perm.completed.1)
                                }
                                
                            } # if (!freq.scale.only)
                            
                        } # if (!all.rarefy | test.omni3)
                        
                        if (all.rarefy | test.omni3) {
                            
                            # pa
                            
                            if (is.null(p.otu.pa.null)) {
                                mat = matrix(otu.pa.perm[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                p.otu.pa.null <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) # pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                            }
                            p.otu.pa.null <- p.otu.pa.null*inv.n.perm.completed
                            
                            med.T.otu.pa = colMaxs(p.otu.pa)
                            
                            med.T.otu.pa.null1 = pmax(p.otu.pa.null[1,,], p.otu.pa[2,]) # dim: J*B
                            med.T.otu.pa.null2 = pmax(p.otu.pa.null[2,,], p.otu.pa[1,]) # dim: J*B
                            med.T.otu.pa.null3 = pmax(p.otu.pa.null[1,,], p.otu.pa.null[2,,]) # dim: J*B
                            med.T.otu.pa.null = pmin(med.T.otu.pa.null1, med.T.otu.pa.null2, med.T.otu.pa.null3) # dim: J*B
                            
                            med.global.pa = sum(1/med.T.otu.pa)
                            med.global.pa.perm = colSums(1/med.T.otu.pa.null) # length: B
                            
                            med.n.global.pa = sum(med.global.pa.perm > med.global.pa - tol.eq) + sum(med.global.pa.perm > med.global.pa + tol.eq)
                            med.p.global.pa = ifelse((med.n.global.pa >= n.rej.stop*2), 0.5*med.n.global.pa*inv.n.perm.completed, (0.5*med.n.global.pa+1)*inv.n.perm.completed.1)
                            
                            if (!is.null(OR)) {
                                if (is.null(p.otu.pa.null.OR)) {
                                    mat = matrix(otu.pa.perm.OR[,,1:n.perm.completed,drop=FALSE], nrow=n.perm.completed, byrow=TRUE)
                                    p.otu.pa.null.OR <- n.perm.completed + 0.5 - array(colRanks(mat), c(n.var1, n.otu.smallp, n.perm.completed)) # pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                                }
                                p.otu.pa.null.OR <- p.otu.pa.null.OR*inv.n.perm.completed
                                
                                med.T.otu.pa.OR = colMaxs(p.otu.pa.OR)
                                
                                med.T.otu.pa.null1.OR = pmax(p.otu.pa.null.OR[1,,], p.otu.pa.OR[2,]) # dim: J*B
                                med.T.otu.pa.null2.OR = pmax(p.otu.pa.null.OR[2,,], p.otu.pa.OR[1,]) # dim: J*B
                                med.T.otu.pa.null3.OR = pmax(p.otu.pa.null.OR[1,,], p.otu.pa.null.OR[2,,]) # dim: J*B
                                med.T.otu.pa.null.OR = pmin(med.T.otu.pa.null1.OR, med.T.otu.pa.null2.OR, med.T.otu.pa.null3.OR) # dim: J*B
                                
                                med.global.pa.OR = sum(1/med.T.otu.pa.OR)
                                med.global.pa.perm.OR = colSums(1/med.T.otu.pa.null.OR) # length: B
                                
                                med.n.global.pa.OR = sum(med.global.pa.perm.OR > med.global.pa.OR - tol.eq) + sum(med.global.pa.perm.OR > med.global.pa.OR + tol.eq)
                                med.p.global.pa.OR = ifelse((med.n.global.pa.OR >= n.rej.stop*2), 0.5*med.n.global.pa.OR*inv.n.perm.completed, (0.5*med.n.global.pa.OR+1)*inv.n.perm.completed.1)
                                
                                # combination test
                                med.p.global.pa.tmp <- 0.5*med.n.global.pa
                                med.p.global.pa.tmp.OR <- 0.5*med.n.global.pa.OR
                                med.p.global.pa.null <- n.perm.completed + 0.5 - rank(med.global.pa.perm[1:n.perm.completed])
                                med.p.global.pa.null.OR <- n.perm.completed + 0.5 - rank(med.global.pa.perm.OR[1:n.perm.completed])
                                
                                med.pmin.global.pa.com <- pmin(med.p.global.pa.tmp, med.p.global.pa.tmp.OR)
                                med.pmin.global.pa.null.com <- pmin(med.p.global.pa.null, med.p.global.pa.null.OR)
                                med.n.global.pa.com <- sum( (med.pmin.global.pa.null.com < c(med.pmin.global.pa.com) - tol.eq) + 0.5 * (abs(med.pmin.global.pa.null.com - c(med.pmin.global.pa.com)) < tol.eq))
                                med.p.global.pa.com = ifelse((med.n.global.pa.com >= n.rej.stop), med.n.global.pa.com*inv.n.perm.completed, (med.n.global.pa.com+1)*inv.n.perm.completed.1)
                            }
                        } # if (all.rarefy | test.omni3)
                        
                        if (test.omni3 | (!all.rarefy & !freq.scale.only)) {
                            
                            # omni
                            
                            if (is.null(med.p.global.freq.tmp)) med.p.global.freq.tmp <- 0.5*med.n.global.freq 
                            if (is.null(med.p.global.tran.tmp)) med.p.global.tran.tmp <- 0.5*med.n.global.tran
                            
                            med.pmin.global.omni <- min(med.p.global.freq.tmp, med.p.global.tran.tmp)
                            if (is.null(med.p.global.freq.null)) med.p.global.freq.null <- n.perm.completed + 0.5 - rank(med.global.freq.perm[1:n.perm.completed])
                            if (is.null(med.p.global.tran.null)) med.p.global.tran.null <- n.perm.completed + 0.5 - rank(med.global.tran.perm[1:n.perm.completed])
                            med.pmin.global.omni.null <- pmin(med.p.global.freq.null, med.p.global.tran.null)
                            med.n.global.omni <- sum(med.pmin.global.omni.null< med.pmin.global.omni + tol.eq) + 0.5*sum(abs(med.pmin.global.omni.null - med.pmin.global.omni) < tol.eq)
                            med.p.global.omni = ifelse((med.n.global.omni >= n.rej.stop), med.n.global.omni*inv.n.perm.completed, (med.n.global.omni+1)*inv.n.perm.completed.1)
                            
                            if (!is.null(OR)) {
                                if (is.null(med.p.global.freq.tmp.OR)) med.p.global.freq.tmp.OR <- 0.5*med.n.global.freq.OR 
                                if (is.null(med.p.global.tran.tmp.OR)) med.p.global.tran.tmp.OR <- 0.5*med.n.global.tran.OR
                                
                                med.pmin.global.omni.OR <- min(med.p.global.freq.tmp.OR, med.p.global.tran.tmp.OR)
                                if (is.null(med.p.global.freq.null.OR)) med.p.global.freq.null.OR <- n.perm.completed + 0.5 - rank(med.global.freq.perm.OR[1:n.perm.completed])
                                if (is.null(med.p.global.tran.null.OR)) med.p.global.tran.null.OR <- n.perm.completed + 0.5 - rank(med.global.tran.perm.OR[1:n.perm.completed])
                                med.pmin.global.omni.null.OR <- pmin(med.p.global.freq.null.OR, med.p.global.tran.null.OR)
                                med.n.global.omni.OR <- sum(med.pmin.global.omni.null.OR< med.pmin.global.omni.OR + tol.eq) + 0.5*sum(abs(med.pmin.global.omni.null.OR - med.pmin.global.omni.OR) < tol.eq)
                                med.p.global.omni.OR = ifelse((med.n.global.omni.OR >= n.rej.stop), med.n.global.omni.OR*inv.n.perm.completed, (med.n.global.omni.OR+1)*inv.n.perm.completed.1)
                                
                                # combination test
                                med.pmin.global.omni.com <- pmin(med.pmin.global.omni, med.pmin.global.omni.OR)
                                med.pmin.global.omni.null.com <- pmin(med.pmin.global.omni.null, med.pmin.global.omni.null.OR)
                                med.n.global.omni.com <- sum( (med.pmin.global.omni.null.com < c(med.pmin.global.omni.com) - tol.eq) + 0.5 * (abs(med.pmin.global.omni.null.com - c(med.pmin.global.omni.com)) < tol.eq))
                                med.p.global.omni.com = ifelse((med.n.global.omni.com >= n.rej.stop), med.n.global.omni.com*inv.n.perm.completed, (med.n.global.omni.com+1)*inv.n.perm.completed.1)
                            }
                            
                            if (test.omni3) {
                                
                                med.pmin.otu.omni3 <- pmin(med.T.otu.freq, med.T.otu.tran, med.T.otu.pa)
                                med.pmin.otu.omni3.null <- pmin(med.T.otu.freq.null, med.T.otu.tran.null, med.T.otu.pa.null)
                                
                                # harmonic mean
                                med.global.harmonic = sum(1/med.pmin.otu.omni3)
                                med.global.harmonic.perm = colSums(1/med.pmin.otu.omni3.null) # B
                                med.n.global.harmonic <- sum((med.global.harmonic.perm > med.global.harmonic + tol.eq) + (med.global.harmonic.perm > med.global.harmonic - tol.eq))
                                med.p.global.harmonic = ifelse((med.n.global.harmonic >= n.rej.stop*2), 0.5*med.n.global.harmonic*inv.n.perm.completed, (0.5*med.n.global.harmonic+1)*inv.n.perm.completed.1)
                                
                                # fisher
                                med.global.fisher = sum(-log(med.pmin.otu.omni3))
                                med.global.fisher.perm = colSums(-log(med.pmin.otu.omni3.null))
                                med.n.global.fisher <- sum((med.global.fisher.perm > med.global.fisher + tol.eq) + (med.global.fisher.perm > med.global.fisher - tol.eq))
                                med.p.global.fisher = ifelse((med.n.global.fisher >= n.rej.stop*2), 0.5*med.n.global.fisher*inv.n.perm.completed, (0.5*med.n.global.fisher+1)*inv.n.perm.completed.1)
                                
                                # omni3
                                if (is.null(med.p.global.pa.tmp)) med.p.global.pa.tmp <- 0.5*med.n.global.pa
                                med.p.global.harmonic.tmp <- 0.5*med.n.global.harmonic
                                med.p.global.fisher.tmp <- 0.5*med.n.global.fisher
                                if (is.null(med.p.global.pa.null)) med.p.global.pa.null <- n.perm.completed + 0.5 - rank(med.global.pa.perm[1:n.perm.completed])
                                med.p.global.harmonic.null <- n.perm.completed + 0.5 - rank(med.global.harmonic.perm)
                                med.p.global.fisher.null <- n.perm.completed + 0.5 - rank(med.global.fisher.perm)
                                
                                med.pmin.global.omni3 <- min(med.p.global.freq.tmp, med.p.global.tran.tmp, med.p.global.pa.tmp,
                                                             med.p.global.harmonic.tmp, med.p.global.fisher.tmp)
                                med.pmin.global.omni3.null <- pmin(med.p.global.freq.null, med.p.global.tran.null, med.p.global.pa.null,
                                                                   med.p.global.harmonic.null, med.p.global.fisher.null)
                                med.n.global.omni3 <- sum( (med.pmin.global.omni3.null < c(med.pmin.global.omni3) - tol.eq) + 0.5 * (abs(med.pmin.global.omni3.null - c(med.pmin.global.omni3)) < tol.eq))
                                med.p.global.omni3 = ifelse((med.n.global.omni3 >= n.rej.stop), med.n.global.omni3*inv.n.perm.completed, (med.n.global.omni3+1)*inv.n.perm.completed.1)
                                
                                if (!is.null(OR)) {
                                    med.pmin.otu.omni3.OR <- pmin(med.T.otu.freq.OR, med.T.otu.tran.OR, med.T.otu.pa.OR)
                                    med.pmin.otu.omni3.null.OR <- pmin(med.T.otu.freq.null.OR, med.T.otu.tran.null.OR, med.T.otu.pa.null.OR)
                                    
                                    # harmonic mean
                                    med.global.harmonic.OR = sum(1/med.pmin.otu.omni3.OR)
                                    med.global.harmonic.perm.OR = colSums(1/med.pmin.otu.omni3.null.OR) # B
                                    med.n.global.harmonic.OR <- sum((med.global.harmonic.perm.OR > med.global.harmonic.OR + tol.eq) + (med.global.harmonic.perm.OR > med.global.harmonic.OR - tol.eq))
                                    med.p.global.harmonic.OR = ifelse((med.n.global.harmonic.OR >= n.rej.stop*2), 0.5*med.n.global.harmonic.OR*inv.n.perm.completed, (0.5*med.n.global.harmonic.OR+1)*inv.n.perm.completed.1)
                                    
                                    # fisher
                                    med.global.fisher.OR = sum(-log(med.pmin.otu.omni3.OR))
                                    med.global.fisher.perm.OR = colSums(-log(med.pmin.otu.omni3.null.OR))
                                    med.n.global.fisher.OR <- sum((med.global.fisher.perm.OR > med.global.fisher.OR + tol.eq) + (med.global.fisher.perm.OR > med.global.fisher.OR - tol.eq))
                                    med.p.global.fisher.OR = ifelse((med.n.global.fisher.OR >= n.rej.stop*2), 0.5*med.n.global.fisher.OR*inv.n.perm.completed, (0.5*med.n.global.fisher.OR+1)*inv.n.perm.completed.1)
                                    
                                    # omni3
                                    if (is.null(med.p.global.pa.tmp.OR)) med.p.global.pa.tmp.OR <- 0.5*med.n.global.pa.OR
                                    med.p.global.harmonic.tmp.OR <- 0.5*med.n.global.harmonic.OR
                                    med.p.global.fisher.tmp.OR <- 0.5*med.n.global.fisher.OR
                                    if (is.null(med.p.global.pa.null.OR)) med.p.global.pa.null.OR <- n.perm.completed + 0.5 - rank(med.global.pa.perm.OR[1:n.perm.completed])
                                    med.p.global.harmonic.null.OR <- n.perm.completed + 0.5 - rank(med.global.harmonic.perm.OR)
                                    med.p.global.fisher.null.OR <- n.perm.completed + 0.5 - rank(med.global.fisher.perm.OR)
                                    
                                    med.pmin.global.omni3.OR <- min(med.p.global.freq.tmp.OR, med.p.global.tran.tmp.OR, med.p.global.pa.tmp.OR,
                                                                    med.p.global.harmonic.tmp.OR, med.p.global.fisher.tmp.OR)
                                    med.pmin.global.omni3.null.OR <- pmin(med.p.global.freq.null.OR, med.p.global.tran.null.OR, med.p.global.pa.null.OR,
                                                                          med.p.global.harmonic.null.OR, med.p.global.fisher.null.OR)
                                    med.n.global.omni3.OR <- sum( (med.pmin.global.omni3.null.OR < c(med.pmin.global.omni3.OR) - tol.eq) + 0.5 * (abs(med.pmin.global.omni3.null.OR - c(med.pmin.global.omni3.OR)) < tol.eq))
                                    med.p.global.omni3.OR = ifelse((med.n.global.omni3.OR >= n.rej.stop), med.n.global.omni3.OR*inv.n.perm.completed, (med.n.global.omni3.OR+1)*inv.n.perm.completed.1)
                                    
                                    # combination test
                                    med.pmin.global.harmonic.com <- pmin(med.p.global.harmonic, med.p.global.harmonic.OR)
                                    med.pmin.global.harmonic.null.com <- pmin(med.p.global.harmonic.null, med.p.global.harmonic.null.OR)
                                    med.n.global.harmonic.com <- sum( (med.pmin.global.harmonic.null.com < c(med.pmin.global.harmonic.com) - tol.eq) + 0.5 * (abs(med.pmin.global.harmonic.null.com - c(med.pmin.global.harmonic.com)) < tol.eq))
                                    med.p.global.harmonic.com = ifelse((med.n.global.harmonic.com >= n.rej.stop), med.n.global.harmonic.com*inv.n.perm.completed, (med.n.global.harmonic.com+1)*inv.n.perm.completed.1)
                                    
                                    med.pmin.global.fisher.com <- pmin(med.p.global.fisher, med.p.global.fisher.OR)
                                    med.pmin.global.fisher.null.com <- pmin(med.p.global.fisher.null, med.p.global.fisher.null.OR)
                                    med.n.global.fisher.com <- sum( (med.pmin.global.fisher.null.com < c(med.pmin.global.fisher.com) - tol.eq) + 0.5 * (abs(med.pmin.global.fisher.null.com - c(med.pmin.global.fisher.com)) < tol.eq))
                                    med.p.global.fisher.com = ifelse((med.n.global.fisher.com >= n.rej.stop), med.n.global.fisher.com*inv.n.perm.completed, (med.n.global.fisher.com+1)*inv.n.perm.completed.1)
                                    
                                    med.pmin.global.omni3.com <- pmin(med.pmin.global.omni3, med.pmin.global.omni3.OR)
                                    med.pmin.global.omni3.null.com <- pmin(med.pmin.global.omni3.null, med.pmin.global.omni3.null.OR)
                                    med.n.global.omni3.com <- sum( (med.pmin.global.omni3.null.com < c(med.pmin.global.omni3.com) - tol.eq) + 0.5 * (abs(med.pmin.global.omni3.null.com - c(med.pmin.global.omni3.com)) < tol.eq))
                                    med.p.global.omni3.com = ifelse((med.n.global.omni3.com >= n.rej.stop), med.n.global.omni3.com*inv.n.perm.completed, (med.n.global.omni3.com+1)*inv.n.perm.completed.1)
                                    
                                }
                            } # if (test.omni3)
                            
                        } # if (test.omni3 | (!all.rarefy & !freq.scale.only))
                        
                        
                        if (test.omni3) {
                            med.meet.all.rej.stop = all(med.n.global.freq >= n.rej.stop*2) & all(med.n.global.tran >= n.rej.stop*2) & all(med.n.global.pa >= n.rej.stop*2) & all(med.n.global.omni3 >= n.rej.stop)
                            if (!is.null(OR)) med.meet.all.rej.stop = med.meet.all.rej.stop & all(med.n.global.freq.OR >= n.rej.stop*2) & all(med.n.global.tran.OR >= n.rej.stop*2) & all(med.n.global.pa.OR >= n.rej.stop*2) & all(med.n.global.omni3.OR >= n.rej.stop)
                        } else if (all.rarefy) {
                            med.meet.all.rej.stop = all(med.n.global.pa >= n.rej.stop*2)
                            if (!is.null(OR)) med.meet.all.rej.stop = med.meet.all.rej.stop & all(med.n.global.pa.OR >= n.rej.stop*2)
                        } else if (freq.scale.only) {
                            med.meet.all.rej.stop = all(med.n.global.freq >= n.rej.stop*2)
                            if (!is.null(OR)) med.meet.all.rej.stop = med.meet.all.rej.stop & all(med.n.global.freq.OR >= n.rej.stop*2)
                        } else {
                            med.meet.all.rej.stop = all(med.n.global.freq >= n.rej.stop*2) & all(med.n.global.tran >= n.rej.stop*2) & all(med.n.global.omni >= n.rej.stop)
                            if (!is.null(OR)) med.meet.all.rej.stop = med.meet.all.rej.stop & all(med.n.global.freq.OR >= n.rej.stop*2) & all(med.n.global.tran.OR >= n.rej.stop*2) & all(med.n.global.omni.OR >= n.rej.stop)
                        }
                        
                    } # test.mediation
                    
                    if (meet.all.rej.stop & med.meet.all.rej.stop) {
                        global.tests.stopped = TRUE
                        cat("global test stopped at permutation", n.perm.completed, "\n")
                        n.global.perm.completed = n.perm.completed
                    }
                    
                    if (global.tests.stopped | n.perm.completed >= n.global.perm.max) {
                        global.tests.done = TRUE
                    }
                    
                } # meet rej stop
            } # !global.tests.done
            
            if (global.tests.done + otu.tests.stopped == 2) break
            
        } # permutation
        
    } # if (n.perm.max > 0)
    
    v.freq = NULL
    d.freq = NULL
    v.tran = NULL
    d.tran = NULL
    if (!all.rarefy) {
        v.freq = fit.ldm$v.freq
        d.freq = fit.ldm$d.freq
        if (!freq.scale.only) {
            v.tran=fit.ldm$v.tran
            d.tran=fit.ldm$d.tran
        }
    }
    
    otu.names <- colnames(otu.table)
    if (!all.rarefy) {
        colnames(ldm.obs.freq$ve.otu) <- otu.names
        if (!freq.scale.only) colnames(ldm.obs.tran$ve.otu) <- otu.names
    }
    if (all.rarefy) {
        colnames(ldm.obs.pa$ve.otu) <- otu.names
    }
    
    # confounder
    
    VE.df.confounders = NULL
    VE.global.freq.confounders = NULL
    VE.otu.freq.confounders = NULL
    VE.global.tran.confounders = NULL
    VE.otu.tran.confounders = NULL
    
    if (!all.rarefy) {
        if (adjust.for.confounders) {
            i.conf <- fit.ldm$low[1]:fit.ldm$up[1]
            VE.df.confounders <- length(i.conf)
            
            VE.global.freq.confounders <- sum((fit.ldm$d.freq[i.conf])^2)
            wt <- fit.ldm$d.freq[i.conf] * t(fit.ldm$v.freq[, i.conf]) 
            if (is.vector(wt)) VE.otu.freq.confounders = wt^2
            else               VE.otu.freq.confounders = colSums(wt^2)
            
            if (!freq.scale.only) {
                VE.global.tran.confounders <- sum((fit.ldm$d.tran[i.conf])^2)
                wt <- fit.ldm$d.tran[i.conf] * t(fit.ldm$v.tran[, i.conf]) 
                if (is.vector(wt)) VE.otu.tran.confounders = wt^2
                else               VE.otu.tran.confounders = colSums(wt^2)
            }
        }
    }
    
    # submodels
    
    VE.df.submodels = NULL
    VE.global.freq.submodels = NULL
    VE.otu.freq.submodels = NULL
    VE.global.tran.submodels = NULL
    VE.otu.tran.submodels = NULL
    VE.global.freq.residuals <- NULL
    VE.global.tran.residuals <- NULL
    i.col = NULL
    beta = NULL
    beta.name = NULL
    
    if (!all.rarefy) {
        fit <- fit.ldm
    } else {
        fit <- fit.ldm.pa
    }
    
    for (k in 1:n.var1) {
        k1 = k + as.numeric(adjust.for.confounders)
        i.m <- fit$low[k1]:fit$up[k1]
        i.col = c(i.col, i.m)
        
        for (im in i.m) {
            corr = t(cor(model[[k1]], fit$x[,im]))
            
            w = which(corr[1,]>0.9999)
            if (length(w)==0) {
                w = which.max(abs(corr[1,]))
            }
            if (corr[1,w] < 0) {
                fit$x[,im] = -fit$x[,im]
                if (!all.rarefy) {
                    fit.ldm$v.freq[,im] = -fit.ldm$v.freq[,im]
                }
            }
            beta.name = c(beta.name, colnames(model[[k1]])[w])
            if (!all.rarefy) {
                beta = rbind(beta, t(fit.ldm$v.freq[,im])*fit.ldm$d.freq[im])
            } else {
                beta = rbind(beta, -t(fit.ldm.pa$x[,im]) %*% fit.ldm.pa$phi)
            }
        }
        
        if (!all.rarefy) {
            VE.df.submodels <- c(VE.df.submodels, length(i.m))
            
            VE.global.freq.submodels <- c(VE.global.freq.submodels, sum((fit.ldm$d.freq[i.m])^2))
            wt <- fit.ldm$d.freq[i.m] * t(fit.ldm$v.freq[, i.m])
            if (is.vector(wt)) VE.otu.freq.submodels = rbind(VE.otu.freq.submodels, wt^2)
            else               VE.otu.freq.submodels = rbind(VE.otu.freq.submodels, colSums(wt^2))
            
            if (!freq.scale.only) {
                VE.global.tran.submodels <- c(VE.global.tran.submodels, sum((fit.ldm$d.tran[i.m])^2))
                wt <- fit.ldm$d.tran[i.m] * t(fit.ldm$v.tran[, i.m])
                if (is.vector(wt)) VE.otu.tran.submodels = rbind(VE.otu.tran.submodels, wt^2)
                else               VE.otu.tran.submodels = rbind(VE.otu.tran.submodels, colSums(wt^2))
            }
        }
    }
    
    beta = data.frame(beta)
    rownames(beta) = beta.name    # Possible error: duplicate 'row.names' are not allowed
    colnames(beta) = colnames(otu.table)
    
    phi = NULL
    if (all.rarefy) phi = fit.ldm.pa$phi
    
    # residuals
    
    if (!all.rarefy) {
        i.all <- fit.ldm$low[1]:fit.ldm$up[n.var]
        
        VE.global.freq.residuals <- fit.ldm$d.freq[-i.all]^2
        VE.global.tran.residuals <- NULL
        if (!freq.scale.only) VE.global.tran.residuals <- fit.ldm$d.tran[-i.all]^2
    }
    
    cov.names = paste("cov", 1:n.var1, sep="")
    rarefy.names = paste("rarefy", 1:n.rarefy, sep="")
    
    p.global.names = list(cov.names, rarefy.names)
    p.otu.names = list(cov.names, otu.names)
    
    if (!is.null(ldm.obs.freq)) dimnames(ldm.obs.freq$ve.global) = p.global.names
    if (!is.null(ldm.obs.tran)) dimnames(ldm.obs.tran$ve.global) = p.global.names
    if (!is.null(ldm.obs.pa)) names(ldm.obs.pa$ve.global) = cov.names
    if (!is.null(ldm.obs.freq)) dimnames(ldm.obs.freq$ve.otu) = p.otu.names
    if (!is.null(ldm.obs.tran)) dimnames(ldm.obs.tran$ve.otu) = p.otu.names
    if (!is.null(ldm.obs.pa)) names(ldm.obs.pa$ve.otu) = p.otu.names
    
    if (!is.null(p.global.freq)) dimnames(p.global.freq) = p.global.names
    if (!is.null(p.global.tran)) dimnames(p.global.tran) = p.global.names 
    if (!is.null(p.global.omni)) dimnames(p.global.omni) = p.global.names
    if (!is.null(p.global.pa)) names(p.global.pa) = cov.names
    if (!is.null(p.global.harmonic)) names(p.global.harmonic) = cov.names
    if (!is.null(p.global.fisher)) names(p.global.fisher) = cov.names
    if (!is.null(p.global.omni3)) names(p.global.omni3) = cov.names 
    
    if (!is.null(p.global.freq.OR)) dimnames(p.global.freq.OR) = p.global.names 
    if (!is.null(p.global.tran.OR)) dimnames(p.global.tran.OR) = p.global.names 
    if (!is.null(p.global.omni.OR)) dimnames(p.global.omni.OR) = p.global.names
    if (!is.null(p.global.pa.OR)) names(p.global.pa.OR) = cov.names
    if (!is.null(p.global.harmonic.OR)) names(p.global.harmonic.OR) = cov.names
    if (!is.null(p.global.fisher.OR)) names(p.global.fisher.OR) = cov.names
    if (!is.null(p.global.omni3.OR)) names(p.global.omni3.OR) = cov.names 
    
    if (!is.null(p.global.freq.com)) names(p.global.freq.com) = cov.names 
    if (!is.null(p.global.tran.com)) names(p.global.tran.com) = cov.names 
    if (!is.null(p.global.omni.com)) names(p.global.omni.com) = cov.names
    if (!is.null(p.global.pa.com)) names(p.global.pa.com) = cov.names
    if (!is.null(p.global.harmonic.com)) names(p.global.harmonic.com) = cov.names
    if (!is.null(p.global.fisher.com)) names(p.global.fisher.com) = cov.names
    if (!is.null(p.global.omni3.com)) names(p.global.omni3.com) = cov.names 
    
    if (!is.null(p.otu.freq)) dimnames(p.otu.freq) = p.otu.names
    if (!is.null(p.otu.tran)) dimnames(p.otu.tran) = p.otu.names
    if (!is.null(p.otu.pa)) dimnames(p.otu.pa) = p.otu.names
    if (!is.null(p.otu.omni)) dimnames(p.otu.omni) = p.otu.names
    if (!is.null(p.otu.omni3)) dimnames(p.otu.omni3) = p.otu.names
    if (!is.null(q.otu.freq)) {
        if (n.otu > 1) q.otu.freq <- t(apply(p.otu.freq, 1, p.adjust, method="BH"))
        dimnames(q.otu.freq) = p.otu.names
    }
    if (!is.null(q.otu.tran)) dimnames(q.otu.tran) = p.otu.names
    if (!is.null(q.otu.pa)) dimnames(q.otu.pa) = p.otu.names
    if (!is.null(q.otu.omni)) dimnames(q.otu.omni) = p.otu.names
    if (!is.null(q.otu.omni3)) dimnames(q.otu.omni3) = p.otu.names
    
    if (!is.null(OR)) {
        if (!is.null(p.otu.freq.OR)) dimnames(p.otu.freq.OR) = p.otu.names
        if (!is.null(p.otu.tran.OR)) dimnames(p.otu.tran.OR) = p.otu.names
        if (!is.null(p.otu.pa.OR)) dimnames(p.otu.pa.OR) = p.otu.names
        if (!is.null(p.otu.omni.OR)) dimnames(p.otu.omni.OR) = p.otu.names
        if (!is.null(p.otu.omni3.OR)) dimnames(p.otu.omni3.OR) = p.otu.names
        if (!is.null(q.otu.freq.OR)) {
            if (n.otu > 1) q.otu.freq.OR <- t(apply(p.otu.freq.OR, 1, p.adjust, method="BH"))
            dimnames(q.otu.freq.OR) = p.otu.names
        }
        if (!is.null(q.otu.tran.OR)) dimnames(q.otu.tran.OR) = p.otu.names
        if (!is.null(q.otu.pa.OR)) dimnames(q.otu.pa.OR) = p.otu.names
        if (!is.null(q.otu.omni.OR)) dimnames(q.otu.omni.OR) = p.otu.names
        if (!is.null(q.otu.omni3.OR)) dimnames(q.otu.omni3.OR) = p.otu.names
        if (!is.null(p.otu.freq.com)) dimnames(p.otu.freq.com) = p.otu.names
        if (!is.null(p.otu.tran.com)) dimnames(p.otu.tran.com) = p.otu.names
        if (!is.null(p.otu.pa.com)) dimnames(p.otu.pa.com) = p.otu.names
        if (!is.null(p.otu.omni.com)) dimnames(p.otu.omni.com) = p.otu.names
        if (!is.null(p.otu.omni3.com)) dimnames(p.otu.omni3.com) = p.otu.names
        if (!is.null(q.otu.freq.com)) {
            if (n.otu > 1) q.otu.freq.com <- t(apply(p.otu.freq.com, 1, p.adjust, method="BH"))
            dimnames(q.otu.freq.com) = p.otu.names
        }
        if (!is.null(q.otu.tran.com)) dimnames(q.otu.tran.com) = p.otu.names
        if (!is.null(q.otu.pa.com)) dimnames(q.otu.pa.com) = p.otu.names
        if (!is.null(q.otu.omni.com)) dimnames(q.otu.omni.com) = p.otu.names
        if (!is.null(q.otu.omni3.com)) dimnames(q.otu.omni3.com) = p.otu.names
    }
    
    
    detected.otu.freq = list()
    detected.otu.tran = list()
    detected.otu.pa = list()
    detected.otu.omni = list()
    detected.otu.omni3 = list()
    detected.otu.freq.OR = list()
    detected.otu.tran.OR = list()
    detected.otu.pa.OR = list()
    detected.otu.omni.OR = list()
    detected.otu.omni3.OR = list()
    detected.otu.freq.com = list()
    detected.otu.tran.com = list()
    detected.otu.pa.com = list()
    detected.otu.omni.com = list()
    detected.otu.omni3.com = list()
    if (!is.null(q.otu.freq)) {for (k in 1:n.var1) detected.otu.freq[[k]] = colnames(q.otu.freq)[which(q.otu.freq[k,]<fdr.nominal)]; names(detected.otu.freq) = cov.names}
    if (!is.null(q.otu.tran)) {for (k in 1:n.var1) detected.otu.tran[[k]] = colnames(q.otu.tran)[which(q.otu.tran[k,]<fdr.nominal)]; names(detected.otu.tran) = cov.names}
    if (!is.null(q.otu.pa)) {for (k in 1:n.var1) detected.otu.pa[[k]] = colnames(q.otu.pa)[which(q.otu.pa[k,]<fdr.nominal)]; names(detected.otu.pa) = cov.names}
    if (!is.null(q.otu.omni)) {for (k in 1:n.var1) detected.otu.omni[[k]] = colnames(q.otu.omni)[which(q.otu.omni[k,]<fdr.nominal)]; names(detected.otu.omni) = cov.names}
    if (!is.null(q.otu.omni3)) {for (k in 1:n.var1) detected.otu.omni3[[k]] = colnames(q.otu.omni3)[which(q.otu.omni3[k,]<fdr.nominal)]; names(detected.otu.omni3) = cov.names}
    if (!is.null(OR)) {
        if (!is.null(q.otu.freq.OR)) {for (k in 1:n.var1) detected.otu.freq.OR[[k]] = colnames(q.otu.freq.OR)[which(q.otu.freq.OR[k,]<fdr.nominal)]; names(detected.otu.freq.OR) = cov.names}
        if (!is.null(q.otu.tran.OR)) {for (k in 1:n.var1) detected.otu.tran.OR[[k]] = colnames(q.otu.tran.OR)[which(q.otu.tran.OR[k,]<fdr.nominal)]; names(detected.otu.tran.OR) = cov.names}
        if (!is.null(q.otu.pa.OR)) {for (k in 1:n.var1) detected.otu.pa.OR[[k]] = colnames(q.otu.pa.OR)[which(q.otu.pa.OR[k,]<fdr.nominal)]; names(detected.otu.pa.OR) = cov.names}
        if (!is.null(q.otu.omni.OR)) {for (k in 1:n.var1) detected.otu.omni.OR[[k]] = colnames(q.otu.omni.OR)[which(q.otu.omni.OR[k,]<fdr.nominal)]; names(detected.otu.omni.OR) = cov.names}
        if (!is.null(q.otu.omni3.OR)) {for (k in 1:n.var1) detected.otu.omni3.OR[[k]] = colnames(q.otu.omni3.OR)[which(q.otu.omni3.OR[k,]<fdr.nominal)]; names(detected.otu.omni3.OR) = cov.names}
        if (!is.null(q.otu.freq.com)) {for (k in 1:n.var1) detected.otu.freq.com[[k]] = colnames(q.otu.freq.com)[which(q.otu.freq.com[k,]<fdr.nominal)]; names(detected.otu.freq.com) = cov.names}
        if (!is.null(q.otu.tran.com)) {for (k in 1:n.var1) detected.otu.tran.com[[k]] = colnames(q.otu.tran.com)[which(q.otu.tran.com[k,]<fdr.nominal)]; names(detected.otu.tran.com) = cov.names}
        if (!is.null(q.otu.pa.com)) {for (k in 1:n.var1) detected.otu.pa.com[[k]] = colnames(q.otu.pa.com)[which(q.otu.pa.com[k,]<fdr.nominal)]; names(detected.otu.pa.com) = cov.names}
        if (!is.null(q.otu.omni.com)) {for (k in 1:n.var1) detected.otu.omni.com[[k]] = colnames(q.otu.omni.com)[which(q.otu.omni.com[k,]<fdr.nominal)]; names(detected.otu.omni.com) = cov.names}
        if (!is.null(q.otu.omni3.com)) {for (k in 1:n.var1) detected.otu.omni3.com[[k]] = colnames(q.otu.omni3.com)[which(q.otu.omni3.com[k,]<fdr.nominal)]; names(detected.otu.omni3.com) = cov.names}
    }
    # mediation
    
    med.detected.otu.freq = NULL
    med.detected.otu.tran = NULL
    med.detected.otu.pa = NULL
    med.detected.otu.omni = NULL
    med.detected.otu.omni3 = NULL
    med.detected.otu.freq.OR = NULL
    med.detected.otu.tran.OR = NULL
    med.detected.otu.pa.OR = NULL
    med.detected.otu.omni.OR = NULL
    med.detected.otu.omni3.OR = NULL
    med.detected.otu.freq.com = NULL
    med.detected.otu.tran.com = NULL
    med.detected.otu.pa.com = NULL
    med.detected.otu.omni.com = NULL
    med.detected.otu.omni3.com = NULL
    if (test.mediation){
        if (!is.null(p.otu.freq)) {
            med.q.otu.freq <- medTest.SBMH(p.otu.freq[1,], p.otu.freq[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
            med.detected.otu.freq = colnames(p.otu.freq)[which(med.q.otu.freq < fdr.nominal)]
        }
        if (!is.null(p.otu.tran)) {
            med.q.otu.tran <- medTest.SBMH(p.otu.tran[1,], p.otu.tran[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
            med.detected.otu.tran = colnames(p.otu.tran)[which(med.q.otu.tran < fdr.nominal)]
        }
        if (!is.null(p.otu.pa)) {
            med.q.otu.pa <- medTest.SBMH(p.otu.pa[1,], p.otu.pa[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
            med.detected.otu.pa = colnames(p.otu.pa)[which(med.q.otu.pa < fdr.nominal)]
        }
        if (!is.null(p.otu.omni)) {
            med.q.otu.omni <- medTest.SBMH(p.otu.omni[1,], p.otu.omni[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
            med.detected.otu.omni = colnames(p.otu.omni)[which(med.q.otu.omni < fdr.nominal)]
        }
        if (!is.null(p.otu.omni3)) {
            med.q.otu.omni3 <- medTest.SBMH(p.otu.omni3[1,], p.otu.omni3[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
            med.detected.otu.omni3 = colnames(p.otu.omni3)[which(med.q.otu.omni3 < fdr.nominal)]
        }
        if (!is.null(OR)) {
            if (!is.null(p.otu.freq.OR)) {
                med.q.otu.freq.OR <- medTest.SBMH(p.otu.freq.OR[1,], p.otu.freq.OR[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.freq.OR = colnames(p.otu.freq.OR)[which(med.q.otu.freq.OR < fdr.nominal)]
            }
            if (!is.null(p.otu.tran.OR)) {
                med.q.otu.tran.OR <- medTest.SBMH(p.otu.tran.OR[1,], p.otu.tran.OR[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.tran.OR = colnames(p.otu.tran.OR)[which(med.q.otu.tran.OR < fdr.nominal)]
            }
            if (!is.null(p.otu.pa.OR)) {
                med.q.otu.pa.OR <- medTest.SBMH(p.otu.pa.OR[1,], p.otu.pa.OR[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.pa.OR = colnames(p.otu.pa.OR)[which(med.q.otu.pa.OR < fdr.nominal)]
            }
            if (!is.null(p.otu.omni.OR)) {
                med.q.otu.omni.OR <- medTest.SBMH(p.otu.omni.OR[1,], p.otu.omni.OR[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.omni.OR = colnames(p.otu.omni.OR)[which(med.q.otu.omni.OR < fdr.nominal)]
            }
            if (!is.null(p.otu.omni3.OR)) {
                med.q.otu.omni3.OR <- medTest.SBMH(p.otu.omni3.OR[1,], p.otu.omni3.OR[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.omni3.OR = colnames(p.otu.omni3.OR)[which(med.q.otu.omni3.OR < fdr.nominal)]
            }
            if (!is.null(p.otu.freq.com)) {
                med.q.otu.freq.com <- medTest.SBMH(p.otu.freq.com[1,], p.otu.freq.com[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.freq.com = colnames(p.otu.freq.com)[which(med.q.otu.freq.com < fdr.nominal)]
            }
            if (!is.null(p.otu.tran.com)) {
                med.q.otu.tran.com <- medTest.SBMH(p.otu.tran.com[1,], p.otu.tran.com[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.tran.com = colnames(p.otu.tran.com)[which(med.q.otu.tran.com < fdr.nominal)]
            }
            if (!is.null(p.otu.pa.com)) {
                med.q.otu.pa.com <- medTest.SBMH(p.otu.pa.com[1,], p.otu.pa.com[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.pa.com = colnames(p.otu.pa.com)[which(med.q.otu.pa.com < fdr.nominal)]
            }
            if (!is.null(p.otu.omni.com)) {
                med.q.otu.omni.com <- medTest.SBMH(p.otu.omni.com[1,], p.otu.omni.com[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.omni.com = colnames(p.otu.omni.com)[which(med.q.otu.omni.com < fdr.nominal)]
            }
            if (!is.null(p.otu.omni3.com)) {
                med.q.otu.omni3.com <- medTest.SBMH(p.otu.omni3.com[1,], p.otu.omni3.com[2,], MCP.type="FDR", t1=fdr.nominal/2, t2=fdr.nominal/2)
                med.detected.otu.omni3.com = colnames(p.otu.omni3.com)[which(med.q.otu.omni3.com < fdr.nominal)]
            }
        }
    }
    
    if (binary & !all.rarefy) {
        ldm.obs.pa = ldm.obs.freq
        p.global.pa = p.global.freq
        p.otu.pa = p.otu.freq
        q.otu.pa = q.otu.freq
        detected.otu.pa = detected.otu.freq
        if (!is.null(OR)) {
            p.global.pa.OR = p.global.freq.OR
            p.global.pa.com = p.global.freq.com
            p.otu.pa.OR = p.otu.freq.OR
            p.otu.pa.com = p.otu.freq.com
            q.otu.pa.OR = q.otu.freq.OR
            q.otu.pa.com = q.otu.freq.com
            detected.otu.pa.OR = detected.otu.freq.OR
            detected.otu.pa.com = detected.otu.freq.com
        }
        if (test.mediation){
            med.p.global.pa = med.p.global.freq
            med.detected.otu.pa = med.detected.otu.freq
            if (!is.null(OR)) {
                med.p.global.pa.OR = med.p.global.freq.OR
                med.p.global.pa.com = med.p.global.freq.com
                med.detected.otu.pa.OR = med.detected.otu.freq.OR
                med.detected.otu.pa.com = med.detected.otu.freq.com
            }
        }
    }
    
    res = list( x=fit$x,
                dist=d.gower,
                mean.freq=mean.freq,
                y.freq=y.freq,
                v.freq=v.freq,
                d.freq=d.freq,
                y.tran=y.tran,
                v.tran=v.tran,
                d.tran=d.tran,
                low=fit$low,
                up=fit$up,
                beta=beta,
                phi=1-phi,
                VE.global.freq.confounders=drop(VE.global.freq.confounders),
                VE.global.freq.submodels=drop(VE.global.freq.submodels),
                VE.global.freq.residuals=drop(VE.global.freq.residuals),
                VE.otu.freq.confounders=drop(VE.otu.freq.confounders),
                VE.otu.freq.submodels=drop(VE.otu.freq.submodels),
                VE.global.tran.confounders=drop(VE.global.tran.confounders),
                VE.global.tran.submodels=drop(VE.global.tran.submodels),
                VE.global.tran.residuals=drop(VE.global.tran.residuals),
                VE.otu.tran.confounders=drop(VE.otu.tran.confounders),
                VE.otu.tran.submodels=drop(VE.otu.tran.submodels),
                VE.df.confounders=drop(VE.df.confounders),
                VE.df.submodels=drop(VE.df.submodels),
                
                F.global.freq=drop(ldm.obs.freq$ve.global),
                F.global.tran=drop(ldm.obs.tran$ve.global),
                F.global.pa=drop(ldm.obs.pa$ve.global),
                F.otu.freq=drop(ldm.obs.freq$ve.otu),
                F.otu.tran=drop(ldm.obs.tran$ve.otu),
                F.otu.pa=drop(ldm.obs.pa$ve.otu),
                
    
                p.otu.freq=drop(p.otu.freq),
                p.otu.tran=drop(p.otu.tran),
                p.otu.pa=drop(p.otu.pa),
                p.otu.omni=drop(p.otu.omni),
                p.otu.omni3=drop(p.otu.omni3),
                q.otu.freq=drop(q.otu.freq),
                q.otu.tran=drop(q.otu.tran),
                q.otu.pa=drop(q.otu.pa),
                q.otu.omni=drop(q.otu.omni),
                q.otu.omni3=drop(q.otu.omni3),
                
                p.otu.freq.OR=drop(p.otu.freq.OR),
                p.otu.tran.OR=drop(p.otu.tran.OR),
                p.otu.pa.OR=drop(p.otu.pa.OR),
                p.otu.omni.OR=drop(p.otu.omni.OR),
                p.otu.omni3.OR=drop(p.otu.omni3.OR), 
                q.otu.freq.OR=drop(q.otu.freq.OR),
                q.otu.tran.OR=drop(q.otu.tran.OR),
                q.otu.pa.OR=drop(q.otu.pa.OR),
                q.otu.omni.OR=drop(q.otu.omni.OR),
                q.otu.omni3.OR=drop(q.otu.omni3.OR), 
                
                p.otu.freq.com=drop(p.otu.freq.com),
                p.otu.tran.com=drop(p.otu.tran.com),
                p.otu.pa.com=drop(p.otu.pa.com),
                p.otu.omni.com=drop(p.otu.omni.com),
                p.otu.omni3.com=drop(p.otu.omni3.com), 
                q.otu.freq.com=drop(q.otu.freq.com),
                q.otu.tran.com=drop(q.otu.tran.com),
                q.otu.pa.com=drop(q.otu.pa.com),
                q.otu.omni.com=drop(q.otu.omni.com),
                q.otu.omni3.com=drop(q.otu.omni3.com), 
                
                
                p.global.freq=drop(p.global.freq), 
                p.global.tran=drop(p.global.tran), 
                p.global.pa=drop(p.global.pa),
                p.global.harmonic=drop(p.global.harmonic),
                p.global.fisher=drop(p.global.fisher),
                p.global.omni=drop(p.global.omni),
                p.global.omni3=drop(p.global.omni3), 
                
                p.global.freq.OR=drop(p.global.freq.OR), 
                p.global.tran.OR=drop(p.global.tran.OR), 
                p.global.pa.OR=drop(p.global.pa.OR),
                p.global.harmonic.OR=drop(p.global.harmonic.OR),
                p.global.fisher.OR=drop(p.global.fisher.OR),
                p.global.omni.OR=drop(p.global.omni.OR),
                p.global.omni3.OR=drop(p.global.omni3.OR), 
                
                p.global.freq.com=drop(p.global.freq.com), 
                p.global.tran.com=drop(p.global.tran.com), 
                p.global.pa.com=drop(p.global.pa.com),
                p.global.harmonic.com=drop(p.global.harmonic.com),
                p.global.fisher.com=drop(p.global.fisher.com),
                p.global.omni.com=drop(p.global.omni.com),
                p.global.omni3.com=drop(p.global.omni3.com), 
                
                
                detected.otu.freq = detected.otu.freq,
                detected.otu.tran = detected.otu.tran,
                detected.otu.pa = detected.otu.pa,
                detected.otu.omni = detected.otu.omni,
                detected.otu.omni3 = detected.otu.omni3,
                
                detected.otu.freq.OR = detected.otu.freq.OR,
                detected.otu.tran.OR = detected.otu.tran.OR,
                detected.otu.pa.OR = detected.otu.pa.OR,
                detected.otu.omni.OR = detected.otu.omni.OR,
                detected.otu.omni3.OR = detected.otu.omni3.OR,
                
                detected.otu.freq.com = detected.otu.freq.com,
                detected.otu.tran.com = detected.otu.tran.com,
                detected.otu.pa.com = detected.otu.pa.com,
                detected.otu.omni.com = detected.otu.omni.com,
                detected.otu.omni3.com = detected.otu.omni3.com,
                
                
                med.p.global.freq = med.p.global.freq, 
                med.p.global.tran = med.p.global.tran, 
                med.p.global.pa = med.p.global.pa,
                med.p.global.harmonic = med.p.global.harmonic,
                med.p.global.fisher = med.p.global.fisher,
                med.p.global.omni = med.p.global.omni,
                med.p.global.omni3 = med.p.global.omni3,
                
                med.p.global.freq.OR = med.p.global.freq.OR, 
                med.p.global.tran.OR = med.p.global.tran.OR, 
                med.p.global.pa.OR = med.p.global.pa.OR,
                med.p.global.harmonic.OR = med.p.global.harmonic.OR,
                med.p.global.fisher.OR = med.p.global.fisher.OR,
                med.p.global.omni.OR = med.p.global.omni.OR,
                med.p.global.omni3.OR = med.p.global.omni3.OR,
                
                med.p.global.freq.com = med.p.global.freq.com, 
                med.p.global.tran.com = med.p.global.tran.com, 
                med.p.global.pa.com = med.p.global.pa.com,
                med.p.global.harmonic.com = med.p.global.harmonic.com,
                med.p.global.fisher.com = med.p.global.fisher.com,
                med.p.global.omni.com = med.p.global.omni.com,
                med.p.global.omni3.com = med.p.global.omni3.com,
                
                
                med.detected.otu.freq = med.detected.otu.freq,
                med.detected.otu.tran = med.detected.otu.tran,
                med.detected.otu.pa = med.detected.otu.pa,
                med.detected.otu.omni = med.detected.otu.omni,
                med.detected.otu.omni3 = med.detected.otu.omni3,
                
                med.detected.otu.freq.OR = med.detected.otu.freq.OR,
                med.detected.otu.tran.OR = med.detected.otu.tran.OR,
                med.detected.otu.pa.OR = med.detected.otu.pa.OR,
                med.detected.otu.omni.OR = med.detected.otu.omni.OR,
                med.detected.otu.omni3.OR = med.detected.otu.omni3.OR,
                
                med.detected.otu.freq.com = med.detected.otu.freq.com,
                med.detected.otu.tran.com = med.detected.otu.tran.com,
                med.detected.otu.pa.com = med.detected.otu.pa.com,
                med.detected.otu.omni.com = med.detected.otu.omni.com,
                med.detected.otu.omni3.com = med.detected.otu.omni3.com,
                
                
                n.perm.completed=n.perm.completed,
                global.tests.stopped=global.tests.stopped,
                otu.tests.stopped=otu.tests.stopped,
                seed=seed)
    
    return(res)
    
} # End of ldm

#' @importFrom modeest mlv
calculate.x.and.resid = function( d.gower, y.freq, y.tran, index, m, adjust.for.confounders) {
    
    n.var = length(index)
    n.otu = ncol(y.freq)
    n.sam = nrow(d.gower)
    ndf.nominal = rep(0, n.var+1)
    
    tol.d = 10^-8
    
    #--------------------------------------------------------------------------
    # construct directions matrix x 
    # from each set of covariates in the list vars
    #--------------------------------------------------------------------------
    
    d.resid = d.gower
    
    for (i in 1:n.var) 
    {
        var = m[,1:index[i]]
        
        svd.var = svd(var)   
        use = (svd.var$d>tol.d)    
        
        hat.matrix = tcrossprod(svd.var$u[, use]) # svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        #---------------------
        # calculate direction
        #---------------------
        
        n.dim = dim( hat.matrix)[1]
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) # es: eigen system in Mathematica
        
        use = ( abs(es.model$values)>tol.d )
        
        ndf.model = sum( use )
        
        x.model = es.model$vectors[, use]
        e.model = es.model$values[use]
        
        hat.matrix.bar = diag(n.dim) - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # end of calculating direction
        #-----------------------------    
        
        if (i==1) {
            x = x.model
            e = e.model
        } else {   
            x = cbind(x, x.model)
            e = c(e, e.model )
        }
        
        ndf.nominal[i] = ndf.model
        
    }
    
    if (!is.null(dim(x))) {
        if (dim(x)[2] > dim(x)[1]) {
            stop("Problems occurred in finding the design matrix with orthogonal columns! More columns than rows!")
        }
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE)
    use = which( abs(es.resid$values)>tol.d )
    
    ndf.nominal[n.var+1] = length(use)
    x = cbind(x, es.resid$vectors[, use])
    e = c(e, es.resid$values[use])
    
    
    #---------------------------------------------------------
    # fit LDM to x
    #---------------------------------------------------------
    
    x = matrix(x, nrow=nrow(d.gower)) 
    
    wt.freq = crossprod(x, y.freq) # t(x) %*% y.freq
    
    d.freq = sqrt(rowSums(wt.freq^2))
    v.freq = t((1/d.freq)*wt.freq)
    d.tran = NULL
    v.tran = NULL
    if (!is.null(y.tran)) {
        wt.tran = crossprod(x, y.tran) # t(x) %*% y.tran    
        d.tran = sqrt(rowSums(wt.tran^2))
        v.tran = t((1/d.tran)*wt.tran)
    }
    
    #-------------------------------------------------
    # low, up
    #-------------------------------------------------
    
    low = rep(NA, n.var)
    up = rep(NA, n.var)
    
    up.prev = 0
    
    for (k in 1:n.var)
    {
        low[k] = up.prev + 1
        up[k] = up.prev + ndf.nominal[k]
        up.prev = up[k]
    }
    
    #-------------------------------------------------
    # calculate resid, ss.tot
    #-------------------------------------------------
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    ss.tot.freq = matrix( rep(NA, n.otu*n.var1), nrow=n.var1)
    resid.freq = array( NA, dim=c( dim(y.freq), n.var1 ) ) 
    ss.tot.tran = NULL
    resid.tran = NULL
    if (!is.null(y.tran)) {
        ss.tot.tran = matrix( rep(NA, n.otu*n.var1), nrow=n.var1)
        resid.tran = array( NA, dim=c( dim(y.tran), n.var1 ) ) 
    }
    
    for (k in 1:n.var1) {
        
        k1 = ifelse(adjust.for.confounders, k+1, k)
        use = setdiff( 1:up[n.var], low[k1]:up[k1] )
        
        resid.freq[,,k] = y.freq - x[,use,drop=FALSE] %*% wt.freq[use,,drop=FALSE]
        ss.tot.freq[k,] = colSums( resid.freq[,,k,drop=FALSE]^2 )
        if (!is.null(y.tran)) {
            resid.tran[,,k] = y.tran - x[,use,drop=FALSE] %*% wt.tran[use,,drop=FALSE]
            ss.tot.tran[k,] = colSums( resid.tran[,,k,drop=FALSE]^2 )
        } 
    }
    
    res = list( x=x,
                d.freq=d.freq,
                v.freq=v.freq,
                d.tran=d.tran,
                v.tran=v.tran,
                low=low,
                up=up,
                resid.freq=resid.freq,
                resid.tran=resid.tran,
                ss.tot.freq=ss.tot.freq,
                ss.tot.tran=ss.tot.tran,
                ndf=ndf.nominal)
    
    return(res)
    
} # calculate.x.and.resid


ldm.stat = function(x, low, up, resid, ss.tot, adjust.for.confounders, comp.anal=FALSE, comp.anal.adjust="median", comp.effect=NULL) {
    
    #---------------------------------------------
    #  calculate FL statistics for each model
    #---------------------------------------------
    
    n.var = length(low)
    n.otu = length(resid[1,,1,1])
    n.rarefy = length(resid[1,1,1,])
    n.sam = length(resid[,1,1,1])
    sqrt.n.sam = sqrt(n.sam)
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    ve.global = matrix(0, n.var1, n.rarefy)
    ve.otu = matrix(0, n.var1, n.otu)

    comp.effect.tmp = NULL
    if (comp.anal) {
        if (is.null(comp.effect)) {
            comp.effect.tmp = rep(0, up[n.var])
        } else {
            comp.effect.tmp = comp.effect
        }
    }
    
    for (r in 1:n.rarefy) {
        for (k in 1:n.var1) {
            
            k1 = k + as.numeric(adjust.for.confounders)
            use = low[k1]:up[k1]
            
            wt = crossprod(x[, use], resid[,,k,r]) # t( x[, use] ) %*% resid[,,k,r]
            
            if (comp.anal) {
                if (is.null(comp.effect)) {
                    if (comp.anal.adjust=="median") {
                        comp.effect.tmp[use] <- rowMedians(wt)
                    } else if (comp.anal.adjust=="meanshift") {
                        comp.effect.tmp[use] <- apply(wt, 1, function(x) modeest::mlv(sqrt.n.sam*x, method = "meanshift", kernel = "gaussian")[1]/sqrt.n.sam)
                    }
                }
                ve.otu.k = colSums( (wt - comp.effect.tmp[use])^2 )
                sigma.k = 1
                
            } else {
            
                ve.otu.k = colSums( wt^2 )
                
                if (n.var==1) {
                    sigma.k = ss.tot[k,,r] - ve.otu.k
                } else {
                    use = 1:up[n.var]
                    wt.cum = crossprod(x[, use], resid[,,k,r]) # t( x[, use] ) %*% resid[,,k,r]
                    sigma.k = ss.tot[k,,r] - colSums( wt.cum^2 )
                }
            }
            
            sigma.tmp = ifelse(sigma.k>1e-16, sigma.k, 1)
            ve.otu[k,] = ve.otu[k,] + ve.otu.k/sigma.tmp
            ve.global[k, r:n.rarefy] = ve.global[k, r:n.rarefy] + sum( ve.otu.k )/sum( sigma.k )
            
        }
    }
    
    out = list( comp.effect = comp.effect.tmp,
                ve.otu=ve.otu, 
                ve.global=ve.global)   
    
    return(out)
    
} # ldm.stat




log.int.seq = function(low, up, log.int) {
    return(ifelse(low==up, 0, sum(log.int[(low+1):up])))
}

calculate.x.and.resid.allrarefy = function( y, index, m, adjust.for.confounders) {
    
    n.var = length(index)
    n.otu = ncol(y)
    n.sam = nrow(y)
    ndf.nominal = rep(0, n.var+1)
    
    tol.d = 10^-8
    
    #--------------------------------------------------------------------------
    # construct directions matrix x 
    # from each set of covariates in the list vars
    #--------------------------------------------------------------------------
    
    # x = gramSchmidt(m)$Q
    
    d.resid = diag(n.sam)
    
    for (i in 1:n.var) 
    {
        var = m[,1:index[i]]
        
        svd.var = svd(var)   
        use = (svd.var$d > tol.d)    
        
        hat.matrix = tcrossprod(svd.var$u[, use]) # svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        #---------------------
        # calculate direction
        #---------------------
        
        n.dim = dim( hat.matrix)[1]
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) # es: eigen system in Mathematica
        
        use = ( abs(es.model$values)>tol.d)
        
        ndf.model = sum( use )
        
        x.model = es.model$vectors[, use]
        e.model = es.model$values[use]
        
        hat.matrix.bar = diag(n.dim)  - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # end of calculating direction
        #-----------------------------    
        
        if (i==1) {
            x = x.model
            e = e.model
        } else {   
            x = cbind(x, x.model)
            e = c(e, e.model )
        }
        
        ndf.nominal[i] = ndf.model
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE)
    use = which( abs(es.resid$values) > tol.d )
    
    ndf.nominal[n.var+1] = length(use)
    x = cbind(x, es.resid$vectors[, use])
    e = c(e, es.resid$values[use])
    
    #-------------------------------------------------
    # low, up
    #-------------------------------------------------
    
    low = rep(NA, n.var)
    up = rep(NA, n.var)
    
    up.prev = 0
    
    for (k in 1:n.var)
    {
        low[k] = up.prev + 1
        up[k] = up.prev + ndf.nominal[k]
        up.prev = up[k]
    }
    
    #-------------------------------------------------
    # calculate phi
    #-------------------------------------------------
    
    lib.size = rowSums(y)
    rarefy.depth = min(rowSums(y))
    log.int = log(1:max(lib.size)) 
    vec.y = as.vector(y)
    vec.tmp = rep(0, length(vec.y))
    phi = rep(NA, length(vec.y))
    
    w = which(lib.size-rarefy.depth-vec.y >= 0)
    phi[-w] = 0
    
    a = mapply(log.int.seq, (lib.size-rarefy.depth-vec.y)[w], (lib.size-rarefy.depth-vec.tmp)[w], MoreArgs=list(log.int=log.int))
    b = mapply(log.int.seq, (lib.size-vec.y)[w], (lib.size-vec.tmp)[w], MoreArgs=list(log.int=log.int))
    phi[w] = exp(a - b)
    phi = matrix(phi, nrow=n.sam)
    phi_1phi = phi*(1-phi)
    
    #-------------------------------------------------
    # calculate resid, ss.tot
    #-------------------------------------------------
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    wt = crossprod(x, phi) # t(x) %*% phi  
    
    resid = array( NA, dim=c( n.sam, n.otu, n.var1 ) ) 
    ss.tot = matrix( NA, nrow=n.var1, ncol=n.otu)
    
    P.resid = array( NA, dim=c(n.sam, n.sam, n.var1) ) 
    ss.tot.1 = matrix(NA, nrow=n.var1, ncol=n.otu)
    
    for (k in 1:n.var1) {
        
        k1 = ifelse(adjust.for.confounders, k+1, k)
        
        use = setdiff( 1:up[n.var], low[k1]:up[k1] )
        
        resid[,,k] = phi - x[,use,drop=FALSE] %*% wt[use,,drop=FALSE]
        ss.tot[k,] = colSums( resid[,,k]^2 )
        
        P.resid[,,k] = diag(n.sam) - tcrossprod(x[,use,drop=FALSE]) # x[,use,drop=FALSE] %*% t(x[,use,drop=FALSE])
        D_k = rowSums(P.resid[,,k]^2)
        ss.tot.1[k,] = colSums(D_k*phi_1phi)
    }
    
    res = list( x=x,
                low=low,
                up=up,
                resid=resid,
                ss.tot=ss.tot,
                P.resid=P.resid,
                ss.tot.1=ss.tot.1,
                phi = phi,
                phi_1phi=phi_1phi)
    
    return(res)
    
} # calculate.x.and.resid.allrarefy


ldm.stat.allrarefy = function(x, low, up, resid, ss.tot, P.resid, ss.tot.1, phi_1phi, adjust.for.confounders) {
    
    #---------------------------------------------
    #  calculate FL statistics for each model
    #---------------------------------------------
    
    n.var = length(low)
    n.sam = dim(resid)[1]
    n.otu = dim(resid)[2]
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    ve.otu = matrix(rep(NA, n.otu*n.var1), nrow=n.var1 )
    ve.global = rep(NA, n.var1)
    
    for (k in 1:n.var1) {
        
        k1 = k + as.numeric(adjust.for.confounders)
        
        use = low[k1]:up[k1]
        
        wt = crossprod(x[, use,drop=FALSE], resid[,,k]) # t( x[, use,drop=FALSE] ) %*% resid[,,k]
        ve.otu.k = colSums( wt^2 )
        
        P.nume = P.resid[,,k] %*% x[, use,drop=FALSE]
        D_k = rowSums(P.nume^2)
        ve.otu.k.1 = colSums(D_k*phi_1phi)
        
        
        use = 1:up[n.var]
        
        wt.cum = crossprod(x[, use,drop=FALSE], resid[,,k]) # t( x[, use,drop=FALSE] ) %*% resid[,,k]
        sigma.k = ss.tot[k,] - colSums( wt.cum^2 )
        
        P.denom = P.resid[,,k] %*% x[, use,drop=FALSE]
        D_k = rowSums(P.denom^2)
        sigma.k.1 = ss.tot.1[k,] - colSums(D_k*phi_1phi)
        
        ve.otu[k,] = ( ve.otu.k + ve.otu.k.1)/( sigma.k + sigma.k.1 )
        ve.global[k] = sum( ve.otu.k + ve.otu.k.1 )/sum( sigma.k + sigma.k.1 )
        
        ve.otu[k, which(is.na(ve.otu[k,]))] = 0  #Debugged: in case the denominator (sigma.k + sigma.k.1) is 0
    }
 
    out = list( ve.otu=ve.otu, 
                ve.global=ve.global)   
    
    return(out)
    
} # ldm.stat.allrarefy



#' PERMANOVA test of association based on the Freedman-Lane permutation scheme
#' 
#' This function performs the PERMANOVA test that can allow adjustment of
#' confounders and control of clustered data. It can also be used for testing 
#' presence-absence associations based on infinite number of rarefaction replicates. 
#' As in \code{ldm}, 
#' \code{permanovaFL} allows multiple sets of covariates to be tested, 
#' in the way that the sets are entered sequentially and the variance 
#' explained by each set is that part that remains after the previous 
#' sets have been fit. It allows testing of a survival outcome, by using the Martingale or deviance residual (from fitting a Cox model to the survival outcome and other covariates) as a covariate in the regression. 
#' It allows multiple distance matrices and provides an omnibus test in such cases. 
#' It also allows testing of the mediation effect of the microbiome in the pathway between the exposure(s) and the outcome(s), 
#' where the exposure(s) and outcomes(s) are specified as the first and second (sets of) covariates. 
#' 
#' @param formula a symbolic description of the model to be fitted in the form
#'   of \code{data.matrix ~ sets of covariates} or \code{data.matrix |
#'   confounders ~ sets of covariates}. The details of model specification are
#'   given in "Details" of \code{ldm}. Additionally, in \code{permanovaFL}, the \code{data.matrix}
#'   can be either an OTU table or a distance matrix. If it is an OTU table,
#'   the distance matrix will be calculated internally using the OTU table, \code{tree} (if required), and 
#'   \code{dist.method}. If \code{data.matrix} is a distance
#'   matrix (having class \code{dist} or \code{matrix}), it can be squared and//or centered by
#'   specifying \code{square.dist} and \code{center.dist} (described below).  Distance matrices are distinguished
#'   from OTU tables by checking for symmetry of \code{as.matrix(data.matrix)}.
#' @param other.surv.resid a vector of data, usually the Martingale or deviance residuals from fitting the Cox model to the survival outcome (if it is the outcome of interest) and other covariates.
#' @param data an optional data frame, list or environment (or object coercible 
#' to a dataframe) containing the covariates of interest and confounding covariates. 
#' If not found in \code{data}, the covariates are taken from environment(formula), 
#' typically the environment from which \code{permanovaFL} is called. The default is .GlobalEnv.
#' @param tree a phylogenetic tree. Only used for calculating a 
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of 
#'   the requested distance does not involve a phylogenetic tree, or if a 
#'   distance matrix is directly imported through \code{formula}.
#' @param dist.method a vector of methods for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package
#'  (i.e., "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", 
#'  "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis")
#'   as well as "hellinger" and "wt-unifrac". 
#'   Not used if a distance matrix is specified in \code{formula} or \code{dist.list}. 
#'   The default is c("bray"). 
#'   For more details, see the \code{dist.method} argument in the \code{ldm} function.
#' @param dist.list a list of pre-calculated distance matrices. 
#' @param test.mediation a logical value indicating whether to perform the mediation analysis. The default is FALSE. 
#' If TRUE, the formula takes the specific form \code{otu.table ~ exposure + outcome} or most generally
#' \code{otu.table or distance matrix | (set of confounders) ~ (set of exposures) + (set of outcomes)}.
#' @param n.cores The number of cores to use in parallel computing, i.e., at most how many child processes will be run simultaneously. 
#' The default is 4.
#' @param cluster.id cluster identifiers. The default is value of NULL should be used if the observations are 
#' not in clusters (i.e., independent).
#' @param strata a factor variable (or, character variable converted into a factor) to define strata (groups), within which to constrain permutations. 
#'   The default is NULL.
#' @param how a permutation control list, for users who want to specify their permutation control list using the \code{how} function 
#'   from the \code{permute} R package.  The default is NULL.
#' @param perm.within.type a character string that takes values "free", "none", "series", or "grid".  
#'   The default is "free" (for random permutations).
#' @param perm.between.type a character string that takes values "free", "none", or "series".  
#'   The default is "none".
#' @param perm.within.nrow a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See the documentation for the R package \code{permute} for further details.
#' @param perm.within.ncol a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See the documentation for the R package \code{permute} for further details.
#' @param n.perm.max the maximum number of permutations.
#'   The default is 5000.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. 
#'   The default is 100.
#' @param seed a user-supplied integer seed for the random number generator in the 
#'   permutation procedure. The default is NULL; with the default value, an integer seed will be 
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case 
#'   the user wants to reproduce the permutation replicates.
#' @param square.dist a logical variable indicating whether to square the
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a vector of logical variables indicating whether to scale the OTU table in calculating the distance matrices in \code{dist.method}.
#'   For count data, this corresponds to dividing by the library size to give
#'   relative abundances. The default is TRUE.
#' @param binary a vector of logical values indicating whether to base the calculation of the distance matrices in \code{dist.method} on presence-absence (binary) data. The default is c(FALSE) (analyzing relative abundance data).
#' @param n.rarefy number of rarefactions. The default is 0 (no rarefaction).
#' @return  a list consisting of 
#' \item{F.statistics}{F statistics for testing each set of covariates}
#' \item{R.squared}{R-squared statistic for each set of covariates}
#' \item{F.statistics.OR, R.squared.OR}{F statistics and R-squared statistic when the last covariate is \code{other.surv.resid}}
#' \item{p.permanova}{p-values for testing each set of covariates} 
#'   \item{p.permanova.omni}{the omnibus p-values (that combines information from multiple distance matrices) for testing each set of covariates} 
#'   \item{med.p.permanova}{p-values for testing mediation}
#'   \item{med.p.permanova.omni}{the omnibus p-values for testing mediation}
#'   \item{p.permanova.OR, p.permanova.omni.OR}{when using \code{other.surv.resid} as the last covariate}
#'   \item{med.p.permanova.OR, med.p.permanova.omni.OR}{when using \code{other.surv.resid} as the outcome in the mediation analysis}
#'   \item{p.permanova.com, p.permanova.omni.com}{the combination test that combines the results from analyzing the Martingale residual and the Deviance residual (one specified in the formula and one specified in \code{other.surv.resid})}
#'   \item{med.p.permanova.com, med.p.permanova.omni.com}{the combination test for the mediation effect}
#' \item{n.perm.completed}{number of permutations completed}
#' \item{permanova.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all tests of covariates}
#' \item{seed}{the seed that is user supplied or internally generated, stored in case 
#'   the user wants to reproduce the permutation replicates}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gsatten@emory.edu>
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom parallel mclapply
#' @importFrom permute how shuffleSet Plots Within
#' @importFrom utils tail
#' @import matrixStats
#' @export
#' @references Hu YJ, Satten GA (2020). Testing hypotheses about the microbiome using the linear decomposition model (LDM) 
#'   Bioinformatics, 36(14), 4106-4115.
#' @references Hu YJ and Satten GA (2021). A rarefaction-without-resampling extension of PERMANOVA for testing presence-absence associations in the microbiome. bioRxiv, https://doi.org/10.1101/2021.04.06.438671.
#' @references Zhu Z, Satten GA, Caroline M, and Hu YJ (2020). Analyzing matched sets of microbiome data using the LDM and PERMANOVA. Microbiome, 9(133), https://doi.org/10.1186/s40168-021-01034-9.
#' @references Hu Y, Li Y, Satten GA, and Hu YJ (2022) Testing microbiome associations with censored survival outcomes at both the community and individual taxon levels. bioRxiv, doi.org/10.1101/2022.03.11.483858.
#' @examples
#'data(throat.otu.tab5)
#'data(throat.meta)
#'res.perm <- permanovaFL(throat.otu.tab5 | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
#'                        data=throat.meta, dist.method="bray", seed=82955)

permanovaFL = function(formula, other.surv.resid=NULL, data=.GlobalEnv, tree=NULL, dist.method=c("bray"), dist.list=NULL, 
                           cluster.id=NULL, strata=NULL, how=NULL,
                           perm.within.type="free", perm.between.type="none",
                           perm.within.ncol=0, perm.within.nrow=0,
                           n.perm.max=5000, n.rej.stop=100, seed=NULL,
                           square.dist=TRUE, center.dist=TRUE, scale.otu.table=c(TRUE), 
                           binary=c(FALSE), n.rarefy=0,
                           test.mediation=FALSE,
                           n.cores=4) {  
    
    #------------------------
    # form.call
    #------------------------
    options(na.action=na.omit) # fixed a bug here
    object=formula
    #
    #   extract cluster.id from dataframe
    #
    cl=match.call()
    mf=match.call(expand.dots=FALSE)
    m=match( x='cluster.id', table=names(mf) )
    mf.string=as.character( mf[c(1L,m)] )
    cluster.name=mf.string[2]
    if (cluster.name=='NULL') {
        cluster.id=NULL
    } else {   
        loc.dollar=utils::tail( gregexpr('\\$', cluster.name)[[1]] , n=1 )
        if (loc.dollar<0)  {
            cluster.id=getElement(data,cluster.name)
            if( is.null(cluster.id) ) cluster.id=get(cluster.name)
        } else {   
            df.name=get( substr(cluster.name, start=1, stop=loc.dollar-1) )
            var.name=substr(cluster.name, start=loc.dollar+1, stop=nchar(cluster.name))            
            cluster.id= getElement(df.name,var.name) 
        }
    }
    #        
    #   extract model from formula    
    #    
    obj=toString(object)
    obj=gsub('\\s','',obj)
    prefix=' ~ + 0 + '
    loc.comma=gregexpr(',',obj)[[1]]
    start.terms=loc.comma[2]
    terms=substr(obj,start=start.terms+1, stop=nchar(obj))
    #
    #   find n.obs and full set of rownames
    #   
    if (class(data)=='data.frame') {
        row.names=rownames(data)
        n.obs=length(row.names)
    } else {   
        df=model.frame( as.formula(paste('~',terms)) , na.action=na.pass )
        row.names=rownames(df)
        n.obs=length(row.names)
    }
    #
    #   check for missing values in cluster.id
    #        
    
    if (is.null(cluster.id)) {
        use.rows=row.names
    } else {   
        use=!is.na(cluster.id)
        use.rows=row.names[use]
    }
    #
    #   check for and extract confounders
    #
    model=list()
    j=1
    loc.bar=regexpr('\\|',obj)[1]
    loc.minus=regexpr('-',obj)[1]
    loc.delim=max( loc.bar, loc.minus)
    if (loc.delim>0) {
        end.confound=loc.comma[2]
        c=substr(obj,start=loc.delim+1, stop=end.confound-1)
        conf=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        model[[j]]=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        #       use.rows=intersect( use.rows, rownames(conf) )
        use.rows=rownames(model[[1]]) 
        j=j+1
    } else {
        conf=NULL
    }     
    #
    #   extract model terms
    #
    #   j=1
    continue=TRUE
    while (continue) {
        if (substr(terms,1,1)=='(') {
            stop=regexpr(')\\+',terms)[1]
        } else {
            stop=regexpr('\\+',terms)[1] - 1
        }          
        
        if (stop<=0) stop=nchar(terms) 
        m=substr(terms, start=1, stop=stop)
        model[[j]]=model.matrix( as.formula( paste(prefix,m) ) , data=data)
        use.rows=intersect( use.rows, rownames(model[[j]]) )
        #        if (j==1) {
        #            use.rows=rownames(model[[1]])
        #            }
        #        else {
        #            use.rows=intersect( use.rows, rownames(model[[j]]) )
        #            }         
        if (stop+2<=nchar(terms)) {
            terms=substr(terms, start=stop+2, stop=nchar(terms))
            j=j+1
        } else {
            continue=FALSE
        }             
    }   
    n.model=j    
    #
    #  extract OTU table
    #      
    if (is.null(conf)) loc.delim=loc.comma[2]
    otu.name=substr(obj, start=loc.comma[1]+1, stop=loc.delim-1)
    #   loc.dollar=regexpr('\\$', otu.name)[1]
    loc.dollar=utils::tail( gregexpr('\\$', otu.name)[[1]] , n=1 )
    if (loc.dollar<0)  {
        if (class(data)=='data.frame') {
            otu.table=getElement(data, otu.name)
            if (is.null(otu.table)) otu.table= get(otu.name) 
            otu.table=as.matrix(otu.table)
        } else {
            otu.table=as.matrix( get(otu.name) )
        }
    } else {
        df.name=get( substr(otu.name, start=1, stop=loc.dollar-1) )
        var.name=substr(otu.name, start=loc.dollar+1, stop=nchar(otu.name))
        otu.table=as.matrix( getElement(df.name,var.name) )
    }        
    
    #    if (is.null(otu.table)) otu.table=as.matrix( getElement(.GlobalEnv,otu.name) )
    if ( nrow(otu.table) != n.obs ) {
        if (ncol(otu.table)==n.obs ) {
            otu.table=t(otu.table)
        } else {   
            stop('OTU table and covariates have different number of observations')
        }
    }   
    
    
    otu.or.dist <- as.matrix(otu.table)
    
    otu.table <- NULL
    dist <- NULL
    if (isSymmetric(otu.or.dist)) {
        dist <- otu.or.dist
    } else {
        otu.table <- otu.or.dist
    }
    
    
    #---------------------------------------
    # checking negative values in otu.table
    #---------------------------------------
    
    if (!isSymmetric(otu.or.dist)) {
        neg.exist = any(otu.table<0)
        if (neg.exist) {
            if (scale.otu.table == TRUE) {
                stop("The OTU table has negative values, so it does not make sense to use 'scale.otu.table=TRUE'")
            }
        }
    }
    
    #
    #   remove rows having NA 
    #    
    for (j in 1:n.model) {
        keep =  rownames( model[[j]] ) %in% use.rows
        model[[j]]=model[[j]][keep,,drop=FALSE]
    }
    if (!is.null(conf)) {
        keep =  rownames(conf) %in% use.rows 
        conf=conf[keep,,drop=FALSE]
    }
    keep=row.names %in% use.rows    
    if (!is.null(cluster.id)) cluster.id=cluster.id[keep]
    
    if (!is.null(dist)) {
        dist = dist[keep, keep]
        if (dim(model[[1]])[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between covariates and the distance matrix' )
    }
    if (!is.null(otu.table)) {
        otu.table = otu.table[keep,,drop=FALSE]
        if (dim(model[[1]])[1] != dim(otu.table)[1]) 
            otu.table <- t(otu.table)
        if (dim(model[[1]])[1] != dim(otu.table)[1]) stop( 'numbers of observations mismatch between covariates and the OTU table' )
    }
    
    #------------------------
    # setup permutation
    #------------------------
    
    if (class(how)=='how') {
        CTRL=how                   # user-provided how list
    }
    else {
        if (is.null(cluster.id)) {
            if (is.null(perm.within.type) & is.null(perm.between.type)) {
                # default when no unclustered data has no type specified is 'free'
                perm.within.type='free'    
            }
            if (is.null(strata)) {
                # setup for unclustered permutation
                CTRL = permute::how( within=permute::Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))  
            }
            else {
                # setup for unclustered, stratified permutation
                strata=as.factor(strata)
                CTRL = permute::how( blocks=strata, within=permute::Within(type=perm.within.type, 
                                                         nrow=perm.within.nrow, 
                                                         ncol=perm.within.ncol))  
            }    
        }
        else {        
            cluster.id=as.factor(cluster.id)
            if (is.null(strata)) {            
                #  clustered but unstratified data
                CTRL = permute::how( plots=permute::Plots(cluster.id, type=perm.between.type ), 
                            within=permute::Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            }
            else {
                #   clustered and stratified data
                strata=as.factor(strata)             
                CTRL = permute::how( blocks=strata, 
                            plots=permute::Plots(cluster.id, type=perm.between.type ), 
                            within=permute::Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            }
        }
    }  
    
    #------------------------
    # setup model
    #------------------------
    
    OR = other.surv.resid # e.g., Deviance residual
    
    adjust.for.confounders = !is.null(conf)
    
    n.var = length(model)
    n.var1 = n.var - as.numeric(adjust.for.confounders)
    n.obs = dim(model[[1]])[1]
    
    center.vars=TRUE
    
    index = rep(0, n.var)
    
    for (i in 1:n.var) {
        m.i = model[[i]]
        if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
        
        if (i==1) {
            m = m.i
            index[i] = dim(m.i)[2] 
        } else {
            m = cbind(m, m.i)   
            index[i] = index[i-1] + dim(m.i)[2]    
        }
        
    }    
    
    if (!is.null(OR)) {
        m.OR = m
        if (center.vars) OR = scale( OR, center=TRUE, scale=FALSE )
        m.OR[,ncol(m.OR)] = OR
    }
    
    #------------------------
    # deciding methods
    #------------------------
    
    if (!is.null(dist.list)) {
        n.dist = length(dist.list)
    } else if (is.null(dist.list)) {
        if (!is.null(dist)) {
            dist.list=list(dist)
            n.dist = 1
        } else {
            dist.list=list(NULL)
            n.dist = length(dist.method)
            if (n.dist!=length(binary)) {
                binary=rep("FALSE", n.dist)
                binary[which(dist.method=="jaccard")] = TRUE
            }
            scale.otu.table=rep(scale.otu.table, n.dist)
            scale.otu.table[which(binary=="TRUE")] = FALSE
            
        }
    }
    
    no_rarefy = (n.rarefy==0)
    if (no_rarefy) n.rarefy=1
    if (!is.null(dist.list[[1]])) n.rarefy=1
    
    #---------------------
    # rarefaction or not?
    #---------------------
    
    resid.dist = array(NA, dim=c(n.dist, n.obs, n.obs, n.var1, n.rarefy))
    resid.dist.OR = array(NA, dim=c(n.dist, n.obs, n.obs, n.var1, n.rarefy))
    
    permanova.obs = array(NA, dim=c(n.dist, n.var1, n.rarefy))
    permanova.perm = array(NA, dim=c(n.dist, n.var1, n.rarefy, n.perm.max))
    permanova.obs.OR = NULL
    permanova.perm.OR = NULL
    if (!is.null(OR)) {
        permanova.obs.OR = array(NA, dim=c(n.dist, n.var1, n.rarefy))
        permanova.perm.OR = array(NA, dim=c(n.dist, n.var1, n.rarefy, n.perm.max))
    }
    
    med.permanova.obs = NULL 
    med.permanova.perm = NULL 
    med.permanova.obs.OR = NULL 
    med.permanova.perm.OR = NULL 
    if (test.mediation) {
        med.permanova.obs = array(NA, dim=c(n.dist, n.rarefy))
        med.permanova.perm = array(NA, dim=c(n.dist, n.rarefy, n.perm.max))
        if (!is.null(OR)) {
            med.permanova.obs.OR = array(NA, dim=c(n.dist, n.rarefy))
            med.permanova.perm.OR = array(NA, dim=c(n.dist, n.rarefy, n.perm.max))
        }
    }
    
    if (is.null(seed)) {
        seed = sample(1:10^6, 1)
    }
    set.seed(seed)
    
    for (r in 1:n.rarefy) {
        
        if (!is.null(otu.table) & is.null(dist.list[[1]])) {
            if (!no_rarefy) {
                otu.rarefy= Rarefy(otu.table)$otu.tab.rff
            } else {
                otu.rarefy = otu.table
            }
        }
        
        for (d in 1:n.dist) {
            
            #------------------------
            # dist matrix
            #------------------------
            
            if (is.null(dist.list[[1]])) {
                if (binary[d]) {
                    otu.rarefy.d = (otu.rarefy>0)*1
                } else {
                    otu.rarefy.d = otu.rarefy
                }
                dist_r_d <- calculate.dist(dist.method=dist.method[d], otu.table=otu.rarefy.d, tree=tree, scale.otu.table=scale.otu.table[d], binary=binary[d])
                d.gower <- gower(d=dist_r_d, square=square.dist, center=center.dist)
            } else {
                d.gower <- gower(d=dist.list[[d]], square=square.dist, center=center.dist)
            }
            
            #---------------------
            # model fitting
            #---------------------
            
            fit.res = fit.permanova( d.gower=d.gower, index=index, m=m, adjust.for.confounders=adjust.for.confounders) 
            resid.dist[d,,,,r] = fit.res$resid.dist
            if (!is.null(OR)) {
                fit.res.OR = fit.permanova( d.gower=d.gower, index=index, m=m.OR, adjust.for.confounders=adjust.for.confounders) 
                resid.dist.OR[d,,,,r] = fit.res.OR$resid.dist
            }
        } # d
        
        if (r==1) {
            x.design = fit.res$x
            if (!is.null(OR)) x.design.OR = fit.res.OR$x
            low = fit.res$low
            up = fit.res$up
            ndf = fit.res$ndf
        }
        
    }# rarefaction
    
    R.squared = array(NA, dim=c(n.dist, n.var1, n.rarefy))
    for (d in 1:n.dist) {
        tmp = permanova.stat(x=x.design, low=low, up=up, resid.dist=resid.dist[d,,,,,drop=FALSE], ndf=ndf, adjust.for.confounders=adjust.for.confounders)
        permanova.obs[d,,] = tmp$permanova
        R.squared[d,,] = tmp$R.squared
    }
    if (test.mediation) {
        med.permanova.obs <- array(permanova.obs[,1,]*permanova.obs[,2,], dim=c(n.dist, n.rarefy))
    }
    
    R.squared.OR = NULL
    if (!is.null(OR)) {
        R.squared.OR = array(NA, dim=c(n.dist, n.var1, n.rarefy))
        for (d in 1:n.dist) {
            tmp.OR = permanova.stat(x=x.design.OR, low=low, up=up, resid.dist=resid.dist.OR[d,,,,,drop=FALSE], ndf=ndf, adjust.for.confounders=adjust.for.confounders)
            permanova.obs.OR[d,,] = tmp.OR$permanova
            R.squared.OR[d,,] = tmp.OR$R.squared
        }
        if (test.mediation) {
            med.permanova.obs.OR <- array(permanova.obs.OR[,1,]*permanova.obs.OR[,2,], dim=c(n.dist, n.rarefy))
        }
    }
    
    p.permanova <- NULL
    p.permanova.OR <- NULL
    p.permanova.com <- NULL
    p.permanova.omni <- NULL
    p.permanova.omni.OR <- NULL
    p.permanova.omni.com <- NULL
    
    med.p.permanova <- NULL
    med.p.permanova.OR <- NULL
    med.p.permanova.com <- NULL
    med.p.permanova.omni <- NULL
    med.p.permanova.omni.OR <- NULL
    med.p.permanova.omni.com <- NULL
    
    n.perm.completed = NULL
    permanova.stopped = FALSE
    
    #---------------------
    # permutation
    #---------------------
    
    if (n.perm.max > 0) {
        
        tol.eq = 10^-8
        n.perm.block = 100
        
        ############################################################################
        
        parallel.perm <- function(i, x, perm, low, up, resid.dist, ndf, adjust.for.confounders) {
            F.stat = array(NA, dim=dim(resid.dist)[c(1,4,5)])
            x.perm = x[perm[i,], ]   
            for (d in 1:dim(F.stat)[1]) {
                F.stat[d,,]=permanova.stat(x=x.perm, low=low, up=up, resid.dist=resid.dist[d,,,,,drop=FALSE], ndf=ndf, adjust.for.confounders=adjust.for.confounders)$permanova
            }
            F.stat
        }
        
        med.permanova.fun <- function(x.perm, x.obs){
            pmax(x.perm[,1,]*x.obs[,2,], x.perm[,2,]*x.obs[,1,], x.perm[,1,]*x.perm[,2,])
        }
        
        ############################################################################
        
        n.perm.completed = 0
        meet.all.rej.stop = TRUE
        med.meet.all.rej.stop = TRUE
        
        n.permanova = array(0, dim=c(n.dist, n.var1, n.rarefy))
        n.permanova.OR = NULL
        if (!is.null(OR)) n.permanova.OR = array(0, dim=c(n.dist, n.var1, n.rarefy))
        
        if (test.mediation) {
            med.n.permanova = array(0, dim=c(n.dist, n.rarefy))
            if (!is.null(OR)) med.n.permanova.OR = array(0, dim=c(n.dist, n.rarefy))
        }
        
        set.seed(seed)
        
        nblock = ceiling(n.perm.max/n.perm.block)
        
        for (i.block in 1:nblock) {
            
            perm = permute::shuffleSet(n.obs, n.perm.block, CTRL)
            
            if (n.cores > 1) { # parallel computing
                
                i.perm = (n.perm.completed+1):(n.perm.completed+n.perm.block)
                
                if (Sys.info()[['sysname']] == 'Windows') {
                    parallel.stat = BiocParallel::bplapply(1:n.perm.block, parallel.perm, BPPARAM = BiocParallel::MulticoreParam(workers=n.cores), x.design, perm, low, up, resid.dist, ndf, adjust.for.confounders)
                } else {
                    parallel.stat = parallel::mclapply(1:n.perm.block, parallel.perm, mc.cores = n.cores, x.design, perm, low, up, resid.dist, ndf, adjust.for.confounders)
                }
                permanova.perm[,,,i.perm] <- array(unlist(parallel.stat), dim=c(dim(permanova.obs), n.perm.block))
                n.rejection <- sapply(parallel.stat, function(x) (x > permanova.obs + tol.eq) + (x > permanova.obs - tol.eq), simplify="array") 
                if (n.dist>1 | n.var1>1 | n.rarefy>1) {
                    n.permanova <- n.permanova + rowSums(n.rejection, dims=3)
                } else {
                    n.permanova <- n.permanova + sum(n.rejection)
                }
                if (test.mediation) {
                    med.permanova.perm[,,i.perm] <- sapply(parallel.stat, med.permanova.fun, x.obs=permanova.obs, simplify="array") 
                    med.n.rejection = array((as.vector(med.permanova.perm[,,i.perm]) > as.vector(med.permanova.obs) + tol.eq) + (as.vector(med.permanova.perm[,,i.perm]) > as.vector(med.permanova.obs) - tol.eq), dim=c(n.dist, n.rarefy, n.perm.block))
                    if (n.dist>1 | n.rarefy>1) {
                        med.n.permanova <- med.n.permanova + rowSums(med.n.rejection, dims=2)
                    } else {
                        med.n.permanova <- med.n.permanova + sum(med.n.rejection)
                    }
                }
                
                if (!is.null(OR)) {
                    if (Sys.info()[['sysname']] == 'Windows') {
                        parallel.stat.OR = BiocParallel::bplapply(1:n.perm.block, parallel.perm, BPPARAM = BiocParallel::MulticoreParam(workers=n.cores), x.design.OR, perm, low, up, resid.dist.OR, ndf, adjust.for.confounders)
                    } else {
                        parallel.stat.OR = parallel::mclapply(1:n.perm.block, parallel.perm, mc.cores = n.cores, x.design.OR, perm, low, up, resid.dist.OR, ndf, adjust.for.confounders)
                    }
                    permanova.perm.OR[,,,i.perm] <- array(unlist(parallel.stat.OR), dim=c(dim(permanova.obs.OR), n.perm.block))
                    n.rejection.OR <- sapply(parallel.stat.OR, function(x) (x > permanova.obs.OR + tol.eq) + (x > permanova.obs.OR - tol.eq), simplify="array") 
                    if (n.dist>1 | n.var1>1 | n.rarefy>1) {
                        n.permanova.OR <- n.permanova.OR + rowSums(n.rejection.OR, dims=3)
                    } else {
                        n.permanova.OR <- n.permanova.OR + sum(n.rejection.OR)
                    }
                    if (test.mediation) {
                        med.permanova.perm.OR[,,i.perm] <- sapply(parallel.stat.OR, med.permanova.fun, x.obs=permanova.obs.OR, simplify="array") 
                        med.n.rejection.OR = array((as.vector(med.permanova.perm.OR[,,i.perm]) > as.vector(med.permanova.obs.OR) + tol.eq) + (as.vector(med.permanova.perm.OR[,,i.perm]) > as.vector(med.permanova.obs.OR) - tol.eq), dim=c(n.dist, n.rarefy, n.perm.block))
                        if (n.dist>1 | n.rarefy>1) {
                            med.n.permanova.OR <- med.n.permanova.OR + rowSums(med.n.rejection.OR, dims=2)
                        } else {
                            med.n.permanova.OR <- med.n.permanova.OR + sum(med.n.rejection.OR)
                        }
                    }
                }
                
                
            } else { # end of parallel computing
                
                for (i in 1:n.perm.block) {
                    
                    i.perm = n.perm.completed + i
                    
                    x.perm = x.design[perm[i,], ]   
                    for (d in 1:n.dist) {
                        permanova.perm[d,,,i.perm] = permanova.stat(x=x.perm, low=low, up=up, resid.dist=resid.dist[d,,,,,drop=FALSE], ndf=ndf, adjust.for.confounders=adjust.for.confounders)$permanova
                    }
                    n.permanova <- n.permanova + (array(permanova.perm[,,,i.perm], dim=dim(permanova.obs)) > permanova.obs + tol.eq) + (array(permanova.perm[,,,i.perm], dim=dim(permanova.obs)) > permanova.obs - tol.eq)
                    
                    if (test.mediation) {
                        med.permanova.null1 <- permanova.perm[,1,,i.perm] * permanova.obs[,2,]
                        med.permanova.null2 <- permanova.perm[,2,,i.perm] * permanova.obs[,1,]
                        med.permanova.null3 <- permanova.perm[,1,,i.perm] * permanova.perm[,2,,i.perm]
                        med.permanova.perm[,,i.perm] <- pmax(med.permanova.null1, med.permanova.null2, med.permanova.null3)
                        med.n.permanova = med.n.permanova + (array(med.permanova.perm[,,i.perm], dim=dim(med.permanova.obs)) > med.permanova.obs + tol.eq) + (array(med.permanova.perm[,,i.perm], dim=dim(med.permanova.obs)) > med.permanova.obs - tol.eq)
                    }
                    
                    if (!is.null(OR)) {
                        x.perm.OR = x.design.OR[perm[i,], ]   
                        for (d in 1:n.dist) {
                            permanova.perm.OR[d,,,i.perm] = permanova.stat(x=x.perm.OR, low=low, up=up, resid.dist=resid.dist.OR[d,,,,,drop=FALSE], ndf=ndf, adjust.for.confounders=adjust.for.confounders)$permanova
                        }
                        n.permanova.OR <- n.permanova.OR + (array(permanova.perm.OR[,,,i.perm], dim=dim(permanova.obs.OR)) > permanova.obs.OR + tol.eq) + (array(permanova.perm.OR[,,,i.perm], dim=dim(permanova.obs.OR)) > permanova.obs.OR - tol.eq)
                        
                        if (test.mediation) {
                            med.permanova.null1.OR <- permanova.perm.OR[,1,,i.perm] * permanova.obs.OR[,2,]
                            med.permanova.null2.OR <- permanova.perm.OR[,2,,i.perm] * permanova.obs.OR[,1,]
                            med.permanova.null3.OR <- permanova.perm.OR[,1,,i.perm] * permanova.perm.OR[,2,,i.perm]
                            med.permanova.perm.OR[,,i.perm] <- pmax(med.permanova.null1.OR, med.permanova.null2.OR, med.permanova.null3.OR)
                            med.n.permanova.OR = med.n.permanova.OR + (array(med.permanova.perm.OR[,,i.perm], dim=dim(med.permanova.obs.OR)) > med.permanova.obs.OR + tol.eq) + (array(med.permanova.perm.OR[,,i.perm], dim=dim(med.permanova.obs.OR)) > med.permanova.obs.OR - tol.eq)
                        }
                    }
                }
            } # non-parallel computing
            
            n.perm.completed = n.perm.completed + n.perm.block
            cat("permutations:", n.perm.completed, "\n")
            
            meet.all.rej.stop = all(n.permanova >= n.rej.stop*2)
            if (!is.null(OR)) meet.all.rej.stop = meet.all.rej.stop & all(n.permanova.OR >= n.rej.stop*2)
            if (test.mediation) {
                med.meet.all.rej.stop = all(med.n.permanova >= n.rej.stop*2)
                if (!is.null(OR)) med.meet.all.rej.stop = med.meet.all.rej.stop & all(med.n.permanova.OR >= n.rej.stop*2)
            }
            
            if (meet.all.rej.stop & med.meet.all.rej.stop) {
                permanova.stopped = TRUE
                cat("PERMANOVA stopped at permutation", n.perm.completed, "\n")
                break
            }
            
        }# permutation
        
        inv.n.perm.completed = 1/n.perm.completed
        inv.n.perm.completed.1 = 1/(n.perm.completed+1)
        
        p.permanova <- ifelse((n.permanova >= n.rej.stop*2), 0.5*n.permanova*inv.n.perm.completed, (0.5*n.permanova+1)*inv.n.perm.completed.1)
        if (!is.null(OR)) p.permanova.OR <- ifelse((n.permanova.OR >= n.rej.stop*2), 0.5*n.permanova.OR*inv.n.perm.completed, (0.5*n.permanova.OR+1)*inv.n.perm.completed.1)

        p.permanova.null <- NULL
        p.permanova.null.OR <- NULL
        
        # combination test
        if (!is.null(OR)) { 
            pmin.permanova.com <- 0.5*pmin(n.permanova, n.permanova.OR)
            p.permanova.null <- n.perm.completed + 0.5 - apply(permanova.perm[,,,1:n.perm.completed, drop=FALSE], c(1,2,3), rank)
            p.permanova.null.OR <- n.perm.completed + 0.5 - apply(permanova.perm.OR[,,,1:n.perm.completed, drop=FALSE], c(1,2,3), rank)
            pmin.permanova.com.null <- pmin(p.permanova.null, p.permanova.null.OR)
            if (length(dim(pmin.permanova.com.null))==4) {
                pmin.permanova.com.null <- aperm(pmin.permanova.com.null, c(2,3,4,1))
            } else {
                pmin.permanova.com.null <- array(pmin.permanova.com.null, c(dim(pmin.permanova.com.null), 1))
            }
            n.permanova.com <- rowSums( (pmin.permanova.com.null < c(pmin.permanova.com) - tol.eq) + 0.5 * (abs(pmin.permanova.com.null - c(pmin.permanova.com)) < tol.eq), dims=3) 
            p.permanova.com = ifelse((n.permanova.com >= n.rej.stop), n.permanova.com*inv.n.perm.completed, (n.permanova.com+1)*inv.n.perm.completed.1)
        }
        
        if (n.dist > 1) {
            pmin.permanova.omni <- apply(0.5*n.permanova , c(2,3), min)
            if (is.null(p.permanova.null)) p.permanova.null <- n.perm.completed + 0.5 - apply(permanova.perm[,,,1:n.perm.completed, drop=FALSE], c(1,2,3), rank)
            pmin.permanova.omni.null <- apply(p.permanova.null, c(1,3,4), min)
            if (length(dim(pmin.permanova.omni.null))==3) {
                pmin.permanova.omni.null <- aperm(pmin.permanova.omni.null, c(2,3,1))
            } else {
                pmin.permanova.omni.null <- array(pmin.permanova.omni.null, c(dim(pmin.permanova.omni.null), 1))
            }
            n.permanova.omni <- rowSums( (pmin.permanova.omni.null < c(pmin.permanova.omni) - tol.eq) + 0.5 * (abs(pmin.permanova.omni.null - c(pmin.permanova.omni)) < tol.eq), dims=2) 
            p.permanova.omni = ifelse((n.permanova.omni >= n.rej.stop), n.permanova.omni*inv.n.perm.completed, (n.permanova.omni+1)*inv.n.perm.completed.1)

            if (!is.null(OR)) {
                pmin.permanova.omni.OR <- apply(0.5*n.permanova.OR , c(2,3), min)
                if (is.null(p.permanova.null.OR)) p.permanova.null.OR <- n.perm.completed + 0.5 - apply(permanova.perm.OR[,,,1:n.perm.completed, drop=FALSE], c(1,2,3), rank)
                pmin.permanova.omni.null.OR <- apply(p.permanova.null.OR, c(1,3,4), min)
                if (length(dim(pmin.permanova.omni.null.OR))==3) {
                    pmin.permanova.omni.null.OR <- aperm(pmin.permanova.omni.null.OR, c(2,3,1))
                } else {
                    pmin.permanova.omni.null.OR <- array(pmin.permanova.omni.null.OR, c(dim(pmin.permanova.omni.null.OR), 1))
                }
                n.permanova.omni.OR <- rowSums( (pmin.permanova.omni.null.OR < c(pmin.permanova.omni.OR) - tol.eq) + 0.5 * (abs(pmin.permanova.omni.null.OR - c(pmin.permanova.omni.OR)) < tol.eq), dims=2) 
                p.permanova.omni.OR = ifelse((n.permanova.omni.OR >= n.rej.stop), n.permanova.omni.OR*inv.n.perm.completed, (n.permanova.omni.OR+1)*inv.n.perm.completed.1)

                # combination test
                pmin.permanova.omni.com <- pmin(pmin.permanova.omni, pmin.permanova.omni.OR)
                pmin.permanova.omni.null.com <- pmin(pmin.permanova.omni.null, pmin.permanova.omni.null.OR)
                n.permanova.omni.com <- rowSums( (pmin.permanova.omni.null.com < c(pmin.permanova.omni.com) - tol.eq) + 0.5 * (abs(pmin.permanova.omni.null.com - c(pmin.permanova.omni.com)) < tol.eq), dims=2) 
                p.permanova.omni.com = ifelse((n.permanova.omni.com >= n.rej.stop), n.permanova.omni.com*inv.n.perm.completed, (n.permanova.omni.com+1)*inv.n.perm.completed.1)
            }
        }
        
        med.p.permanova.null <- NULL
        med.p.permanova.null.OR <- NULL
        
        if (test.mediation) {
            med.p.permanova <- ifelse((med.n.permanova >= n.rej.stop*2), 0.5*med.n.permanova*inv.n.perm.completed, (0.5*med.n.permanova+1)*inv.n.perm.completed.1)
            if (!is.null(OR)) med.p.permanova.OR <- ifelse((med.n.permanova.OR >= n.rej.stop*2), 0.5*med.n.permanova.OR*inv.n.perm.completed, (0.5*med.n.permanova.OR+1)*inv.n.perm.completed.1)
            
            # combination test
            if (!is.null(OR)) {
                med.pmin.permanova.com <- 0.5*pmin(med.n.permanova, med.n.permanova.OR)
                med.p.permanova.null <- n.perm.completed + 0.5 - apply(med.permanova.perm[,,1:n.perm.completed, drop=FALSE], c(1,2), rank)
                med.p.permanova.null.OR <- n.perm.completed + 0.5 - apply(med.permanova.perm.OR[,,1:n.perm.completed, drop=FALSE], c(1,2), rank)
                med.pmin.permanova.com.null <- pmin(med.p.permanova.null, med.p.permanova.null.OR)
                if (length(dim(med.pmin.permanova.com.null))==3) {
                    med.pmin.permanova.com.null <- aperm(med.pmin.permanova.com.null, c(2,3,1))
                } else {
                    med.pmin.permanova.com.null <- array(med.pmin.permanova.com.null, c(dim(med.pmin.permanova.com.null), 1))
                }
                med.n.permanova.com <- rowSums( (med.pmin.permanova.com.null < c(med.pmin.permanova.com) - tol.eq) + 0.5 * (abs(med.pmin.permanova.com.null - c(med.pmin.permanova.com)) < tol.eq), dims=2)
                med.p.permanova.com = ifelse((med.n.permanova.com >= n.rej.stop), med.n.permanova.com*inv.n.perm.completed, (med.n.permanova.com+1)*inv.n.perm.completed.1)
            }
            
            # omnibus
            if (n.dist > 1) {
                med.pmin.permanova.omni <- apply(0.5*med.n.permanova, 2, min)
                if (is.null(med.p.permanova.null)) med.p.permanova.null <- n.perm.completed + 0.5 - apply(med.permanova.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                med.pmin.permanova.omni.null <- apply(med.p.permanova.null, c(1,3), min)
                if (length(dim(med.pmin.permanova.omni.null))==2) {
                    med.pmin.permanova.omni.null <- aperm(med.pmin.permanova.omni.null, c(2,1))
                } else {
                    med.pmin.permanova.omni.null <- array(med.pmin.permanova.omni.null, c(dim(med.pmin.permanova.omni.null), 1))
                }
                med.n.permanova.omni <- rowSums( (med.pmin.permanova.omni.null < c(med.pmin.permanova.omni) - tol.eq) + 0.5 * (abs(med.pmin.permanova.omni.null - c(med.pmin.permanova.omni)) < tol.eq)) 
                med.p.permanova.omni = ifelse((med.n.permanova.omni >= n.rej.stop), med.n.permanova.omni*inv.n.perm.completed, (med.n.permanova.omni+1)*inv.n.perm.completed.1)

                if (!is.null(OR)) {
                    med.pmin.permanova.omni.OR <- apply(0.5*med.n.permanova.OR , 2, min)
                    if (is.null(med.p.permanova.null.OR)) med.p.permanova.null.OR <- n.perm.completed + 0.5 - apply(med.permanova.perm.OR[,,1:n.perm.completed, drop=FALSE], c(1,2), rank)
                    med.pmin.permanova.omni.null.OR <- apply(med.p.permanova.null.OR, c(1,3), min)
                    if (length(dim(med.pmin.permanova.omni.null.OR))==2) {
                        med.pmin.permanova.omni.null.OR <- aperm(med.pmin.permanova.omni.null.OR, c(2,1))
                    } else {
                        med.pmin.permanova.omni.null.OR <- array(med.pmin.permanova.omni.null.OR, c(dim(med.pmin.permanova.omni.null.OR), 1))
                    }
                    med.n.permanova.omni.OR <- rowSums( (med.pmin.permanova.omni.null.OR < c(med.pmin.permanova.omni.OR) - tol.eq) + 0.5 * (abs(med.pmin.permanova.omni.null.OR - c(med.pmin.permanova.omni.OR)) < tol.eq))
                    med.p.permanova.omni.OR = ifelse((med.n.permanova.omni.OR >= n.rej.stop), med.n.permanova.omni.OR*inv.n.perm.completed, (med.n.permanova.omni.OR+1)*inv.n.perm.completed.1)

                    # combination test
                    med.pmin.permanova.omni.com <- pmin(med.pmin.permanova.omni, med.pmin.permanova.omni.OR)
                    med.pmin.permanova.omni.null.com <- pmin(med.pmin.permanova.omni.null, med.pmin.permanova.omni.null.OR)
                    med.n.permanova.omni.com <- rowSums( (med.pmin.permanova.omni.null.com < c(med.pmin.permanova.omni.com) - tol.eq) + 0.5 * (abs(med.pmin.permanova.omni.null.com - c(med.pmin.permanova.omni.com)) < tol.eq))
                    med.p.permanova.omni.com = ifelse((med.n.permanova.omni.com >= n.rej.stop), med.n.permanova.omni.com*inv.n.perm.completed, (med.n.permanova.omni.com+1)*inv.n.perm.completed.1)
                }
            }
        }
    }# if (n.perm.max > 0)
    
    name.list <- list(paste("dist", 1:n.dist, sep=""), 
                      paste("cov", 1:n.var1, sep=""), 
                      paste("rarefy", 1:n.rarefy, sep=""))
    dimnames(permanova.obs) <- name.list
    dimnames(R.squared) <- name.list
    if (!is.null(p.permanova) & length(p.permanova)!=1) dimnames(p.permanova) <- name.list
    if (!is.null(p.permanova.OR) & length(p.permanova.OR)!=1) dimnames(p.permanova.OR) <- name.list
    if (!is.null(p.permanova.com) & length(p.permanova.com)!=1) dimnames(p.permanova.com) <- name.list
    if (!is.null(p.permanova.omni) & length(p.permanova.omni)!=1) dimnames(p.permanova.omni) <- name.list[-1]
    if (!is.null(p.permanova.omni.OR) & length(p.permanova.omni.OR)!=1) dimnames(p.permanova.omni.OR) <- name.list[-1]
    if (!is.null(p.permanova.omni.com) & length(p.permanova.omni.com)!=1) dimnames(p.permanova.omni.com) <- name.list[-1]
    
    if (test.mediation) {
        if (!is.null(med.p.permanova) & length(med.p.permanova)!=1)dimnames(med.p.permanova) <- name.list[-2]
        if (!is.null(med.p.permanova.OR) & length(med.p.permanova.OR)!=1) dimnames(med.p.permanova.OR) <- name.list[-2]
        if (!is.null(med.p.permanova.com) & length(med.p.permanova.com)!=1) dimnames(med.p.permanova.com) <- name.list[-2]
        if (!is.null(med.p.permanova.omni) & length(med.p.permanova.omni)!=1) names(med.p.permanova.omni) <- name.list[[3]]
        if (!is.null(med.p.permanova.omni.OR) & length(med.p.permanova.omni.OR)!=1) names(med.p.permanova.omni.OR) <- name.list[[3]]
        if (!is.null(med.p.permanova.omni.com) & length(med.p.permanova.omni.com)!=1) names(med.p.permanova.omni.com) <- name.list[[3]]
    }
    
    res = list( F.statistics=drop(permanova.obs),
                R.squared=drop(R.squared),
                
                F.statistics.OR=drop(permanova.obs.OR),
                R.squared.OR=drop(R.squared.OR),
                
                p.permanova=drop(p.permanova), 
                p.permanova.omni=drop(p.permanova.omni),
                med.p.permanova=drop(med.p.permanova), 
                med.p.permanova.omni=drop(med.p.permanova.omni),
                
                p.permanova.OR=drop(p.permanova.OR), 
                p.permanova.omni.OR=drop(p.permanova.omni.OR),
                med.p.permanova.OR=drop(med.p.permanova.OR), 
                med.p.permanova.omni.OR=drop(med.p.permanova.omni.OR),
                
                p.permanova.com=drop(p.permanova.com), 
                p.permanova.omni.com=drop(p.permanova.omni.com),
                med.p.permanova.com=drop(med.p.permanova.com), 
                med.p.permanova.omni.com=drop(med.p.permanova.omni.com),
                
                n.perm.completed=n.perm.completed, 
                permanova.stopped=permanova.stopped,
                seed=seed)
    return(res)
    
} # End of permanovaFL


permanova.stat = function(x, low, up, resid.dist, ndf, adjust.for.confounders) {
    
    #---------------------------------------------
    #  calculate FL statistics for each model
    #---------------------------------------------
    n.var = length(low)
    n.sam = nrow(x)
    n.rarefy = dim(resid.dist)[5]
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    permanova = matrix(0, n.var1, n.rarefy)
    ve = matrix(0, n.var1+1, n.rarefy)
    
    for (r in 1:n.rarefy) {
        
        use = 1:up[n.var]
        Hcum = tcrossprod(x[, use]) # x[, use] %*% t( x[, use] )
        I_Hcum = diag(n.sam) - Hcum
        
        for (k in 1:n.var1) {
            
            k1 = k + as.numeric(adjust.for.confounders)
            
            use = low[k1]:up[k1]
            Hk = tcrossprod(x[, use])
            
            numerator =  sum(Hk * resid.dist[1,,,k,r]) # sum( diag(Hk %*% resid.dist[1,,,k,r]) )   # the other Hk and I_Hcum is not needed
            denominator = sum(I_Hcum * resid.dist[1,,,k,r]) # sum( diag(I_Hcum %*% resid.dist[1,,,k,r]) ) 

            ve[k, r:n.rarefy] = ve[k, r:n.rarefy] + numerator
            permanova[k, r:n.rarefy] = permanova[k, r:n.rarefy] + numerator/ denominator
        }
        ve[n.var1+1, r:n.rarefy] = ve[n.var1+1, r:n.rarefy] + denominator
    }
    
    R.squared = sweep(ve, 2, colSums(ve), FUN = "/")
    
    var = ifelse(adjust.for.confounders, 2:n.var, 1:n.var)
    out = list( R.squared=R.squared[1:n.var1,],
                permanova=permanova * ndf[n.var+1] / ndf[var]) 
    
    return(out)
    
} # permanova.stat


fit.permanova = function( d.gower, index, m, adjust.for.confounders) {
    
    n.var = length(index)
    n.otu = ncol(d.gower)
    n.sam = nrow(d.gower)
    ndf.nominal = rep(0, n.var+1)
    
    tol.d = 10^-8
    
    #--------------------------------------------------------------------------
    # construct directions matrix x 
    # from each set of covariates in the list vars
    #--------------------------------------------------------------------------
    
    d.resid = d.gower
    
    for (i in 1:n.var) 
    {
        
        var = m[,1:index[i]]
        
        svd.var = svd(var)   
        use = (svd.var$d > tol.d)    
        
        hat.matrix = tcrossprod(svd.var$u[, use]) # svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        #---------------------
        # calculate direction
        #---------------------
        
        n.dim = dim( hat.matrix)[1]
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) # es: eigen system in Mathematica
        
        use = ( abs(es.model$values) > tol.d )
        ndf.model = sum( use )
        
        x.model = es.model$vectors[, use]
        e.model = es.model$values[use]
        
        hat.matrix.bar = diag(n.dim)  - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # end of calculating direction
        #-----------------------------    
        
        if (i==1) {
            x = x.model
            e = e.model
        } else {   
            x = cbind(x, x.model)
            e = c(e, e.model )
        }
        
        ndf.nominal[i] = ndf.model
        
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE)
    use = which( abs(es.resid$values) > tol.d )
    
    ndf.nominal[n.var+1] = length(use)
    x = cbind(x, es.resid$vectors[, use])
    e = c(e, es.resid$values[use])
    
    #-------------------------------------------------
    # low, up
    #-------------------------------------------------
    
    low = rep(NA, n.var)
    up = rep(NA, n.var)
    
    up.prev = 0
    
    for (k in 1:n.var)
    {
        low[k] = up.prev + 1
        up[k] = up.prev + ndf.nominal[k]
        up.prev = up[k]
    }
    
    #---------------------
    # permanova: resid.dist
    #---------------------
    
    if (n.var==1) {
        resid.dist = array( NA, dim=c( dim(d.gower), n.var ) ) 
        resid.dist[,,1] = d.gower
    }
    else {
        n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
        
        resid.dist = array( NA, dim=c( dim(d.gower), n.var1 ) ) 
        
        for (k in 1:n.var1) {
            k1 = ifelse(adjust.for.confounders, k+1, k)
            use = setdiff( 1:up[n.var], low[k1]:up[k1] )
            
            hat.matrix = tcrossprod(x[,use,drop=FALSE]) # x[,use,drop=FALSE] %*% t( x[,use,drop=FALSE] )
            hat.matrix.bar = diag(n.sam) - hat.matrix
            resid.dist[,,k] = hat.matrix.bar %*% d.gower 
            resid.dist[,,k] = resid.dist[,,k] %*% hat.matrix.bar
        }
    }
    
    
    res = list( x=x,
                low=low,
                up=up,
                resid.dist=resid.dist,
                ndf = ndf.nominal)
    return(res)
    
} # fit.permanova


#################################
# Expectation of distance matrix
#################################

#' Expected value of the Jaccard distance matrix
#' 
#' This function computes the expected value of the Jaccard distance matrix over rarefaction replicates.
#' 
#' @param otu.table the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' @param rarefy.depth rarefaction depth. The default is the minimum library size observed in the OTU table.
#' @param first.order.approx.only a logical value indicating whether to calculate the expected value 
#' using the first order approixmation by the delta method. 
#' The default is FALSE, using the second order approixmation.
#' @return a list consisting of 
#'   \item{jac.mean.o1}{Expected Jaccard distance matrix by the first order approixmation.}
#'   \item{jac.mean.o2}{Expected Jaccard distance matrix by the second order approixmation.} 
#'   \item{jac.mean.sq.o1}{Expected squared Jaccard distance matrix by the first order approixmation.} 
#'   \item{jac.mean.sq.o2}{Expected squared Jaccard distance matrix by the second order approixmation.} 
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gsatten@emory.edu>
#' @export
#' @examples
#' data(throat.otu.tab5)
#' res.jaccard <- jaccard.mean( throat.otu.tab5 )

jaccard.mean = function( otu.table, rarefy.depth=min(rowSums(otu.table)), first.order.approx.only=FALSE) {
    if (!first.order.approx.only) {
        out = jaccard.mean.o1o2(otu.table, rarefy.depth)
    }
    else {
        res = jaccard.mean.o1(otu.table, rarefy.depth)
        out=list( jac.mean.o1=res$jac.mean.1, 
                  jac.mean.o2=NA, 
                  jac.mean.sq.o1=res$jac.mean.1^2, 
                  jac.mean.sq.o2=NA)
    }
    return(out)
}

# (old) jaccard.mean.fast
jaccard.mean.o1o2 = function( otu.table, rarefy.depth=min(rowSums(otu.table)) ) {
    otu.table=as.matrix(otu.table)
    n.obs=dim(otu.table)[1]
    n.taxa=dim(otu.table)[2]
    n.elem=n.taxa*(n.taxa-1)/2
    m.elem=max( rowSums( ifelse(otu.table>0,1,0) ) )
    index.mu=matrix(0,nrow=n.obs, ncol=m.elem)
    n.index.mu=rep(0,n.obs)
    m.elem=m.elem*(m.elem-1)/2
    #
    #	calculate mu, and mu.dot
    #
    d=matrix(0, nrow=n.obs, ncol=n.elem)
    n.i=rowSums( otu.table )
    l.denom=lchoose(n.i, rarefy.depth)
    jac.sq=diag(n.obs)
    mu=matrix(0,nrow=n.obs, ncol=n.taxa)
    mu.dot=rep(0,n.obs)
    for (i in 1:n.obs) {
        use=which( otu.table[i,]>0 )
        l.pr=lchoose( n.i[i] - otu.table[i,use], rarefy.depth ) - l.denom[i]
        mu[i,use]=1 - exp(l.pr)
        #		mu[i,use]=1 - choose( n.i[i] - otu.table[i,use], rarefy.depth)/denom[i]
        mu.dot[i]=sum( mu[i,use] )
        n.index.mu[i]=length(use)
        index.mu[i, 1:n.index.mu[i] ]=use
    }
    #
    #	calculate d, psi.11 and psi.12
    #		
    psi.11=psi.12=psi.21=matrix(0, nrow=n.obs, ncol=n.obs)
    psi.22=matrix(0, nrow=n.obs, ncol=n.obs)
    index.d=matrix(0, nrow=n.obs, ncol=m.elem)
    n.index.d=rep(0,n.obs)
    for (i in 1:n.obs) {
        use=which( otu.table[i,]>0 )
        n.use=length(use)
        k.1=use[rep( 1:(n.use-1), times=(n.use-1):1 )]
        k.2=use[unlist(sapply(2:n.use, FUN = function(x){x:n.use}))]
        k=n.taxa*(k.1-1) + k.2 - k.1*(k.1 + 1)/2
        ln.pr=lchoose( n.i[i] - otu.table[i,k.1] - otu.table[i,k.2], rarefy.depth ) - l.denom[i]
        d[i,k]=exp(ln.pr) + mu[i,k.1] + mu[i,k.2] - 1
        #		
        cs.1=colSums( t(mu[,k.1])*d[i,k] )
        cs.2=colSums( t(mu[,k.2])*d[i,k] )
        cs.d=sum( d[i,k] )
        #		
        psi.11[i,]=cs.d -   cs.1 -   cs.2
        psi.12[i,]=cs.d -   cs.1 - 2*cs.2
        psi.21[i,]=cs.d - 2*cs.1 -   cs.2
        psi.22[i,]=cs.d - 2*cs.1 - 2*cs.2
        #		
        n.index.d[i]=n.use*(n.use-1)/2
        index.d[i,1:n.index.d[i]]=k
        
        # psi.11[i,]=colSums( t(1-mu[,k.1]-mu[,k.2]) * d[i,k] )
        # psi.12[i,]=colSums( t(1-mu[,k.1]-2*mu[,k.2]) * d[i,k] )
        # psi.21[i,]=colSums( t(1-2*mu[,k.1]-mu[,k.2]) * d[i,k] )
        # psi.22[i,]=colSums( t(1-2*mu[,k.1]-2*mu[,k.2]) * d[i,k] )
    }
    #
    #	calculate 2nd order approximation to average jaccard 
    #	
    jac.mean.1=jac.mean.2=matrix(0, nrow=n.obs, ncol=n.obs)
    jac.mean.sq.1=jac.mean.sq.2=matrix(0, nrow=n.obs, ncol=n.obs)
    jac.sq.ave=u22=u11=u12=matrix(0, nrow=n.obs, ncol=n.obs)
    for (i in 1:(n.obs-1)) {
        for (j in (i+1):n.obs) {
            #			use=which( otu.table[i,]*otu.table[j,]>0 )
            use=intersect( index.mu[i,1:n.index.mu[i]], index.mu[j,1:n.index.mu[j]] )
            dot.product.mu=sum( mu[i,use]*mu[j,use] )
            use=intersect( index.d[i,1:n.index.d[i]], index.d[j, 1:n.index.d[j]] )
            dot.product.d=sum( d[i,use]*d[j,use] )
            e.1=mu.dot[i]+mu.dot[j]-dot.product.mu
            e.2=e.1 - dot.product.mu
            u.11=mu.dot[i] + mu.dot[j] + 2*mu.dot[i]*mu.dot[j] - 3*dot.product.mu + 2*(psi.11[i,j]+psi.11[j,i]) + 2*dot.product.d
            u.12=mu.dot[i] + mu.dot[j] + 2*mu.dot[i]*mu.dot[j] - 4*dot.product.mu + (psi.12[i,j]+psi.12[j,i]+psi.21[i,j]+psi.21[j,i]) + 4*dot.product.d
            u.22=mu.dot[i] + mu.dot[j] + 2*mu.dot[i]*mu.dot[j] - 4*dot.product.mu + 2*(psi.22[i,j]+psi.22[j,i]) + 8*dot.product.d
            jac.mean.1[i,j]=e.2/e.1
            jac.mean.2[i,j]=jac.mean.1[i,j] + (u.11*e.2 - u.12*e.1)/(e.1^3)
            jac.mean.sq.1[i,j]=jac.mean.1[i,j]^2
            jac.mean.sq.2[i,j]=jac.mean.sq.1[i,j] + (u.22*e.1^2-4*u.12*e.1*e.2+3*u.11*e.2^2)/(e.1^4)
            jac.sq.ave[i,j]=u.22/u.11
            u11[i,j]=u.11
            u22[i,j]=u.22
            u12[i,j]=u.12
        }
    }	
    jac.mean.1=jac.mean.1 + t(jac.mean.1)
    jac.mean.2=jac.mean.2 + t(jac.mean.2)
    jac.mean.sq.1=jac.mean.sq.1 + t(jac.mean.sq.1)
    jac.mean.sq.2=jac.mean.sq.2 + t(jac.mean.sq.2)
    jac.sq.ave=jac.sq.ave + t(jac.sq.ave)	
    u11=u11 + t(u11)
    u22=u22 + t(u22)
    u12=u12 + t(u12)
    
    res=list( jac.mean.o1=jac.mean.1, 
              jac.mean.o2=jac.mean.2, 
              jac.mean.sq.o1=jac.mean.sq.1, 
              jac.mean.sq.o2=jac.mean.sq.2)
    return(res)
    
} # jaccard.mean.fast


#(old) jaccard.ave.1
jaccard.mean.o1 = function( otu.table, rarefy.depth=min(rowSums(otu.table)) ) {
    #
    #	calculates only the zeroth order (first term) of the jaccard mean 
    #
    otu.table=as.matrix(otu.table)
    n.obs=dim(otu.table)[1]
    n.taxa=dim(otu.table)[2]
    n.elem=n.taxa*(n.taxa-1)/2
    m.elem=max( rowSums( ifelse(otu.table>0,1,0) ) )
    index.mu=matrix(0,nrow=n.obs, ncol=m.elem)
    n.index.mu=rep(0,n.obs)
    m.elem=m.elem*(m.elem-1)/2
    #
    #	calculate mu, and mu.dot
    #
    n.i=rowSums( otu.table )
    l.denom=lchoose(n.i, rarefy.depth)
    mu=matrix(0,nrow=n.obs, ncol=n.taxa)
    mu.dot=rep(0,n.obs)
    for (i in 1:n.obs) {
        use=which( otu.table[i,]>0 )
        l.pr=lchoose( n.i[i] - otu.table[i,use], rarefy.depth ) - l.denom[i]
        mu[i,use]=1 - exp(l.pr)
        #		mu[i,use]=1 - choose( n.i[i] - otu.table[i,use], rarefy.depth)/denom[i]
        mu.dot[i]=sum( mu[i,use] )
        n.index.mu[i]=length(use)
        index.mu[i, 1:n.index.mu[i] ]=use
    }
    #
    #	calculate 2nd order approximation to average jaccard 
    #	
    jac.mean.1=e1=e2=matrix(0, nrow=n.obs, ncol=n.obs)
    for (i in 1:(n.obs-1)) {
        for (j in (i+1):n.obs) {
            use=intersect( index.mu[i,1:n.index.mu[i]], index.mu[j,1:n.index.mu[j]] )
            dot.product.mu=sum( mu[i,use]*mu[j,use] )
            e.1=mu.dot[i]+mu.dot[j]-dot.product.mu
            e.2=e.1 - dot.product.mu
            jac.mean.1[i,j]=e.2/e.1
            e1[i,j]=e.1
            e2[i,j]=e.2
        }
    }	
    jac.mean.1=jac.mean.1 + t(jac.mean.1)
    e1=e1 + t(e1)
    e2=e2 + t(e2)
    
    res=list( jac.mean.o1=jac.mean.1)
    return(res)
    
} # jaccard.ave.1



#' Expected value of the unweighted UniFrac distance matrix
#' 
#' This function computes the expected value of the unweighted UniFrac distance matrix over rarefaction replicates.
#' 
#' @param otu.table the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' @param tree the phylogeneic tree.
#' @param rarefy.depth rarefaction depth. The default is the minimum library size observed in the OTU table.
#' @param first.order.approx.only a logical value indicating whether to calculate the expected value 
#' using the first order approixmation by the delta method. The default is FALSE, 
#' using the second order approixmation.
#' @return a list consisting of 
#'   \item{unifrac.mean.o1}{Expected unweighted UniFrac distance matrix by the first order approixmation.}
#'   \item{unifrac.mean.o2}{Expected unweighted UniFrac distance matrix by the second order approixmation.} 
#'   \item{unifrac.mean.sq.o1}{Expected squared unweighted UniFrac distance matrix by the first order approixmation.} 
#'   \item{unifrac.mean.sq.o2}{Expected squared unweighted UniFrac distance matrix by the second order approixmation.} 
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gsatten@emory.edu>
#' @importFrom castor get_subtree_with_tips
#' @importFrom phangorn Ancestors Descendants
#' @export
#' @examples
#' data(throat.otu.tab5)
#' data(throat.tree)
#' res.unifrac <- unifrac.mean( throat.otu.tab5, throat.tree)

unifrac.mean = function( otu.table, tree, rarefy.depth=min(rowSums(otu.table)), first.order.approx.only=FALSE) {
    if (!first.order.approx.only) {
        out = unifrac.mean.o1o2(otu.table, tree, rarefy.depth)
    }
    else {
        res = unifrac.mean.o1(otu.table, tree, rarefy.depth)
        out=list( unifrac.mean.o1=res$unifrac.mean.1, 
                  unifrac.mean.o2=NA, 
                  unifrac.mean.sq.o1=res$unifrac.mean.1^2, 
                  unifrac.mean.sq.o2=NA)
    }
    return(out)
}

`%notin%` <- Negate(`%in%`)


#(old) unifrac.ave.sq.fast
unifrac.mean.o1o2 = function( otu.table, tree, rarefy.depth=min(rowSums(otu.table)), trim=TRUE, keep.root=TRUE) {
    
    # requireNamespace('ips')
    # requireNamespace('castor')
    # requireNamespace('phangorn')
    
    otu.table=as.matrix(otu.table)
    n.obs=dim(otu.table)[1]
    n.taxa=dim(otu.table)[2]                # note for later:  remove cols with zero colSums?
    #	
    #
    #	sort OTU table to agree with ordering of tips in tree
    #
    tree.names=tree$tip.label
    otu.names=colnames(otu.table)
    #
    #	if trim=TRUE remove any taxa that do not appear in OTU table and then prune tree accordingly
    #
    if (trim==TRUE) {
        keep=which(colSums( otu.table )>0)
        if (sum(keep)<n.taxa) {
            print('number of empty taxa dropped from OTU table:')
            print( n.taxa-sum(keep) )
            otu.table=otu.table[,keep]
            n.taxa=dim(otu.table)[2]
            otu.names=colnames(otu.table)
        }
        if (length(tree.names)<length(otu.names)) {
            print('error: there are more taxa in OTU table than tree')
            return
        }
        else if (length(tree.names)>length(otu.names)) {
            drop.names=setdiff( tree.names, otu.names )
            drop.list=which( tree.names %in% drop.names)
            tree=castor::get_subtree_with_tips(tree,omit_tips=drop.list, force_keep_root=keep.root)$subtree
            #			tree=drop.tip(tree,drop.list,collapse.singles=!keep.root)
            tree.names=tree$tip.label
            otu.table=otu.table[,tree.names]
            print('number of tips dropped from tree:')
            print( length(drop.list) )
        }
    }	
    if( !all(sort(tree.names)==sort(otu.names)) ) {
        print('taxa names in tree and OTU table to not agree')
        return
    }
    otu.table=otu.table[,tree.names]
    #
    #	set up branch lengths
    #	
    edge.2=tree$edge[,2]
    n.internal=max(edge.2) - n.taxa - 1
    n.branch=length(edge.2)
    n.node=n.taxa + n.internal
    n.elem=n.node*(n.node-1)/2
    edge.2=tree$edge[,2]
    tips=which( edge.2<=n.taxa )
    nodes=which( edge.2>n.taxa+1)
    #	edge.2[tips]=edge.2[tips]+n.taxa-2
    edge.2[tips]=edge.2[tips]+n.internal
    edge.2[nodes]=edge.2[nodes]-n.taxa-1
    ord=order(edge.2)
    edge.2=edge.2[ord]
    br.length=tree$edge.length[ord]
    #
    #	set up list of descendants for each node (nodes  1:(n.internal) are internal nodes; nodes (n.internal):(n.internal+n.taxa-1) are tips)
    #	
    #	descendants=allDescendants(tree)
    descendants=phangorn::Descendants(tree,node=c((n.taxa+2):(n.taxa+1+n.internal), 1:n.taxa), type='tips')
    # k.1=rep( 1:(n.node-1), times=(n.node-1):1 )
    # k.2=unlist(sapply(2:n.node, FUN = function(x){x:n.node}))
    # k=n.node*(k.1-1) + k.2 - k.1*(k.1+1)/2
    ancestry=matrix(FALSE,nrow=n.taxa, ncol=n.node)
    for (k in 1:n.node) {
        ancestry[ descendants[[k]], k]=TRUE
    }
    
    #
    #	set up list of tips and internal nodes of ancestors for each observation; calculate mu, mu.bar and mu.bar.2
    #	
    tips.i=list()
    #	nodes.i=list()
    all.nodes.i=list()
    mu=mu.bar=mu.bar.2=matrix(0,nrow=n.obs, ncol=n.node)
    mu.dot.1=rep(0,n.obs)
    mu.dot.2=rep(0,n.obs)
    n.i=rowSums( otu.table )
    l.denom=lchoose(n.i, rarefy.depth)
    for (i in 1:n.obs) {
        tips=which( otu.table[i,]>0 )
        nodes=sort( unique( unlist( phangorn::Ancestors(x=tree, node=tips, type='all') ) ) ) - n.taxa - 1 
        nodes=nodes[nodes>0]
        all.nodes=union(nodes,tips+n.internal)
        tips.i=append(tips.i, list(tips) )
        all.nodes.i=append(all.nodes.i, list( all.nodes ))
        n.nodes=length(all.nodes)
        #		use.k=lapply( descendants[all.nodes], FUN=function(x) intersect(x,tips) )
        #		c.sum.k=unlist( lapply( use.k, FUN=function(x) sum( otu.table[i,x] ) ) )
        c.sum.k=colSums( otu.table[i,tips]*ancestry[tips,all.nodes] )
        mu.bar[i,all.nodes]=exp( lchoose( n.i[i]-c.sum.k, rarefy.depth ) - l.denom[i] )
        mu[i,all.nodes]=1-mu.bar[i,all.nodes]
        mu.bar.2[i,all.nodes]=1-mu.bar[i,all.nodes]
        mu.dot.1[i]=sum( mu[i,all.nodes]*br.length[all.nodes] )
        mu.dot.2[i]=sum( mu[i,all.nodes]*br.length[all.nodes]^2 )
    }
    #
    #   calculate d, psi.11, psi.12 and psi.22
    #	
    d=matrix(0,nrow=n.obs, ncol=n.elem)
    psi.11=psi.12=psi.21=matrix(0, nrow=n.obs, ncol=n.obs)
    psi.22=matrix(0, nrow=n.obs, ncol=n.obs)
    
    
    for (i in 1:n.obs) {
        all.nodes=all.nodes.i[[i]]
        tips=tips.i[[i]]
        n.nodes=length(all.nodes)
        k.1=all.nodes[rep( 1:(n.nodes-1), times=(n.nodes-1):1 )]
        k.2=all.nodes[unlist(sapply(2:n.nodes, FUN = function(x){x:n.nodes}))]
        k=n.node*(k.1-1) + k.2 - k.1*(k.1+1)/2
        #		c.sum.k1.k2=mapply( descendants[k.1], descendants[k.2], FUN=function(x,y) sum( otu.table[i,intersect( union(x,y),tips )]), SIMPLIFY=TRUE)  #this is slower than loop that follows
        #		n.k=length(k)
        #		c.sum.k1.k2=rep(0, n.k)
        #		for (j in 1:n.k) {
        #			c.sum.k1.k2[j]=sum( otu.table[i, intersect( union( descendants[[k.1[j]]], descendants[[k.2[j]]] ), tips) ] )
        #			}
        
        c.sum.k1.k2=colSums( otu.table[i,tips]*(ancestry[tips,k.1] | ancestry[tips,k.2]) )
        
        
        d[i,k]=exp( lchoose( n.i[i]-c.sum.k1.k2, rarefy.depth ) - l.denom[i] ) + mu[i,k.1] + mu[i,k.2] - 1
        #		
        cs.1=colSums( t(mu[,k.1])*d[i,k]*br.length[k.1]*br.length[k.2] )
        cs.2=colSums( t(mu[,k.2])*d[i,k]*br.length[k.1]*br.length[k.2] )
        cs.d=sum( d[i,k]*br.length[k.1]*br.length[k.2] )
        #		
        psi.11[i,]=cs.d -   cs.1 -   cs.2
        psi.12[i,]=cs.d -   cs.1 - 2*cs.2
        psi.21[i,]=cs.d - 2*cs.1 -   cs.2
        psi.22[i,]=cs.d - 2*cs.1 - 2*cs.2
        #		
    }	
    #
    #	calculate 2nd order approximation to average unifrac
    #	
    unifrac.mean.1=unifrac.mean.2=matrix(0, nrow=n.obs, ncol=n.obs)
    unifrac.mean.sq.1=unifrac.mean.sq.2=unifrac.sq.ave=matrix(0, nrow=n.obs, ncol=n.obs)
    u11=u12=u22=e1=e2=matrix(0,nrow=n.obs, ncol=n.obs)
    #
    k.index=rep( 1:(n.node-1), times=(n.node-1):1 )
    k1.index=unlist(sapply(2:n.node, FUN = function(x){x:n.node}))
    #
    for (i in 1:(n.obs-1)) {
        for (j in (i+1):n.obs) {
            use=which( mu[i,]*mu[j,]>0 )
            dot.product.mu.1=sum( mu[i,use]*mu[j,use]*br.length[use] )
            dot.product.mu.2=sum( mu[i,use]*mu[j,use]*br.length[use]^2 )
            use=which( d[i,]*d[j,]>0 )
            k.node=k.index[use]
            k1.node=k1.index[use]
            dot.product.d=sum( d[i,use]*d[j,use]*br.length[k.node]*br.length[k1.node] )
            e.1=mu.dot.1[i]+mu.dot.1[j]-dot.product.mu.1
            e.2=e.1 - dot.product.mu.1
            u.11=mu.dot.2[i] + mu.dot.2[j] + 2*mu.dot.1[i]*mu.dot.1[j] - 3*dot.product.mu.2 + 2*(psi.11[i,j]+psi.11[j,i]) + 2*dot.product.d
            u.12=mu.dot.2[i] + mu.dot.2[j] + 2*mu.dot.1[i]*mu.dot.1[j] - 4*dot.product.mu.2 + (psi.12[i,j]+psi.12[j,i]+psi.21[i,j]+psi.21[j,i]) + 4*dot.product.d
            u.22=mu.dot.2[i] + mu.dot.2[j] + 2*mu.dot.1[i]*mu.dot.1[j] - 4*dot.product.mu.2 + 2*(psi.22[i,j]+psi.22[j,i]) + 8*dot.product.d
            unifrac.mean.1[i,j]=e.2/e.1
            unifrac.mean.2[i,j]=unifrac.mean.1[i,j] + (u.11*e.2 - u.12*e.1)/(e.1^3)
            unifrac.mean.sq.1[i,j]=unifrac.mean.1[i,j]^2
            unifrac.mean.sq.2[i,j]=unifrac.mean.sq.1[i,j] + (u.22*e.1^2-4*u.12*e.1*e.2+3*u.11*e.2^2)/(e.1^4)
            unifrac.sq.ave[i,j]=u.22/u.11
            u11[i,j]=u.11
            u12[i,j]=u.12
            u22[i,j]=u.22
            e1[i,j]=e.1
            e2[i,j]=e.2
        }
    }	
    unifrac.mean.1=unifrac.mean.1 + t(unifrac.mean.1)
    unifrac.mean.2=unifrac.mean.2 + t(unifrac.mean.2)
    unifrac.mean.sq.1=unifrac.mean.sq.1 + t(unifrac.mean.sq.1)
    unifrac.mean.sq.2=unifrac.mean.sq.2 + t(unifrac.mean.sq.2)
    unifrac.sq.ave=unifrac.sq.ave + t(unifrac.sq.ave)	
    
    u11=u11 + t(u11)
    u12=u12 + t(u12)
    u22=u22 + t(u22)
    e1=e1 + t(e1)
    e2=e2 + t(e2)
    
    res=list( unifrac.mean.o1=unifrac.mean.1, 
              unifrac.mean.o2=unifrac.mean.2, 
              unifrac.mean.sq.o1=unifrac.mean.sq.1, 
              unifrac.mean.sq.o2=unifrac.mean.sq.2)
    return(res)
    
} #unifrac.ave.sq.fast

# (old) unifrac.ave.1
unifrac.mean.o1 = function( otu.table, tree, rarefy.depth=min(rowSums(otu.table)), trim=TRUE, keep.root=TRUE) {
    #
    #	calculates only the zeroth order (first) term for the unifrac distance
    #
    # requireNamespace('ips')
    # requireNamespace('castor')
    # requireNamespace('phangorn')
    
    otu.table=as.matrix(otu.table)
    n.obs=dim(otu.table)[1]
    n.taxa=dim(otu.table)[2]                
    if (trim==TRUE) {	
        #
        #		remove zero columns from OTU table
        #		
        keep.list=which( colSums(otu.table)>0 )
        otu.table=otu.table[,keep.list]
        print( c('note: empty taxa removed, number of taxa changed from',n.taxa,'to',length(keep.list) ) )
        n.taxa=dim(otu.table)[2]
    }
    #	
    #
    #	harmonize lists of taxa from tree, OTU table and then sort OTU table to agree with ordering of tips in tree
    #
    tree.names=tree$tip.label
    otu.names=colnames(otu.table)
    n.tree.names=length(tree.names)
    n.otu.names=length(otu.names)
    if (n.otu.names<n.tree.names) {
        #
        #		remove tree nodes that do not correspond to taxa in otu.table
        #	
        drop.list=which( tree.names %notin% otu.names )
        tree=castor::get_subtree_with_tips(tree,omit_tips=drop.list, force_keep_root=keep.root)$subtree
        tree.names=tree$tip.label
        print( c('note: number of nodes in tree reduced from', n.tree.names, 'to', length(tree.names) ), quote=FALSE)
    }
    if( !all(sort(tree.names)==sort(otu.names)) ) {
        #
        #		remove otus that are not in tree
        #	
        keep.list=which( otu.names %in% tree.names )
        otu.table=otu.table[,keep.list]
        n.taxa=dim(otu.table)[2]
        print( c('note: number of taxa in OTU table reduced from', n.otu.names, 'to', n.taxa), quote=FALSE )
    }
    otu.table=otu.table[,tree.names]
    #
    #	set up branch lengths
    #	
    edge.2=tree$edge[,2]
    n.internal=max(edge.2) - n.taxa - 1
    n.branch=length(edge.2)
    n.node=n.taxa + n.internal
    n.elem=n.node*(n.node-1)/2
    edge.2=tree$edge[,2]
    tips=which( edge.2<=n.taxa )
    nodes=which( edge.2>n.taxa+1)
    #	edge.2[tips]=edge.2[tips]+n.taxa-2
    edge.2[tips]=edge.2[tips]+n.internal
    edge.2[nodes]=edge.2[nodes]-n.taxa-1
    ord=order(edge.2)
    edge.2=edge.2[ord]
    br.length=tree$edge.length[ord]
    #
    #	set up list of descendants for each node (nodes  1:(n.internal) are internal nodes; nodes (n.internal):(n.internal+n.taxa-1) are tips)
    #	
    #	descendants=allDescendants(tree)
    descendants=phangorn::Descendants(tree,node=c((n.taxa+2):(n.taxa+1+n.internal), 1:n.taxa), type='tips')
    # k.1=rep( 1:(n.node-1), times=(n.node-1):1 )
    # k.2=unlist(sapply(2:n.node, FUN = function(x){x:n.node}))
    # k=n.node*(k.1-1) + k.2 - k.1*(k.1+1)/2
    ancestry=matrix(FALSE,nrow=n.taxa, ncol=n.node)
    for (k in 1:n.node) {
        ancestry[ descendants[[k]], k]=TRUE
    }
    
    #
    #	set up list of tips and internal nodes of ancestors for each observation; calculate mu, mu.bar and mu.bar.2
    #	
    mu=matrix(0,nrow=n.obs, ncol=n.node)
    n.i=rowSums( otu.table )
    l.denom=lchoose(n.i, rarefy.depth)
    mu.dot.1=rep(0,n.obs)
    for (i in 1:n.obs) {
        tips=which( otu.table[i,]>0 )
        nodes=sort( unique( unlist( phangorn::Ancestors(x=tree, node=tips, type='all') ) ) ) - n.taxa - 1 
        nodes=nodes[nodes>0]
        all.nodes=union(nodes,tips+n.internal)
        #		use.k=lapply( descendants[all.nodes], FUN=function(x) intersect(x,tips) )
        #		c.sum.k=unlist( lapply( use.k, FUN=function(x) sum( otu.table[i,x] ) ) )
        c.sum.k=colSums( otu.table[i,tips]*ancestry[tips,all.nodes] )
        mu[i,all.nodes]=1- exp( lchoose( n.i[i]-c.sum.k, rarefy.depth ) - l.denom[i] )
        mu.dot.1[i]=sum( mu[i,all.nodes]*br.length[all.nodes] )
    }
    
    #		
    unifrac.mean.1=e1=e2=matrix(0, nrow=n.obs, ncol=n.obs)
    #
    for (i in 1:(n.obs-1)) {
        for (j in (i+1):n.obs) {
            use=which( mu[i,]*mu[j,]>0 )
            dot.product.mu.1=sum( mu[i,use]*mu[j,use]*br.length[use] )
            e.1=mu.dot.1[i]+mu.dot.1[j]-dot.product.mu.1
            e.2=e.1 - dot.product.mu.1
            unifrac.mean.1[i,j]=e.2/e.1
            e1[i,j]=e.1
            e2[i,j]=e.2
        }
    }	
    unifrac.mean.1=unifrac.mean.1 + t(unifrac.mean.1)
    
    e1=e1 + t(e1)
    e2=e2 + t(e2)
    
    res=list( unifrac.mean.o1=unifrac.mean.1)
    return(res)
    
} # unifrac.ave.1


#####################################################################
# no use (memory efficient version)
#####################################################################

jaccard.mean.fast.small = function( otu.table, rarefy.depth=min(rowSums(otu.table)) ) {
    otu.table=as.matrix(otu.table)
    n.obs=dim(otu.table)[1]
    n.taxa=dim(otu.table)[2]
    n.elem=n.taxa*(n.taxa-1)/2
    m.elem=max( rowSums( ifelse(otu.table>0,1,0) ) )
    index.mu=matrix(0,nrow=n.obs, ncol=m.elem)
    n.index.mu=rep(0,n.obs)
    m.elem=m.elem*(m.elem-1)/2
    #
    #	calculate mu, and mu.dot
    #
    n.i=rowSums( otu.table )
    l.denom=lchoose(n.i, rarefy.depth)
    jac.sq=diag(n.obs)
    mu=matrix(0,nrow=n.obs, ncol=n.taxa)
    mu.dot=rep(0,n.obs)
    for (i in 1:n.obs) {
        use=which( otu.table[i,]>0 )
        l.pr=lchoose( n.i[i] - otu.table[i,use], rarefy.depth ) - l.denom[i]
        mu[i,use]=1 - exp(l.pr)
        #		mu[i,use]=1 - choose( n.i[i] - otu.table[i,use], rarefy.depth)/denom[i]
        mu.dot[i]=sum( mu[i,use] )
        n.index.mu[i]=length(use)
        index.mu[i, 1:n.index.mu[i] ]=use
    }
    #
    #	calculate d, psi.11 and psi.12
    #
    max.taxa=max( rowSums( ifelse(otu.table>0,1,0) ) )
    n.elem.d=max.taxa*(max.taxa-1)/2
    psi.11=psi.12=psi.21=matrix(0, nrow=n.obs, ncol=n.obs)
    psi.22=matrix(0, nrow=n.obs, ncol=n.obs)
    d=matrix(0, nrow=n.obs, ncol=n.elem.d)
    index.d=matrix(0, nrow=n.obs, ncol=n.elem.d)
    n.index.d=rep(0,n.obs)
    for (i in 1:n.obs) {
        use=which( otu.table[i,]>0 )
        n.use=length(use)
        k.1=use[rep( 1:(n.use-1), times=(n.use-1):1 )]
        k.2=use[unlist(sapply(2:n.use, FUN = function(x){x:n.use}))]
        k=n.taxa*(k.1-1) + k.2 - k.1*(k.1 + 1)/2
        ln.pr=lchoose( n.i[i] - otu.table[i,k.1] - otu.table[i,k.2], rarefy.depth ) - l.denom[i]
        n.index.d[i]=n.use*(n.use-1)/2
        index.d[i,1:n.index.d[i]]=k
        d[i,1:n.index.d[i]]=exp(ln.pr) + mu[i,k.1] + mu[i,k.2] - 1
        #		
        cs.1=colSums( t(mu[,k.1])*d[i,1:n.index.d[i]] )
        cs.2=colSums( t(mu[,k.2])*d[i,1:n.index.d[i]] )
        cs.d=sum( d[i,1:n.index.d[i]] )
        #		
        psi.11[i,]=cs.d -   cs.1 -   cs.2
        psi.12[i,]=cs.d -   cs.1 - 2*cs.2
        psi.21[i,]=cs.d - 2*cs.1 -   cs.2
        psi.22[i,]=cs.d - 2*cs.1 - 2*cs.2
        #		
        
        # psi.11[i,]=colSums( t(1-mu[,k.1]-mu[,k.2]) * d[i,k] )
        # psi.12[i,]=colSums( t(1-mu[,k.1]-2*mu[,k.2]) * d[i,k] )
        # psi.21[i,]=colSums( t(1-2*mu[,k.1]-mu[,k.2]) * d[i,k] )
        # psi.22[i,]=colSums( t(1-2*mu[,k.1]-2*mu[,k.2]) * d[i,k] )
    }
    #
    #	calculate 2nd order approximation to average jaccard 
    #	
    jac.ave.1=jac.ave.2=matrix(0, nrow=n.obs, ncol=n.obs)
    jac.ave.sq.1=jac.ave.sq.2=matrix(0, nrow=n.obs, ncol=n.obs)
    jac.sq.ave=u22=u11=u12=matrix(0, nrow=n.obs, ncol=n.obs)
    for (i in 1:(n.obs-1)) {
        for (j in (i+1):n.obs) {
            #			use=which( otu.table[i,]*otu.table[j,]>0 )
            use=intersect( index.mu[i,1:n.index.mu[i]], index.mu[j,1:n.index.mu[j]] )
            dot.product.mu=sum( mu[i,use]*mu[j,use] )
            # 			use=intersect( index.d[i,1:n.index.d[i]], index.d[j, 1:n.index.d[j]] )
            use.i=which( index.d[i,1:n.index.d[i]] %in% index.d[j,1:n.index.d[j]] )
            use.j=which( index.d[j,1:n.index.d[j]] %in% index.d[i,1:n.index.d[i]] )
            dot.product.d=sum( d[i,use.i]*d[j,use.j] )
            e.1=mu.dot[i]+mu.dot[j]-dot.product.mu
            e.2=e.1 - dot.product.mu
            u.11=mu.dot[i] + mu.dot[j] + 2*mu.dot[i]*mu.dot[j] - 3*dot.product.mu + 2*(psi.11[i,j]+psi.11[j,i]) + 2*dot.product.d
            u.12=mu.dot[i] + mu.dot[j] + 2*mu.dot[i]*mu.dot[j] - 4*dot.product.mu + (psi.12[i,j]+psi.12[j,i]+psi.21[i,j]+psi.21[j,i]) + 4*dot.product.d
            u.22=mu.dot[i] + mu.dot[j] + 2*mu.dot[i]*mu.dot[j] - 4*dot.product.mu + 2*(psi.22[i,j]+psi.22[j,i]) + 8*dot.product.d
            jac.ave.1[i,j]=e.2/e.1
            jac.ave.2[i,j]=jac.ave.1[i,j] + (u.11*e.2 - u.12*e.1)/(e.1^3)
            jac.ave.sq.1[i,j]=jac.ave.1[i,j]^2
            jac.ave.sq.2[i,j]=jac.ave.sq.1[i,j] + (u.22*e.1^2-4*u.12*e.1*e.2+3*u.11*e.2^2)/(e.1^4)
            jac.sq.ave[i,j]=u.22/u.11
            u11[i,j]=u.11
            u22[i,j]=u.22
            u12[i,j]=u.12
        }
    }	
    jac.ave.1=jac.ave.1 + t(jac.ave.1)
    jac.ave.2=jac.ave.2 + t(jac.ave.2)
    jac.ave.sq.1=jac.ave.sq.1 + t(jac.ave.sq.1)
    jac.ave.sq.2=jac.ave.sq.2 + t(jac.ave.sq.2)
    jac.sq.ave=jac.sq.ave + t(jac.sq.ave)	
    u11=u11 + t(u11)
    u22=u22 + t(u22)
    u12=u12 + t(u12)
    
    res=list( jac.ave.1=jac.ave.1, 
              jac.ave.2=jac.ave.2, 
              jac.ave.sq.1=jac.ave.sq.1, 
              jac.ave.sq.2=jac.ave.sq.2)
    return(res)
    
} #jaccard.mean.fast.small

unifrac.ave.sq.fast.small = function( otu.table, tree, rarefy.depth=min(rowSums(otu.table)), trim=FALSE, keep.root=TRUE, n.batch=1) {

    # requireNamespace('ips')
    # requireNamespace('castor')
    # requireNamespace('phangorn')
    
    otu.table=as.matrix(otu.table)
    n.obs=dim(otu.table)[1]
    n.taxa=dim(otu.table)[2]                
    if (trim==TRUE) {	
        #
        #		remove zero columns from OTU table
        #		
        keep.list=which( colSums(otu.table)>0 )
        otu.table=otu.table[,keep.list]
        print( c('note: empty taxa removed, number of taxa changed from',n.taxa,'to',length(keep.list) ) )
        n.taxa=dim(otu.table)[2]
    }
    #	
    #
    #	harmonize lists of taxa from tree, OTU table and then sort OTU table to agree with ordering of tips in tree
    #
    tree.names=tree$tip.label
    otu.names=colnames(otu.table)
    n.tree.names=length(tree.names)
    n.otu.names=length(otu.names)
    if (n.otu.names<n.tree.names) {
        #
        #		remove tree nodes that do not correspond to taxa in otu.table
        #	
        drop.list=which( tree.names %notin% otu.names )
        tree=castor::get_subtree_with_tips(tree,omit_tips=drop.list, force_keep_root=keep.root)$subtree
        tree.names=tree$tip.label
        print( c('note: number of nodes in tree reduced from', n.tree.names, 'to', length(tree.names) ), quote=FALSE)
    }
    if( !all(sort(tree.names)==sort(otu.names)) ) {
        #
        #		remove otus that are not in tree
        #	
        keep.list=which( otu.names %in% tree.names )
        otu.table=otu.table[,keep.list]
        n.taxa=dim(otu.table)[2]
        print( c('note: number of taxa in OTU table reduced from', n.otu.names, 'to', n.taxa), quote=FALSE )
    }
    otu.table=otu.table[,tree.names]
    #
    #	set up branch lengths
    #	
    edge.2=tree$edge[,2]
    n.internal=max(edge.2) - n.taxa - 1
    n.branch=length(edge.2)
    n.node=n.taxa + n.internal
    n.elem=n.node*(n.node-1)/2
    edge.2=tree$edge[,2]
    tips=which( edge.2<=n.taxa )
    nodes=which( edge.2>n.taxa+1)
    #	edge.2[tips]=edge.2[tips]+n.taxa-2
    edge.2[tips]=edge.2[tips]+n.internal
    edge.2[nodes]=edge.2[nodes]-n.taxa-1
    ord=order(edge.2)
    edge.2=edge.2[ord]
    br.length=tree$edge.length[ord]
    #
    #	set up list of descendants for each node (nodes  1:(n.internal) are internal nodes; nodes (n.internal):(n.internal+n.taxa-1) are tips)
    #	
    #	descendants=allDescendants(tree)
    descendants=phangorn::Descendants(tree,node=c((n.taxa+2):(n.taxa+1+n.internal), 1:n.taxa), type='tips')
    # k.1=rep( 1:(n.node-1), times=(n.node-1):1 )
    # k.2=unlist(sapply(2:n.node, FUN = function(x){x:n.node}))
    # k=n.node*(k.1-1) + k.2 - k.1*(k.1+1)/2
    ancestry=matrix(FALSE,nrow=n.taxa, ncol=n.node)
    for (k in 1:n.node) {
        ancestry[ descendants[[k]], k]=TRUE
    }
    
    #
    #	set up list of tips and internal nodes of ancestors for each observation; calculate mu, mu.bar and mu.bar.2
    #	
    tips.i=list()
    #	nodes.i=list()
    all.nodes.i=list()
    mu=mu.bar=mu.bar.2=matrix(0,nrow=n.obs, ncol=n.node)
    mu.dot.1=rep(0,n.obs)
    mu.dot.2=rep(0,n.obs)
    n.i=rowSums( otu.table )
    l.denom=lchoose(n.i, rarefy.depth)
    n.nodes.max=0
    for (i in 1:n.obs) {
        tips=which( otu.table[i,]>0 )
        nodes=sort( unique( unlist( phangorn::Ancestors(x=tree, node=tips, type='all') ) ) ) - n.taxa - 1 
        nodes=nodes[nodes>0]
        all.nodes=union(nodes,tips+n.internal)
        tips.i=append(tips.i, list(tips) )
        all.nodes.i=append(all.nodes.i, list( all.nodes ))
        n.nodes=length(all.nodes)
        n.nodes.max=max(n.nodes.max, n.nodes)
        #		use.k=lapply( descendants[all.nodes], FUN=function(x) intersect(x,tips) )
        #		c.sum.k=unlist( lapply( use.k, FUN=function(x) sum( otu.table[i,x] ) ) )
        c.sum.k=colSums( otu.table[i,tips]*ancestry[tips,all.nodes] )
        mu.bar[i,all.nodes]=exp( lchoose( n.i[i]-c.sum.k, rarefy.depth ) - l.denom[i] )
        mu[i,all.nodes]=1-mu.bar[i,all.nodes]
        mu.bar.2[i,all.nodes]=1-mu.bar[i,all.nodes]
        mu.dot.1[i]=sum( mu[i,all.nodes]*br.length[all.nodes] )
        mu.dot.2[i]=sum( mu[i,all.nodes]*br.length[all.nodes]^2 )
    }
    #
    #   calculate d, psi.11, psi.12 and psi.22
    #	
    n.max.elem=n.nodes.max*(n.nodes.max-1)/2
    #	d=matrix(0,nrow=n.obs, ncol=n.max.elem)
    d=list()
    #	index.d=matrix(0, nrow=n.obs, ncol=n.max.elem)
    index.d=list()
    n.index.d=rep(0,n.obs)
    psi.11=psi.12=psi.21=matrix(0, nrow=n.obs, ncol=n.obs)
    psi.22=matrix(0, nrow=n.obs, ncol=n.obs)
    
    
    for (i in 1:n.obs) {
        all.nodes=all.nodes.i[[i]]
        tips=tips.i[[i]]
        n.nodes=length(all.nodes)
        k.1=all.nodes[rep( 1:(n.nodes-1), times=(n.nodes-1):1 )]
        k.2=all.nodes[unlist(sapply(2:n.nodes, FUN = function(x){x:n.nodes}))]
        k=n.node*(k.1-1) + k.2 - k.1*(k.1+1)/2
        
        #		c.sum.k1.k2=mapply( descendants[k.1], descendants[k.2], FUN=function(x,y) sum( otu.table[i,intersect( union(x,y),tips )]), SIMPLIFY=TRUE)  #this is slower than loop that follows
        #		n.k=length(k)
        #		c.sum.k1.k2=rep(0, n.k)
        #		for (j in 1:n.k) {
        #			c.sum.k1.k2[j]=sum( otu.table[i, intersect( union( descendants[[k.1[j]]], descendants[[k.2[j]]] ), tips) ] )
        #			}
        
        if (n.batch>1) {
            n.tips=length(tips)
            use.batch=batch( n.tips, n.batch )
            n.batch.use=length(use.batch)
            tips.use=tips[ use.batch[[1]] ]
            c.sum.k1.k2=colSums( otu.table[i, tips.use]*(ancestry[tips.use, k.1, drop=FALSE] | ancestry[tips.use, k.2, drop=FALSE]) )
            for (nb in 2:n.batch.use) {
                tips.use=tips[ use.batch[[nb]] ]
                c.sum.k1.k2=c.sum.k1.k2+colSums( otu.table[i, tips.use]*(ancestry[tips.use, k.1, drop=FALSE] | ancestry[tips.use, k.2, drop=FALSE]) )
            }
        }	
        else {	
            c.sum.k1.k2=colSums( otu.table[i,tips]*(ancestry[tips,k.1] | ancestry[tips,k.2]) )
        }
        
        n.index.d[i]=n.nodes*(n.nodes-1)/2
        #		index.d[i,1:n.index.d[i]]=k
        index.d=c(index.d,list(k))
        
        #		d[i,1:n.index.d[i]]=exp( lchoose( n.i[i]-c.sum.k1.k2, rarefy.depth ) - l.denom[i] ) + mu[i,k.1] + mu[i,k.2] - 1
        d=c(d, list(exp( lchoose( n.i[i]-c.sum.k1.k2, rarefy.depth ) - l.denom[i] ) + mu[i,k.1] + mu[i,k.2] - 1 ) )
        #		
        #		cs.1=colSums( t(mu[,k.1])*d[i,1:n.index.d[i]]*br.length[k.1]*br.length[k.2] )
        #		cs.2=colSums( t(mu[,k.2])*d[i,1:n.index.d[i]]*br.length[k.1]*br.length[k.2] )
        #		cs.d=sum( d[i,1:n.index.d[i]]*br.length[k.1]*br.length[k.2] )
        cs.1=colSums( t(mu[,k.1])*d[[i]]*br.length[k.1]*br.length[k.2] )
        cs.2=colSums( t(mu[,k.2])*d[[i]]*br.length[k.1]*br.length[k.2] )
        cs.d=sum( d[[i]]*br.length[k.1]*br.length[k.2] )
        #		
        psi.11[i,]=cs.d -   cs.1 -   cs.2
        psi.12[i,]=cs.d -   cs.1 - 2*cs.2
        psi.21[i,]=cs.d - 2*cs.1 -   cs.2
        psi.22[i,]=cs.d - 2*cs.1 - 2*cs.2
        #		
    }	
    #
    #	calculate 2nd order approximation to average unifrac
    #	
    unifrac.ave.1=unifrac.ave.2=matrix(0, nrow=n.obs, ncol=n.obs)
    unifrac.ave.sq.1=unifrac.ave.sq.2=unifrac.sq.ave=matrix(0, nrow=n.obs, ncol=n.obs)
    u11=u12=u22=e1=e2=matrix(0,nrow=n.obs, ncol=n.obs)
    #
    k.index=rep( 1:(n.node-1), times=(n.node-1):1 )
    k1.index=unlist(sapply(2:n.node, FUN = function(x){x:n.node}))
    #
    for (i in 1:(n.obs-1)) {
        for (j in (i+1):n.obs) {
            use=which( mu[i,]*mu[j,]>0 )
            dot.product.mu.1=sum( mu[i,use]*mu[j,use]*br.length[use] )
            dot.product.mu.2=sum( mu[i,use]*mu[j,use]*br.length[use]^2 )
            #			use=intersect( index.d[i,1:n.index.d[i]], index.d[j,1:n.index.d[j]] )
            #			use.i=which( index.d[i,1:n.index.d[i]] %in% use )
            #			use.j=which( index.d[j,1:n.index.d[j]] %in% use )
            use=intersect( index.d[[i]], index.d[[j]] )
            use.i=which( index.d[[i]] %in% use )
            use.j=which( index.d[[j]] %in% use )
            k.node=k.index[use]
            k1.node=k1.index[use]
            # 			dot.product.d=sum( d[i,use.i]*d[j,use.j]*br.length[k.node]*br.length[k1.node] )
            dot.product.d=sum( (d[[i]][use.i])*(d[[j]][use.j])*br.length[k.node]*br.length[k1.node] )
            e.1=mu.dot.1[i]+mu.dot.1[j]-dot.product.mu.1
            e.2=e.1 - dot.product.mu.1
            u.11=mu.dot.2[i] + mu.dot.2[j] + 2*mu.dot.1[i]*mu.dot.1[j] - 3*dot.product.mu.2 + 2*(psi.11[i,j]+psi.11[j,i]) + 2*dot.product.d
            u.12=mu.dot.2[i] + mu.dot.2[j] + 2*mu.dot.1[i]*mu.dot.1[j] - 4*dot.product.mu.2 + (psi.12[i,j]+psi.12[j,i]+psi.21[i,j]+psi.21[j,i]) + 4*dot.product.d
            u.22=mu.dot.2[i] + mu.dot.2[j] + 2*mu.dot.1[i]*mu.dot.1[j] - 4*dot.product.mu.2 + 2*(psi.22[i,j]+psi.22[j,i]) + 8*dot.product.d
            unifrac.ave.1[i,j]=e.2/e.1
            unifrac.ave.2[i,j]=unifrac.ave.1[i,j] + (u.11*e.2 - u.12*e.1)/(e.1^3)
            unifrac.ave.sq.1[i,j]=unifrac.ave.1[i,j]^2
            unifrac.ave.sq.2[i,j]=unifrac.ave.sq.1[i,j] + (u.22*e.1^2-4*u.12*e.1*e.2+3*u.11*e.2^2)/(e.1^4)
            unifrac.sq.ave[i,j]=u.22/u.11
            u11[i,j]=u.11
            u12[i,j]=u.12
            u22[i,j]=u.22
            e1[i,j]=e.1
            e2[i,j]=e.2
        }
    }	
    unifrac.ave.1=unifrac.ave.1 + t(unifrac.ave.1)
    unifrac.ave.2=unifrac.ave.2 + t(unifrac.ave.2)
    unifrac.ave.sq.1=unifrac.ave.sq.1 + t(unifrac.ave.sq.1)
    unifrac.ave.sq.2=unifrac.ave.sq.2 + t(unifrac.ave.sq.2)
    unifrac.sq.ave=unifrac.sq.ave + t(unifrac.sq.ave)	
    
    u11=u11 + t(u11)
    u12=u12 + t(u12)
    u22=u22 + t(u22)
    e1=e1 + t(e1)
    e2=e2 + t(e2)
    
    res=list( unifrac.ave.1=unifrac.ave.1, 
              unifrac.ave.2=unifrac.ave.2, 
              unifrac.ave.sq.1=unifrac.ave.sq.1, 
              unifrac.ave.sq.2=unifrac.ave.sq.2, 
              unifrac.sq.ave=unifrac.sq.ave)
    return(res)
    
} # unifrac.ave.sq.fast.small



batch = function( n, n.batch ) {
    max.num=floor( n/n.batch )
    x=seq(n)
    levels=floor(x/max.num)
    x.split=split(x, factor(levels))
    return(x.split)
}


#####################################
# Multi-Med
####################################

#############################################################################
###
###  medTest.SBMH (Mediator Test based on Bogomolov & Heller)
###
#############################################################################

###INPUT: 
###pEM:     a vector of size m (where m = number of mediators). Entries are the p-values for the E,M_j relationship 
###pMY:     a vector of size m (where m = number of mediators). Entries are the p-values for the M_j,Y|E relationship
###MCP.type:    multiple comparison procedure - either "FWER" or "FDR"
###t1:   threshold for determining the cutoff to be one of the top S_1 E/M_j relationships 
###t2:   threshold for determining the cutoff to be one of the top S_2 M_j/Y relationships 
###adaptive:  FALSE/TRUE depending on whether an adaptive threshold should be used
###
###OUTPUT:
###m x 1 matrix - either p-values (if MCP.type = "FWER") or q-values (if MCP.type = "FDR")

medTest.SBMH <- function(pEM,pMY,MCP.type="FDR",t1=0.05,t2=0.05,lambda=0){
    if (MCP.type=="FDR")   possVal <- proc.intersection.adaptiveFDR( pEM,pMY, t1=t1, t2=t2,lambda=lambda)$r.value
    possVal <- ifelse(is.na(possVal),1,possVal)
    ##threshold values at 1
    ifelse(possVal > 1, 1, possVal)
}


##FDR approach for medTest.SBMH
proc.intersection.adaptiveFDR = function(pv1,pv2, t1=0.025, t2=0.025,lambda=0.0){
    #if adaptive=TRUE, estimates the fraction of zeros of study i in selection j, for (i,j) = (1,2) or (2,1)
    #incorporate it into procedure by multiplying the p-value of study i by this fraction.   
    #print(c(t1,t2))
    #print(summary(pv1))
    #print(summary(pv2))
    #print(lambda)
    if (lambda>0){
        t1 = min(t1,lambda)
        t2= min(t2,lambda)
    }
    selected =which((pv1<=t1) & (pv2<=t2))
    R1 = sum(pv1<=t1)
    R2 = sum(pv2<=t2)
    if(length(selected)==0){
        return(list(r.value = rep(NA,length(pv1)), R1=R1, R2=R2, selected = selected))
    }
    
    S1 = (pv1<=t1) #the selected from study 1
    pi2 = (1+sum(pv2[S1]>lambda))/(R1*(1-lambda)) #the fraction of nulls in study 2 among S1
    
    S2 = (pv2<=t2) #the selected from study 2
    pi1 = (1+sum(pv1[S2]>lambda))/(R2*(1-lambda)) #the fraction of nulls in study 1 among S2
    
    if (lambda==0){ #ie nonadaptive
        pi2=1
        pi1=1
    }
    
    if((pi1 > 1) | (pi2 > 1))
    {
        warning("At least one of the estimated fraction of nulls is > 1: Using lambda=0 to set them to 1.")
        pi2=1
        pi1=1
    }
    
    Z.selected <- pmax(pi1*R2*pv1/0.5, pi2*R1*pv2/0.5)[selected]
    oz <- order(Z.selected, decreasing =TRUE)
    ozr <- order(oz)
    r.value.selected <- cummin( (Z.selected/rank(Z.selected, ties.method= "max"))[oz] )[ozr]
    r.value.selected <- pmin(r.value.selected,1)
    r.value  = rep(NA,length(pv1))
    r.value[selected] = r.value.selected
    return(list(r.value = r.value, R1=R1, R2=R2, selected = selected, pi1=pi1, pi2=pi2))
}





