---
name: Distributions
topic: Probability Distributions
maintainer: Christophe Dutang, Patrice Kiener, Bruce J. Swihart
email: dutangc@gmail.com
version: 2023-08-09
source: https://github.com/cran-task-views/Distributions/
---


For most of the classical distributions, base R provides probability
distribution functions (p), density functions (d), quantile functions
(q), and random number generation (r). Beyond this basic functionality,
many CRAN packages provide additional useful distributions. In
particular, multivariate distributions as well as copulas are available
in contributed packages.

Ultimate bibles on probability distributions are:

-   Different volumes of N. L. Johnson, S. Kotz and N. Balakrishnan
    books, e.g. Continuous Univariate Distributions, Vol. 1,
-   Thesaurus of univariate discrete probability distributions by G.
    Wimmer and G. Altmann.
-   Statistical Distributions by M. Evans, N. Hastings, B. Peacock.
-   Distributional Analysis with L-moment Statistics using the R
    Environment for Statistical Computing, Asquith (2011).

The maintainers gratefully acknowledge Achim Zeileis, David Luethi,
Tobias Verbeke, Robin Hankin, Mathias Kohl, G. Jay Kerns, Kjetil
Halvorsen, William Asquith for their useful comments/suggestions. If you
think information is not accurate or not complete, please send an e-mail
to the maintainer or submit an issue or pull request in the GitHub
repository linked above.

### Table of contents
- [Base functionality](#Base) 

- [Discrete distributions](#Discrete) 
  - [Univariate Discrete](#UnivariateDiscrete) 
  - [Multivariate Discrete](#MultivariateDiscrete) 
  
- [Continuous distributions](#Continuous) 
  - [Univariate Continuous](#UnivariateContinuous) 
  - [Multivariate Continuous](#MultivariateContinuous) 

- [Other distributions](#Other)
  - [Mixed-type distributions](#MixedType)
  - [Mixture of probability laws](#Mixture)
  - [Random matrices](#Matrix) 
  - [Copulas](#Copulas) 
  - [Compound, composite, discretized, exponentiated and transformation of distributions](#Transform)

- [Moments, skewness, kurtosis and etc](#Moments) 
- [Random number generators (RNG)](#Random) 
- [Miscellaneous](#Misc) 


# [Base functionality:]{#Base}


-   Base R provides probability distribution functions `p`*foo*`()`
    density functions `d`*foo*`()`, quantile functions `q`*foo*`()`,
    and random number generation `r`*foo*`()` where *foo* indicates
    the type of distribution: beta (*foo* = `beta`), binomial `binom`,
    Cauchy `cauchy`, chi-squared `chisq`, exponential `exp`, Fisher F
    `f`, gamma `gamma`, geometric `geom`, hypergeometric `hyper`,
    logistic `logis`, lognormal `lnorm`, negative binomial `nbinom`,
    normal `norm`, Poisson `pois`, Student t `t`, uniform `unif`,
    Weibull `weibull`. Following the same naming scheme, but somewhat
    less standard are the following distributions in base R:
    probabilities of coincidences (also known as "birthday paradox")
    `birthday` (only p and q), studentized range distribution `tukey`
    (only p and q), Wilcoxon signed rank distribution `signrank`,
    Wilcoxon rank sum distribution `wilcox`.
-   Base R provides various one-sample or two-sample tests for univariate
    distributions, e.g., `ks.test`, `shapiro.test`, `ansari.test`, `chisq.test`,
    `poisson.test`. `r pkg("Ecume")` provides non-parametric two-sample (or k-sample) 
    distribution comparisons in the univariate or multivariate case allowing
    observation weights and thresholds.

-   Probability generating function: no longer implemented.

Some packages may optionally provide the symbolic derivatives with respect 
to the parameters for the probability functions. 
For instance, the first and second derivatives of the log-density can be 
of some help in estimation and inference tasks, and the derivatives of 
the quantile function can help when inferring on a given quantile.
For that purpose, the following base R functions can be used
`stats::D()` for derivatives w.r.t. a single parameter, 
or `stats::deriv()` for (partial) derivatives w.r.t. multiple parameters.
The `r pkg("Deriv")` package provides a much more flexible 
symbolic differentiation interface.
One can also use Stan Math library through `r pkg("StanHeaders")` package,
see e.g.
[this blog](https://www.jchau.org/2022/01/24/automatic-differentiation-in-r-with-stan-math/).
The `r pkg("nieve")` package provides symbolic differentiation for
two probability distribution (Generalized Pareto and Generalized
Extreme Value) in order to compute the log-likelihood for example.


# [Discrete distributions:]{#Discrete}

## [Discrete univariate distributions:]{#UnivariateDiscrete}


-   *Beta-binomial distribution:* provided in
    `r pkg("VGAM", priority = "core")`,
    `r pkg("extraDistr")`, `r pkg("rmutil")`,
    `r pkg("emdbook")`. ZI/ZM beta binomial distributions
    are implemented in
    `r pkg("gamlss.dist", priority = "core")`.
-   *Beta-geometric distribution:* provided in
    `r pkg("VGAM")`.
-   *Binomial (including Bernoulli) distribution:* provided in **stats**
    . Zero-modified, zero-inflated, truncated versions are provided in
    `r pkg("gamlss.dist")`,
    `r pkg("extraDistr")`,
    `r pkg("actuar", priority = "core")` and in
    `r pkg("VGAM")`. `r pkg("LaplacesDemon")`
    provides dedicated functions for the Bernoulli distribution.
    `r pkg("rmutil")` provides the double binomial and the
    multiplicative binomial distributions.
    
      ---------------------- ------------- ------------------ -----------------------
      *Distribution name*    *Packages*    *Functions*          *Distribution suffix*
      binomial               stats         `d`, `p`, `q`, `r`   `binom`
      zero-infl. binomial    extraDistr    `d`, `p`, `q`, `r`   `zib`
      zero-infl. binomial    VGAM          `d`, `p`, `q`, `r`   `zibinom`
      zero-infl. binomial    gamlss.dist   `d`, `p`, `q`, `r`   `ZIBI`
      zero mod. binomial     VGAM          `d`, `p`, `q`, `r`   `zabinom`
      zero mod. binomial     actuar        `d`, `p`, `q`, `r`   `zmbinom`
      zero mod. binomial     gamlss.dist   `d`, `p`, `q`, `r`   `ZABI`
      zero trunc. binomial   actuar        `d`, `p`, `q`, `r`   `ztbinom`
      trunc. binomial        extraDistr    `d`, `p`, `q`, `r`   `tbinom`
      ---------------------- ------------- ------------------ -----------------------

      :  Summary for Binomial-related distributions

-   *Bell Touchard distribution:* standard and zero-inflated
    provided in `r pkg("countDM")`.
-   *Benford distribution:* provided in `r pkg("VGAM")` and
    `r pkg("BenfordTests")`.
-   *Bernoulli distribution:* provided in
    `r pkg("extraDistr")`.
-   *Borel-Tanner distribution:* provided in
    `r pkg("VGAM")`.
-   *Delaporte distribution:* provided in
    `r pkg("gamlss.dist")` and
    `r pkg("Delaporte")`.
-   *Dirac distribution:* provided in
    `r pkg("distr", priority = "core")`.
-   *Discrete categorical distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Discrete exponential distribution:* provided in
    `r pkg("poweRlaw")`.
-   *Discrete gamma distribution:* provided in
    `r pkg("extraDistr")`.
-   *Discrete inverse Weibull distribution:*
    `r pkg("DiscreteInverseWeibull")` provides d, p, q, r
    functions for the inverse Weibull as well as hazard rate function
    and moments.
-   *Discrete Laplace distribution:* The discrete Laplace distribution
    is provided in `r pkg("extraDistr")` (d, p, r). The skew
    discrete Laplace distribution has two parametrization (DSL and
    ADSL), both provided in `r pkg("DiscreteLaplace")` and
    DSL in `r pkg("disclap")`.
    `r pkg("LaplacesDemon")` also provides the DSL
    parametrization only.
-   *Discrete lognormal distribution:* provided in
    `r pkg("poweRlaw")`.
-   *Discrete normal distribution:* provided in
    `r pkg("extraDistr")`.
-   *Discrete power law distribution:* provided in
    `r pkg("poweRlaw")`.
-   *Discrete uniform distribution:* can be easily obtained with the
    functions `sum,cumsum,sample` and is provided in
    `r pkg("extraDistr")`.
-   *Discrete Weibull distribution:* provided in
    `r pkg("DiscreteWeibull")`: d, p, q, r, m for disc.
    Weib. type 1, d, p, q, r, m, h for disc. Weib. type 3.
    `r pkg("extraDistr")` provides d, p, q, r for Type 1.
-   *Felix distribution:* provided in `r pkg("VGAM")`.
-   *gamma count distribution:* provided in
    `r pkg("rmutil")`.
-   *Geometric distribution:* provided in **stats** . Zero-modified,
    zero-inflated, truncated versions are provided in
    `r pkg("gamlss.dist")`, `r pkg("actuar")`
    and in `r pkg("VGAM")`. The time-varying geometric is
    provided in `r pkg("tvgeom")`.
-   *Geometric (compound) Poisson distribution (also known Polya-Aeppli
    distribution):* provided in `r pkg("polyaAeppli")`.
-   *Generalized/fractional binomial distribution:* 
    `r pkg("GenBinomApps")` provides the generalized binomial distribution.
    `r pkg("frbinom")` provides the fractional binomial distribution
    where trials are from a generlized Bernoulli process.
-   *Generalized Hermite distribution:* provided in
    `r pkg("hermite")`.
-   *Good distribution:* provided in
    `r pkg("good")`.
-   *Hypergeometric distribution:* provided in **stats** . Non-central
    hypergeometric distribution is provided in
    `r pkg("MCMCpack")` (d,r). Extended hypergeometric
    distribution can be found in `r pkg("BiasedUrn")`
    package, which provides not only p, d, q, r functions but also mean,
    variance, mode functions. Generalized hypergeometric distribution is
    implemented in `r pkg("SuppDists")`. Negative
    hypergeometric distribution is provided in
    `r pkg("tolerance")`, `r pkg("extraDistr")`.
-   *Lagrangian Poisson distribution:*
    `r pkg("RMKdiscrete")` provides d, p, q, r functions for
    the univariate and the bivariate Lagrangian Poisson distribution.
-   *Lindley's power series distribution:* provided in
    `r pkg("LindleyPowerSeries")`.
-   *Logarithmic distribution:* This can be found in
    `r pkg("extraDistr")`, `r pkg("VGAM")`,
    `r pkg("actuar")`, and
    `r pkg("gamlss.dist")`. Zero-modified and zero-truncated
    versions is provided in `r pkg("actuar")`. A fast random
    generator is available for the logarithmic distribution is
    implemented in `r pkg("Runuran")` as well as the
    'density' function.
-   *Poisson distribution:* provided in **stats** and in
    `r pkg("poweRlaw")`. Zero-modified, zero-inflated,
    truncated versions are provided in
    `r pkg("extraDistr")`,
    `r pkg("gamlss.dist")`, `r pkg("actuar")`
    and in `r pkg("VGAM")`.
    `r pkg("extraDistr")` provides the truncated Poisson
    distribution. `r pkg("LaplacesDemon")` provides the
    generalized Poisson distribution. `r pkg("rmutil")`
    provides the double Poisson, the multiplicative Poisson and the
    Power variance function Poisson distributions.
    `r pkg("poibin")` and
    `r pkg("PoissonBinomial")` provide the Poisson binomial
    distribution. See the mixture section such as the Poisson-lognormal
    mixture.
-   *Poisson-Lindley distribution:* provided in
    `r pkg("tolerance")`.
-   *Power law distribution:* provided in
    `r pkg("poweRlaw")`.
-   *Mana Clash distribution:* provided in
    `r pkg("RMKdiscrete")`.
-   *Negative binomial distribution:* provided in **stats** .
    Zero-modified, zero-inflated, truncated versions are provided in
    `r pkg("gamlss.dist")`,
    `r pkg("extraDistr")`, `r pkg("emdbook")`,
    `r pkg("actuar")` and in `r pkg("VGAM")`.
    New parametrization of the negative binomial distribution is
    available in `r pkg("RMKdiscrete")`.
    `r pkg("nbconv")` provides p, q, r functions for convolutions
    of negative binomial distributions.
-   *Sichel distribution:* provided in
    `r pkg("gamlss.dist")`.
-   *Skellam distribution:* provided in
    `r pkg("extraDistr")`, `r pkg("VGAM")` and
    `r pkg("skellam")`.
-   *Waring distribution:* sampling in `r pkg("degreenet")`.
-   *Yule-Simon distribution:* provided in `r pkg("VGAM")`
    and sampling in `r pkg("degreenet")`.
-   *Zeta and Haight's Zeta distribution:* provided in
    `r pkg("VGAM")`, `r pkg("tolerance")`.
-   *Zipf distribution and extensions:* d, p, q, r functions of the Zipf
    and the Zipf-Mandelbrot distributions are provided in
    `r pkg("tolerance")`, `r pkg("VGAM")`.
    Package `r pkg("zipfR")` provides tools for distribution
    of word frequency, such as the Zipf distribution.
    `r pkg("zipfextR")` provides three extensions of the
    Zipf distribution: the Marshall-Olkin Extended Zipf, the
    Zipf-Poisson Extreme and the Zipf-Poisson Stopped Sum distributions.

## [Discrete multivariate distributions:]{#MultivariateDiscrete}


-   *Hyper Dirichlet distribution:* provided in
    `r pkg("hyper2")` package.
-   *Multinomial distribution:* stats, `r pkg("mc2d")`,
    `r pkg("extraDistr")` packages provide d, r functions. r
    is provided in `r pkg("MultiRNG")` and
    `r pkg("compositions")`. p function is provided by
    `r pkg("pmultinom")`.
-   *Multinomial Dirichlet distribution:* functions d, r are provided in
    `r pkg("MCMCpack")`, `r pkg("mc2d")`,
    `r pkg("dirmult")`, `r pkg("extraDistr")`
    and `r pkg("bayesm")`. r is provided in
    `r pkg("MultiRNG")`.
-   *Multivariate Ewens distribution:* not yet implemented?
-   *Multivariate geometric:* d, r functions provided in
    `r pkg("bivgeom")` for the bivariate geometric distribution. 
    `r pkg("BivGeo")` provides the Basu-Dhar bivariate geometric distribution.
-   *Multivariate hypergeometric distribution:* provided in
    `r pkg("extraDistr")`. The conditional hypergeometric distribution is
    provided (d, p, q, r) in `r pkg("chyper")`.
-   *Multivariate logarithmic distribution:* the bivariate
    logarithmic distribution is provided in
    `r pkg("trawl")`.    
-   *Multiplicative multinomial distribution:* The multiplicative multinomial
    distribution is implemented in `r pkg("MM")`.
-   *Multivariate negative distribution:* A bivariate distribution with
    negative-binomial marginals is available in
    `r pkg("RMKdiscrete")` and `r pkg("trawl")`. 
    `r pkg("MNB")` provides a generator and 
    diagnostic tool for multivariate negative binomial distribution.  
    `r pkg("bzinb")` provides a random generator for the bivariate negative
    binomial (classic and zero-inflated) distribution.
-   *Multivariate Poisson distribution:*
    `r pkg("compositions")` provides a random generator.
    `r pkg("bzinb")` provides a random generator for the bivariate Poisson 
    (classic and zero-inflated) distribution.
-   *Multivariate Poisson-lognormal:* the bivariate
    Poisson-lognormal distribution is provided in
    `r pkg("poilog")`.
-   *Multivariate Dirichlet (also known as Polya) distribution:* 
    functions d, r of the Dirichlet distribution are
    provided in `r pkg("extraDistr")`,
    `r pkg("LaplacesDemon")`, `r pkg("DirichletReg")` and
    `r pkg("Compositional")`.
-   *Truncated Stick-Breaking distribution:* provided in
    `r pkg("LaplacesDemon")`.
    
# [Continuous distributions:]{#Continuous}     

## [Continuous univariate distributions:]{#UnivariateContinuous}


-   *Arcsine distribution:* implemented in package
    `r pkg("distr")`.
-   *Beta distribution and its extensions:* Base R provides the d, p, q,
    r functions for this distribution (see above).
    `r pkg("extraDistr")` provides the beta distribution
    parametrized by the mean and the precision.
    `r pkg("actuar")` provides moments and limited expected
    values. `r pkg("sadists")` implements Gram Charlier,
    Edgeworth and Cornish-Fisher approximations for doubly non central
    beta distribution for computing d, p, q, r functions.
    `r pkg("extraDistr")` provides the four-parameter beta
    with lower and upper bounds. The generalized beta of the first kind
    (GB1) (exponentiation of beta 1) is provided in
    `r pkg("gamlss.dist")`, `r pkg("mbbefd")`,
    `r pkg("actuar")`. `r pkg("betafunctions")`
    provides the four-parameter beta (that is with location and scale
    parameters), the beta parametrized by the mean and the variance as
    well as the beta compound beta distribution. The beta prime (or beta
    of the second kind), which is the distribution of X/(1-X) when X
    follows a beta distribution of the first kind, is provided in
    `r pkg("VGAM")`, `r pkg("extraDistr")`,
    `r pkg("LaplacesDemon")` and
    `r pkg("mc2d")`. The zero and one inflated beta
    distribution can be found in `r pkg("gamlss.dist")`. The
    generalized beta of the second kind (GB2) is provided in
    `r pkg("gamlss.dist")`, `r pkg("GB2")`.
    Several special cases of the generalized beta distribution are also
    implemented in `r pkg("VGAM")`,
    `r pkg("mc2d")`: Lomax, inverse Lomax, Dagum,
    Singh-Maddala, Pert distributions. `r pkg("actuar")`
    provides the Feller-Pareto distribution as special cases Burr,
    loglogistic, paralogistic, generalized Pareto, Pareto, see also the
    Pareto subsection. `r pkg("llogistic")` provides the
    log-logistic parametrized by the median.
    
      ------------------------- -------------- -------------------  -----------------------
      *Distribution name*       *Packages*     *Functions*           *Distribution suffix*
      Beta (1st kind)           stats           d, p, q, r           `beta`
      Beta                      actuar          m, mgf, lev          `beta`
      Beta                      betafunctions   d, p, q, r           `Beta.4P`
      Doubly non central beta   sadists         d, p, q, r           `nbeta`
      4-param beta              extraDistr      d, p, q, r           `nsbeta`
      zero-infl beta            gamlss.dist     d, p, q, r           `BEZI`
      one-infl beta             gamlss.dist     d, p, q, r           `BEOI`
      one-infl beta             mbbefd          d, p, q, r, m, ec    `oibeta`
      GB1                       gamlss.dist     d, p, q, r           `GB1`
      GB1                       mbbefd          d, p, q, r, m, ec    `gbeta`
      GB1                       actuar          d, p, q, r, m, lev   `genbeta`
      one-infl GB1              mbbefd          d, p, q, r, m, ec    `oigbeta`
      ------------------------- --------------- -------------------- -----------------------

      :  Summary for Beta-related distributions

    
    
      --------------------- ------------------------  -------------------- ---------------------
      *Distribution name*   *Packages*                *Functions*          *Distribution suffix*
      Beta (2nd kind)       VGAM                      d, p, q, r           `beta`
      Beta (2nd kind)       extraDistr                d, p, q, r           `invbeta`
      Beta (2nd kind)       LaplacesDemon             d, r                 `betapr`
      GB2                   VGAM                      d, p, q, r           `genbetaII`
      GB2                   gamlss.dist               d, p, q, r           `GB2`
      GB2                   GB2                       d, p, q, r           `gb2`
      Trans beta 2          actuar                    d, p, q, r, m, lev   `trbeta`
      --------------------- ------------------------  -------------------- ---------------------

      :  Summary for Beta-2-related distributions

    
-   *Bell-G distribution:* `r pkg("BGFD")` provides d, p, q, r
    functions for Bell exponential, Bell extended exponential, Bell
    Weibull, Bell extended Weibull, Bell-Fisk, Bell-Lomax, Bell
    Burr-XII, Bell Burr-X, complementary Bell exponential,
    complementary Bell extended exponential, complementary Bell
    Weibull, complementary Bell extended Weibull, complementary
    Bell-Fisk, complementary Bell-Lomax, complementary Bell 
    Burr-XII and complementary Bell Burr-X distribution.  
    The package also provides hazard function and an estimation procedure.
-   *Benini distribution:* provided in `r pkg("VGAM")`.
-   *Bezier-Montenegro-Torres distribution:* provided in
    `r pkg("BMT")`.
-   *Bhattacharjee (normal+uniform) distribution:* provided in package
    `r pkg("extraDistr")`.
-   *Birnbaum-Saunders distribution:* provided in package
    `r pkg("VGAM")` and `r pkg("extraDistr")`.
-   *Bridge distribution:* provided in
    `r pkg("bridgedist")`, as detailed in Wang and Louis
    (2003). The distribution of random intercept that allows a
    marginalized random intercept logistic regression to also be
    logistic regression.
-   *Box Cox distribution:* `r pkg("gamlss.dist")` provides
    the Box-Cox normal, the Box-Cox power exponential and the Box-Cox t
    distributions. `r pkg("rmutil")` provides the Box-Cox
    normal.
-   *Burr distribution:* see Pareto.
-   *Cardioid distribution:* provided in `r pkg("VGAM")`
    (d,p,q,r) and `r pkg("CircStats")`,
    `r pkg("circular")` (d,r).
-   *Carthwrite's Power-of-Cosine distribution:* provided in
    `r pkg("circular")` (d,r).
-   *Cauchy distribution:* Base R provides the d, p, q, r functions for
    this distribution (see above). Other implementations are available
    in `r pkg("lmomco", priority = "core")` and
    `r pkg("sgt")`. The skew Cauchy distribution is provided
    in `r pkg("sn")`. `r pkg("LaplacesDemon")`
    provides d, p, q, r functions for the Half-Cauchy distribution. The
    wrapped Cauchy distribution is provided in
    `r pkg("CircStats")`.
-   *Chen distribution:* no longer implemented.
-   *Chernoff distribution:*  `r pkg("ChernoffDist")` provides
    d, p, q functions of the distribution of the maximizer of the 
    two-sided Brownian motion minus quadratic drift, known as Chernoff's
    distribution.
-   *Chi(-squared or not) distribution:* Base R provides the d, p, q, r
    functions for the chi-squared distribution, both central and
    non-central (see above). Moments, limited expected values and the
    moment generating function are provided in
    `r pkg("actuar")`. `r pkg("extraDistr")`
    provides d, p, q, r functions for inverse chi-squared distribution
    (standard and scaled). Only d,r functions are available for the
    inverse chi-squared distribution in package
    `r pkg("LaplacesDemon")`. A fast random generator is
    available for the Chi distribution is implemented in
    `r pkg("Runuran")` as well as the density function. The
    non-central Chi distribution is not yet implemented. The
    chi-bar-squared distribution is implemented in
    `r pkg("emdbook")`. `r pkg("sadists")`
    implements Gram Charlier, Edgeworth and Cornish-Fisher
    approximations for sums of non central chi-squared raised to powers
    distribution and sums of log of non central chi-squared for
    computing d, p, q, r functions.
    
      ---------------------------- -----------------------  ------------- -----------------------
      *Distribution name*          *Packages*               *Functions*   *Distribution suffix*
      Chi-squared                  stats                    d, p, q, r    `chisq`
      Chi-squared                  actuar                   m, mgf, lev   `chisq`
      Chi-squared                  Runuran                  d, r          `chisq`
      Chi-bar-squared              emdbook                  d, p, q, r    `chibarsq`
      Chi                          Runuran                  d, r          `chi`
      Inverse Chi-squared          extraDistr               d, p, q, r    `invchisq`
      Scaled Inverse Chi-squared   extraDistr               d, p, q, r    `invchisq`
      Sum of power Chi-squared     sadists                  d, p, q, r    `sumchisqpow`
      Sum of log Chi-squared       sadists                  d, p, q, r    `sumlogchisq`
      ---------------------------- -----------------------  ------------- -----------------------

      :  Summary for Chi-related distributions

    
-   *Circular distribution:* uniform circular provided in
    `r pkg("circular")` (d,r); Generalized von Mises
    circular provided in `r pkg("circular")` (d).
-   *Consul distribution:* see `r pkg("rmutil")`.
-   *Continuous binomial distribution:* `r pkg("cbinom")`
    provides the d/p/q/r functions for a continuous analog to the
    standard discrete binomial with continuous size parameter and
    continuous support with x in \[0, size + 1\].
-   *Dagum distribution:* see beta.
-   *Davies distribution:* The Davies distribution is provided in
    `r pkg("Davies")` package.
-   *(non-central) Dunnett's test distribution:* provided in
    `r pkg("nCDunnett")`.
-   *Eta-mu distribution:* provided in `r pkg("lmomco")`.
    `r pkg("sadists")` implements Gram Charlier, Edgeworth
    and Cornish-Fisher approximations for doubly non central eta
    distribution for computing d, p, q, r functions.
-   *Exponential distribution and its extensions:* Base R provides the
    d, p, q, r functions for this distribution (see above).
    `r pkg("actuar")` provides additional functions such as
    the moment generating function, moments and limited expected values.
    It also has the d, p, q, r for the inverse exponential distribution.
    The shifted (or two-parameter exponential) and the truncated
    exponential distributions are implemented in
    `r pkg("lmomco")` and `r pkg("tolerance")`
    packages with d, p, q, r functions. Exponential Power distribution
    is also known as General Error Distribution: d, p, q, r functions
    for the power and the skew power exponential type 1-4 distributions
    are implemented in `r pkg("gamlss.dist")` and
    `r pkg("lmomco")`. The power exponential distribution is
    also provided in `r pkg("normalp")`,
    `r pkg("rmutil")`, `r pkg("LaplacesDemon")`. 
    The skew power exponential is
    provided `r pkg("mixSPE")`. A fast
    random generator is available for the power Exponential distribution
    is implemented in `r pkg("Runuran")` as well as the
    density function.
    `r pkg("AEP")` implements the Asymmetric Exponential Power Distribution.
    
      --------------------------------------------------  ---------------------   ---------------------- ----------------------------
      *Distribution name*                                 *Packages*             *Functions*             *Distribution suffix*
      Exponential                                         stats                   d, p, q, r             `exp`
      Exponential                                         actuar                  m, mgf, lev            `exp`
      Exponential                                         gamlss.dist             d, p, q, r             `EXP`
      Exponential                                         poweRlaw                d, p, q, r             `exp`
      Inverse exponential                                 actuar                  d, p, q, r, m, lev     `invexp`
      Shifted exponential                                 lmomco                  d, p, q, r, lm, tlmr   `exp`
      Shifted exponential                                 tolerance               d, p, q, r             `2exp`
      Truncated exponential                               lmomco                  d, p, q, r, lm, tlmr   `texp`
      Truncated exponential                               ReIns                   d, p, q, r             `texp`
      Power exponential                                   normalp                 d, p, q, r             `normp`
      Power exponential                                   Runuran                 d, r                   `exp`
      Power exponential                                   rmutil                  d, r                   `powexp`
      Power exponential                                   LaplacesDemon           d, p, q, r             `pe`
      Skew power exp.                                     lmomco                  d, p, q, r, lm, tlmr   `aep4`
      Power and skew power exp.                           mixSPE                  r                      `pe, spe`
      Power and skew power exp.                           gamlss.dist             d, p, q, r             `PE, SEP`
      -------------------------------------------------- ----------------------   ---------------------- ----------------------------

      :  Summary for exponential-related distributions

    
-   *Externally studentized midrange distribution:* Package
    `r pkg("SMR")` computes the studentized midrange
    distribution (d, p, q, r).
-   *Fisher-Snedecor (or F) distribution:* Base R provides the d, p, q,
    r functions for the F distribution, possibly with a non-central
    parameter. `r pkg("sadists")` implements Gram Charlier,
    Edgeworth and Cornish-Fisher approximations for doubly non central
    Fisher distribution (and product of multiple doubly non central
    Fisher distribution) for computing d, p, q, r functions.
    `r pkg("flexsurv")` provides d, p, q, r functions as
    well as hazard (h) and integrated hazard rate (i) functions for the
    generalized F distribution. `r pkg("fpow")` returns the
    noncentrality parameter of the noncentral F distribution if
    probability of type I and type II error, degrees of freedom of the
    numerator and the denominator are given.
-   *Frechet distribution:* provided in `r pkg("VGAM")`,
    `r pkg("RTDE")`, `r pkg("ReIns")`,
    `r pkg("extraDistr")`,
    `r pkg("distributionsrd")` and
    `r pkg("evd")`. A fast random generator is available for
    the Frechet distribution is implemented in
    `r pkg("Runuran")` as well as the density function. The
    truncated Frechet distribution is provided in
    `r pkg("ReIns")`.
-   *Friedman's Chi distribution:* provided in
    `r pkg("SuppDists")`.
-   *Gamma distribution and its extensions:* Base R provides the d, p,
    q, r functions for this distribution (see above).
    `r pkg("EnvStats")` provides d, p, q, r functions of the
    gamma parametrized by the mean and the coefficient of variation.
    `r pkg("actuar")` provides d, p, q, r functions of the
    inverse, the inverse transformed and the log gamma distributions
    while `r pkg("ghyp")` provides those functions for the
    variance gamma distribution. `r pkg("extraDistr")` and
    `r pkg("LaplacesDemon")` provide the inverse gamma
    distribution. `r pkg("CaDENCE")` provides the
    zero-inflated gamma distribution.
    `r pkg("VarianceGamma")` provides d, p, q, r functions
    for the variance gamma distribution as well as moments (skewness,
    kurtosis, \...). `r pkg("VGAM")`,
    `r pkg("ggamma")` provide d, p, q, r functions of the
    log gamma and the generalized gamma distribution. The generalized
    gamma distribution can also be found in
    `r pkg("gamlss.dist")`. See Pearson III for a
    three-parameter gamma distribution with a location parameter.
    `r pkg("flexsurv")` provides d, p, q, r functions as
    well as hazard (h) and integrated hazard rate (i) functions for the
    generalized gamma distribution. `r pkg("coga")` provides
    d, p, r functions for a sum of independent but not identically
    distributed gamma distributions. `r pkg("MCMCpack")`
    provides d, r functions of the Inverse Gamma.
    `r pkg("rmutil")` provides the generalized Gamma.
    `r pkg("distTails")` provides the full-tail gamma
    distribution
    `r pkg("sglg")` provides the generalized log-Gamma along with
    various functions to fit semi-parametric regression models.
    `r pkg("ollggamma")` provides d, p, q, r for the Odd Log-Logistic Generalized Gamma. 
    
      ---------------------- ------------------   ------------------------- -----------------------
      *Distribution name*    *Packages*           *Functions*               *Distribution suffix*
      Gamma                  stats                d, p, q, r                `gamma`
      Gamma                  actuar               m, mgf, lev               `gamma`
      Gamma                  EnvStats             d, p, q, r                `gammaAlt`
      zero-inflated Gamma    CaDENCE              d, p, q, r                `bgamma`
      Inverse gamma          actuar               d, p, q, r, m, lev, mgf   `invgamma`
      Inverse gamma          extraDistr           d, p, q, r                `invgamma`
      Inverse gamma          LaplacesDemon        d, r                      `invgamma`
      Inverse gamma          MCMCpack             d, r                      `invgamma`
      Log-gamma              actuar               d, p, q, r, m, lev        `lgamma`
      Log-gamma              VGAM                 d, p, q, r                `lgamma`
      Variance gamma         ghyp                 d, p, q, r                `VG`
      Variance gamma         VarianceGamma        d, p, q, r, m             `vg`
      Generalized gamma      flexsurv             d, p, q, r, h, i          `gengamma`
      Generalized gamma      gamlss.dist          d, p, q, r                `GG`
      Generalized gamma      VGAM                 d, p, q, r                `gengamma.stacy`
      Generalized gamma      rmutil               d, p, q, r                `ggamma`
      Generalized gamma      ggamma               d, p, q, r                `ggamma`
      convolution of gamma   coga                 d, p, r                   `coga`
      Full-taill gamma       distTails            d, p, r                   `dFTG`
      Generalized log-gamma  sglg                 d, p, q, r                `glg`
      ---------------------- ------------------   ------------------------- -----------------------

      :  Summary for gamma-related distributions

    
-   *Gaussian (or normal) distribution and its extensions:* Base R
    provides the d, p, q, r functions for this distribution (see above).
    `r pkg("actuar")` provides the moment generating
    function and moments. The `r pkg("truncnorm")` package
    provides d, p, q, r functions for the truncated gaussian
    distribution as well as functions for the first two moments.
    `r pkg("EnvStats")`
    provides d, p, q, r functions for the truncated normal distribution
    and the zero-modified distribution.
    `r pkg("extraDistr")` provides the truncated normal.
    `r pkg("LaplacesDemon")` provides d, p, q, r functions
    for the Half-normal distribution. The wrapped normal distribution is
    provided in `r pkg("CircStats")`.
    `r pkg("lmomco")` implements the generalized normal
    distribution. The Exponentially modified Gaussian is available in
    `r pkg("emg")` and `r pkg("gamlss.dist")`
    `r pkg("sn")` implements the skew normal distribution.
    `r pkg("greybox")` implements the folded normal
    distribution. `r pkg("VGAM")` implements the folded and
    the skewed normal distribution, and `r pkg("csn")`
    provides d, r functions for the closed skew normal distribution.
    `r pkg("CompQuadForm")` provides the distribution
    function of quadratic forms in normal variates.
    `r pkg("NormalLaplace")` provides d, p, q, r functions
    for the sum of a normal and a Laplace random variables, while
    `r pkg("LaplacesDemon")` provides d, r function of the
    sum of a normal and a Laplace random variables.
    
      --------------------------------- -----------------   --------------- -----------------------
      *Distribution name*               *Packages*          *Functions*     *Distribution suffix*
      Normal                            stats               d, p, q, r      `norm`
      Normal                            actuar              m, mgf          `norm`
      Truncated normal                  truncnorm           d, p, q, r, m   `truncnorm`
      Truncated normal                  EnvStats            d, p, q, r      `normTrunc`
      Truncated normal                  extraDistr          d, p, q, r      `tnorm`
      Truncated normal                  crch                d, p, q, r      `cnorm`
      Generalized normal                lmomco              d, p, q, r      `gno`
      Zero modified Gaussian            EnvStats            d, p, q, r      `zmnorm`
      Exponentially modified Gaussian   emg                 d, p, q, r      `emg`
      Exponentially modified Gaussian   gamlss.dist         d, p, q, r      `exGAUSS`
      Folded and skew normal            gamlss.dist         d, p, q, r      `SN1, SN2`
      Folded normal                     greybox             d, p, q, r      `fnorm`
      Closed skew normal                csn                 d, p, q, r      `csn`
      Skew normal                       sn                  d, p, q, r      `sn`
      --------------------------------- -----------------   --------------- -----------------------

      :  Summary for Gaussian-related distributions

    
-   *General error distribution (also known as exponential power
    distribution):* see *exponential* item.
-   *Generalized extreme value distribution:* d, p, q provided in
    `r pkg("lmomco")`; d, p, q, r, provided in
    `r pkg("VGAM")`, `r pkg("evd")`,
    `r pkg("evir")`, `r pkg("FAdist")`,
    `r pkg("extraDistr")`, `r pkg("EnvStats")`,
    `r pkg("TLMoments")`, `r pkg("rmutil")`,
    `r pkg("QRM")`, `r pkg("ROOPSD")` and
    `r pkg("fExtremes")`. 
    `r pkg("revdbayes")` provide d,p,q,r functions of the
    GEV distribution in a Bayesian setting.
-   *Gompertz distribution:* provided in 
    `r pkg("flexsurv")`, `r pkg("extraDistr")`.
    `r pkg("flexsurv")` also provides hazard (h) and
    integrated hazard rate (i) functions. The shifted Gompertz
    distribution is implemented in `r pkg("extraDistr")`.
    The unit-Gompertz is provided in `r pkg("ugomquantreg")`.
-   *Govindarajulu distribution:* provided in
    `r pkg("lmomco")`.
-   *Gumbel distribution:* provided in packages
    `r pkg("lmomco")`, `r pkg("VGAM")`,
    `r pkg("gamlss.dist")`, `r pkg("FAdist")`,
    `r pkg("extraDistr")`, 
    `r pkg("QRM")`, `r pkg("TLMoments")`,
    `r pkg("dgumbel")`, `r pkg("EnvStats")` and
    `r pkg("evd")`. `r pkg("actuar")` provides
    the raw moments and the moment generating function (mgf) in addition
    to the d, p, q, r functions. A fast random generator is available
    for the Gumbel distribution is implemented in
    `r pkg("Runuran")` as well as the density function. The
    reverse Gumbel distribution is implemented in
    `r pkg("lmomco")` and `r pkg("gamlss.dist")`.
    `r pkg("bgumbel")` provides the bimodel Gumbel distribution.
-   *Hjorth distribution:* provided in `r pkg("rmutil")`.
-   *Huber distribution:* Huber's least favourable distribution
    provided in package `r pkg("smoothmest")` (d, r), and in
    `r pkg("VGAM")`, `r pkg("marg")`,
    `r pkg("extraDistr")` (d, p, q, r).
-   *(generalized) G-and-K, G-and-H distributions:*
    `r pkg("gk")` provides d, p, q, r functions for the
    g-and-k and generalized g-and-h distributions which are nonlinear
    transforms of the Gaussian variables.
-   *(generalized) Hyperbolic distribution:*
    `r pkg("fBasics")`, `r pkg("ghyp")`,
    `r pkg("GeneralizedHyperbolic")` and
    `r pkg("HyperbolicDist")` packages provide d, p, q, r
    functions for the generalized hyperbolic distribution.
    `r pkg("QRM")` provides d, r functions for the
    generalized hyperbolic distribution.
    `r pkg("SkewHyperbolic")` provides the skewed Hyperbolic
    Student t-Distribution. `r pkg("fBasics")` also
    implements the standardized generalized Hyperbolic distribution. A
    fast random generator is available for the hyperbolic distribution
    is implemented in `r pkg("Runuran")` as well as the
    density function.
-   *Hyperbolic sine distribution and extension:*
    `r pkg("gamlss.dist")` provides the sinh and the asinh
    distributions. Generalized Power Hyperbolic sine distributions are
    provided in `r pkg("FatTailsR")`.
-   *Inverse Gaussian (also known Wald) distribution:* d, p, q, and r
    functions of the inverse Gaussian are provided in
    `r pkg("statmod")`, `r pkg("extraDistr")`,
    `r pkg("SuppDists")`, `r pkg("rmutil")`. `r pkg("LaplacesDemon")`
    provides d, r functions for the inverse Gaussian distribution.
    `r pkg("actuar")` provides d, p, q, r, m, lev, mgf
    functions for the Inverse Gaussian distribution.
    `r pkg("SuppDists")` also provides a function that
    returns moments, skewness, kurtosis. `r pkg("fBasics")`
    the normal inverse Gaussian and standardized normal inverse Gaussian
    distributions. The generalized inverse gaussian distribution can be
    found in `r pkg("gamlss.dist")`,
    `r pkg("QRM")`, `r pkg("rmutil")`, and
    `r pkg("HyperbolicDist")`. A random generator is
    available for the (generalized) Inverse Gaussian distribution is
    implemented in `r pkg("Runuran")` as well as the density
    function. `r pkg("GIGrvg")` generates random variables
    from the generalized inverse Gaussian distribution.
-   *Johnson distribution:* provided in
    `r pkg("SuppDists")`. `r pkg("ForestFit")`
    provides d, p of Johnson SB distribution.
-   *Jones and Pewsey distribution:* provided in
    `r pkg("circular")` (d).
-   *K-prime distribution:* `r pkg("sadists")` implements
    Gram Charlier, Edgeworth and Cornish-Fisher approximations for
    K-prime distribution for computing d, p, q, r functions.
-   *Kappa distribution:* A 4-parameter Kappa distribution is provided
    in `r pkg("lmomco")` and `r pkg("FAdist")`.
-   *Kappa-mu distribution:* provided in `r pkg("lmomco")`.
-   *Kato-Jones distribution:* provided in
    `r pkg("circular")` (d, r).
-   *Kendall's tau distribution:* provided in
    `r pkg("SuppDists")`.
-   *Kiener distribution:* a family of distributions generalizing
    hyperbolic sine distributions (see hyperbolic sine section), d, p,
    q, r, m provided in `r pkg("FatTailsR")`.
-   *Kruskal Wallis distribution:* provided in
    `r pkg("SuppDists")`.
-   *Kumaraswamy distribution:* provided in packages
    `r pkg("VGAM")`, `r pkg("extraDistr")` and
    `r pkg("lmomco")`. `r pkg("elfDistr")`
    provides the Kumaraswamy Complementary Weibull Geometric Probability
    Distribution.
-   *(Tukey) Lambda distribution and its extensions:* The generalized
    Lambda distribution (GLD) is well known for its wide range of
    shapes. The original Tukey Lambda distribution can be obtained as a
    special case of the generalized Lambda distribution. There exists
    different parametrization of GLD in the literature: RS
    (Ramberg-Schmeiser or tail-index param), FMKL
    (Freimer-Mudholkar-Kollia-Lin), FM5 (Five-parameter version of FKML
    by Gilchrist), GPD (gen. Pareto dist.) and AS (Asymmetry-steepness).
    The following packages implement such distributions (with d, p, q, r
    functions): `r pkg("gld")` (RS, FKML, FM5, GPD),
    `r pkg("Davies")` (RS), `r pkg("gb")` (RS),
    `r pkg("lmomco")` (FMKL),
    `r pkg("extraDistr")` (original Tukey).
    `r pkg("ecd")` provides the elliptic lambda distribution
    and its use for financial pricing.
-   *Tukey's G/H distribution:* 
    no longer provided directly, but
    Tukey's H distribution is provided as a special case of Lambert W x
    F distribution.
-   *Lambda-prime distribution:* `r pkg("sadists")`
    implements Gram Charlier, Edgeworth and Cornish-Fisher
    approximations for K-prime distribution for computing d, p, q, r
    functions.
-   *Lambert W x F distribution:* `r pkg("LambertW")`
    package provides d, p, q, r functions as well as the first 4 central
    moments and a qqplot.
-   *Laplace (also called double exponential distribution) and asymmetric Laplace distribution:* provided in
    `r pkg("distr")`, `r pkg("lmomco")`,
    `r pkg("LaplacesDemon")`, `r pkg("L1pack")`, `r pkg("VGAM")`,
    `r pkg("sgt")`, `r pkg("extraDistr")`,
    `r pkg("greybox")`, `r pkg("rmutil")` and
    `r pkg("HyperbolicDist")` packages.
    `r pkg("LaplacesDemon")` provides the Laplace
    distribution parametrized by the precision parameter as well as the
    skew Laplace distribution. Asymetric Laplace distribution is
    implemented in `r pkg("ald")`,
    `r pkg("greybox")`. A fast random generator is available
    for the Laplace distribution is implemented in
    `r pkg("Runuran")` as well as the density function.
    `r pkg("smoothmest")` implements the density and the
    random generator. The skew Laplace distribution is available in
    `r pkg("sgt")`. `r pkg("LaplacesDemon")`
    provides the log-Laplace distribution.
-   *LASSO distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Lévy distribution:* provided in `r pkg("rmutil")`.
-   *Lindley distribution:* provided in `r pkg("VGAM")` and
    `r pkg("gambin")`.
-   *Linear failure rate distribution:* no longer implemented.
-   *Loglog distribution:* no longer implemented.
-   *Lomax distribution:* see beta.
-   *Logistic distribution and its extensions:* Base R provides the d,
    p, q, r functions for this distribution (see above).
    `r pkg("actuar")` and `r pkg("VGAM")`
    provide d, p, q, r functions for the log logistic (also called
    Fisk), the paralogistic and the inverse paralogistic distributions.
    `r pkg("FAdist")` the log-logistic distribution with two
    and three parameters. The generalized logistic distribution (Type I,
    also known as skew-logistic distribution) is provided in
    `r pkg("lmomco")`, `r pkg("sld")`,
    `r pkg("rmutil")`, `r pkg("SCI")` and
    `r pkg("glogis")`.
    `r pkg("GTDL")` implements generalized Time-Dependent Logistic distribution.
    
      ---------------------- --------------   -------------------- -----------------------
      *Distribution name*    *Packages*       *Functions*          *Distribution suffix*
      Logistic               stats            d, p, q, r           `logis`
      Logistic               actuar           m, mgf               `logis`
      Log logistic           actuar           d, p, q, r, m, lev   `llogis`
      Log logistic           VGAM             d, p, q, r           `fisk`
      Log logistic           FAdist           d, p, q, r           `llog, llog3`
      Paralogistic           actuar           d, p, q, r, m, lev   `paralogis`
      Paralogistic           VGAM             d, p, q, r           `paralogistic`
      Inv. paralogistic      actuar           d, p, q, r, m, lev   `invparalogis`
      Inv. paralogistic      VGAM             d, p, q, r           `inv.paralogistic`
      Truncated logistic     crch             d, p, q, r           `tlogis`
      Generalized logistic   glogis           d, p, q, r           `glogis`
      Generalized logistic   SCI              d, p, q              `genlog`
      Generalized logistic   lmomco           d, p, q, r           `glo`
      Generalized logistic   sld              d, p, q, r           `sl`
      Generalized logistic   rmutil           d, p, q, r           `glogis`
      ---------------------- --------------   -------------------- -----------------------

      :  Summary for Logistic-related distributions

    
-   *Logit-normal distribution:* provided in
    `r pkg("logitnorm")`.
-   *Log-normal distribution and its extensions:* The log normal
    distribution is implemented in Base R (see above) and
    `r pkg("poweRlaw")`. The log normal distribution
    parametrized by its mean and its coefficient of variation is also
    provided in `r pkg("EnvStats")`.
    `r pkg("LaplacesDemon")` provides the lognormal
    parametrized by the precision parameter. The truncated lognormal
    distribution is provided in `r pkg("EnvStats")` with two
    possible parametrizations as well as in
    `r pkg("ReIns")`. The 3-parameter lognormal distribution
    is available in `r pkg("lmomco")`,
    `r pkg("greybox")`, `r pkg("TLMoments")`,
    `r pkg("EnvStats")` and `r pkg("FAdist")`.
    The package `r pkg("loglognorm")` implements d, p, q, r
    functions for the double lognormal distribution, as well as the raw
    moment, the expected value and the variance functions.
    `r pkg("EnvStats")` provides d, p, q, r functions for
    the zero-modified lognormal distribution with two possible
    parametrizations. `r pkg("distributionsrd")` provides
    the double Pareto-lognormal distribution, the left Pareto-lognormal
    distribution, the truncated lognormal distribution.
-   *Makeham distribution:* provided in `r pkg("VGAM")` and
-   *Maxwell distribution:* provided in `r pkg("VGAM")`.
-   *Minimax distribution:* provided in `r pkg("minimax")`.
-   *Mittag-Leffler distribution:* d, p, q, r functions provided in
    `r pkg("MittagLeffleR")`.
-   *Nakagami distribution:* provided in `r pkg("VGAM")`.
-   *Pareto distribution:* d, p, q, r functions are implemented in
    `r pkg("VGAM")` for the Pareto distribution type IV
    (which includes Burr's distribution, Pareto type III, Pareto type
    II (also called the lomax distribution) and Pareto type I) and the
    (upper/lower) truncated Pareto distribution. In an actuarial
    context, `r pkg("actuar")` provides d, p, q, r functions
    as well as moments and limited expected values for the Pareto I and
    II, the inverse Pareto, the 'generalized pareto' distributions,
    the Burr and the inverse Burr distributions, all special cases of
    the transformed beta II distribution. A fast random generator for
    the Burr and the Pareto II distribution is implemented in
    `r pkg("Runuran")` as well as the density.
    `r pkg("EnvStats")` and
    `r pkg("LaplacesDemon")` provides d, p, q, r functions
    for Pareto I distribution. `r pkg("extremefit")`
    provides the Burr, the Pareto II, mixture of Pareto I distributions
    and a composite distribution of two Pareto I distributions.
    `r pkg("lmomco")`, `r pkg("evd")`,
    `r pkg("fExtremes")`, `r pkg("extraDistr")`,
    `r pkg("QRM")`, `r pkg("Renext")`,
    `r pkg("revdbayes")`, `r pkg("FAdist")`,
    `r pkg("LaplacesDemon")`,
    `r pkg("TLMoments")` `r pkg("qrmtools")` and
    `r pkg("evir")` packages implement the Generalized
    Pareto Distribution (from Extreme Value Theory), which is depending
    the shape parameter's value a Pareto II distribution, a shifted
    exponential distribution or a generalized beta I distribution.
    `r pkg("ParetoPosStable")` implements the Pareto
    positive stable distribution. The extended Pareto distribution is
    implemented in `r pkg("RTDE")` and the shifted truncated
    (to unit interval) Pareto is implemented in
    `r pkg("mbbefd")`. `r pkg("ReIns")` provides
    Burr, extended Pareto, generalized Pareto, Pareto 1 distributions
    and their truncated version. `r pkg("CaDENCE")` provides
    the Pareto 2 and the zero-inflated Pareto 2 distribution.
    `r pkg("Pareto")` provides the Pareto 1, piecewise Pareto 
    and the generalized Pareto (from actuarial theory).
    
      -------------------------- ----------------- -------------------- -----------------------
      *Distribution name*        *Packages*        *Functions*          *Distribution suffix*
      Pareto I                   VGAM              d, p, q, r           `paretoI`
      Pareto I                   actuar            d, p, q, r, m, lev   `pareto1`
      Pareto I                   EnvStats          d, p, q, r           `pareto`
      Pareto I                   extraDistr        d, p, q, r           `pareto`
      Pareto I                   ReIns             d, p, q, r           `pareto`
      Pareto I                   LaplacesDemon     d, p, q, r           `pareto`
      Pareto I                   distributionsrd   d, p, q, r           `pareto`
      Pareto I                   Pareto            d, p, q, r           `Pareto`
      Trunc. Pareto I            ReIns             d, p, q, r           `tpareto`
      Pareto II                  VGAM              d, p, q, r           `paretoII`
      Pareto II                  actuar            d, p, q, r, m, lev   `pareto, pareto2`
      Pareto II                  Runuran           d, r                 `pareto`
      Pareto II                  extraDistr        d, p, q, h           `lomax`
      Pareto II                  extremefit        d, p, q, h           `pareto`
      Pareto II                  Renext            d, p, q, r           `lomax`
      Pareto II                  rmutil            d, p, q, r           `pareto`
      Pareto II                  CaDENCE           d, p, q, r           `pareto2`
      zero-inflated Pareto II    CaDENCE           d, p, q, r           `bpareto2`
      Pareto III                 VGAM              d, p, q, r           `paretoIII`
      Pareto III                 actuar            d, p, q, r           `pareto3`
      Pareto IV                  VGAM              d, p, q, r           `paretoIV`
      Pareto IV                  actuar            d, p, q, r           `pareto4`
      Inverse Pareto             actuar            d, p, q, r, m, lev   `invpareto`
      Inverse Pareto             distributionsrd   d, p, q, r, m, lev   `invpareto`
      Extended Pareto            RTDE              d, p, q, r           `EPD`
      Extended Pareto            ReIns             d, p, q, r           `epd`
      Shift. trunc. Pareto       mbbefd            d, p, q, r, m, ec    `stpareto`
      Gen. Pareto (actuarial)    actuar            d, p, q, r, m, lev   `genpareto`
      Gen. Pareto (actuarial)    Pareto            d, p, q, r           `GenPareto`
      Gen. Pareto (EVT)          lmomco            d, p, q, r           `gpa`
      Gen. Pareto (EVT)          evd               d, p, q, r           `gpd`
      Gen. Pareto (EVT)          fExtremes         d, p, q, r           `gpd`
      Gen. Pareto (EVT)          evir              d, p, q, r           `gpd`
      Gen. Pareto (EVT)          extraDistr        d, p, q, r           `gpd`
      Gen. Pareto (EVT)          QRM               d, p, q, r           `GPD`
      Gen. Pareto (EVT)          ReIns             d, p, q, r           `gpd`
      Gen. Pareto (EVT)          LaplacesDemon     d, r                 `gpd`
      Gen. Pareto (EVT)          TLMoments         d, p, q, r           `gpd`
      Trunc. Gen. Pareto (EVT)   ReIns             d, p, q, r           `tgpd`
      Gen. Pareto (EVT)          revdbayes         d, p, q, r           `gp`
      Gen. Pareto (EVT)          Renext            d, p, q, r           `GPD`
      Gen. Pareto (EVT)          qrmtools          d, p, q, r           `GPD`
      Gen. Pareto (EVT)          ROOPSD            d, p, q, r           `gpd`
      Feller-Pareto              actuar            d, p, q, r, m, lev   `fpareto`
      Burr                       actuar            d, p, q, r, m, lev   `burr`
      Burr                       extremefit        d, p, q, r           `burr`
      Burr                       ReIns             d, p, q, r           `burr`
      Burr                       rmutil            d, p, q, r           `burr`
      Trunc. Burr                ReIns             d, p, q, r           `tburr`
      Inverse Burr               actuar            d, p, q, r, m, lev   `invburr`
      -------------------------- ----------------  -------------------- -----------------------

      :  Summary for Pareto-related distributions

    
-   *Pearson's distribution:* Pearson type III available in
    `r pkg("lmomco")` and `r pkg("FAdist")`. A
    log-Pearson type III distribution is also available in
    `r pkg("FAdist")`.
    `r pkg("PearsonDS", priority = "core")` provides the d,
    p, q, r functions as well as the first four moments for the Pearson
    distributions: types I, II, III, IV, V, VI, VII.
-   *Pearson's Rho distribution:* provided in
    `r pkg("SuppDists")`.
-   *Perks distribution:* provided in `r pkg("VGAM")`.
-   *Planck's distribution:* a random generator is available in
    `r pkg("Runuran")`.
-   *Phase-type distribution:* provided in `r pkg("actuar")`,
    `r pkg("mapfit")`, `r pkg("matrixdist")`, `r pkg("PhaseTypeR")`.
-   *Power distribution:* `r 
    `r pkg("poweRlaw")` implement the exponential power
    distribution. Two-sided power distribution provided in
    `r pkg("rmutil")`.
-   *Proportion distribution:* this is the distribution for the
    difference between two independent beta distributions. d, p, q, r
    functions in `r pkg("tolerance")`.
-   *Rayleigh distribution:* provided in packages
    `r pkg("VGAM")`, `r pkg("extraDistr")` and
    `r pkg("lmomco")`. 
-   *Response time distribution:* `r pkg("rtdists")`
    provides d, p, q, r functions for the (Ratcliff) diffusion
    distribution and for the linear ballistic accumulator (LBA) with
    different underlying drift-distributions (Normal, Gamma, Frechet,
    and log-normal).
-   *Rice distribution:* provided in `r pkg("VGAM")` and
    `r pkg("lmomco")`.
-   *Simplex distribution:* provided in `r pkg("rmutil")`.
-   *Singh-Maddala distribution:* see beta.
-   *Slash distribution:* provided in `r pkg("lmomco")`,
    `r pkg("extraDistr")` and `r pkg("VGAM")`.
-   *Spearman's Rho distribution:* provided in
    `r pkg("SuppDists")`.
-   *Stable distribution:* d, p, q, r functions are available in
    `r pkg("fBasics")` and `r pkg("stabledist")`, the functions use the approach of
    J.P. Nolan for general stable distributions.
    `r pkg("stable")` (d, p, q, r, h) is also used for general stable and uses
    a modified Buck parametrization.
    `r pkg("MixedTS")` provides mixed tempered stable
    distribution (d, p, q, r). `r pkg("FMStable")` provides
    (d, p, q) the extremal or maximally skew stable and the finite
    moment log stable distributions.
    `r pkg("SymTS")` provides (d, p, q, r) functions for symmetric stable, symmetric 
    classical tempered stable, and symmetric power tempered stable distributions.
-   *Student distribution and its extensions:* Base R provides the d, p,
    q, r functions for Student and non central Student distribution (see
    above). `r pkg("extraDistr")` and
    `r pkg("LaplacesDemon")` provides the Student
    distribution with location and scale parameters.
    `r pkg("LaplacesDemon")` provides d, p, q, r functions
    for the Half-Student distribution. `r pkg("sadists")`
    implements Gram Charlier, Edgeworth and Cornish-Fisher
    approximations for doubly non central Student distribution for
    computing d, p, q, r functions. The skewed Student distribution is
    provided in `r pkg("skewt")`, `r pkg("sn")`
    and `r pkg("gamlss.dist")` packages. The generalized
    skew distribution is provided in `r pkg("sgt")`. d, p,
    q, r functions for the generalized t-distribution can be found in
    `r pkg("gamlss.dist")`. `r pkg("fBasics")`
    provides d, p, q, r functions for the skew and the generalized
    hyperbolic t-distribution. The L-moments of the Student t
    (3-parameter) are provided in `r pkg("lmomco")`.
    `r pkg("crch")` provides d, p, q, r functions for the 
    truncated student distribution.
    
      ----------------------------- -----------------   ------------- ---------------------------
      *Distribution name*           *Packages*          *Functions*   *Distribution suffix*
      Student                       stats               d, p, q, r    `t`
      Student with loc. and scal.   extraDistr          d, p, q, r    `lst`
      Student with loc. and scal.   LaplacesDemon       d, p, q, r    `st`
      Doubly non central St.        sadists             d, p, q, r    `dnt`
      Skew Student                  skewt               d, p, q, r    `skt`
      Skew Student                  sn                  d, p, q, r    `st`
      Skew St. Type 1-5             gamlss.dist         d, p, q, r    `ST1, ST2, ST3, ST4, ST5`
      Gen. Student                  gamlss.dist         d, p, q, r    `GT`
      Gen. Hyp. Student             fBasics             d, p, q, r    `ght`
      Skew Gen. Student             sgt                 d, p, q, r    `sgt`
      ----------------------------- -----------------   ------------- ---------------------------

      :  Summary for Student-related distributions

    
-   *Triangle/trapezoidal distribution:* packages
    `r pkg("triangle")`, `r pkg("extraDistr")`,
    `r pkg("mc2d")`, `r pkg("EnvStats")` and
    `r pkg("VGAM")` provide d, p, q, r functions for the
    triangle or triangular distribution, while the package
    `r pkg("trapezoid")` provides d, p, q, r functions for
    the Generalized Trapezoidal Distribution.
    `r pkg("CircStats")`, `r pkg("circular")`
    provide d, r functions for triangular distribution. A fast random
    generator is available for the triangle distribution is implemented
    in `r pkg("Runuran")` as well as the density function.
-   *Tsallis or q-Exponential distribution:*
    `r pkg("tsallisqexp")` provides d, p, q, r functions for
    two parametrizations of the Tsallis distribution and also implements
    a left-censored version.
-   *Tweedie distribution:* the Tweedie distribution is implemented in
    package `r pkg("tweedie")`. Let us note that the Tweedie
    distribution is not necessarily continuous, a special case of it is
    the Poisson distribution.
-   *Uniform distribution:* d, p, q, r functions are of course provided
    in R. See section RNG for random number generation topics.
    `r pkg("KScorrect")` provides d, p, q, r functions for
    the log-uniform distribution.
-   *Upsilon distribution:* `r pkg("sadists")` implements
    Gram Charlier, Edgeworth and Cornish-Fisher approximations for
    Upsilon distribution for computing d, p, q, r functions.
-   *Vasicek distribution:* 
    `r pkg("vasicek")` implements d, p, r functions.
    `r pkg("vasicekreg")` implements d, p, q, r functions.    
-   *von Mises distribution:* The `r pkg("CircStats")`
    package provides d, p, r functions; the
    `r pkg("circular")` package provides d, p, q, r
    functions.
    `r pkg("rvMF")` package provides a fast random generator for von Mises
    Fisher distribution.
    
-   *Wakeby distribution:* A 5-parameter Wakeby is provided in
    `r pkg("lmomco")`.
-   *Weibull distribution and its extensions:* Base R provides the d, p,
    q, r functions for this distribution (see above). The inverse
    Weibull is provided in `r pkg("actuar")` package and
    also the moments and the limited expected value for both the raw and
    the inverse Weibull distribution. `r pkg("FAdist")`
    implements the three-parameter Weibull distribution. 
    Furthermore, `r pkg("lmomco")` implements
    the Weibull distribution while `r pkg("evd")` implements
    the reverse Weibull distribution. The reverse generalized extreme
    value distribution are provided in
    `r pkg("gamlss.dist")` (d, p, q, r) and the shifted left
    truncated Weibull distribution is provided in
    `r pkg("Renext")`. The right truncated Weibull is
    provided in `r pkg("ReIns")`. The generalized Weibull is
    provided in `r pkg("rmutil")`. The tail Weibull is
    provided in `r pkg("distTails")`.
    `r pkg("CaDENCE")` provides the zero-inflated Weibull
    distribution.
-   *(first-passage time) of a Wiener process:*    
    `r pkg("WienR")` provides d, p functions of the 
    first-passage time of a diffusion model.
mode

## [Continuous multivariate distributions:]{#MultivariateContinuous}


-   *Bivariate Pareto:* `r pkg("Bivariate.Pareto")` provides
    a random generator for the bivariate Pareto distribution.
-   *Multivariate beta distribution:*
    `r pkg("NonNorMvtDist")` provides d, p, q, r, s
    functions for inverted beta distribution.
-   *Multivariate Burr distribution:*
    `r pkg("NonNorMvtDist")` provides d, p, q, r, s
    functions.
-   *Multivariate Cauchy distribution:* `r pkg("sn")`
    provide d, p, r functions for the multivariate skew Cauchy
    distribution, while `r pkg("LaplacesDemon")` provides d,
    r functions for the multivariate Cauchy distribution parametrized
    either by sigma, by the Cholesky decomposition of sigma, by the
    precision matrix omega or by the Cholesky decomposition of omega.
-   *Cook-Johnson's Multivariate Uniform Distribution:*
    `r pkg("NonNorMvtDist")` provides d, p, q, r, s
    functions.
-   *Multivariate Dirichlet distribution:*
    `r pkg("Compositional")`,
    `r pkg("LaplacesDemon")`,
    `r pkg("MCMCpack")` packages provide d, r functions as
    well as a fitting function for `r pkg("Compositional")`.
    `r pkg("compositions")`, `r pkg("bayesm")`
    provide r function.
    `r pkg("SGB")` provides a generalization of the Dirichlet 
    distribution called Simplicial Generalized Beta distribution.
-   *Multivariate exponential distribution:* while
    `r pkg("LaplacesDemon")` provides d, r functions for the
    multivariate power exponential distribution parametrized either by
    sigma, or by the Cholesky decomposition of sigma.
-   *Multivariate F distribution:* `r pkg("NonNorMvtDist")`
    provides d, p, q, r, s functions.
-   *Multivariate Gaussian (or normal) distribution:* The multivariate
    Gaussian distribution is provided in the packages
    `r pkg("mvtnorm", priority = "core")` (d, p, r),
    `r pkg("mnormt", priority = "core")` (d, p, r),
    `r pkg("mniw")` (d, r),
    `r pkg("Compositional")` (r),
    `r pkg("compositions")` (r). `r pkg("pbv")`
    provides d, p functions for bivariate normal distributions.
    `r pkg("symmoments")` computes central and non-central
    moments of the multivariate Gaussian distribution.
    `r pkg("LaplacesDemon")` provides d, r functions for the
    multivariate normal distribution parametrized either by sigma, by
    the Cholesky decomposition of sigma, by the precision matrix omega
    or by the Cholesky decomposition of omega. Futhermore, the
    multivariate truncated normal is implemented in
    `r pkg("TruncatedNormal")` for d, p, r functions;
    `r pkg("tmvtnorm")` for p, q, r, m(oments) functions;
    `r pkg("tmvmixnorm")` for a fast RNG.
    `r pkg("sparseMVN")` implements very fast algorithms to
    compute the density and generate random variates of a multivariate
    normal distribution for which the covariance matrix or precision
    matrix is sparse. `r pkg("cmvnorm")` implements the
    complex multivariate normal distribution (d, r). Finally,
    `r pkg("condMVNorm")` implements d, p, r functions for
    the conditional multivariate normal distribution. Furthermore,
    `r pkg("sn")` besides providing facilities for their
    distribution functions, `r pkg("sn")` allows the
    creation of S4 objects which encapsulate these distributions and
    provide facilities for plotting, summary, marginalization,
    conditioning, affine transformations of these S4 objects.
    `r pkg("Compositional")` provides random generator for
    the multivariate normal distribution on the simplex and multivariate
    skew normal distribution on the simplex. A random generator of the
    multivariate normal is provided in `r pkg("MultiRNG")`.
-   *Multivariate generalized hyperbolic distribution:*
    `r pkg("QRM")` provides d, r functions of the standard
    and the symmetric multivariate generalized hyperbolic distribution.
    `r pkg("ghyp")` provides d, p, r functions of the
    standard multivariate generalized hyperbolic distribution.
-   *Multivariate generalized extreme value distribution:* Both
    bivariate and multivariate Extreme Value distributions as well as
    order/maxima/minima distributions are implemented in
    `r pkg("evd")` (d, p, r).
-   *Multivariate Laplace distribution:*
    `r pkg("LaplacesDemon")` provides d, r functions for the
    multivariate Laplace distribution parametrized either by sigma, or
    by the Cholesky decomposition of sigma. r is provided in
    `r pkg("MultiRNG")`.
    `r pkg("L1pack")` provides d, r functions of the multivariate
    Laplace distribution.
-   *Multivariate logistic distribution:* `r pkg("VGAM")`
    package implements the bivariate logistic distribution, while
    `r pkg("NonNorMvtDist")` implements the multivariate
    logistic distribution.
-   *Multivariate lognormal distribution:*
    `r pkg("compositions")` provides r function.
-   *Multivariate Pareto distribution:* `r pkg("evd")` provides
    the density for the multivariate generalized Pareto type I. 
    `r pkg("NonNorMvtDist")`
    provides d, p, q, r, s functions for multivariate Lomax (type II)
    distributions and its generalized version.
    `r pkg("NonNorMvtDist")` provides d, p, q, r, s
    functions for Mardia's Multivariate Pareto Type I Distribution
-   *Multivariate Stable distribution:* For elliptically contoured (subgaussian 
     stable), `r pkg("alphastable")` provides d, r functions as well as a 
     fitting function, `r pkg("mvgb")` provides p function.
-   *Multivariate Student distribution:* The multivariate Student
    distribution is provided in the packages
    `r pkg("mvtnorm")` (d, r), `r pkg("mnormt")`
    (d, p, r), `r pkg("Compositional")` (r),
    `r pkg("tmvmixnorm")` (r), `r pkg("QRM")`
    (d, r), `r pkg("bayesm")` (r), `r pkg("MVT")` (r). 
    `r pkg("TruncatedNormal")` for d, p, r functions;
    `r pkg("tmvtnorm")` for d, p, q, r functions.
    `r pkg("sn")` provides d, p, r functions for the
    multivariate skew t distribution.
    `r pkg("LaplacesDemon")` provides d, r functions for the
    multivariate Student distribution parametrized either by sigma, by
    the Cholesky decomposition of sigma, by the precision matrix omega
    or by the Cholesky decomposition of omega. Random generator r is
    provided in `r pkg("MultiRNG")`. A special case of a
    bivariate noncentral t-distribution called Owen distribution is
    provided in `r pkg("OwenQ")`.
-   *Multivariate Uniform distribution:* r is provided in
    `r pkg("MultiRNG")`. `r pkg("compositions")`
    provides a random generator on the simplex.

# [Other distributions:]{#Other}

## [Mixed-type distributions:]{#MixedType}


-   *Maxwell-Boltzmann-Bose-Einstein-Fermi-Dirac (MBBEFD) distribution
    :* provided in `r pkg("mbbefd")`.
-   *Mixed ordinal and normal distribution:* provided in
    `r pkg("OrdNor")`.
-   *One-inflated distributions:* a generic distribution as well as
    special cases (OI-beta, OI-uniform, OI-GB1, OI-Pareto) are provided
    in `r pkg("mbbefd")`. The zero and one inflated beta
    distribution can be found in `r pkg("gamlss.dist")`.
-   *Zero-modified distributions:* `r pkg("EnvStats")`
    provides the zero-modified normal distribution and the zero-modified
    lognormal distribution.

## [Mixture of probability laws:]{#Mixture}


-   *Bernoulli-dist mixture:* d, p, q, r functions for
    Bernoulli-exponential, Bernoulli-Gamma, Bernoulli-lognormal,
    Bernoulli-Weibull distributions are provided in
    `r pkg("qmap")`.
-   *Cauchy-polynomial quantile mixture:* d, p, q, r functions are
    provided in `r pkg("Lmoments")`.
-   *Chi-square mixture:* d, p, q, r functions are provided in
    `r pkg("emdbook")`.
-   *Gaussian mixture:* Functions d, r are provided in
    `r pkg("mixtools")`, `r pkg("bmixture")`
    package when dealing with finite mixture models.
    `r pkg("nor1mix")`, `r pkg("extraDistr")`,
    `r pkg("mclust")`, `r pkg("LaplacesDemon")`,
    `r pkg("KScorrect")` provides d, p, r functions for
    Gaussian mixture. `r pkg("EnvStats")` provides d, p, q,
    r functions for mixture of two normal distributions.
    `r pkg("bayesm")` provides d function for the mixture of
    multivariate normals.
-   *Gamma Poisson:* provided in `r pkg("extraDistr")`.
-   *Gamma mixture:* Ga `r pkg("GSM")` package provides d,
    p, r, `r pkg("bmixture")` provides d, r,
    `r pkg("evmix")` provides d, p, q, r.
-   *Generic mixtures:* there is an implementation via S4-class
    UnivarMixingDistribution in package `r pkg("distr")`.
    `r pkg("gendist")` provides d, p, q, r functions for
    two-distribution mixture models working with any distribution
    defined by its d, p, q, r functions.
-   *Horseshoe distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Laplace mixture distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Log normal mixture:* d, p, q, r functions are provided in
    `r pkg("EnvStats")` with two possible parametrizations.
-   *Normal-polynomial quantile mixture:* d, p, q, r functions are
    provided in `r pkg("Lmoments")`.
-   *Pareto distribution:* `r pkg("extremefit")` implements
    the mixture of two Pareto I distributions.
-   *Poisson beta distribution:* provided in
    `r pkg("scModels")`.
-   *Poisson Binomial distribution:* `r pkg("poibin")`
    implements the Poisson Binomial distribution.
-   *Poisson lognormal distribution:* `r pkg("poilog")`
    implements the Poisson lognormal distribution.
-   *Poisson mixture:* provided in `r pkg("extraDistr")`.
-   *Poisson-Tweedie exponential family models:* provided in
    `r pkg("poistweedie")`.
-   *Student mixture:* The `r pkg("AdMit")` package provides
    d, r functions for Student mixtures in the context of Adaptive
    Mixture of Student-t distributions. `r pkg("bmixture")`
    package also provide d, r functions for mixture of Student-t
    distributions.
-   *von Mises Fisher (or Langevin) mixture:* The
    `r pkg("movMF")` and `r pkg("CircStats")`
    packages provide d, r functions for finite von Mises Fisher
    mixtures.
    
    
## [Random matrices:]{#Matrix}


-   *Huang-Wan distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Inverse matrix gamma distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Inverse Wishart distribution:* `r pkg("LaplacesDemon")`
    provides inverse Wishart distribution parametrized either by Sigma
    or by its Cholesky decomposition.
    `r pkg("LaplacesDemon")` provides the scaled inverse
    Wishart distribution. `r pkg("MCMCpack")` and
    `r pkg("mniw")` provides the inverse Wishart
    distribution.
-   *Marcenko-Pastur distribution:* provided in
    `r pkg("RMTstat")`, `r pkg("MCMCpack")` and
    `r pkg("bayesm")`.
-   *Matrix gamma distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Matrix normal distribution:* `r pkg("MBSP")` (r)
    provides a random generator using a Cholesky decomposition;
    `r pkg("matrixsampling")` (r) provides a random
    generator using a spectral decomposition;
    `r pkg("LaplacesDemon")` and `r pkg("mniw")`
    (d, r); `r pkg("matrixNormal")` (d, p, r) collects these
    forms in one place and allows users to be flexible in simulating
    random variates (Cholesky, spectral, SVD).
    
-   *Matrix student distribution:* provided in
    `r pkg("mniw")`.
-   *Normal Inverse Wishart distribution:* provided in
    `r pkg("LaplacesDemon")`, `r pkg("mniw")`.
-   *Normal Wishart distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Tracy-Widom distribution:* provided in
    `r pkg("RMTstat")`, `r pkg("MCMCpack")` and
    `r pkg("bayesm")`: supported beta values are 1 (Gaussian
    Orthogonal Ensemble), 2 (Gaussian Unitary Ensemble), and 4 (Gaussian
    Symplectic Ensemble).
-   *Sparse matrix:* `r pkg("spam")` provides
    functionalities to draw random numbers from a user-supplied RNG
    (e.g. `rexp`) or from a multivariate normal distribution for large
    sparse matrices: typically for sparse covariance matrices.
-   *Spiked Wishart Maximum Eigenvalue Distribution:* provided in
    `r pkg("RMTstat")`, `r pkg("MCMCpack")` and
    `r pkg("bayesm")`.
-   *Wishart distributions:* Base R provides the r function for the
    Wishart distribution. `r pkg("MCMCpack")`,
    `r pkg("RMTstat")`, `r pkg("bayesm")`,
    `r pkg("mniw")` provides d, r functions,
    `r pkg("bayesm")` provides r function.
    `r pkg("LaplacesDemon")` provides Wishart distribution
    parametrized either by Sigma or by its Cholesky decomposition.
-   *White Wishart Maximum Eigenvalue Distribution:* provided in
    `r pkg("RMTstat")`, `r pkg("MCMCpack")` and
    `r pkg("bayesm")`.
-   *Yang-Berger distribution:* provided in
    `r pkg("LaplacesDemon")`.
-   *Zellner distribution:* provided in
    `r pkg("LaplacesDemon")`.

## [Copulas:]{#Copulas}


-   *Unified approaches:* The packages
    `r pkg("fCopulae", priority = "core")`,
    `r pkg("copula", priority = "core")`, and
    `r pkg("copBasic")` provide a lot of general
    functionality for copulas. Although lacking support for many
    existing copulas themselves, `r pkg("copBasic")` is
    primarily oriented around utility functions for the general
    mathematics of copulas as described in the well known introduction
    to copulas by Nelsen.
-   *Archimedean copulas:* `r pkg("gumbel")` is a standalone
    package for the Gumbel copula `r pkg("fCopulae")`
    implements the 22 Archimedean copulas of Nelsen (1998, *Introduction
    to Copulas* , Springer-Verlag) including Gumbel, Frank, Clayton, and
    Ali-Mikhail-Haq. `r pkg("VGAM")` provides
    Ali-Mikhail-Haq, Clayton, Frank, Frechet copulas.
    `r pkg("copula")` provides Ali-Mikhail-Haq, Clayton,
    Frank, Gumbel and Joe copulas. The Frank bivariate distribution is
    available in `r pkg("RTDE")`.
    `r pkg("VineCopula")` provides Clayton, Gumbel, Frank,
    Joe, BB1, BB6, BB7 and BB8 copulas. Nested Archimedean copulas are
    available in the `r pkg("HAC")` package. 
    `r pkg("copBasic")` provides functions for
    Ali-Mikhail-Haq, Clayton, Frechet copulas.
    `r pkg("QRM")` provides pdf and random generator for
    Clayton, Gumbel, Frank, BB9 copula.
    `r pkg("Bivariate.Pareto")` provides a random generator
    for the Frank copula with Pareto margins.
    `r pkg("nCopula")`, `r pkg("HAC")` provide
    hierarchical archimedean copulas. `r pkg("lcopula")`
    provides the Liouville copula.
    `r pkg("CopulaGAMM")` provides the bivariate version of Frank, FGM,
    Galambos, Gumbel, Huesler-Reiss, Joe, MTCJ, Plackett copulas.
-   *Blomqvist copula:* provided in `r pkg("copBasic")`.
-   *Composition of copula:* `r pkg("copBasic")` provides
    functions for composition of a single symmetric copula and
    composition of two copulas.
-   *Cubic copula:* Not yet implemented?
-   *Dirichlet copula:* Not yet implemented?
-   *Empirical copula:* provided in `r pkg("copBasic")`,
    `r pkg("HAC")`.
    `r pkg("GenOrd")` provides sampling function for
    multivariate discrete random vectors with a specified correlation
    matrix.
-   *Elliptical copulas:* Gaussian, Student and Cauchy copulas are
    implemented in `r pkg("fCopulae")` for the bivariate
    cases. `r pkg("copula")`, `r pkg("VGAM")`,
    `r pkg("VineCopula")` provide the Gaussian and the
    Student copulas. `r pkg("QRM")` provides pdf and random
    generator for Gaussian, Student copulas.
    `r pkg("relliptical")` provides a random generator for multivariate
    truncated Normal, Student-t, Power Exponential,
    Pearson VII, Slash and Contaminated Normal distributions.
    `r pkg("CopulaGAMM")` provides the bivariate Gaussian and 
    student copula.
-   *Extreme value copulas:* `r pkg("fCopulae")` provides
    the following copulas Gumbel, Galambos, Husler-Reiss, Tawn, or BB5.
    `r pkg("copula")` implements Gumbel, Galambos and
    Husler-Reiss.
-   *Eyraud-Farlie-Gumbel-Morgenstern copula:* provided in
    `r pkg("VGAM")`, `r pkg("RTDE")`, and
    `r pkg("copula")`.
-   *Integrated gamma copula:* provided in `r pkg("igcop")`.
-   *Mardia copula:* Not yet implemented?
-   *Nested copulas:* arbitrary nested versions of copulas can be
    implemented in `r pkg("copula")`.
-   *Plackett:* provided in `r pkg("VGAM")`,
    `r pkg("copBasic")` and `r pkg("copula")`.
-   *Vine copulas:* Package `r pkg("vines")` provides
    functions for C- and D-vine copulas and
    `r pkg("VineCopula")` for general R-vine copulas.

## [Compound, composite, discretized, exponentiated and transformation of distributions:]{#Transform}


-   *Absolute value or half distribution:* Half-Cauchy, half normal and
    half-student are implemented both in
    `r pkg("extraDistr")` and in
    `r pkg("LaplacesDemon")`.
-   *Composite distribution also known as spliced distribution:*
    Split-normal (also known as
    the two-piece normal distribution) not yet implemented.
    Split-student provided in package `r pkg("dng")`.
    `r pkg("evmix")` provides d, p, q, r of the following
    composite distributions: gamma-GPD, lognormal GPD, normal-GPD,
    Weibull-GPD as well as bulk models such as GPD-normal-GPD
    distribution. `r pkg("gendist")` provides d, p, q, r
    functions for composite models working with any distribution defined
    by its d, p, q, r functions.
-   *Compound distribution:* 
    `r pkg("kdist")` provides d, p, q, r functions of the K
    distribution.
-   *Discretized distribution:* `r pkg("distcrete")` allows
    discretised versions of continuous distribution by mapping
    continuous values to an underlying discrete grid, based on a
    (uniform) frequency of discretisation, a valid discretisation point,
    and an integration range.
-   *Quantile-based asymmetric (QBA) family of distributions:*
    no longer implemented.
-   *Transformed distribution:* `r pkg("Newdistns")`
    provides G-transformed distributions for a selected number of
    distributions which includes Marshall Olkin G distribution,
    exponentiated G distribution, beta G distribution, gamma G
    distribution, Kumaraswamy G distribution, generalized beta G
    distribution, beta extended G distribution, gamma G distribution,
    gamma uniform G distribution, beta exponential G distribution,
    Weibull G distribution, log gamma G1/G2 distribution, exponentiated
    generalized G distribution, exponentiated Kumaraswamy G
    distributions, geometric exponential Poisson G distribution,
    truncated-exponential skew-symmetric G distribution, modified beta G
    distribution, and exponentiated exponential Poisson G distribution.
    `r pkg("MPS")` provides also G-transformed
    distributions, such as beta exponential G distribution, beta G
    distribution, exponentiated exponential Poisson G distribution,
    exponentiated G distribution, exponentiated generalized G
    distribution, exponentiated Kumaraswamy G distribution, gamma
    uniform G distribution, gamma uniform type I/II G distribution,
    generalized beta G distribution, geometric exponential Poisson G
    distribution, gamma-X family of modified beta exponential G
    distribution, exponentiated exponential Poisson G distribution,
    gamma-X generated of log-logistic-X familiy of G distribution,
    Kumaraswamy G distribution, log gamma G type I/II distribution,
    modified beta G distribution, Marshall-Olkin Kumaraswamy G
    distribution, odd log-logistic G distribution, truncated-exponential
    skew-symmetric G distribution, T-X{log-logistic}G distribution,
    Weibull G distribution. `r pkg("gendist")` provides d,
    p, q, r functions for composite models, folded models, skewed
    symmetric models and arctan models working with any distribution
    defined by its d, p, q, r functions.
    `r pkg("ComRiskModel")` provides also G-transformed such as
     binomial-G, complementary negative binomial-G and complementary 
     geometric-G families of distributions taking baseline models 
     such as exponential, extended exponential, Weibull, extended 
     Weibull, Fisk, Lomax, Burr-XII and Burr-X. 
    
-   *Truncated distribution:* A generic code snippet is available [in
    the JSS](https://www.jstatsoft.org/article/view/v016c02) . This code
    is now available in two packages: `r pkg("truncdist")`
    is a dedicated package providing d, p, q, r, m(oments) functions for
    a univariate truncated distribution given a user-supplied
    distribution; `r pkg("LaplacesDemon")` provides a
    generic function in a Bayesian environment.

# [Moments, skewness, kurtosis and etc:]{#Moments}


-   *Empirical mean, standard deviation and variance:* base R provides
    `mean()`, `sd()`, `var()` functions to compute the mean, standard
    deviation and variance, respectively.
-   *Empirical skewness:* available in `r pkg("agricolae")`,
    `r pkg("e1071")`, `r pkg("GLDEX")`,
    `r pkg("HyperbolicDist")`,
    `r pkg("modeest")`, `r pkg("moments")`,
    `r pkg("s20x")`, `r pkg("fromo")`,
    `r pkg("DistributionUtils")`,
    `r pkg("EnvStats")`, `r pkg("parameters")`
    packages.
-   *Empirical kurtosis:* available in `r pkg("agricolae")`,
    `r pkg("DistributionUtils")`,
    `r pkg("e1071")`, `r pkg("EnvStats")`,
    `r pkg("GLDEX")`, `r pkg("HyperbolicDist")`,
    `r pkg("fromo")`, `r pkg("moments")`,
    `r pkg("parameters")` packages. The raw or centered
    moments are provided in `r pkg("e1071")`,
    `r pkg("moments")`.
-   *Empirical L-moments:* L-moments are available in
    `r pkg("lmom")`, `r pkg("lmomco")`,
    `r pkg("Lmoments")`, `r pkg("GLDEX")`,
    `r pkg("EnvStats")`, trimmed L-moments are available in
    `r pkg("lmomco")`, `r pkg("TLMoments")` and
    `r pkg("Lmoments")`, right-censored L-moments are
    available in `r pkg("lmomco")`, and cumulants in
    `r pkg("GLDEX")`. `r pkg("TLMoments")`
    provides a function to convert them to some distribution parameters.
-   *Empirical probability weighted moments:* Probability weighted
    moments are available in `r pkg("EnvStats")` and
    `r pkg("fromo")`.
-   *Empirical cumulants:* `r pkg("fromo")` provides
    centered and standardized cumulants.
-   *Mode estimation:* Package `r pkg("modeest")` provides
    mode computation of known distributions and mode estimation 
    on datasets in the unimodal case.
    Package `r pkg("ModEstM")` provides mode estimation in unimodal 
    and multimodal cases.
    Package `r pkg("multimode")` provides for testing and exploring
    the number of modes on data using non-parametric procedures.
-   *Order statistics:* Distribution function of the jth order statistic
    can be obtained with base R functions.
    `r pkg("orders")` allows to generate samples of k-th order statistics and 
    others quantities of interest for the following distributions:
    Burr, Feller-Pareto, Generalized Pareto, 
    The Inverse Paralogistic, Marshall-Olkin G, 
    exponentiated G, beta G, gamma G, 
    Kumaraswamy G, generalized beta G, beta extended G, 
    gamma G, gamma uniform G, beta exponential G, Weibull G, 
    log gamma G I/II, exponentiated generalized G, 
    exponentiated Kumaraswamy G, geometric exponential Poisson G,
    truncated-exponential skew-symmetric G, modified beta G, 
    exponentiated exponential Poisson G, Poisson-inverse gaussian, 
    Skew normal type 1, Skew student t, Sinh-Arcsinh, 
    Sichel, Zero inflated Poisson. 
-   *Empirical characteristic function:* `r pkg("empichar")`
    evaluates the empirical characteristic function of univariate and
    multivariate samples.
-   *Dispersion index:* Package `r pkg("GWI")` provides
    univariate dispersion index against a particular distribution.    
-   *Theoretical moments:*
    -   *common distributions:* The `r pkg("actuar")`
        package implements raw moments, limited expected values and
        moment generating function for base R distributions.
        `r pkg("lmomco")` provides L-moments (L), trimmed
        L-moments (TL), and right-censored \[RC\] for the following
        distributions: Asymmetric Exponential Power (L), Cauchy (TL),
        Eta-Mu (L), Exponential (L), Gamma (L), Generalized Extreme
        Value (L), Generalized Lambda (L and TL), Generalized Logistic
        (L), Generalized Normal (L), Generalized Pareto (L\[RC\] and
        TL), Govindarajulu (L), Gumbel (L), Kappa (L), Kappa-Mu (L),
        Kumaraswamy (L), Laplace (L), Normal (L), 3-parameter log-Normal
        (L), Pearson Type III (L), Rayleigh (L), Reverse Gumbel
        (L\[RC\]), Rice/Rician (L), Slash (TL), 3-parameter Student T
        (L), Truncated Exponential (L), Wakeby (L), and Weibull (L).
        Multivariate L-moments (L-comoments).
    -   *hyperbolic distributions:*
        `r pkg("HyperbolicDist")` provides the mean,
        variance, skewness, kurtosis, mode, raw and centered moments for
        the hyperbolic, the generalized hyperbolic and the generalized
        inverse Gaussian distributions.
    -   *Lambda distribution:* `r pkg("GLDEX")` also
        provides the mean, variance, skewness, kurtosis of generalized
        Lambda distribution.
    -   *multivariate distributions:* `r pkg("MomTrunc")`
        provides mean vector, covariance matrices and raw moments for
        truncated or folded of the following multivariate distributions:
        normal, skew normal, extended skew normal and student.



# [Random number generators (RNG):]{#Random}


-   *Basic functionality:* R provides several random number generators
    (RNGs). The random seed can be provided via `set.seed` and the kind
    of RNG can be specified using `RNGkind`. The default RNG is the
    Mersenne-Twister algorithm. Other generators include Wichmann-Hill,
    Marsaglia-Multicarry, Super-Duper, Knuth-TAOCP, Knuth-TAOCP-2002, as
    well as user-supplied RNGs. For normal random numbers, the following
    algorithms are available: Kinderman-Ramage, Ahrens-Dieter,
    Box-Muller, Inversion (default). In addition to the tools above,
    `r pkg("setRNG")` provides an easy way to set, retain
    information about the setting, and reset the RNG.
-   *Pseudo-randomness:* `r pkg("RDieHarder")` offers
    several dozen new RNGs from the GNU GSL.
    `r pkg("randtoolbox")` provides more recent RNGs such as
    SF Mersenne-Twister and WELL, which are generators of Mersenne
    Twister type, but with improved quality parameters.
    `r pkg("SuppDists")` implements two RNGs of G.
    Marsaglia. `r pkg("dqrng")` provides PCG family by
    O'Neill (2014) as well as Xoroshiro128+ and Xoshiro256+ by Blackman
    and Vigna (2018).
    -   Support for several independent streams:
        `r pkg("rstream")` focuses on multiple independent
        streams of random numbers from different sources (in an object
        oriented approach). `r pkg("dqrng")` provides RNG
        for parallel computation either in R or in C++.
    -   For non-uniform generation, the `r pkg("Runuran")`
        package interfaces to the UNU.RAN library for universal
        non-uniform generation as well as customised distributions based
        on polynomial interpolation of the inverse cumulative
        distribution function. `r pkg("rust")` performs
        non-uniform random variate generation from unimodal
        (low-dimensional) multivariate continuous distributions, using
        the generalized ratio-of-uniforms method.
        `r pkg("UnivRNG")` provides 17 non-uniform
        generators either using an acceptance/rejection algorithm or the
        inverse CDF method. `r pkg("MultiRNG")` provides 11
        multivariate generators, see each distribution.
        `r pkg("Tinflex")` provides a non-uniform random number generator 
        for quite arbitrary distributions with piecewise twice differentiable densities.
    -   `r pkg("kernelboot")` provides functions for random
        generation from univariate and multivariate kernel densities (in
        particular multivariate Gaussian kernels).
-   *Quasi-randomness:* The `r pkg("randtoolbox")` provides
    the following quasi random sequences: the Sobol sequence, the Halton
    (hence Van Der Corput) sequence and the Torus sequence (also known
    as Kronecker sequence). `r pkg("lhs")` and
    `r pkg("mc2d")` packages implement the latin hypercube
    sampling, an hybrid quasi/pseudo random method.
    `r pkg("sfsmisc")` also provides the Halton sequence.
    `r pkg("qrng")` provides Korobov, generalize Halton and
    Sobol quasi-random sequences.
    `r pkg("spacefillr")` provides Halton and Sobol sequences.
-   *True randomness:* The `r pkg("random")` package
    provides several functions that access the true random number
    service at [random.org](https://www.random.org/) .
-   *RNG tests:* `r pkg("RDieHarder")` offers numerous tests
    of RNGs based on a reimplementation and extension of Marsaglia's
    DieHarder battery. `r pkg("randtoolbox")` provides basic
    RNG tests.
-   *Parallel computing:* Random-number generators for parallel
    computing are available via the `r pkg("rlecuyer")`
    package. See the `r view("HighPerformanceComputing")`
    task view for more details.
-   *Multivariate random vectors:* for parametric multivariate distributions,
    we refer to [Multivariate Continuous](#MultivariateContinuous)
    and [Multivariate Discrete](#MultivariateDiscrete).
    For non-parametric distributions, `r pkg("SimJoint")` offers 
    various to simulate multivariate distributions with non-parametric marginals
    given a Pearson or Spearman correlation matrix.
-   *Unit sphere and other:* `r pkg("simdd")` provides a generator for the Fisher Bingham 
    distribution on the unit sphere, the matrix Bingham distribution on
    a Grassmann manifold, the matrix Fisher distribution on SO(3), and the bivariate 
    von Mises sin model on the torus.
    `r pkg("uniformly")` provides sampling on various geometric shapes, 
    such as spheres, ellipsoids, simplices.
    `r pkg("watson")` allows simulating mixtures of Watson distributions.
    

# [Miscellaneous:]{#Misc}


-   *Computation:*
    -   *Approximation of d, p, q, r functions:*
        `r pkg("PDQutils")` provides tools for computing the
        density, cumulative distribution, and quantile functions of a
        distribution when the cumulants or moments are given, using the
        classical Gram Charlier, Edgeworth and Cornish-Fisher
        approximations. `r pkg("sadists")` is a showcase for
        PDQutils, providing density, cumulative distribution, quantile,
        and random generation for the doubly non-central t, doubly
        non-central F, K-prime, Lambda-prime, Upsilon, and sum of
        (non-central) chi-squares to powers distributions. Various
        approximations and alternative computations for d, p, q
        functions of probability distributions in R are given
        `r pkg("DPQ")`.
    -   For non-uniform generation, see the
        `r pkg("Runuran")` above.

-   *Non parametric models:*
    -   *Binned Empirical distributions:* The
        `r pkg("HistogramTools")` package provides a number
        of methods for manipulating empirical data that has been binned
        into histogram form, including: (1) the empirical cumulative
        distribution function, (2) the empirical quantile, and (3)
        information loss metrics associated with binning.
    -   *Empirical distribution:* Base R provides functions for
        univariate analysis: (1) the empirical density (see
        density()), (2) the empirical cumulative distribution function
        (see ecdf()), (3) the empirical quantile (see quantile())
        and (4) random sampling (see sample()).
        `r pkg("distributionsrd")` provides d, p, q, r
        user-friendly functions for the empirical distributions as well
        as moments. `r pkg("mded")` provides a function for
        measuring the difference between two independent or
        non-independent empirical distributions and returning a
        significance level of the difference.
        `r pkg("MEPDF")` provides functions to compute and
        visualize empirical density functions for multivariate data.
    -   *Non Parametric distributions :* `r pkg("spd")`
        provides the Semi Parametric Piecewise Distribution, while
        `r pkg("fBasics")` implements spline smoothed
        distributions.
-   *Hierarchical models:* Distributions whose some parameters are no
    longer constant but random according to a particular distribution.
    `r pkg("VGAM")` provides a lot of hierarchical models:
    beta/binomial, beta/geometric and beta/normal distributions.
    `r pkg("bayesm")` implements: binary logit, linear,
    multivariate logit and negative binomial models. Furthermore
    `r pkg("LearnBayes")` and
    `r pkg("MCMCpack")` provides poisson/gamma,
    beta/binomial, normal/normal and multinomial/Dirichlet models.
-   *Unified interface to handle distributions:*
    -   *S3 Object-orientation:* `r pkg("distributions3")`
        provides tools to create and to manipulate probability
        distributions using S3, that is
        `r pkg("distributions3")`, generics `random()`,
        `pdf()`, `cdf()` and `quantile()` provide replacements for base
        R's `r/d/p/q` style functions.
        `r pkg("distributional")` also provides tools to
        create and to manipulate probability distributions using S3,
        with `cdf()`, `density()`, `hdr()`, `mean()`, `median()`,
        `quantile()`,\...
    -   *S4 Object-orientation:* General discrete and continuous
        distributions are implemented in package
        `r pkg("distr")` respectively via S4-class
        DiscreteDistribution and AbscontDistribution providing the
        classic d, p, q and r functions. `r pkg("distrEx")`
        extends available distributions to multivariate and conditional
        distributions as well as methods to compute useful statistics
        (expectation, variance,\...) and distances between distributions
        (Hellinger, Kolmogorov,\... distance). Finally package
        `r pkg("distrMod")` provides functions for the
        computation of minimum criterion estimators (maximum likelihood
        and minimum distance estimators). See other packages of the
        distr-family (`r pkg("distrSim")`,
        `r pkg("distrTEst")`,
        `r pkg("distrTeach")`,
        `r pkg("distrDoc")`,
        `r pkg("distrEllipse")`).
    -   *R6 Object-orientation:* 
        `r pkg("ROOPSD")` provides a R6 class interface to
        classic statistical distribution.
    -   *Transformation:* Lebesgue decomposition are implemented in
        `r pkg("distr")`, as well as Convolution, Truncation
        and Huberization of distributions. Furthermore,
        `r pkg("distr")` provides distribution of the
        maximum or minimum of two distributions. See Object-orientation
        above.
-   *Transversal functions:*
    -   *Histogram, tail plots, distance estimation:*
        `r pkg("DistributionUtils")` provides log-histogram,
        tail plots, functions for testing distributions using inversion
        tests and the Massart inequality.
        `r pkg("visualize")` provides functions to plot the
        pdf or pmf with highlights on area or when probability is
        present in user defined locations, as well as the graph is the
        mean and variance of the distribution.
        `r pkg("visualize")` provides lower tail, bounded,
        upper tail, and two tail calculations.
        `r pkg("visualize")` contains convenience functions
        for constructing and plotting bivariate probability
        distributions (probability mass functions, probability density
        functions and cumulative distribution functions).
        `r pkg("vistributions")` provides visualization
        tools for a selected number of distributions.
    -   *Parameter estimation:* `r pkg("lmomco")` and
        `r pkg("Lmoments")` focus on univariate/multivariate
        (L-)moments estimation. `r pkg("VGAM")` provides a
        lot of parameter estimation for usual and "exotic"
        distributions. `r pkg("gaussDiff")` provides a
        collection difference measures for multivariate Gaussian
        probability density functions Package
        `r pkg("MASS")` implements the flexible `fitdistr`
        function for parameter estimations.
        `r pkg("fitdistrplus")` greatly enlarges `fitdistr` and
        enhances the tools to fit a user-supplied probability distribution.
        `r pkg("OneStep")` is based upon `r pkg("fitdistrplus")` to provide
        one-step estimation procedures.
        `r pkg("EnvStats")`, `r pkg("fitteR")`, `r pkg("ExtDist")` also provide
        tools to fit and select a set of probability distributions. `r pkg("flexsurv")` and
        `r pkg("msm")` provides a quantile function for a
        generic distribution based on numerical computation based on a
        dichotomic search.



### Links

-   [Clickable diagram of distribution relationships](http://www.johndcook.com/distribution_chart.html)
-   [Diagram of discrete distribution relationships](http://www.stat.rice.edu/~dobelman/courses/texts/Distributions.Discrete.Kendall.jpg)
-   [Diagram of continuous distribution relationships](http://www.stat.rice.edu/~dobelman/courses/texts/Distributions.Chart.C&B.pdf)
-   [List and diagram of distribution relationship.](http://www.stat.rice.edu/~dobelman/courses/texts/leemis.distributions.2008amstat.pdf)
-   [Compendium of distributions.](http://www.stat.rice.edu/~dobelman/courses/DistributionCompendium.pdf)
-   [A comprehensive list of data types](https://en.wikipedia.org/wiki/Statistical_data_type)
-   [Journal of Statistical Software: R programs for truncated distributions](https://www.jstatsoft.org/v16/c02/)


