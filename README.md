# BayesRegDTR

Optimal Dynamic Treatment Regime Inference using Bayesian Likelihood-Based 
Regression Method

# About this package
This package provides methods to estimate optimal dynamic treatment regimes 
using Bayesian likelihood-based regression approach as described in  Yu, W., & 
Bondell, H. D. (2023). Uses backward induction and dynamic programming theory 
for computing expected values. Offers options for future parallel computing.

------------------------------------------------------------------------

# Pre-installation instructions (Mac Users Only)

To install this package in Mac requires a Fortran compiler (through its
RcppArmadillo dependency). Chances are, your current Fortran compiler is not
up-to-date. To update your Fortran compiler, simply follow the steps
here: <br />  

1.  In your Mac App Store, search "Xcode" and install. <br />
2.  Open Terminal application. Type in

``` {eval="FALSE"}
xcode-select --install
```

     and follow the instructions.<br />      3. Click on the link
[here](https://github.com/fxcoudert/gfortran-for-macOS/releases).
Download the gfortan dmg file according to your MacOS version. <br />  
   4. Open the dmg file, run the gfortran installer, follow all the
instructions.

An alternative recommended method is to use the packet manager
[Homebrew](https://docs.brew.sh/Installation):

  1. Check if you have homebrew with

``` {eval="FALSE"}
$ brew doctor
```

     If you don't have it installed, use the following code from the
Homebrew webiste. Check the website that it hasn't changed since. It
will ask for your user password (you won't see characters as you type).
Follow the instructions.

``` {eval="FALSE"}
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

     2. Install GFortran using gcc (contains GFortran).

``` {eval="FALSE"}
brew install gcc
```

------------------------------------------------------------------------

# Installation

```{r}
# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("jlimrasc/BayesRegDTR")
```

------------------------------------------------------------------------

# Toy example

```{r}
library(BayesRegDTR)
# -----------------------------
# Set Up Parallelism & Progress Bar
# -----------------------------
progressr::handlers("cli")          # Set handler to something with title/text
numCores <- parallel::detectCores() # Detect number of cores, use max
future::plan(future::multisession,  # Or plan(multicore, workers) on Unix
         workers = numCores)        # Set number of cores to use
doFuture::registerDoFuture()        # Or doParallel::registerDoParallel()
                                    # if no progress bar is needed and future
                                    # is unwanted

## UVT
# -----------------------------
# Initialise Inputs
# -----------------------------
num_stages  <- 5
t           <- 3
p_list      <- rep(1, num_stages)
num_treats  <- rep(2, num_stages)
n.train     <- 5000
n.pred      <- 10
# -----------------------------
# Generate Dataset
# -----------------------------
Dat.train  <- generate_dataset(n.train,  num_stages, p_list, num_treats)
Dat.pred  <- generate_dataset(n.pred,  num_stages, p_list, num_treats)
Dat.pred  <- Dat.pred[-1]
Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]

# -----------------------------
# Main
# -----------------------------
gcv_uvt <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred,
                                 num_stages, num_treats,
                                 p_list, t, R = 30,
                                 tau = 0.01, B = 500, nu0 = NULL,
                                 V0 = NULL, alph = 3, gam = 4)

## MVT
# -----------------------------
# Initialise Inputs
# -----------------------------
num_stages  <- 3
t           <- 2
p_list      <- rep(2, num_stages)
num_treats  <- rep(2, num_stages)
n.train     <- 5000
n.pred      <- 10

# -----------------------------
# Generate Dataset
# -----------------------------
Dat.train <- generate_dataset_mvt(n.train, num_stages, p_list, num_treats)
Dat.pred  <- generate_dataset_mvt(n.pred,  num_stages, p_list, num_treats)
Dat.pred  <- Dat.pred[-1]
Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]

# -----------------------------
# Main
# -----------------------------
gcv_res <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred,
                                 num_stages, num_treats,
                                 p_list, t, R = 30,
                                 tau = 0.01, B = 500, nu0 = 3,
                                 V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                                 alph = 3, gam = 4)
```

------------------------------------------------------------------------
