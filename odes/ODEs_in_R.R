#' \documentclass[pdftex,12pt,final,reqno]{article}
#' \usepackage{times}
#' \usepackage[utf8]{inputenc}
#' \usepackage[round]{natbib}
#' \usepackage{paralist}
#' \usepackage{longtable}
#' \usepackage[colorlinks=true]{hyperref}
#' \hypersetup{pdftitle={ODEs in R}}
#' \usepackage{amsmath}
#' \usepackage{amsthm}
#' \usepackage{amsfonts}
#' 
#' \setlength{\textwidth}{6.25in}
#' \setlength{\textheight}{8.75in}
#' \setlength{\evensidemargin}{0in}
#' \setlength{\oddsidemargin}{0in}
#' \setlength{\topmargin}{-0.35in}
#' \setlength{\parskip}{0.1in}
#' \setlength{\parindent}{0.0in}
#' \setcounter{secnumdepth}{1}
#' \setcounter{tocdepth}{2}
#' 
#' \theoremstyle{definition}
#' \newtheorem{exercise}{Exercise}
#' \newtheorem{challenge}[exercise]{*Exercise}
#' \theoremstyle{remark}
#' \newtheorem*{solution}{Solution}
#' 
#' \newcommand{\code}[1]{\texttt{#1}}
#' \newcommand{\pkg}[1]{\textsf{#1}}
#' \newcommand{\Prob}[1]{\mathbb{P}\left[#1\right]}
#' \newcommand{\expect}[1]{\mathbb{E}\left[#1\right]}
#' \newcommand{\var}[1]{\mathrm{Var}\left[#1\right]}
#' 
#' %%\newcommand<>{\emphcol}[1]{\textcolor#2{red!50!blue}{#1}}
#' \newcommand\scinot[2]{$#1 \times 10^{#2}$}
#' 
#' %%%%%%%%%%%%%%%%%%
#' 
#' \title{Integrating ordinary differential equations in \pkg{R}}
#' 
#' \author{Aaron A. King\\
#' with contributions from\\
#' Ben Bolker, John Drake, Pej Rohani, and Dave Smith}
#' 
#' \date{\today\\
#' \vspace{1em}
#' Licensed under the Creative Commons attribution-noncommercial license, \texttt{http://creativecommons.org/licenses/by-nc/3.0/}.
#' Please share and remix noncommercially, mentioning its origin. \parbox[bottom]{57pt}{\includegraphics[height=20pt]{cc-by-nc}}}
#' 
#' 
#' \begin{document}
#' 
#' \maketitle
#' \thispagestyle{empty}
#' 
#' \section{Introduction}
#' 
#' Here we begin our study of computational techniques for studying epidemiological models.
#' In this session we introduce the numerical solution (or integration) of nonlinear differential equations using the sophisticated solvers found in the package \pkg{deSolve}.
#' Numerical integration is one of the most important tools we have for the analysis of epidemiological models.
#' 
#' \section{The SIR model}
#' 
#' As we saw in the lecture, the classical SIR compartmental model divides a population of hosts into three classes: susceptible, infected, recovered.
#' The model describes how the fraction of a population in each of these classes changes with time.
#' Alternatively, the model can track the number of individuals in each class.
#' Births are modeled as flows from ``nowhere'' into the susceptible class; deaths are modeled as flows from the S, I, or R compartment into ``nowhere''.
#' If $S$, $I$, and $R$ refer to the fractions of indivduals in each compartment, then these \emph{state variables} change according to the following system of differential equations:
#' \begin{equation*}
#'   \begin{aligned}
#'     \frac{dS}{dt} &= B-\lambda\,S-\mu\,S\\
#'     \frac{dI}{dt} &= \lambda\,S-\gamma\,I-\mu\,I\\
#'     \frac{dR}{dt} &= \gamma\,I-\mu\,R\\
#'   \end{aligned}
#' \end{equation*}
#' Here, $B$ is the crude birth rate (births per unit time), $\mu$ is the death rate and $\gamma$ is the recovery rate.
#' We'll assume that the force of infection, $\lambda$, has the form
#' \begin{equation*}
#'   \lambda = \beta\,I
#' \end{equation*}
#' so that the risk of infection a susceptible faces is proportional to the \emph{prevalence} (the fraction of the population that is infected).
#' This is known as the assumption of frequency-dependent transmission.
#' %  Notice that we allow for the possibility of a contact rate, $\beta(t)$, that varies in time.
#' 
#' \section{Solving ODEs in R}
#' 
#' Like almost all epidemiological models, one can't solve these equations analytically.
#' However, we can compute the \emph{trajectories} of a continuous-time model such as this one by integrating the equations numerically.
#' Doing this accurately involves a lot of calculation, and there are smart ways and not-so-smart ways of going about it.
#' This very common problem has been very thoroughly studied by numerical analysts for generations so that, when the equations are smooth, well-behaved functions, excellent numerical integration algorithms are readily available to compute approximate solutions to high precision.
#' In particular, \pkg{R} has several sophisticated ODE solvers which (for many problems) will give highly accurate solutions.
#' These algorithms are flexible, automatically perform checks, and give informative errors and warnings.
#' To use the numerical differential equation solver package, we load the \pkg{deSolve} package
## ------------------------------------------------------------------------

require(deSolve)


#' [Here, you may get an error saying that the \pkg{deSolve} package is not installed.
#'   If you do, run the following command:
## ----eval=F--------------------------------------------------------------
## install.packages("deSolve")

#' ]
#' The ODE solver we'll use is called \code{ode}.
#' Let's have a look at the help page for this function.
## ----eval=F--------------------------------------------------------------
## ?ode

#' \code{ode} needs to know the \emph{initial values} of the state variables (\code{y}), the \code{times} at which we want solutions, the right-hand side of the ODE \code{func}.
#' The latter can optionally depend on some parameters (\code{parms}).
#' 
#' \section{SIR for a closed epidemic}
#' 
#' Let's study the SIR model for a closed population, i.e., one in which we can neglect births and deaths.
#' Recall that the differential equations for the closed epidemic are
#' \begin{equation*}
#'   \begin{aligned}
#'     \frac{dS}{dt} &= -\beta\,S\,I\\
#'     \frac{dI}{dt} &= \beta\,S\,I-\gamma\,I\\
#'     \frac{dR}{dt} &= \gamma\,I\\
#'   \end{aligned}
#' \end{equation*}
#' To encode these equations in a form suitable for use as the \code{func} argument to \code{ode}, we'll need to write a function.
#' For example:
## ----closed-sir-model-defn-----------------------------------------------

closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  ## now code the model equations
  dSdt <- -beta*S*I
  dIdt <- beta*S*I-gamma*I
  dRdt <- gamma*I
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}


#' Note that the order and type of the arguments and output of this function must exactly match \code{ode}'s expectations.
#' Thus, for instance, the time variable \code{t} must be the first argument even if, as is the case here, nothing in the function depends on time.
#' [When the RHS of the ODE are independent of time, we say the ODE are \emph{autonomous}.]
#' Note also, that \code{ode} expects the values of the ODE RHS to be the first element of a \code{list}.
#' 
#' Now we can call \code{ode} to compute trajectories of the model.
#' To do this, we'll need some values of the parameters.
#' If we're thinking of a disease something like measles, and measuring time in years, we might use something like:
## ----set-closed-params---------------------------------------------------
params <- c(beta=400,gamma=365/13)

#' \begin{exercise}
#'   What is the infectious period of this disease?
#'   What is $R_0$ in this case?
#' \end{exercise}
#' 
#' We now state the times at which we want solutions and specify the \emph{initial conditions},
#' i.e., the starting values of the state variables $S$, $I$, and $R$:
## ----set-times-ics-------------------------------------------------------
times <- seq(from=0,to=60/365,by=1/365/4) # returns a sequence
xstart <- c(S=0.999,I=0.001,R=0.000)     # initial conditions

#' 
#' Next, we compute a model trajectory with the \code{ode} command and store the result in a data-frame:
## ----solve-closed-sir----------------------------------------------------
out <- as.data.frame(
                     ode(
                         func=closed.sir.model,
                         y=xstart,
                         times=times,
                         parms=params
                         )
                     )

#' and plot the results in Figure~\ref{fig:closed-epi} using the commands:
## ----epi-curve-plot,eval=F-----------------------------------------------
## plot(I~time,data=out,type='l')

#' 
#' \begin{figure}
#'   \begin{center}
## ----fig.height=3,fig.width=6,echo=F-------------------------------------
op <- par(mgp=c(2,1,0),mar=c(3,3,1,1))
plot(I~time,data=out,type='l')
par(op)

#'   \end{center}
#'   \caption{
#'     Trajectory of the SIR model of an epidemic in a closed population.
#'     $I$ is the fraction of the population infected.
#'   }
#'   \label{fig:closed-epi}
#' \end{figure}
#' 
#' \begin{exercise}
#'   Suppose that you'd rather measure time in days.
#'   Modify the parameters accordingly and verify your modifications.
#' \end{exercise}
#' 
#' \clearpage
#' 
#' Let's study how the epidemic curve depends on the transmission rate, $\beta$, and the infectious period.
#' In particular, we'll investigate how the epidemic curve changes as we vary $\beta$ from 20 to 500 and the infectious period from 5 to 30 days.
#' 
## ----nine-curves---------------------------------------------------------

betavals <- c(20,50,500)
ips <- c(5,10,30)
gammavals <- 365/ips

## set some plot parameters
op <- par(mgp=c(2,1,0),mar=c(3,3,1,1),mfrow=c(3,3))

for (beta in betavals) {
  for (gamma in gammavals) {
    params <- c(beta=beta,gamma=gamma)
    out <- as.data.frame(
                         ode(
                             func=closed.sir.model,
                             y=xstart,
                             times=times,
                             parms=params
                             )
                         )
    title <- bquote(list(beta==.(beta),"IP"==.(365/gamma)~"da"))
    plot(I~time,data=out,type='l',main=title)
  }
}
par(op)      # restore old settings


#' 
#' Simulation is a useful tool, but its power is limited.
#' The next exercise demonstrates the importance of being able to analyze the equations as well.
#' 
#' \begin{exercise}
#'   For each of the above parameter combinations, describe the system's behavior.
#'   Compute $R_0$ for each parameter combination and relate it to the behavior of the system.
#' \end{exercise}
#' 
#' \begin{exercise}
#'   Use the ODE solver to study the dependence of the epidemic's \emph{final size} on $R_0$.
#'   Compare your results with the predictions of the final size equation
#'   \begin{equation*}
#'     1-R(\infty)=S(0)\,e^{-R(\infty)\,R_0}=e^{-R(\infty)\,R_0}
#'   \end{equation*}
#'   solutions of which are plotted in Fig.~\ref{fig:finalsize}.
#' \end{exercise}
#' 
#' \begin{figure}
## ----finalsize,echo=F----------------------------------------------------
## don't worry about executing the following
## it just generates the figure in the document,
## which is all you need to do the exercise
require(emdbook)
finalsize <- function (R0) {1+lambertW(-R0*exp(-R0))/R0}
curve(finalsize,from=1,to=5,xlab=expression(R[0]),
      ylab=expression(R(infinity)),bty='l')

#'   \caption{The final size, $R(\infty)$, of an SIR epidemic depends only on $R_0$.}
#'   \label{fig:finalsize}
#' \end{figure}
#' 
#' \begin{solution}
#' Here's one code that will solve the exercise.
#' Others are certainly possible.
## ----finalsize-soln------------------------------------------------------

times <- seq(from=0,to=20,by=1/12) # you might think 20 yr would be enough!
xstart <- c(S=0.999,I=0.001,R=0.000)
R0vals <- seq(0,5,by=0.1)
gamma <- 365/13
betavals <- R0vals*gamma                # R0 = beta/gamma
finalsize <- numeric(length(R0vals))    # a vector to hold the solutions
for (k in seq_along(betavals)) {
    params <- c(beta=betavals[k],gamma=gamma)
    out <- as.data.frame(
                         ode(
                             func=closed.sir.model,
                             y=xstart,
                             times=times,
                             parms=params
                             )
                         )
    finalsize[k] <- tail(out$R,1)       # the final value of R
}
plot(finalsize~R0vals,type='o')


#' \end{solution}
#' 
#' 
#' \clearpage
#' \section{SIR dynamics in an open population}
#' 
#' Over a sufficiently short time scale, the assumption that the population is closed is reasonable.
#' To capture the dynamics over the longer term, we'll need to account for births and deaths, i.e., allow the population to be an \emph{open} one.
#' As we've seen, if we further assume that the birth rate equals the death rate, then the SIR equations become
#' \begin{equation*}
#'     \begin{aligned}
#'       \frac{dS}{dt} &= \mu -\beta\,S\,I-\mu\,S\\
#'       \frac{dI}{dt} &= \beta\,S\,I-\gamma\,I-\mu\,I\\
#'       \frac{dR}{dt} &= \gamma\,I-\mu\,R\\
#'     \end{aligned}
#' \end{equation*}
#' 
#' We must modify the ODE function accordingly:
## ----open-sir-model-defn-------------------------------------------------

open.sir.model <- function (t, x, params) {
  beta <- params["beta"]
  mu <- params["mu"]
  gamma <- params["gamma"]
  dSdt <- mu*(1-x[1])-beta*x[1]*x[2]
  dIdt <- beta*x[1]*x[2]-(mu+gamma)*x[2]
  dRdt <- gamma*x[2]-mu*x[3]
  list(c(dSdt,dIdt,dRdt))
}


#' 
#' \clearpage
#' We'll need to specify a birth/death rate in addition to the two parameters we specified before:
## ----set-open-params-----------------------------------------------------

params <- c(mu=1/50,beta=400,gamma=365/13)


#' We integrate the equations as before:
## ----solve-open-sir------------------------------------------------------

times <- seq(from=0,to=25,by=1/365)
out <- as.data.frame(
                     ode(
                         func=open.sir.model,
                         y=xstart,
                         times=times,
                         parms=params
                         )
                     )


#' 
#' We can plot each of the state variables against time, and $I$ against $S$, as shown in Fig.~\ref{fig:open-epi}, using the commands:
## ----open-epi-plot,eval=F------------------------------------------------
## 
## op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
##           mar=c(3,3,1,1),mgp=c(2,1,0))
## plot(S~time,data=out,type='l',log='y')
## plot(I~time,data=out,type='l',log='y')
## plot(R~time,data=out,type='l',log='y')
## plot(I~S,data=out,log='xy',pch='.',cex=0.5)
## par(op)
## 

#' 
#' \begin{figure}
#'   \begin{center}
## ----fig.height=6,fig.width=6,echo=F-------------------------------------

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,1,1),mgp=c(2,1,0))
plot(S~time,data=out,type='l',log='y')
plot(I~time,data=out,type='l',log='y')
plot(R~time,data=out,type='l',log='y')
plot(I~S,data=out,log='xy',pch='.',cex=0.5)
par(op)


#'   \end{center}
#'   \caption{
#'     Trajectory of the SIR model.
#'     Each state variable is plotted against time, and $I$ is also plotted against $S$.
#'   }
#'   \label{fig:open-epi}
#' \end{figure}
#' 
#' \clearpage
#' 
#' \begin{exercise}
#'   Explore the dynamics of the system for different values of the $\beta$ and $\gamma$ parameters by simulating and plotting trajectories as time series and in phase space (e.g., $I$ vs. $S$).
#'   Use the same values of $\beta$ and $\gamma$ we looked at above.
#'   How does the value of $R_0$ affect the results?
#' \end{exercise}
#' 
#' \begin{exercise}
#'   Under the assumptions of this model, the average host lifespan is $1/\mu$.
#'   Explore how host lifespan affects the dynamics by integrating the differential equations for lifespans of 20 and 200 years.
#' \end{exercise}
#' 
#' The compartmental modeling strategy can be put to use in modeling a tremendous range of infections.
#' The following exercises make some first steps in this direction.
#' 
#' \begin{challenge}
#'   The SIR model assumes lifelong sterilizing immunity following infection.
#'   For many infections, immunity is not permanent.
#'   Make a compartment diagram for an SIRS model, in which individuals lose their immunity after some time.
#'   Write the corresponding differential equations and modify the above codes to study its dynamics.
#'   Compare the SIR and SIRS dynamics for the parameters $\mu=1/50$, $\gamma=365/13$, $\beta=400$ and assuming that, in the SIRS model, immunity lasts for 10 years.
#' \end{challenge}
#' 
#' \begin{challenge}
#'   Make a diagram, write the equations, and study the dynamics of the SEIR model for the dynamics of an infection with a latent period.
#'   Compare the dynamics of SIR and SEIR models for the parameters $\mu=1/50$, $\gamma=365/5$, $\beta=1000$ and assuming that, in the SEIR model, the latent period has duration 8 days.
#' \end{challenge}
#' 
#' \section{Nonautonomous equations}
#' 
#' \subsection*{SIR with seasonal transmission}
#' 
#' The simple SIR model always predicts damped oscillations towards an equilibrium (or pathogen extinction if $R_0$ is too small).
#' This is at odds with the recurrent outbreaks seen in many real pathogens.
#' Sustained oscillations require some additional drivers in the model.
#' An important driver in childhood infections of humans (e.g., measles) is seasonality in contact rates because of aggregation of children the during school term.
#' We can analyze the consequences of this by assuming sinusoidal forcing on $\beta$ according to $\beta(t)=\beta_0\,(1+\beta_1\cos(2\,\pi\,t))$.
#' We can modify the code presented above to solve the equations for a seasonally forced epidemic.
## ----seas-sir,cache=T,fig.height=10--------------------------------------

seas.sir.model <- function (t, x, params) {
  beta0 <- params["beta0"]
  beta1 <- params["beta1"]
  mu <- params["mu"]
  gamma <- params["gamma"]
  beta <- beta0*(1+beta1*cos(2*pi*t))
  dSdt <- mu*(1-x[1])-beta*x[1]*x[2]
  dIdt <- beta*x[1]*x[2]-(mu+gamma)*x[2]
  dRdt <- gamma*x[2]-mu*x[3]
  list(c(dSdt,dIdt,dRdt))
}

params <- c(mu=1/50,beta0=400,beta1=0.15,gamma=365/13)
xstart <- c(S=0.07,I=0.00039,R=0.92961)
times <- seq(from=0,to=30,by=7/365)
out <- as.data.frame(
                     ode(
                         func=seas.sir.model,
                         y=xstart,
                         times=times,
                         parms=params
                         )
                     )

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,1,1),mgp=c(2,1,0))
plot(S~time,data=out,type='l',log='y')
plot(I~time,data=out,type='l',log='y')
plot(R~time,data=out,type='l',log='y')
plot(I~S,data=out,log='xy',pch='.',cex=0.5)
par(op)


#' 
#' \begin{exercise}
#'   Explore the dynamics of the seasonally forced SIR model for increasing amplitude $\beta_1$.
#'   Be sure to distinguish between transient and asymptotic dynamics.
#' \end{exercise}
#' 
#' \subsection*{Forcing with a covariate}
#' 
#' When a covariate forces the equations, we must interpolate the covariate.
#' To give an example, let's suppose that the transmission rate depends on rainfall, $R(t)$, and that we have data on rainfall (in mm/mo).
## ----dacca-rain,cache=T--------------------------------------------------
rain <- read.csv("https://kingaa.github.io/thid/data/dacca_rainfall.csv")
rain$time <- with(rain,year+(month-1)/12)

plot(rainfall~time,data=rain,type='l')

#' 
#' Let's assume that transmission depends on rainfall according to
#' \begin{equation*}
#'   \log\beta(t) = \frac{a\,R(t)}{b+R(t)}
#' \end{equation*}
#' Since the data are accumulated monthly rainfall figures but the ODE integrator will need to evaluate $R(t)$ at arbitrary times, we'll need some way of interpolating the rainfall data.
#' \pkg{R} affords us numerous ways of doing this.
#' The \pkg{R} functions \code{approxfun} and \code{splinefun} construct interpolating functions in different ways;
#' see the documentation on these functions.
## ----spline-force--------------------------------------------------------

interpol <- with(rain,approxfun(time,rainfall,rule=2,method='constant'))

data.frame(time=seq(from=1920,to=1930,by=1/365)) -> smoothed.rain
smoothed.rain$rainfall <- sapply(smoothed.rain$time,interpol)

plot(rainfall~time,data=rain,col='black',cex=2,pch=16,log='')
lines(rainfall~time,data=smoothed.rain,col='red')


#' 
## ----solve-rain-forced,cache=T,fig.height=10-----------------------------

rain.sir.model <- function (t, x, params) {
  a <- params["a"]
  b <- params["b"]
  mu <- params["mu"]
  gamma <- params["gamma"]
  R <- interpol(t)
  beta <- a*R/(b+R)
  dSdt <- mu*(1-x[1])-beta*x[1]*x[2]
  dIdt <- beta*x[1]*x[2]-(mu+gamma)*x[2]
  dRdt <- gamma*x[2]-mu*x[3]
  list(c(dSdt,dIdt,dRdt))
}

params <- c(a=500,b=50,mu=1/50,gamma=365/13)
xstart <- c(S=0.07,I=0.00039,R=0.92961)
times <- seq(from=1920,to=1930,by=7/365)
out <- as.data.frame(
                     ode(
                         func=rain.sir.model,
                         y=xstart,
                         times=times,
                         parms=params
                         )
                     )

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,1,1),mgp=c(2,1,0))
plot(S~time,data=out,type='l',log='y')
plot(I~time,data=out,type='l',log='y')
plot(R~time,data=out,type='l',log='y')
plot(I~S,data=out,log='xy',pch='.',cex=1)
par(op)


#' 
#' More generally, any fitting method that has an associated \code{predict} method can be used.
#' For example, we can use local polynomial regression, \code{loess}, to smooth the rainfall data.
#' Do \verb+?loess+ and \verb+?predict.loess+ for more details.
## ----cov-force-----------------------------------------------------------

fit <- loess(log1p(rainfall)~time,data=rain,span=0.05)

data.frame(time=seq(from=1920,to=1930,by=1/365)) -> smoothed.rain
smoothed.rain$rainfall <- expm1(predict(fit,newdata=smoothed.rain))

plot(rainfall~time,data=rain,col='black',cex=2,pch=16,log='')
lines(rainfall~time,data=smoothed.rain,col='red')


#' 
#' 
#' \end{document}
#' 
#' 
#' 
#' \clearpage
#' \bibliographystyle{ecology}
#' \bibliography{biblios}
#' 
#' 
#' Second lab:
#' * linear regression on initial segment of flu data to estimate R0
#' * exercise: try this for 4, 5, 6 points, look at both R0hat and se(R0)
#' * least squares fit of closed epidemic model to data
#' * one parameter at a time for fixed values of the others
#' * again on log scale
#' *
