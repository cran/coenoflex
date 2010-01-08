\name{plot.coenoflex}
\alias{plot.coenoflex}
\title{Plotting Coenoflex Results}
\description{Plots the various components of a coenoflex simulation results}

\usage{
\method{plot}{coenoflex}(x,which='all',\dots)
}

\arguments{
\item{x}{an object of class \sQuote{coenoflex}}
\item{which}{a switch to control plotiing just \sQuote{plots}, \sQuote{species},
            \sQuote{amplitude}, or \sQuote{all} attributes}
\item{\dots}{ancillary arguments for the plot generic function}
}

\details{coenoflex generates a complex object of class \sQuote{coenoflex} with
simulated vegetation samples and species response curves.  \code{plot.coenoflex}
graphs the results in a series of panels depending on the argument
\sQuote{which}.
}

\value{plot.coenoflex produces only graphical out and does not return anything.}

\author{
David W. Roberts
\email{droberts@montana.edu}
}

\examples{
example <- coenoflex(numgrd=3,numplt=100,numspc=100,grdtyp=c('e','e','e'),
              grdlen=c(300,200,100),width=c(100,75,50),variab=c(20,20,20),
              grdprd=c(0,0,0),alphad=c(1,1,1),pdist='r',sdist='r',
              skew=3.0,aacorr=0.0,cmpasy=3.0,cmpphy=0.0,maxtot=100,
              noise=20,slack=0.2,autlin='irm(1,2,3)')
\dontrun{plot(example,which='all')}
}

\keyword{hplot}

