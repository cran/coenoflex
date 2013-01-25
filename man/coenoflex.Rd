\name{coenoflex}
\alias{coenoflex}
\title{Coenoflex - Flexible Coenospace Simulation}
\description{
Coenoflex simulates complex vegetation data sets.  Species have modal
physiological response curves along gradients, which can be either complex or
direct environmental gradients, or resource gradients sensu Austin and Smith
1989.  The integrated response of species can be a nested response of specific
gradient response functions.  Simulated samples can be placed at random or on a
grid on up to 10 dimensions; species modes can be placed at random or on a grid.
}
\usage{coenoflex(numgrd,numplt,numspc,grdtyp,grdlen,width,variab,grdprd,
                 alphad,pdist,sdist,skew,aacorr,cmpasy,cmpphy,
                 maxtot,noise,slack,autlin)}
\arguments{
    \item{numgrd}{number of gradients to simulate}
    \item{numplt}{number of sample plots}
    \item{numspc}{number of species to simuulate}
    \item{grdtyp}{a vector of gradient type equal in length to the number of
gradients}
    \item{grdlen}{a vector of gradient length equal in length to the number of 
gradients}
    \item{width}{a vector of species ecological amplitude width equal in length to the number of
gradients}
    \item{variab}{a vector of gradient variability equal in length to the number of
gradients}
    \item{grdprd}{a vector of gradient productivity equal in length to the
number of gradients}
    \item{alphad}{a vector of alpha-diversities equal in length to the
number of gradients}
    \item{pdist}{sample plot distribution, \sQuote{r} for random, \sQuote{g} for
grid}
    \item{sdist}{species distribution, \sQuote{r} for random, \sQuote{g} for
grid}
    \item{skew}{the skewness in the pseudo-log normal distribution of species
maximum abundances}
    \item{aacorr}{the correlation between maximum abundance and amplitude width}
    \item{cmpasy}{competition asymmetry}
    \item{cmpphy}{physiological competition coefficient, currently not used}
    \item{maxtot}{the maximum total species cover in sample plots}
    \item{noise}{simulated sample error on species abundance}
    \item{slack}{the probability of species absence in suitable sites}
    \item{autlin}{the autecological function specifying how to combine
gradient-specific response functions}
}

\details{Coenoflex is a complex function to simulate vegetation and environment
data sets.  It is somewhat unusual as an R function, in that it can be run in
classical command line prompt and response mode by calling the coenoflex function with no
arguments.  Any prompt can be replied to with \sQuote{?} to invoke user help for
that parameter.  Alternatively, it can be called in classical R mode by
specifying all parameters in the function call.  Surprisingly, there are no 
defaults for parameters, so if you fail to specify a value for a parameter 
you will still be prompted for a value, even if you specified most values
in the function call.

Coenoflex has many parameters, many of which are vector-valued with lengths
equal to the number of gradients simulated.  Accordingly, calling coenoflex can
be a long-winded affair.  In more detail, the parameters are:

\sQuote{numgrd} The number of gradients to simulate.  \sQuote{Gradients} are
central to modern vegetation theory, and are used to describe environmental
variability in a stylized manner, similar to the role of fields in physics or
vector calculus.

\sQuote{numplt} The number of sample plots to simulate.  Sample plots can be
visualized as points in a high dimensional space defined by the gradients.
Accordingly, as the number of gradients simulated goes up, the number of sample
points should also go up exponentially to maintain a reasonable density of
points.

\sQuote{numspc} The number of species to simulate.  The sample plot species
richness is a function of the number of species simulated, their amplitude
compared to gradient length, and the number of gradients simulated.  As the mean
species amplitude (or width) gets smaller relative to gradient length, the
species richness goes down.  As the number of gradients simulated goes up, 
fewer species will occcur in any given plot due to the curse of dimensionality.
Again, as the number of gradients goes up the number of species should go up 
exponentially.

\sQuote{grdtyp} Gradient type, with one value specified for each gradient.
Gradient types can be either \sQuote{e} for
\sQuote{environmental} or \sQuote{r} for \sQuote{resource}.  On environmental
gradients, species have modal distributions defined by a minimum, optimum, and
maximum gradient values.  None of these values need occur within the range of
the gradient values simulated, but in general species exhibit such a modal curve.  On
resource gradients, all species have their maximim value at the end of the
gradient, with grid-based or random minimum values.  See Austin and Smith (1989)
for details. 

\sQuote{grdlen} Gradient length, with one value specified for each gradient.
Coenoflex calculates gradient lengths in arbitrary units.
These units are used to compare gradient lengths in
relative values, so that for example a gradient of 200
units is twice as long as a gradient of 100 units.
In addition, species amplitudes are measured in the same
units so that a species with an amplitude of 50 would
extend about one half the length of a gradient of 100
units.  

\sQuote{width} Species mean ecological amplitude width, with one value specified for each gradient.
Coenoflex calculates the ecological amplitude of species in
gradient units.  Thus, if the gradient length is 100, and
a species amplitude is 50, it will extend about one half
the length of the gradient.  Species amplitude, along with
the number of species, determines the number of species
per plot and the rate of species turnover.

\sQuote{variab} Species amplitude variability in percent, with one value for each gradient.
Coenoflex varies individual species amplitudes according to
the variance parameter.  Previously you specified a mean
species amplitude along each gradient.  The variance
parameter modifies that value for each species by +/- \%.
For example, if you specified a mean amplitude of 50 and a
variance of 50\%, species amplitudes would vary randomly
from 25 units to 75 units.

\sQuote{grdprd} Gradient productivity, with one value for each gradient.
Coenoflex can vary the productivity along gradients
systematically.  To increase productivity along the
gradient, enter the percent increase in productivity
from the low to high end of the gradient.  Resource
gradients should have a significant increase associated
with them; environmental gradients may or may not.

\sQuote{alphad} Alpha diversity, with one value for each gradient.
Coenoflex can systematically vary the alpha-diversity along
gradients by modifying the location of species centroids.
Values from 0.5 to 1.0 increase the alpha-diversity at the
high end of the gradient, and values from 1.0 to 2.0
increase alpha-diversity at the low end of the gradient.
Values near 1.0 have small effects; the farther from 1.0
the more dramatic the effect.

\sQuote{pdist} Sample plot distribution.
Coenoflex can establish plot centers randomly or along
a grid.  If you select randomly, plot centers will have
random coordinates along each gradient with approximately
constant density throughout the sample space.
If you select a grid, Coenoflex will attempt to space plots
along the gradients in proportion to their length, again
attempting to maintain approximately constant density.
If necessary, Coenoflex will change the number of plots
specified to achieve equal spacing.

\sQuote{sdist} Species centroid distribution.
Coenoflex can establish species modes randomly or along
a grid.  If you select randomly, species modes will have
random coordinates along each gradient with approximately
constant density throughout the sample space.  Many species 
will have their centroid outside the sample space, depending
on the ratio of species amplitude to gradient length.
If you select a grid, Coenoflex will attempt to space
species modes along the gradients in proportion to their
length, again attempting to maintain approximately
constant density.  If necessary, Coenoflex will change
the number of species specified to achieve equal spacing.

\sQuote{skew} Species maximum abundance distribution skew.
Coenoflex simulates the maximum abundance of organisms
according an equation similar to a log-random
distribution.  The skew coefficient controls how skewed
the distribution is toward low values.  A skew of 1.0
results in a normal distribution with a mean of 50\% cover.
As the value increases, the mode moves left and the tail to
the right grows longer.

\sQuote{aacorr} Correlation of species amplitude and maximum abundance.  
Coenoflex allows species abundances to be correlated
with their ecological amplitude as proposed in the Core-Satelite 
hypothesis and described in the 
Hierarchical Continuum Model (Collins et al. 1993).
Values greater than 0.0 result in positive correlation
between amplitude and abundance while 0.0 results
in independence.

\sQuote{cmpasy} Competition asymmetry.
This is the primary mechanism to simulate
competition.  The competition asymmetry coefficient
simulates a competitive advantage to larger plants
intended to simulate asymmetric competition such as
for light.  The algorithm scales species abundance
proportional to raw abundance raised to the competition
asymmetry coefficient.  Accordingly, a value of 1.0
confers no advantage, and increasingly larger values
give large species more advantage.  

\sQuote{maxtot} Maximum total abundance in sample plots.
The typical total abundance of sample plots is normally a simple
sum of the abundances of species in that sample.  It's possible 
to standardize the value to any specified total, including the 
commonly used 100\%.  The \sQuote{maxtot} parameter interacts with
the competition asymmetry parameter by specifying the total abundance
to be partitioned among species.  If no total is specified, or a total
of zero is specified, no competition occurs.

\sQuote{noise} Noise from sample error.
Coenoflex allows you to simulate noise in the species
response model.  Values are +/- in percent.  For example,
a response of 25 randomly adds plus or minus 25\% of the actual
value.

\sQuote{slack} Probability of absence in suitable sites.
Coenoflex includes a form of noise known as slack.  In
contrast to noise, slack causes a species to be absent
in a favorable environment, due perhaps to poor dispersal
(or incorrect identification in the field).  The values
are probabilities of absence in favorable conditions, for
example a value of 0.1 means a species will be absent in
10\% of favorble plots.

\sQuote{autlin} Autecological response model.  
Coenoflex combines the physiological response of a species
along each gradient into an autecological response using:

\tabular{ll}{
min(x,y,...) \tab  = minimum response of gradients x,y,... \cr
max(x,y,...) \tab  = maximum response of gradients x,y,... \cr
irm(x,y,...) \tab  = integrated rate methodology response 
of gradients x,y,... \cr
geo(x,y,...) \tab  = geometric mean response of gradients 
x,y,... \cr
ave(x,y,...) \tab  = arithmetic mean response of gradients 
of gradients x,y,... \cr
}

Functions can be nested as e.g.  min(1,irm(2,3),geo(4,5))
takes the integrated rate methodology of gradients 2 and 3, the geometric
mean of gradients 4 and 5, and then finally the minimum value of gradients 1 and
those combinations. 
}

\value{
a list object with class \sQuote{coenoflex} with components:
\item{params}{a copy of all the parameters used to specify the model}
\item{site}{the sample plot coordinates on each gradient}
\item{veg}{the vegetation data matrix}
\item{spcamp}{the parameters that describe the species amplitude response
functions}
\item{maxabu}{the list of maximum abundance for each species}
\item{pltprd}{the plot productivity for each sample plot} 
}

\references{
Austin, M.P. and T.M. Smith.  1989.  A new model for the continuum concept.
Vegetatio 83:35-47.

Collins, S.L., S. Glenn, and D.W. Roberts.  1993.  The hierarchical
continuum concept.  J. Veg. Sci. 4:149-156.

Roberts, D.W.  2008.  Statistical analysis of multidimensional
fuzzy set ordination.  Ecology 89:1246-1260.}

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
}
\keyword{data}

