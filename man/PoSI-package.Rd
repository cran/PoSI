\name{PoSI-package}
\alias{PoSI-package}
\docType{package}
\title{

Valid Post-Selection Inference for Linear LS Regression
%% What the package does (short line)
%% ~~ package title ~~
}
\description{

  In linear LS regression, calculate for a given regressor matrix the
  multiplier K of coefficient standard errors such that the confidence
  intervals [b - K*SE(b), b + K*SE(b)] have a guaranteed coverage
  probability for all coefficient estimates \code{b} in any submodels
  after performing arbitrary model selection.

%% More about what it does (maybe more than one line)
%% ~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab PoSI\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2020-10-24\cr
License: \tab GPL-3 \cr
}
%% ~~ An overview of how to use the package, including the most important ~~
%% ~~ functions ~~
}
\author{Andreas Buja and Kai Zhang
%% Who wrote it

Maintainers: Andreas Buja <buja.at.wharton@gmail.com>,
            Kai Zhang <zhangk@email.unc.edu> and
            Wan Zhang <wanz63@live.unc.edu>
%% Who to complain to <yourfault@somewhere.net>
%% ~~ The author and/or maintainer of the package ~~
}
\references{
%% ~~ Literature or other references for background information ~~
``Valid Post-Selection Inference,''
Berk, R., Brown, L., Buja, A., Zhang, K., Zhao, L.,
The Annals of Statistics, 41 (2), 802--837~(2013).
}
%% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%% ~~ the R documentation directory ~~
\keyword{ post-selection inference }
\keyword{ family-wise error }
\keyword{ model selection }
\keyword{ LS regression }
\seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
  \code{\link{lm}}, \code{link{model.matrix}}
}
\examples{
%% ~~ simple examples of the most important functions ~~
data(Boston, package="MASS")
summary(PoSI(Boston[,-14]))
}
