\name{gsa.read.squidmetrics}
\alias{gsa.read.squidmetrics}
\title{
gsa.read.squidmetrics
}
\description{
Reads metrics for a specified SQUID project into a dataframe.
}
\usage{
gsa.read.squidmetrics("C315")
}
\arguments{
  \item{project}{
The project for which metrics should be obtained.
}
  \item{bylane}{
If TRUE, obtains per-lane metrics rather than the default per-sample metrics.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a data frame with samples (or lanes) as the row and the metric as the column.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kiran Garimella
}
\note{
This method will only work within the Broad Institute internal network.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Obtain metrics for project C315.
d = gsa.read.squidmetrics("C315");
}
\keyword{ ~kwd1 }
