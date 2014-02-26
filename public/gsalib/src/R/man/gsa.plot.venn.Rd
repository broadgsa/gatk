\name{gsa.plot.venn}
\alias{gsa.plot.venn}
\title{
Plot a proportional venn diagram
}
\description{
Plot a proportional venn diagram (two or three-way venns allowed)
}
\usage{
gsa.plot.venn(a, b, c = 0, a_and_b, a_and_c = 0, b_and_c = 0, col = c("#FF6342", "#63C6DE", "#ADDE63"), pos = c(0.2, 0.2, 0.8, 0.82), debug = 0)
}
\arguments{
  \item{a}{
size of 'a' circle
}
  \item{b}{
size of 'b' circle
}
  \item{c}{
size of 'c' circle
}
  \item{a_and_b}{
size of a and b overlap
}
  \item{a_and_c}{
size of a and c overlap
}
  \item{b_and_c}{
size of b and c overlap
}
  \item{col}{
vector of colors for each venn piece
}
  \item{pos}{
vector of positional elements
}
  \item{debug}{
if 1, set debug mode and print useful information
}
}
\details{
Plots a two-way or three-way proportional Venn diagram.  Internally, this method uses the Google Chart API to generate the diagram, then renders it into the plot window where it can be annotated in interesting ways.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
}
\author{
Kiran Garimella
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Plot a two-way Venn diagram
gsa.plot.venn(1000, 750, 0, 400);

## Plot a three-way Venn diagram
gsa.plot.venn(1000, 750, 900, 400, 650, 500);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
