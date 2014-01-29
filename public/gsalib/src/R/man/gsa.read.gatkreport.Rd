\name{gsa.read.gatkreport}
\alias{gsa.read.gatkreport}
\title{
gsa.read.gatkreport
}
\description{
Reads a GATKReport file - a multi-table document - and loads each table as a separate data.frame object in a list.
}
\usage{
gsa.read.gatkreport(filename)
}
\arguments{
  \item{filename}{
The path to the GATKReport file.
}
}
\details{
The GATKReport format replaces the multi-file output format used by many GATK tools and provides a single, consolidated file format.  This format accomodates multiple tables and is still R-loadable - through this function.

The file format looks like this:
\preformatted{##:GATKReport.v0.1 TableName : The description of the table
col1   col2                      col3
0      0.007451835696110506      25.474613284804366
1      0.002362777171937477      29.844949954504095
2      9.087604507451836E-4      32.87590975254731
3      5.452562704471102E-4      34.498999090081895
4      9.087604507451836E-4      35.14831665150137
}

}
\value{
Returns a list object, where each key is the TableName and the value is the data.frame object with the contents of the table.  If multiple tables with the same name exist, each one after the first will be given names of "TableName.v1", "TableName.v2", ..., "TableName.vN".
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kiran Garimella
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
report = gsa.read.gatkreport("/path/to/my/output.gatkreport");
}
\keyword{ ~kwd1 }
