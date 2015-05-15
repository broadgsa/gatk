\name{gsa.read.gatkreport}
\alias{gsa.read.gatkreport}
\title{
Function to read in a GATKReport
}
\description{
This function reads in data from a GATKReport. A GATKReport is a document containing multiple tables produced by the GATK. Each table is loaded as a separate data.frame object in a list.
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
The GATKReport format replaces the multi-file output format used previously by many GATK tools and provides a single, consolidated file format.  This format accommodates multiple tables and is still R-loadable through this function. 
}
\value{
Returns a LIST object, where each key is the TableName and the value is the data.frame object with the contents of the table.  If multiple tables with the same name exist, each one after the first will be given names of TableName.v1, TableName.v2, ..., TableName.vN.
}
\references{
http://www.broadinstitute.org/gatk/guide/article?id=1244
}
\author{
Kiran Garimella
}
\note{
This function accepts different versions of the GATKReport format by making internal calls to gsa.read.gatkreportv0() or gsa.read.gatkreportv1() as appropriate.
}
\examples{
test_file = system.file("extdata", "test_gatkreport.table", package = "gsalib");
report = gsa.read.gatkreport(test_file);
}
\keyword{ manip }
