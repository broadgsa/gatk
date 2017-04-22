\name{gsa.read.gatkreportv0}
\alias{gsa.read.gatkreportv0}
\title{
Function to read in an old-style GATKReport
}
\description{
This function reads in data from a version 0.x GATKReport. It should not be called directly; instead, use gsa.read.gatkreport()
}
\usage{
gsa.read.gatkreportv0(lines)
}
\arguments{
  \item{lines}{
The lines read in from the input file.
}
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
\keyword{ manip }
