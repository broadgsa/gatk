\name{gsa.read.vcf}
\alias{gsa.read.vcf}
\title{
gsa.read.vcf
}
\description{
Reads a VCF file into a table.  Optionally expands genotype columns into separate columns containing the genotype, separate from the other fields specified in the FORMAT field.
}
\usage{
gsa.read.vcf(vcffile, skip=0, nrows=-1, expandGenotypeFields = FALSE)
}
\arguments{
  \item{vcffile}{
The path to the vcf file.
}
  \item{skip}{
The number of lines of the data file to skip before beginning to read data.
}
  \item{nrows}{
The maximum number of rows to read in.  Negative and other invalid values are ignored.
}
  \item{expandGenotypeFields}{
If TRUE, adds an additional column per sample containing just the genotype.
}
}
\details{
The VCF format is the standard variant call file format used in the GATK.  This function reads that data in as a table for easy analysis.
}
\value{
Returns a data.frame object, where each column corresponds to the columns in the VCF file.
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
vcf = gsa.read.vcf("/path/to/my/output.vcf");
}
\keyword{ ~kwd1 }
