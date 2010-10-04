\name{gsalib-package}
\alias{gsalib-package}
\alias{gsalib}
\docType{package}
\title{
GATK utility analysis functions
}
\description{
Utility functions for analyzing GATK-processed NGS data
}
\details{
This package contains functions for working with GATK-processed NGS data.  These functions include a command-line parser that also allows a script to be used in interactive mode (good for developing scripts that will eventually be automated), a proportional Venn diagram generator, convenience methods for parsing VariantEval output, and more.
}
\author{
Genome Sequencing and Analysis Group

Medical and Population Genetics Program

Maintainer: Kiran Garimella
}
\references{
GSA wiki page: http://www.broadinstitute.org/gsa/wiki

GATK help forum: http://www.getsatisfaction.com/gsa
}
\examples{
## get script arguments in interactive and non-interactive mode
cmdargs = getargs( list(
    requiredArg1 = list(
        value = NA,
        doc   = "Documentation for requiredArg1"
    ),

    optionalArg1 = list(
        value = 3e9,
        doc   = "Documentation for optionalArg1"
    )
) );

## plot a proportional Venn diagram
plot.venn(500, 250, 0, 100);

## read a GATKReport file
report = gatk.report("/path/to/my/output.gatkreport");

## emit a message
gsa.message("This is a message");

## emit a warning message
gsa.message("This is a warning message");

## emit an error message
gsa.message("This is an error message");
}
\keyword{ package }
