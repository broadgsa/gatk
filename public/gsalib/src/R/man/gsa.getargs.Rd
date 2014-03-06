\name{gsa.getargs}
\alias{gsa.getargs}
\title{
Get script arguments
}
\description{
Get script arguments given a list object specifying arguments and documentation.  Can be used in command-line or interactive mode.  This is helpful when developing scripts in interactive mode that will eventually become command-line programs.  If no arguments are specified or help is requested in command-line mode, the script will print out a usage statement with available arguments and exit.
}
\usage{
gsa.getargs(argspec, doc = NA)
}
\arguments{
  \item{argspec}{
A list object.  Each key is an argument name.  The value is another list object with a 'value' and 'doc' keys.  For example:
\preformatted{argspec = list(
    arg1 = list(value=10, doc="Info for optional arg1"),
    arg2 = list(value=NA, doc="Info for required arg2")
);
}

If the value provided is NA, the argument is considered required and must be specified when the script is invoked.  For command-line mode, this means the argument must be specified on the command-line.  In interactive mode, there are two ways of specifying these arguments.  First, if a properly formatted list argument called 'cmdargs' is present in the current environment (i.e. the object returned by gsa.getargs() from a previous invocation), the value is taken from this object.  Otherwise, the argument is prompted for.
}

  \item{doc}{
An optional string succinctly documenting the purpose of the script.
}
}
\details{
Interactive scripts typically make use of hardcoded filepaths and parameter settings.  This makes testing easy, but generalization to non-interactive mode more difficult.  This utility provides a mechanism for writing scripts that work properly in both interactive and command-line modes.

To use this method, specify a list with key-value pairs representing the arguments as specified above.  In command-line mode, if no arguments are specified or the user specifies '-h' or '-help' anywhere on the command string, a help message indicating available arguments, their default values, and some documentation about the argument are provided.
}
\value{
Returns a list with keys matching the argspec and values representing the specified arguments.

\item{arg1 }{Value for argument 1}
\item{arg2 }{Value for argument 2}
...etc.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kiran Garimella
}
\examples{
argspec = list(
    file    = list(value="/my/test.vcf", doc="VCF file"),
    verbose = list(value=0,              doc="If 1, set verbose mode"),
    test2   = list(value=2.3e9,          doc="Another argument that does stuff")
);
        
cmdargs = gsa.getargs(argspec, doc="My test program");

print(cmdargs$file);  # will print '[1] "/my/test.vcf"'
}
\keyword{ ~kwd1 }
