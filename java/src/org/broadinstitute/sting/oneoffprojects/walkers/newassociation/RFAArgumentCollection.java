package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;

import java.io.File;
import java.util.List;

/**
 * Argument collection for the read feature associator and related walkers
 */
public class RFAArgumentCollection {

    @Argument(doc="ReadFeatures you want to test. None specified = all will be tested.",required=false,shortName="f",fullName="Feature")
    public List<String> inputFeatures = null;

    @Argument(doc="Size of window on which to perform tests",required=false,fullName="windowSize")
    public int windowSize = 50;

    @Argument(doc="Size of the jump between tested windows",required=false,fullName="windowJump")
    public int windowJump = 10;

    @Argument(doc="File containing a list of case samples",required=false,shortName="case",fullName="case")
    public File caseFile = null;

    @Argument(doc="File containing a list of control samples",required=false,shortName="control",fullName="control")
    public File controlFile = null;

    @Argument(doc="Fixed significance level, as a Z-score",required=false,shortName="z",fullName="fixedZ")
    public double fixedZ = 6.0;

    @Argument(doc="Insert size below which to flag a read as aberrant",required=false,shortName="LIS",fullName="LowInsertSize")
    public int lowInsertSize = 50;

    @Argument(doc="Insert size above which to flag a read as aberrant",required=false,shortName="HIS",fullName="HighInsertSize")
    public int highInsertSize = 450;

    @Argument(doc="Significance level for determining whether a sample contains a significant proportion of affected reads",required=false,shortName="sz",fullName="perSampleZ")
    public double sampleZThresh = 3.0;

    @Argument(doc="Lower bound for significant proportion of affected reads",shortName="se",fullName="sampleEpsilon",required=false)
    public double EPSILON = 0.05;

    @Argument(doc="Number of clipped bases for a read to be considered \"split\"",required=false,shortName="cb",fullName="clippedBases")
    public short clippedBases = 4;

    /* todo -- these for a possible "split read" binarification of clipped bases -- if there's a fast way to traverse clipped bases
    @Argument(doc="Minimum base quality for a clipped base to be indicative of a split read",required=false,shortName="CLBQ",fullName="CLBQ")
    public short clbq = 10;

    @Argument(doc="Minimum base quality sum for clipped bases (as a string) to be indicative of a split read",required=false,shortName="CLBS",fullName="CLBS")
    public short clbs=50;
    */
}
