package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * samtools pileup [-f in.ref.fasta] [-t in.ref_list] [-l in.site_list] [-iscg] [-T theta] [-N nHap] [-r pairDiffRate] <in.alignment>
 *
 * Print the alignment in the pileup format. In the pileup format, each line represents a genomic position,
 * consisting of chromosome name, coordinate, reference base, read bases, read qualities and alignment mapping
 * qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all
 * encoded at the read base column. At this column, a dot stands for a match to the reference base on the forward strand,
 * a comma for a match on the reverse strand, ‘ACGTN’ for a mismatch on the forward strand and ‘acgtn’ for a mismatch on the
 * reverse strand.
 *
 * A pattern ‘\+[0-9]+[ACGTNacgtn]+’ indicates there is an insertion between this reference position and the next
 * reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence.
 * Similarly, a pattern ‘-[0-9]+[ACGTNacgtn]+’ represents a deletion from the reference.
 * Also at the read base column, a symbol ‘^’ marks the start of a read segment which is a contiguous subsequence on the read
 * separated by ‘N/S/H’ CIGAR operations. The ASCII of the character following ‘^’ minus 33 gives the mapping quality.
 * A symbol ‘$’ marks the end of a read segment.
 */
public class PileupWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    @Argument(fullName="alwaysShowSecondBase",doc="If true, prints dummy bases for the second bases in the BAM file where they are missing",required=false)
    public boolean alwaysShowSecondBase = false;

    @Argument(fullName="showSecondBaseQuals",doc="If true, prints out second base qualities in the pileup",required=false)
    public boolean showSecondBaseQuals = false;


    @Argument(fullName="extended",shortName="ext",doc="extended",required=false)
    public boolean EXTENDED = false;

    public boolean FLAG_UNCOVERED_BASES = true;     // todo: how do I make this a command line argument?
    
    public void initialize() {
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        
        if ( bases.equals("") && FLAG_UNCOVERED_BASES ) {
            bases = "***UNCOVERED_SITE***";
        }

        StringBuilder extras = new StringBuilder();

        String secondBasePileup = pileup.getSecondaryBasePileup();
        if ( secondBasePileup == null && alwaysShowSecondBase ) {
            secondBasePileup = Utils.dupString('N', bases.length());
        }
        if ( secondBasePileup != null ) extras.append(" ").append(secondBasePileup);

        if ( showSecondBaseQuals ) {
            String secondQualPileup = pileup.getSecondaryQualPileup();
            if ( secondQualPileup == null )
             secondQualPileup = Utils.dupString((char)(33), bases.length());
            extras.append(" ").append(secondQualPileup);
        }

        String rodString = "";
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if ( datum != null && ! (datum instanceof rodDbSNP)) {
                //System.out.printf("rod = %s%n", datum.toSimpleString());
                rodString += datum.toSimpleString();
                //System.out.printf("Rod string %s%n", rodString);
            }
        }
        
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp != null )
            rodString += dbsnp.toMediumString();

        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        //if ( context.getLocation().getStart() % 1 == 0 ) {
        out.printf("%s%s %s%n", pileup.getPileupString(), extras, rodString);
        //}

        if ( EXTENDED ) {
            String probDists = pileup.getProbDistPileup();
            out.println(probDists);
        }

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return treeReduce(sum,value);
    }
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }
}
