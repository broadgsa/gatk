package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.picard.reference.ReferenceSequence;

import java.util.List;
import java.io.File;
import java.io.IOException;

/**
 * samtools pileup [-f in.ref.fasta] [-t in.ref_list] [-l in.site_list] [-iscg] [-T theta] [-N nHap] [-r pairDiffRate] <in.alignment>
 *
 * Print the alignment in the pileup format. In the pileup format, each line represents a genomic position,
 * consisting of chromosome name, coordinate, reference base, read bases, read qualities and alignment mapping
 * qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all
 * encoded at the read base column. At this column, a dot stands for a match to the reference base on the forward strand,
 * a comma for a match on the reverse strand, 'ACGTN' for a mismatch on the forward strand and 'acgtn' for a mismatch on the
 * reverse strand.
 *
 * A pattern '\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next
 * reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence.
 * Similarly, a pattern '-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference.
 * Also at the read base column, a symbol '^' marks the start of a read segment which is a contiguous subsequence on the read
 * separated by 'N/S/H' CIGAR operations. The ASCII of the character following '^' minus 33 gives the mapping quality.
 * A symbol '$' marks the end of a read segment.
 */
public class PileupWithContextWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    @Argument(fullName="alwaysShowSecondBase",doc="If true, prints dummy bases for the second bases in the BAM file where they are missing",required=false)
    public boolean alwaysShowSecondBase = false;

    //@Argument(fullName="showSecondBaseQuals",doc="If true, prints out second base qualities in the pileup",required=false)
    //public boolean showSecondBaseQuals = false;

    @Argument(fullName="qualsAsInts",doc="If true, prints out qualities in the pileup as comma-separated integers",required=false)
    public boolean qualsAsInts = false;

    @Argument(fullName="extended",shortName="ext",doc="extended",required=false)
    public boolean EXTENDED = false;

    @Argument(fullName="flagUncoveredBases",shortName="fub",doc="Flag bases with zero coverage",required=false)
    public boolean FLAG_UNCOVERED_BASES = true;

    @Argument(fullName="contextBases",shortName="cb",doc="How much context around the locus should we show?",required=false)
    public int contextBases = 0;

    public IndexedFastaSequenceFile refSeq;
    public ReferenceSequence contigSeq = null;
    public String contig = null;

    public void initialize() {
        File refFile = this.getToolkit().getArguments().referenceFile;

        try {
            refSeq = new IndexedFastaSequenceFile(refFile);
        } catch (IOException e) {
            refSeq = null;
        }
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

        /*
        if ( showSecondBaseQuals ) {
            String secondQualPileup = pileup.getSecondaryQualPileup();
            if ( secondQualPileup == null )
             secondQualPileup = Utils.dupString((char)(33), bases.length());
            extras.append(" ").append(secondQualPileup);
        }
        */
        if (contig == null || !context.getContig().equals(contig)) {
            contig = context.getContig();
            contigSeq = refSeq.getSequence(contig);
        }

        if (contextBases > 0 && refSeq != null) {
            long startPos = context.getPosition() - contextBases <= 0 ? 1 : context.getPosition() - contextBases;
            long stopPos = context.getPosition() + contextBases > contigSeq.length() ? contigSeq.length() : context.getPosition() + contextBases;

            ReferenceSequence prevRefSequence = refSeq.getSubsequenceAt(context.getContig(), startPos, context.getPosition());
            ReferenceSequence nextRefSequence = refSeq.getSubsequenceAt(context.getContig(), context.getPosition(), stopPos);

            extras.append(" prev=").append(new String(prevRefSequence.getBases()));
            extras.append(" next=").append(new String(nextRefSequence.getBases()));
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
        //System.out.printf("quals as ints %b%n", qualsAsInts);
        out.printf("%s%s %s%n", pileup.getPileupString(qualsAsInts), extras, rodString);
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