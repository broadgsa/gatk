package org.broadinstitute.sting.playground.gatk.walkers.graphalign;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.fasta.FastaReferenceWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.StringUtil;

@WalkerName("GraphReferenceBuilder")
@Requires(value={DataSource.REFERENCE})
public class GraphReferenceBuilder extends RefWalker<Integer, Integer> {
    @Argument(fullName="graphFile", shortName="GF", doc="", required=true)
    String graphFile = null;

    @Argument(fullName="DEBUG", shortName="DB", doc="", required=false)
    boolean DEBUG = false;

    @Argument(fullName="VALIDATE", shortName="VD", doc="", required=false)
    boolean VALIDATE_GRAPH = false;

    @Argument(fullName="printFrequency", shortName="F", doc="", required=false)
    int printFrequency = 10000;

    ObjectOutputStream graphSerialStream = null;

    ReferenceGraph graphRef = null;
    ReferenceSequenceFile flatReferenceFile = null;
    
    public void initialize() {
        super.initialize();

        graphRef = new ReferenceGraph(DEBUG);

        try {
            graphSerialStream = new ObjectOutputStream( new FileOutputStream( graphFile ) );
        } catch ( FileNotFoundException e ) {
            throw new StingException("Couldn't open file " + graphFile, e);
        } catch ( IOException e ) {
            throw new StingException("Couldn't write to file " + graphFile, e);
        }

        flatReferenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.getToolkit().getArguments().referenceFile);

        ReferenceSequence refSeq = flatReferenceFile.nextSequence();
        do {
            //logger.info("Read " + refSeq);
            graphRef.bindRefenceSequence(refSeq);
            logger.info(String.format("contig %s has length %d", refSeq.getName(), refSeq.length()));
            refSeq = flatReferenceFile.nextSequence();
        } while ( refSeq != null );

        System.out.println(graphRef.toBriefString());
    }

    int counter = printFrequency;
    public Integer map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
//        if ( context.getLocation().getStart() == 59384 ) {
//            try {
//                Thread.currentThread().sleep(5000);
//            } catch (InterruptedException e) {
//                ;
//            }
//        }

        boolean alreadyAddedAtThisLoc = false;
        for ( ReferenceOrderedDatum rod : rodData.getAllRods() ) {
            if ( rod instanceof Variation && ! alreadyAddedAtThisLoc ) {
                // if we have multiple variants at a locus, just take the first damn one we see for now
                Variation variant = (Variation) rod;
                // todo -- getAlternativeBases should be getAlleles()
                GenomeLoc loc = variant.getLocation();
                String[] allAllelesList = null; // variant.getAlternateBases().split(""); // todo fixme
                if ( allAllelesList.length >= 3 ) { // bad dbSNP format :-(
                    List<String> alleles = Arrays.asList(allAllelesList).subList(1,3);
                    //logger.info(String.format("Adding %s %s", loc, alleles));
                    graphRef.addVariation(variant, loc, alleles);
                    //logger.info(String.format("  Added %s %s", loc, alleles));
                    alreadyAddedAtThisLoc = true;
                    if ( counter-- == 0 ) {
                        logger.info(String.format("Added %s %s %s", loc, alleles, graphRef.toBriefString()));
                        counter = printFrequency;
                        if ( VALIDATE_GRAPH )
                            graphRef.validateGraph();
                    }
                }
            }
        }

        return null;
    }

    // todo -- graph should be the reduce result
    public Integer reduceInit() {
        return null;
    }

	public Integer reduce(Integer value, Integer sum) {
		return sum;
	}

    public void onTraversalDone(Integer sum) {
        super.onTraversalDone(sum);
        try {
            graphSerialStream.writeObject(graphRef);
            graphSerialStream.close();
        } catch ( IOException e ) {
            throw new StingException("Couldn't write to file " + graphFile, e);
        }
    }
}

