/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.graphalign;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.List;

/**
 * A completely experimental walker that constructs a graphical reference that incorporates variation from provided
 * RODs [Not for public use and will change drastically in the future].
 */
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
        for ( VariantContext vc : rodData.getAllVariantContexts(ref)) {
            if ( ! alreadyAddedAtThisLoc ) {
                // if we have multiple variants at a locus, just take the first damn one we see for now
                // todo -- getAlternativeBases should be getAlleles()
                GenomeLoc loc = vc.getLocation();
                String[] allAllelesList = null; // variant.getAlternateBases().split(""); // todo fixme
                if ( allAllelesList.length >= 3 ) { // bad dbSNP format :-(
                    List<String> alleles = Arrays.asList(allAllelesList).subList(1,3);
                    //logger.info(String.format("Adding %s %s", loc, alleles));
                    graphRef.addVariation(vc, loc, alleles);
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

