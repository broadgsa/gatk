/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: April 7, 2011
 */
public class ReduceReadsWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Output
    protected StingSAMFileWriter out;

    @Argument(fullName = "SNPContextSize", shortName = "SCS", doc = "", required = true)
    protected int SNPContextSize;

    @Argument(fullName = "IndelContextSize", shortName = "ICS", doc = "", required = true)
    protected int IndelContextSize;

    protected ReducingSAMFileWriter reducingOut;
    protected int totalReads = 0;

    @Override
    public void initialize() {
        reducingOut = new ReducingSAMFileWriter(out, SNPContextSize, IndelContextSize);
    }

    @Override
    public SAMRecord map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        for ( GATKFeature feature : metaDataTracker.getAllCoveringRods() ) {
            if ( feature.getUnderlyingObject() instanceof VariantContext ) {
                VariantContext vc = (VariantContext)feature.getUnderlyingObject();
                reducingOut.addVariant(vc);
            }
        }

        totalReads++;
        return read; // all the work is done in the reduce step for this walker
    }


    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    @Override
    public SAMFileWriter reduceInit() {
        return reducingOut;
    }

    /**
     * given a read and a output location, reduce by emitting the read
     * @param read the read itself
     * @param output the output source
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public SAMFileWriter reduce( SAMRecord read, SAMFileWriter output ) {
        output.addAlignment(read);
        return output;
    }

    public void onTraversalDone( SAMFileWriter reduceResult ) {
        logger.info("Compressed reads: " + reducingOut.getNCompressedReads());
        logger.info("Total reads     : " + totalReads);
        // todo -- fixme
        //reducingOut.close();
    }
    
}
