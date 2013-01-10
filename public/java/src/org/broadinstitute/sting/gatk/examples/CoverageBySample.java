/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.examples;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Computes the coverage per sample for every position (use with -L argument!).
 */
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
public class CoverageBySample extends LocusWalker<Integer, Integer> {
    @Output
    protected PrintStream out;    

    private HashSet<String> sampleNames = new HashSet<String>();

    public boolean requiresReads() { return true; }

    public void initialize() {

        List<SAMReadGroupRecord> read_groups = this.getToolkit().getSAMFileHeader().getReadGroups();

        for ( SAMReadGroupRecord record : read_groups ) {
            String sample = record.getSample();
            if ( sample != null )
                sampleNames.add(sample);
        } 
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        HashMap<String, Integer> depthBySample = new HashMap<String, Integer>();
        for ( String sample : sampleNames )
            depthBySample.put(sample, 0);

        ReadBackedPileup pileup = context.getPileup();
        for ( PileupElement p : pileup ) {

            SAMReadGroupRecord readGroup = p.getRead().getReadGroup();
            if ( readGroup == null )
                continue;

            String sample = readGroup.getSample();
            if ( sample != null ) {
                int oldDepth = depthBySample.get(sample);
                depthBySample.put(sample, oldDepth + 1);
            }
        }

        for ( Map.Entry<String, Integer> sample : depthBySample.entrySet() ) {
            out.printf("  %s %8d%n", sample.getKey(), sample.getValue());
        }

        return 1;
    }


    public void onTraversalDone(Integer result) {
        out.println("Processed " + result + " loci.");
    }

    public Integer reduceInit() {
		return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }
}
