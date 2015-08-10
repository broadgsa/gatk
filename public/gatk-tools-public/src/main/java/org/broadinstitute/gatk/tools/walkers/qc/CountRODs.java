/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.qc;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.engine.walkers.NanoSchedulable;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.ExpandingArrayList;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;
import java.util.*;

/**
 * Count the number of ROD objects encountered
 *
 * <p>CountRods is a RODWalker, and so traverses the data by ROD (reference ordered data). For example if the ROD passed to it is a VCF file,
 * it will count the variants in the file.</p>
 *
 * <p>Note that this tool is different from CountRodsByRef which is a RefWalker, and so traverses the data by
 * position along the reference. CountRodsByRef can count ROD elements (such as, but not limited to, variants) found
 * at each position or within specific intervals if you use the -L argument (see CommandLineGATK).</p>
 *
 * <p>Both these tools are different from CountVariants in that they are more generic (they can also count RODs that
 * are not variants) and CountVariants is more detailed, in that it computes additional statistics (type of variants
 * being indels vs. SNPs etc). </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more ROD files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Number of RODs seen.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CountRODs \
 *   -R reference.fasta \
 *   -o output.txt \
 *   --rod input.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class CountRODs extends RodWalker<CountRODs.Datum, Pair<ExpandingArrayList<Long>, Long>> implements TreeReducible<Pair<ExpandingArrayList<Long>, Long>>, NanoSchedulable {
    @Output
    public PrintStream out;

    /**
     * One or more input rod files
     */
    @Input(fullName="rod", shortName = "rod", doc="Input VCF file(s)", required=true)
    public List<RodBinding<Feature>> rods = Collections.emptyList();

    @Argument(fullName = "verbose", shortName = "v", doc="If true, this tool will print out detailed information about the rods it finds and locations", required = false)
    public boolean verbose = false;

    @Argument(fullName = "showSkipped", shortName = "s", doc="If true, this tool will print out the skipped locations", required = false)
    public boolean showSkipped = false;

    @Override
    public Pair<ExpandingArrayList<Long>, Long> treeReduce(Pair<ExpandingArrayList<Long>, Long> lhs, Pair<ExpandingArrayList<Long>, Long> rhs) {
        ExpandingArrayList<Long> nt = new ExpandingArrayList<Long>();
        nt.addAll(lhs.first);
        int index = 0;
        for (Long l : rhs.first) {
            if (nt.get(index) == null)
                nt.add(l);
            else
                nt.set(index,nt.get(index) + l);
            index++;
        }
        return new Pair<ExpandingArrayList<Long>, Long>(nt, lhs.second + rhs.second);
    }

    public class Datum {
        public long nRodsAtThisLocation = 0;
        public long nSkippedBases =0;
        public long nTotalBases = 0;

        public Datum( long nRodsAtThisLocation, long nSkippedBases, long nTotalBases ) {
            this.nRodsAtThisLocation = nRodsAtThisLocation;
            this.nSkippedBases = nSkippedBases;
            this.nTotalBases = nTotalBases;
        }

        public String toString() {
            return String.format("<%d %d %d>", nRodsAtThisLocation, nSkippedBases, nTotalBases);
        }
    }

    public Datum map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc cur = context.getLocation();

        if ( verbose && showSkipped ) {
            for(long i = context.getSkippedBases(); i >= 0; i--) {
                SAMSequenceDictionary dictionary = getToolkit().getReferenceDataSource().getReference().getSequenceDictionary();
                SAMSequenceRecord contig = dictionary.getSequence(cur.getContig());
                if(cur.getStop() < contig.getSequenceLength())
                    cur = getToolkit().getGenomeLocParser().incPos(cur,1);
                else
                    cur = getToolkit().getGenomeLocParser().createGenomeLoc(dictionary.getSequence(contig.getSequenceIndex()+1).getSequenceName(),1,1);
                out.printf("%s: skipped%n", cur);

            }
        }

        long nRodsHere = 0;
        long nTotalBases = 0;

        if ( ref == null ) {
            // we're getting the last skipped update
            if ( verbose )
                out.printf("Last position was %s: skipping %d bases%n",
                        context.getLocation(), context.getSkippedBases() );
            nRodsHere = -1; // don't update this
            nTotalBases = context.getSkippedBases();
        } else {
            Collection<RODRecordList> rods = new LinkedList<RODRecordList>();
            for ( RODRecordList rod : tracker.getBoundRodTracks() ) {
                //System.out.printf("Considering rod %s%n", rod);
                if ( rod.getLocation().getStart() == context.getLocation().getStart() && ! rod.getName().equals("interval") ) {
                    // only consider the first element
                    //System.out.printf("adding it%n");
                    rods.add(rod);
                }
            }

            nRodsHere = rods.size();

            if ( nRodsHere > 0 ) {
                if ( verbose ) {
                    List<String> names = new ArrayList<String>();
                    for ( RODRecordList rod : rods ) {
                        names.add(rod.getName());
                    }

                    //System.out.printf("context is %s", context.getSkippedBases());
                    out.printf("At %s: found %d rod(s) [%s] after skipping %d bases%n",
                            context.getLocation(), nRodsHere, Utils.join(",", names), context.getSkippedBases() );
                }
            }

            nTotalBases = context.getSkippedBases() + 1;
        }

        return new Datum(nRodsHere, context.getSkippedBases(), nTotalBases);
    }

    public Pair<ExpandingArrayList<Long>, Long> reduceInit() {
        return new Pair<ExpandingArrayList<Long>, Long>(new ExpandingArrayList<Long>(), 0l);
    }

    private void updateCounts(ExpandingArrayList<Long> counts, long nRods, long nObs) {
        if ( nRods >= 0 ) {
            long prev = counts.get((int)nRods) == null ? 0l : counts.get((int)nRods);
            counts.set((int)nRods, nObs + prev);
        }
    }

    public Pair<ExpandingArrayList<Long>, Long> reduce(Datum point, Pair<ExpandingArrayList<Long>, Long> sum) {
        ExpandingArrayList<Long> counts = sum.getFirst();
        updateCounts(counts, point.nRodsAtThisLocation, 1);
        updateCounts(counts, 0, point.nSkippedBases);

        Pair<ExpandingArrayList<Long>, Long> r = new Pair<ExpandingArrayList<Long>, Long>(counts, point.nTotalBases + sum.getSecond());

        //System.out.printf("Reduce: %s %s => %s%n", point, sum, r);
        return r;
    }
}