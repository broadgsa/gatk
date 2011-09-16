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

package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;
import java.util.*;

/**
 * Prints out counts of the number of reference ordered data objects encountered.
 *
 *
 * <h2>Input</h2>
 * <p>
 * One or more rod files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Number of rods seen.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountRODs \
 *   -o output.txt \
 *   --rod input.vcf
 * </pre>
 *
 */
public class CountRODsWalker extends RodWalker<CountRODsWalker.Datum, Pair<ExpandingArrayList<Long>, Long>> implements TreeReducible<Pair<ExpandingArrayList<Long>, Long>> {
    @Output
    public PrintStream out;

    /**
     * One or more input rod files
     */
    @Input(fullName="rod", shortName = "rod", doc="Input VCF file(s)", required=false)
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