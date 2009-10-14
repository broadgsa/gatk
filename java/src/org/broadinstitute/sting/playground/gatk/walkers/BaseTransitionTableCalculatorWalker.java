package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.HashSet;
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;
import java.io.PrintStream;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;


class ReferenceContextWindow {

    protected int windowSize;
    protected int nPrevBases;
    protected LinkedList<AlignmentContext> prevAlignments;
    protected LinkedList<ReferenceContext> prevRefs;
    protected LinkedList<Boolean> usePrevious;
    protected boolean initialized;

    public ReferenceContextWindow( int nPrevBases ) {
        windowSize = 2*nPrevBases + 1;
        this.nPrevBases = nPrevBases;
        prevAlignments = new LinkedList<AlignmentContext>();
        prevRefs = new LinkedList<ReferenceContext>();
        usePrevious = new LinkedList<Boolean>();
        initialized = false;
    }

    public void update( ReferenceContext ref, AlignmentContext context, boolean useLocus ) {
        if ( ! initialized ) {
            prevAlignments.add(context);
            prevRefs.add(ref);
            usePrevious.add(useLocus);
            if ( prevAlignments.size() == windowSize ) {
                initialized = true;
            }
        } else {
            prevAlignments.removeFirst();
            prevRefs.removeFirst();
            usePrevious.removeFirst();
            prevAlignments.add(context);
            prevRefs.add(ref);
            usePrevious.add(useLocus);
        }
    }

    public String getReferenceString() {
        String ref = "";
        for ( ReferenceContext c : prevRefs ) {
            ref = ref + c.getBase();
        }

        return ref;
    }

    public String getForwardRefString() {
        String ref = "";
        for ( ReferenceContext c : prevRefs.subList(0,nPrevBases+1) ) {
            ref = ref + c.getBase();
        }

        return ref;
    }

    public String getReverseRefString() { // todo -- make sure we want to flip this done (yes we do)
        String ref = "";
        for ( int base = prevRefs.size()-1; base >= nPrevBases; base -- ) {
            ref = ref + prevRefs.get(base).getBase();
        }

        return ref;
    }

    public AlignmentContext getContext() {
        // because lists are 0-indexed, this returns the alignments
        // to the middle base in the window.
        return prevAlignments.get(nPrevBases);
    }

    public ReferenceContext getMiddleReferenceContext() {
        return prevRefs.get(nPrevBases);
    }

    public boolean isValidWindow() {
        boolean valid;
        if ( ! initialized ) {
            valid = false;
        } else {
            valid = true;
            // check if everything is confident ref
            for ( Boolean b : usePrevious ) {
                if ( !b ) {
                    valid = false;
                    break;
                }
            }
            // if still valid, check distances
            if ( valid ) {
                ListIterator<ReferenceContext> iter = prevRefs.listIterator();
                ReferenceContext prev = iter.next();
                while ( iter.hasNext() ) {
                    ReferenceContext cur = iter.next();

                    if ( cur.getLocus().distance(prev.getLocus()) > 1 ) {
                        valid = false;
                        break;
                    }
                
                    prev = cur;
                }
            }
        }

        return valid;
    }

    public int getWindowSize() {
        return windowSize;
    }

}