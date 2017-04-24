package org.labkey.discvrseq;

import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.commandline.Output;

import java.io.PrintStream;

/**
 * Created by bimber on 1/1/2015.
 */
public class Walker1 extends LocusWalker {
    @Output
    PrintStream out;

    @Override
    public Object map(RefMetaDataTracker refMetaDataTracker, ReferenceContext referenceContext, AlignmentContext alignmentContext) {
        return null;
    }

    @Override
    public Object reduceInit() {
        return null;
    }

    @Override
    public Object reduce(Object o, Object o2) {
        return null;
    }
}

