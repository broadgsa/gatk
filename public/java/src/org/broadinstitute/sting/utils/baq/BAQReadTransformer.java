package org.broadinstitute.sting.utils.baq;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Applies Heng's BAQ calculation to a stream of incoming reads
 */
public class BAQReadTransformer extends ReadTransformer {
    private BAQ baqHMM;
    private IndexedFastaSequenceFile refReader;
    private BAQ.CalculationMode cmode;
    private BAQ.QualityMode qmode;

    @Override
    public ApplicationTime initializeSub(final GenomeAnalysisEngine engine, final Walker walker) {
        final BAQMode mode = WalkerManager.getWalkerAnnotation(walker, BAQMode.class);
        this.refReader = engine.getReferenceDataSource().getReference();
        this.cmode = engine.getArguments().BAQMode;
        this.qmode = mode.QualityMode();
        baqHMM = new BAQ(engine.getArguments().BAQGOP);

        if ( qmode == BAQ.QualityMode.DONT_MODIFY )
            throw new ReviewedStingException("BUG: shouldn't create BAQ transformer with quality mode DONT_MODIFY");

        if ( mode.ApplicationTime() == ReadTransformer.ApplicationTime.FORBIDDEN && enabled() )
            throw new UserException.BadArgumentValue("baq", "Walker cannot accept BAQ'd base qualities, and yet BAQ mode " + cmode + " was requested.");

        return mode.ApplicationTime();
    }

    @Override
    public boolean enabled() {
        return cmode != BAQ.CalculationMode.OFF;
    }

    @Override
    public GATKSAMRecord apply(final GATKSAMRecord read) {
        baqHMM.baqRead(read, refReader, cmode, qmode);
        return read;
    }
}
