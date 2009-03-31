package org.broadinstitute.sting.gatk.dataSources.datum;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
/**
 *
 * User: aaron
 * Date: Mar 30, 2009
 * Time: 2:53:37 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Mar 30, 2009
 * <p/>
 * Class ReadDatum
 * <p/>
 * The base read datum class.
 */
public class ReadDatum implements Datum {

    // our SAM record
    final private SAMRecord sam;

    // our locus context
    final private LocusContext locus;

    // the constructor, taking a sam read and a locus
    public ReadDatum(SAMRecord r, LocusContext locus) {
        this.sam = r;
        this.locus = locus;
    }

    // get the SAMRecord
    public SAMRecord getRead() {
        return this.sam;
    }

    // get the locus context
    public LocusContext getLocus() {
        return this.locus;
    }

    /**
     * gets the region that our read spans
     *
     * @return a genome loc that details the region that our read spans.
     */
    public GenomeLoc getSequenceLocation() {
        return Utils.genomicLocationOf(sam);
    }
}
