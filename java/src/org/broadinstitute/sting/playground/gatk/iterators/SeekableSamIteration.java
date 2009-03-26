package org.broadinstitute.sting.playground.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordComparator;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 26, 2009
 * Time: 9:27:59 AM
 * To change this template use File | Settings | File Templates.
 */
public interface SeekableSamIteration {
    public boolean supportsSeeking();
    public void queryOverlapping( final String contig, final int start, final int stop );
    public void query(final String contig, final int start, final int stop, final boolean contained);
    public void queryContained(final String contig, final int start, final int stop);
}
