package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 30, 2009
 */
public interface Covariate {
	public Comparable<?> getValue(SAMRecord read, int offset, char[] refBases); // used to pick out the value from the read and etc
	public Comparable<?> getValue(String str); // used to get value from input file
}
