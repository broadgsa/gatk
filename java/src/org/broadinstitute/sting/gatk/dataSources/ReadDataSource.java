package org.broadinstitute.sting.gatk.dataSources;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.util.Iterator;
import java.io.File;

import org.broadinstitute.sting.gatk.iterators.VerifyingSamIterator;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 10:35:40 AM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class ReadDataSource {

    /**
     * our SAM data files
     */
    // our SAM reader
    private SAMFileReader samReader = null;
    // iterator over the sam records in the readsFile
    private Iterator<SAMRecord> samReadIter = null;

    // The verifying iterator, it does checking
    VerifyingSamIterator verifyingSamReadIter = null;


    /**
     * our reference data source
     */
    // The reference data -- filename, refSeqFile, and iterator
    private File refFileName = null;                        // the name of the reference file
    //private ReferenceSequenceFile refFile = null;
    private FastaSequenceFile2 refFile = null;              // todo: merge FastaSequenceFile2 into picard!
    private ReferenceIterator refIter = null;
}
