package org.broadinstitute.sting.utils.codecs.vcf;

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.*;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * tests out the various functions in the index factory class
 */
public class IndexFactoryUnitTest extends BaseTest {

    File inputFile = new File("public/testdata/HiSeq.10000.vcf");
    File outputFile = new File("public/testdata/onTheFlyOutputTest.vcf");
    File outputFileIndex = Tribble.indexFile(outputFile);

    private SAMSequenceDictionary dict;

    @BeforeTest
    public void setup() {
        try {
            dict = new CachingIndexedFastaSequenceFile(new File(b37KGReference)).getSequenceDictionary();
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(b37KGReference,ex);
        }
    }

    //
    // test out scoring the indexes
    //
    @Test
    public void testOnTheFlyIndexing1() throws IOException {
        Index indexFromInputFile = IndexFactory.createIndex(inputFile, new VCFCodec());
        if ( outputFileIndex.exists() ) {
            System.err.println("Deleting " + outputFileIndex);
            outputFileIndex.delete();
        }

        for ( int maxRecords : Arrays.asList(0, 1, 10, 100, 1000, -1)) {
            BasicFeatureSource<VariantContext> source = new BasicFeatureSource<VariantContext>(inputFile.getAbsolutePath(), indexFromInputFile, new VCFCodec());

            int counter = 0;
            VCFWriter writer = new StandardVCFWriter(outputFile, dict);
            writer.writeHeader((VCFHeader)source.getHeader());
            CloseableTribbleIterator<VariantContext> it = source.iterator();
            while (it.hasNext() && (counter++ < maxRecords || maxRecords == -1) ) {
                VariantContext vc = it.next();
                writer.add(vc);
            }
            writer.close();

            // test that the input index is the same as the one created from the identical input file
            // test that the dynamic index is the same as the output index, which is equal to the input index
            WalkerTest.assertOnDiskIndexEqualToNewlyCreatedIndex(outputFileIndex, "unittest", outputFile);
        }
    }
}
