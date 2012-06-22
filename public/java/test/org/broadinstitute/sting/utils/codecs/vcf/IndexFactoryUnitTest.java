package org.broadinstitute.sting.utils.codecs.vcf;

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriterFactory;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.EnumSet;

/**
 * tests out the various functions in the index factory class
 */
public class IndexFactoryUnitTest extends BaseTest {

    File inputFile = new File(privateTestDir + "HiSeq.10000.vcf");
    File outputFile = new File(privateTestDir + "onTheFlyOutputTest.vcf");
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
        Index indexFromInputFile = IndexFactory.createDynamicIndex(inputFile, new VCFCodec());
        if ( outputFileIndex.exists() ) {
            System.err.println("Deleting " + outputFileIndex);
            outputFileIndex.delete();
        }

        for ( int maxRecords : Arrays.asList(0, 1, 10, 100, 1000, -1)) {
            AbstractFeatureReader<VariantContext> source = AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), new VCFCodec(), indexFromInputFile);

            int counter = 0;
            final EnumSet<Options> options = EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
            VariantContextWriter writer = VariantContextWriterFactory.create(outputFile, dict, options);
            writer.writeHeader((VCFHeader)source.getHeader());
            CloseableTribbleIterator<VariantContext> it = source.iterator();
            while (it.hasNext() && (counter++ < maxRecords || maxRecords == -1) ) {
                VariantContext vc = it.next();
                writer.add(vc);
            }
            writer.close();

            // test that the input index is the same as the one created from the identical input file
            // test that the dynamic index is the same as the output index, which is equal to the input index
            //WalkerTest.assertOnDiskIndexEqualToNewlyCreatedIndex(outputFileIndex, "unittest", outputFile);
        }
    }
}
