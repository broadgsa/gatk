package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.TribbleRMDTrackBuilder;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.Assert;
import org.junit.Test;
import org.junit.BeforeClass;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFWriterUnitTest
 *         <p/>
 *         This class tests out the ability of the VCF writer to correctly write VCF files
 */
public class VCFWriterUnitTest extends BaseTest {
    private Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
    private Set<String> additionalColumns = new HashSet<String>();
    private File fakeVCFFile = new File("FAKEVCFFILEFORTESTING.vcf");

    @BeforeClass
    public static void beforeTests() {
        try {
            IndexedFastaSequenceFile seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
            GenomeLocParser.setupRefContigOrdering(seq);
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
    }

    /** test, using the writer and reader, that we can output and input a VCF file without problems */
    @Test
    public void testBasicWriteAndRead() {
        VCFHeader header = createFakeHeader(metaData,additionalColumns);
        VCFWriter writer = new VCFWriter(fakeVCFFile);
        writer.writeHeader(header);
        writer.addRecord(createVCFRecord(header));
        writer.addRecord(createVCFRecord(header));
        writer.close();
        VCFReader reader = new VCFReader(fakeVCFFile);
        int counter = 0;
        // validate what we're reading in
        validateHeader(reader.getHeader());
        for (VCFRecord rec : reader) {
            counter++;
        }
        Assert.assertEquals(2,counter);
        reader.close();
        new File(fakeVCFFile + TribbleRMDTrackBuilder.linearIndexExtension).delete();
        fakeVCFFile.delete();
    }

    /**
     * create a fake header of known quantity
     * @param metaData           the header lines
     * @param additionalColumns  the additional column names
     * @return a fake VCF header
     */
    public static VCFHeader createFakeHeader(Set<VCFHeaderLine> metaData, Set<String> additionalColumns) {
        metaData.add(new VCFHeaderLine(VCFHeader.FILE_FORMAT_KEY, VCFHeader.VCF_VERSION));
        metaData.add(new VCFHeaderLine("two", "2"));
        additionalColumns.add("FORMAT");
        additionalColumns.add("extra1");
        additionalColumns.add("extra2");
        return new VCFHeader(metaData, additionalColumns);
    }

    /**
     * create a fake VCF record
     * @param header the VCF header
     * @return a VCFRecord
     */
    private VCFRecord createVCFRecord(VCFHeader header) {
        List<VCFGenotypeEncoding> altBases = new ArrayList<VCFGenotypeEncoding>();
        altBases.add(new VCFGenotypeEncoding("C"));
        altBases.add(new VCFGenotypeEncoding("D1"));
        Map<String,String> infoFields = new HashMap<String,String>();
        infoFields.put("DP","50");

        List<VCFGenotypeRecord> gt = new ArrayList<VCFGenotypeRecord>();
        for (String name : header.getGenotypeSamples()) {
            List<VCFGenotypeEncoding> myAlleles = new ArrayList<VCFGenotypeEncoding>();
            myAlleles.add(new VCFGenotypeEncoding("C"));
            myAlleles.add(new VCFGenotypeEncoding("D1"));
            VCFGenotypeRecord rec = new VCFGenotypeRecord(name, myAlleles, VCFGenotypeRecord.PHASE.PHASED);
            rec.setField("bb", "0");
            gt.add(rec);
        }
        return new VCFRecord("A","chr1",1,"RANDOM",altBases,0,".",infoFields, "GT:AA",gt);

    }


    /**
     * validate a VCF header
     * @param header the header to validate
     */
    public void validateHeader(VCFHeader header) {
        // check the fields
        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
            Assert.assertEquals(VCFHeader.HEADER_FIELDS.values()[index], field);
            index++;
        }
        Assert.assertEquals(metaData.size(), header.getMetaData().size());
        index = 0;
        for (String key : header.getGenotypeSamples()) {
            Assert.assertTrue(additionalColumns.contains(key));
            index++;
        }
        Assert.assertEquals(additionalColumns.size(), index+1 /* for the header field we don't see */);
    }
}
