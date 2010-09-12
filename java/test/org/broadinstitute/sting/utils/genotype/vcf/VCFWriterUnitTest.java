package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.TribbleRMDTrackBuilder;
import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.Assert;
import org.junit.Test;
import org.junit.BeforeClass;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

import net.sf.picard.reference.IndexedFastaSequenceFile;


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
        IndexedFastaSequenceFile seq = new IndexedFastaSequenceFile(new File(hg18Reference));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    /** test, using the writer and reader, that we can output and input a VCF file without problems */
    @Test
    public void testBasicWriteAndRead() {
        VCFHeader header = createFakeHeader(metaData,additionalColumns);
        VCFWriter writer = new StandardVCFWriter(fakeVCFFile);
        writer.writeHeader(header);
        writer.add(createVC(header),"A".getBytes()[0]);
        writer.add(createVC(header),"A".getBytes()[0]);
        writer.close();
        VCFCodec reader = new VCFCodec();
        AsciiLineReader lineReader;
        VCFHeader headerFromFile = null;
        try {
             lineReader = new AsciiLineReader(new FileInputStream(fakeVCFFile));
             headerFromFile = (VCFHeader)reader.readHeader(lineReader);
        }
        catch (FileNotFoundException e ) {
            throw new GATKException(e.getMessage());
        }

        int counter = 0;

        // validate what we're reading in
        validateHeader(headerFromFile);
        
        try {
            while(true) {
                String line = lineReader.readLine();
                if (line == null)
                    break;

                VariantContext vc = (VariantContext)reader.decode(line);
                counter++;
            }
            Assert.assertEquals(2,counter);
            new File(fakeVCFFile + TribbleRMDTrackBuilder.indexExtension).delete();
            fakeVCFFile.delete();
        }
        catch (IOException e ) {
            throw new GATKException(e.getMessage());
        }

    }

    /**
     * create a fake header of known quantity
     * @param metaData           the header lines
     * @param additionalColumns  the additional column names
     * @return a fake VCF header
     */
    public static VCFHeader createFakeHeader(Set<VCFHeaderLine> metaData, Set<String> additionalColumns) {
        metaData.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_0.getFormatString(), VCFHeaderVersion.VCF4_0.getVersionString()));
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
    private VariantContext createVC(VCFHeader header) {

        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1",1);
        List<Allele> alleles = new ArrayList<Allele>();
        Set<String> filters = null;
        Map<String, String> attributes = new HashMap<String,String>();
        Map<String, Genotype> genotypes = new HashMap<String,Genotype>();

        alleles.add(Allele.create("-",true));
        alleles.add(Allele.create("CC",false));

        attributes.put("DP","50");
        for (String name : header.getGenotypeSamples()) {
            Map<String, String> gtattributes = new HashMap<String,String>();
            gtattributes.put("BB","1");
            Genotype gt = new Genotype(name,alleles.subList(1,2),0,null,gtattributes,true);

            genotypes.put(name,gt);
            
        }
        return new VariantContext("RANDOM",loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes, 0, filters, attributes);


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
