/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.importer.genotype;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.genotype.GeliFileWriter;
import edu.mit.broad.picard.genotype.GenotypeLikelihoods;
import edu.mit.broad.picard.genotype.GenotypeLikelihoodsCodec;
import edu.mit.broad.picard.genotype.GenotypeLikelihoods.GenotypeLikelihoodsComparator;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.SAMTextHeaderCodec;
import edu.mit.broad.sam.util.AsciiLineReader;
import edu.mit.broad.sam.util.SortingCollection;

/**
 * Converts a BED/BIM/FAM file trio to a number of GELI files (1 per individual).
 * BED files come in 2 formats, individual-major and snp-major. The former lists all SNPs for the 
 * first individual then all SNPs for the second individual, etc. The latter list all individuals 
 * for first SNP then all individuals for second SNP, etc. The order for snps is dictated by
 * the bim file and the order for individuals is dictated by the fam file.
 * <p>
 * See <a href="http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml">this page</a> for details
 * of the format.
 *
 * @author Doug Voet
 */
public class BedToGeli extends CommandLineProgram {
    static final float LIKELIHOOD = 500;
    private static final Log log = Log.getInstance(BedToGeli.class);
    
    @Usage(programVersion="1.0")
    public final String USAGE = "";
    
    @Option(doc="The bed file name.", mutex="BFILE")
    public File BED;

    @Option(doc="The bim file name.", mutex="BFILE")
    public File BIM;

    @Option(doc="The fam file name.", mutex="BFILE")
    public File FAM;

    @Option(doc="The root file name of the bed, bim & fam files.", mutex={"BED", "BIM", "FAM"})
    public String BFILE;
    
    @Option(doc="The directory to write the output GELI files", shortName="D")
    public File OUTPUT_DIR;
    
    @Option(doc="Set to 'true' if the family name should be included in the output file names, default false", 
            shortName="F",
            optional=true)
    public Boolean USE_FAMILY = Boolean.FALSE;
    
    @Option(doc="Name of file containing sequence dictionary to embed in new GELI files",
            shortName="DICT")
    public File SEQUENCE_DICTIONARY;
    
    private List<SNP> snpCache;
    private List<String> geliFileNames;
    private List<SAMSequenceRecord> sequenceDictionary;
    private Map<String, Byte> referenceIndexes;

    public static void main(String[] argv) {
        System.exit(new BedToGeli().instanceMain(argv));
    }

    @Override
    protected int doWork() {
        populateFileNames();
        IoUtil.assertFileIsReadable(this.BED);
        IoUtil.assertFileIsReadable(this.BIM);
        IoUtil.assertFileIsReadable(this.FAM);
        IoUtil.assertFileIsReadable(this.SEQUENCE_DICTIONARY);
        IoUtil.assertDirectoryIsWritable(this.OUTPUT_DIR);
        
        populateSequenceDictionary();
        
        BedFileReader bedReader = new BedFileReader(this.BED);
        if (bedReader.getMode() == BedFileReader.MODE_INDIVIDUAL_MAJOR) {
            log.debug("Detected BED file in individual-major mode");
            parseIndividualMajor(bedReader);
        } else {
            log.debug("Detected BED file in snp-major mode");
            parseSnpMajor(bedReader);
        }
        
        return 0;
    }

    /**
     * loads the SEQUENCE_DICTIONARY file
     */
    private void populateSequenceDictionary() {
        try {
            final SAMFileHeader header = new SAMTextHeaderCodec().decode(new AsciiLineReader(new FileInputStream(this.SEQUENCE_DICTIONARY)), null);
            this.sequenceDictionary = header.getSequences();
            
            this.referenceIndexes = new HashMap<String, Byte>();
            for (byte i = 0; i < sequenceDictionary.size(); i++) {
                this.referenceIndexes.put(sequenceDictionary.get(i).getSequenceName().intern(), i);
            }
        } catch (FileNotFoundException e) {
            throw new PicardException("Unexpected exception", e);
        }
    }

    private void parseIndividualMajor(BedFileReader bedReader) {
        cacheSnps();
        BasicTextFileParser famReader = new BasicTextFileParser(true, this.FAM);
        for (String[] famFields : famReader) {
            GeliFileWriter geliWriter = getGeliFileWriter(getGeliFileName(famFields[0], famFields[1]), false);
            for (SNP snp : this.snpCache) {
                GenotypeLikelihoods genotypeLikelihoods = constructGenotypeLikelihoods(
                        bedReader, snp);
                if (genotypeLikelihoods != null) {
                    geliWriter.addGenotypeLikelihoods(genotypeLikelihoods);
                }
            }
            bedReader.dropRemainingBlock();
            geliWriter.close();
        }
        famReader.close();
    }

    /**
     * @return null if for a no-call or the snp has no position on the genome
     */
    private char[] getNextGenotype(BedFileReader bedReader, SNP snp) {
        char[] genotype = new char[2];
        byte genotypeCode = bedReader.nextGenotype();
        if (snp == null) {
            // unplaced marker... we need to read the genotype off the reader so we don't lose
            // our place, but we cannot put the marker in the geli file.
            return null;
        }
        switch (genotypeCode) {
        case BedFileReader.GENOTYPE_AA:
            genotype[0] = (char) snp.getAllele1();
            genotype[1] = (char) snp.getAllele1();
            break;
        case BedFileReader.GENOTYPE_AB:
            genotype[0] = (char) snp.getAllele1();
            genotype[1] = (char) snp.getAllele2();
            break;
        case BedFileReader.GENOTYPE_BB:
            genotype[0] = (char) snp.getAllele2();
            genotype[1] = (char) snp.getAllele2();
            break;
        case BedFileReader.GENOTYPE_NO_CALL:
            // don't record a genotype likelihood for a no call
            return null;
        default:
            throw new PicardException("Unknown genotype code: " + Integer.toBinaryString(genotypeCode));
        }
        return genotype;
    }

    private void cacheSnps() {
        BasicTextFileParser bimReader = null;
        try {
            bimReader = new BasicTextFileParser(true, this.BIM);
            this.snpCache = new LinkedList<SNP>();
            for (String[] bimFields : bimReader) {
                SNP snp = constructSnp(bimFields);
                snpCache.add(snp);
            }
        } finally {
            try {
                bimReader.close();
            } catch (Exception e) {
            }
        }
    }

    private SNP constructSnp(String[] bimFields) {
        byte referenceIndex = getReferenceIndex(bimFields[0]);
        if (referenceIndex == -1) {
            return null;
        }
        SNP snp = new SNP(
                referenceIndex,
                Integer.parseInt(bimFields[3]),
                bimFields[4].toUpperCase().getBytes()[0],
                bimFields[5].toUpperCase().getBytes()[0]);
        return snp;
    }

    /**
     * determines the index in the sequence dictionary for the given chromosome
     */
    private byte getReferenceIndex(String chromosome) {
        final String referenceName;
        int chromosomeNumber;
        try {
            chromosomeNumber = Integer.parseInt(chromosome);
        } catch (NumberFormatException e) {
            chromosomeNumber = -1;
        }
        
        if (chromosomeNumber >= 1 && chromosomeNumber <= 22) {
            referenceName = ("chr" + chromosome).intern();
        } else if (chromosomeNumber == 26 || chromosome.equalsIgnoreCase("MT")) {
            referenceName = "chrM";
        } else if (chromosomeNumber == 23 || chromosomeNumber == 25 ||
                chromosome.equalsIgnoreCase("XY") || chromosome.equalsIgnoreCase("X")) {
            referenceName = "chrX";
        } else if (chromosomeNumber == 24 || chromosome.equalsIgnoreCase("Y")) {
            referenceName = "chrY";
        } else {
            // unplaced marker
            return -1;
        }

        Byte referenceIndex = this.referenceIndexes.get(referenceName);
        if (referenceIndex == null) {
            throw new PicardException("Reference sequence [" + referenceName + "] not found in sequence dictionary");
        }
        return referenceIndex;
    }

    private void cacheGELIFileNames() {
        BasicTextFileParser famReader = null;
        try {
            famReader = new BasicTextFileParser(true, this.FAM);
            this.geliFileNames = new LinkedList<String>();
            for (String[] fields : famReader) {
                this.geliFileNames.add(getGeliFileName(fields[0], fields[1]));
            }
        } finally {
            try {
                famReader.close();
            } catch (Exception e) {
            }
        }
    }
    
    private void parseSnpMajor(BedFileReader bedReader) {
        cacheGELIFileNames();        
        BasicTextFileParser bimReader = new BasicTextFileParser(true, this.BIM);
        Map<String, SortingCollection<GenotypeLikelihoods>> likelihoodsByFile =
            new HashMap<String, SortingCollection<GenotypeLikelihoods>>(
                    (int) Math.ceil(this.geliFileNames.size() * 1.34));
        
        int maxRecordsInRam = calculateMaxRecordsInRam();
        for (String geliFileName : this.geliFileNames) {
            likelihoodsByFile.put(geliFileName, SortingCollection.newInstance(
                    GenotypeLikelihoods.class, 
                    new GenotypeLikelihoodsCodec(), 
                    new GenotypeLikelihoodsComparator(), 
                    maxRecordsInRam));
        }
        
        for (String[] bimFields : bimReader) {
            for (String fileName : this.geliFileNames) {
                SNP snp = constructSnp(bimFields);
                GenotypeLikelihoods genotypeLikelihoods = constructGenotypeLikelihoods(
                        bedReader, snp);
                if (genotypeLikelihoods != null) {
                    likelihoodsByFile.get(fileName).add(genotypeLikelihoods);
                }
            }
            bedReader.dropRemainingBlock();
        }
        bimReader.close();
        
        writeGeliFiles(likelihoodsByFile);
    }

    /**
     * @return
     */
    private int calculateMaxRecordsInRam() {
        Runtime.getRuntime().gc();
        double memoryToUse = Runtime.getRuntime().maxMemory() * .8; // use up to 80%
        int objectCountLimit = (int) (memoryToUse / GenotypeLikelihoods.OBJECT_SIZE_BYTES);
        return objectCountLimit / this.geliFileNames.size();
    }

    /**
     * @param likelihoodsByFile
     */
    private void writeGeliFiles(
            Map<String, SortingCollection<GenotypeLikelihoods>> likelihoodsByFile) {

        for (Map.Entry<String, SortingCollection<GenotypeLikelihoods>> entry : likelihoodsByFile.entrySet()) {
            GeliFileWriter fileWriter = getGeliFileWriter(entry.getKey(), true);
            for (GenotypeLikelihoods likelihoods : entry.getValue()) {
                fileWriter.addGenotypeLikelihoods(likelihoods);
            }
            fileWriter.close();
        }
    }

    private GeliFileWriter getGeliFileWriter(
            String fileName, boolean presorted) {
        File geliFile = new File(this.OUTPUT_DIR, fileName);
        GeliFileWriter fileWriter = new GeliFileWriter(geliFile, presorted);
        SAMFileHeader header = new SAMFileHeader();
        header.setAttribute(SAMFileHeader.VERSION_TAG, "1.0");
        header.setSequences(this.sequenceDictionary);
        fileWriter.setHeader(header);
        return fileWriter;
    }

    /**
     * @param bedReader
     * @param snp
     * @return
     */
    private GenotypeLikelihoods constructGenotypeLikelihoods(
            BedFileReader bedReader, SNP snp) {
        char[] genotype = getNextGenotype(bedReader, snp);
        if (genotype == null) {
            // no call or unplaced marker
            return null;
        }
        
        GenotypeLikelihoods genotypeLikelihoods = new GenotypeLikelihoods();
        genotypeLikelihoods.setLikelihood(
                GenotypeLikelihoods.getLikelihoodIndex(genotype), 
                LIKELIHOOD);
        genotypeLikelihoods.setReferenceIndex(snp.getReferenceIndex());
        genotypeLikelihoods.setPosition(snp.getPosition());
        return genotypeLikelihoods;
    }

    /**
     * populates bed/bim/fam if bfile option is used
     */
    private void populateFileNames() {
        if (this.BFILE != null) {
            this.BED = new File(this.BFILE + ".bed");
            this.BIM = new File(this.BFILE + ".bim");
            this.FAM = new File(this.BFILE + ".fam");
        }
    }

    /**
     * @return the appropriate name taking into account this.USE_FAMILY
     */
    private String getGeliFileName(String family, String individual) {
        StringBuilder fileName = new StringBuilder(individual).append(".geli");
        if (this.USE_FAMILY) {
            fileName.insert(0, "_").insert(0, family);
        }
        return fileName.toString();
    }
}
