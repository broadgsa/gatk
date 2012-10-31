package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.alignment.bwa.java.AlignmentMatchSequence;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.*;
import java.util.*;

/**
 * Converts a VCF file to a binary plink Ped file (.bed/.bim/.fam)
 */
@DocumentedGATKFeature( groupName = "Variant Evaluation and Manipulation Tools", extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=100))
public class VariantsToBinaryPed extends RodWalker<Integer,Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * The metaData file can take two formats, the first of which is the first 6 lines of the standard ped file. This
     * is what Plink describes as a fam file. An example fam file is (note that there is no header):
     * <p><p>
     * CEUTrio NA12878 NA12891 NA12892 2 -9</p><p>
     * CEUTrio NA12891 UNKN1 UNKN2 2 -9</p><p>
     * CEUTrio NA12892 UNKN3 UNKN4 1 -9</p><p>
     * </p>
     * where the entries are (FamilyID IndividualID DadID MomID Phenotype Sex)
     * <p>
     * An alternate format is a two-column key-value file
     * </p><p><p>
     * NA12878        fid=CEUTrio;dad=NA12891;mom=NA12892;sex=2;phenotype=-9</p><p>
     * NA12891        fid=CEUTrio;sex=2;phenotype=-9</p><p>
     * NA12892        fid=CEUTrio;sex=1;phenotype=-9</p><p>
     * </p><p>
     * wherein unknown parents needn't be specified. The columns are the individual ID, and a list of key-value pairs.
     * </p><p>
     * Regardless of which file is specified, the walker will output a .fam file alongside the bed file. If the
     * command line has "-md [name].fam", the fam file will simply be copied. However, if a metadata file of the
     * alternate format is passed by "-md [name].txt", the walker will construct a formatted .fam file from the data.
     * </p>
     */
    @Input(shortName="m",fullName = "metaData",required=true,doc="Sample metadata file. You may specify a .fam file " +
            "(in which case it will be copied to the file you provide as fam output).")
    File metaDataFile;

    @Input(shortName="mode",fullName="outputMode",required=false,doc="The output file mode (SNP major or individual major)")
    OutputMode mode = OutputMode.INDIVIDUAL_MAJOR;

    @Output(shortName="bed",fullName = "bed",required=true,doc="output ped file")
    PrintStream outBed;

    @Output(shortName="bim",fullName="bim",required=true,doc="output map file")
    PrintStream outBim;

    @Output(shortName="fam",fullName="fam",required=true,doc="output fam file")
    PrintStream outFam;

    @Argument(shortName="mgq",fullName="minGenotypeQuality",required=true,doc="If genotype quality is lower than this value, output NO_CALL")
    int minGenotypeQuality = 0;

    @Argument(fullName="majorAlleleFirst",required=false,doc="Sets the major allele to be 'reference' for the bim file, rather than the ref allele")
    boolean majorAlleleFirst = false;

    enum OutputMode { INDIVIDUAL_MAJOR,SNP_MAJOR }

    private static double APPROX_CM_PER_BP = 1000000.0/750000.0;

    private static final byte HOM_REF = 0x0;
    private static final byte HOM_VAR = 0x3;
    private static final byte HET =     0x2;
    private static final byte NO_CALL = 0x1;

    private static final int BUFFER_SIZE = 1000; //4k genotypes per sample = Nmb for N*1000 samples

    private static final String PLINK_DELETION_MARKER = "-";

    // note that HET and NO_CALL are flipped from the documentation: that's because
    // plink actually reads these in backwards; and we want to use a shift operator
    // to put these in the appropriate location

    private Map<String,OutputStream> printMap = new HashMap<String,OutputStream>();
    private Map<String,File> tempFiles = new HashMap<String,File>();
    private Map<String,byte[]> genotypeBuffer = new HashMap<String,byte[]>();
    private int genotypeCount = 0;
    private int byteCount = 0;
    private List<String> famOrder = new ArrayList<String>();
    private long totalByteCount = 0l;
    private long totalGenotypeCount = 0l;

    public void initialize() {
        writeBedHeader();
        Map<String,Map<String,String>> sampleMetaValues = parseMetaData();
        // create temporary output streams and buffers

        // family ID, individual ID, Paternal ID, Maternal ID, Sex, Phenotype
        int dummyID = 0; // increments for dummy parental and family IDs used
        // want to be especially careful to maintain order here
        Map<String,VCFHeader> headers = VCFUtils.getVCFHeadersFromRods(getToolkit());
        for ( Map.Entry<String,VCFHeader> header : headers.entrySet() ) {
            if ( ! header.getKey().equals(variantCollection.variants.getName()) && ! metaDataFile.getAbsolutePath().endsWith(".fam") ) {
                continue;
            }
            for ( String sample : header.getValue().getGenotypeSamples() ) {
                if ( ! metaDataFile.getAbsolutePath().endsWith(".fam") ) {
                    Map<String,String> mVals = sampleMetaValues.get(sample);
                    if ( mVals == null ) {
                        throw new UserException("No metadata provided for sample "+sample);
                    }
                    if ( ! mVals.containsKey("phenotype") ) {
                        throw new UserException("No phenotype data provided for sample "+sample);
                    }
                    String fid = mVals.containsKey("fid") ? mVals.get("fid") : String.format("dummy_%d",++dummyID);
                    String pid = mVals.containsKey("dad") ? mVals.get("dad") : String.format("dummy_%d",++dummyID);
                    String mid = mVals.containsKey("mom") ? mVals.get("mom") : String.format("dummy_%d",++dummyID);
                    String sex = mVals.containsKey("sex") ? mVals.get("sex") : "3";
                    String pheno = mVals.get("phenotype");
                    outFam.printf("%s\t%s\t%s\t%s\t%s\t%s%n",fid,sample,pid,mid,sex,pheno);
                } else {
                    // even if a fam file is input, we can't diverge the bed file from the fam file, which
                    // could lead to a malformed plink trio. Fail fast if there's any extra sample in the VCF.
                    if ( ! sampleMetaValues.containsKey(sample) ) {
                        throw new UserException("No metadata provided for sample "+sample);
                    }
                }
                if ( mode == OutputMode.INDIVIDUAL_MAJOR ) {
                    // only need to instantiate the files and buffers if in individual major.
                    // Cut down on memory.
                    try {
                        File temp = File.createTempFile("VariantsToBPed_"+sample, ".tmp");
                        printMap.put(sample,new PrintStream(temp));
                        tempFiles.put(sample,temp);
                    } catch (IOException e) {
                        throw new ReviewedStingException("Error creating temporary file",e);
                    }
                    genotypeBuffer.put(sample,new byte[BUFFER_SIZE]);
                }
                famOrder.add(sample);
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return 0;
        }

        VariantContext vc = tracker.getFirstValue(variantCollection.variants,context.getLocation());
        if ( vc == null || vc.isFiltered() || ! vc.isBiallelic() ) {
            return 0;
        }
        try {
            validateVariantSite(vc,ref,context);
        } catch (TribbleException e) {
            throw new UserException("Input VCF file is invalid; we cannot guarantee the resulting ped file. "+
            "Please run ValidateVariants for more detailed information. This error is: "+e.getMessage());
        }

        String refOut;
        String altOut;
        String vcRef = getReferenceAllele(vc);
        String vcAlt = getAlternateAllele(vc);
        boolean altMajor;
        if ( majorAlleleFirst ) {
            // want to use the major allele as ref
            HashMap<String,Object> ats = new HashMap<String,Object>(vc.getAttributes());
            if ( ! vc.hasAttribute("AF") ) {
                VariantContextUtils.calculateChromosomeCounts(vc,ats,true);
            }
            if ( getAF(ats.get("AF")) > 0.5 ) {
                refOut = vcAlt;
                altOut = vcRef;
                altMajor = true;
            } else {
                refOut = vcRef;
                altOut = vcAlt;
                altMajor = false;
            }
        } else {
            refOut = vcRef;
            altOut = vcAlt;
            altMajor = false;
        }
        // write an entry into the map file
        outBim.printf("%s\t%s\t%.2f\t%d\t%s\t%s%n",vc.getChr(),getID(vc),APPROX_CM_PER_BP*vc.getStart(),vc.getStart(),
                refOut,altOut);
        if ( mode == OutputMode.INDIVIDUAL_MAJOR ) {
            writeIndividualMajor(vc,altMajor);
        } else {
            writeSNPMajor(vc,altMajor);
        }


        return 1;
    }

    public void writeIndividualMajor(VariantContext vc, boolean altMajor) {
        // store genotypes per sample into the buffer
        for ( Genotype g : vc.getGenotypes() ) {
            ++totalGenotypeCount;
            String sample = g.getSampleName();
            byte[] samBuf = genotypeBuffer.get(sample);
            byte enc = getEncoding(g,genotypeCount,altMajor);
            samBuf[byteCount] |= enc;
        }

        genotypeCount++;
        if ( genotypeCount % 4 == 0 ) {
            byteCount++;
            if ( byteCount >= BUFFER_SIZE ) {
                // dump the buffer to the print streams
                for ( String sample : printMap.keySet() ) {
                    OutputStream samOut = printMap.get(sample);
                    // print the buffer for this sample
                    try {
                        samOut.write(genotypeBuffer.get(sample));
                    } catch ( IOException e ) {
                        throw new ReviewedStingException("Error writing to temporary bed file.",e);
                    }
                    // reset the buffer for this sample
                    genotypeBuffer.put(sample,new byte[BUFFER_SIZE]);
                }
                byteCount = 0;
            }
            genotypeCount = 0;
        }
    }

    public void writeSNPMajor(VariantContext vc, boolean altMajor) {
        // for each sample, write the genotype into the bed file, in the
        // order of the fam file
        genotypeCount = 0;
        byteCount = 0;
        byte[] bytes = new byte[(3+famOrder.size())/4]; // this exploits java integer fractions, which round down by default (1-4) -> 1, (5-8) -> 2
        for ( Genotype g : vc.getGenotypesOrderedBy(famOrder) ) {
            byte enc = getEncoding(g,genotypeCount,altMajor);
            bytes[byteCount] |= enc;
            genotypeCount++;
            if ( genotypeCount % 4 == 0 ) {
                byteCount++;
                genotypeCount = 0;
            }
        }
        totalGenotypeCount += famOrder.size();
        totalByteCount += bytes.length;
        try {
            outBed.write(bytes);
        } catch (IOException e) {
            throw new ReviewedStingException("Error writing to output bed file",e);
        }
    }

    public Integer reduce(Integer m, Integer r) {
        return r + m;
    }

    public Integer reduceInit() {
        return 0;
    }

    public void onTraversalDone(Integer numSites) {
        logger.info(String.format("%d sites processed for a total of %d genotypes encoded in %d bytes",numSites,totalGenotypeCount,totalByteCount));

        if ( mode == OutputMode.INDIVIDUAL_MAJOR ) {
            mergeGenotypeTempFiles(numSites);
        }

    }

    private void mergeGenotypeTempFiles(int numSites) {
        // push out the remaining genotypes and close stream
        for ( String sample : printMap.keySet() ) {
            try {
                int lim = byteCount + (genotypeCount > 0 ? 1 : 0);
                printMap.get(sample).write(genotypeBuffer.get(sample),0,lim);
            } catch (IOException e) {
                throw new ReviewedStingException("Error closing temporary file.",e);
            }

            try {
               printMap.get(sample).close();
            } catch (IOException e) {
                throw new ReviewedStingException("Error closing temporary file.",e);
            }
        }
        for ( String sample : famOrder ) {
            logger.info("Merging genotypes for "+sample);
            FileInputStream inStream;
            try {
                inStream = new FileInputStream(tempFiles.get(sample));
            } catch (IOException e) {
                throw new ReviewedStingException("Error opening temp file for input.",e);
            }


            try {
                int ttr = numSites/4 + (genotypeCount > 0 ? 1 : 0);
                for ( ; ttr > BUFFER_SIZE ; ttr -= BUFFER_SIZE ) {
                    byte[] readGenotypes = new byte[BUFFER_SIZE];
                    inStream.read(readGenotypes);
                    outBed.write(readGenotypes);
                    totalByteCount += BUFFER_SIZE;
                }
                if ( ttr > 0 ) {
                    byte[] readGenotypes = new byte[ttr];
                    inStream.read(readGenotypes);
                    outBed.write(readGenotypes);
                    totalByteCount += ttr;
                }
                inStream.close();
            } catch (IOException e) {
                throw new ReviewedStingException("Error reading form temp file for input.",e);
            }
        }
    }

    private byte getEncoding(Genotype g, int offset, boolean altMajor) {
        if ( ! altMajor ) {
            return getStandardEncoding(g,offset);
        }

        return getFlippedEncoding(g,offset);
    }

    private byte getStandardEncoding(Genotype g, int offset) {
        byte b;
        if ( ! checkGQIsGood(g) ) {
            b = NO_CALL;
        } else if ( g.isHomRef() ) {
            b = HOM_REF;
        } else if ( g.isHomVar() ) {
            b = HOM_VAR;
        } else if ( g.isHet() ) {
            b = HET;
        } else {
            b = NO_CALL;
        }

        return (byte) (b << (2*offset));
    }

    private byte getFlippedEncoding(Genotype g, int offset) {
        byte b;
        if ( ! checkGQIsGood(g) ) {
            b = NO_CALL;
        } else if ( g.isHomRef() ) {
            b = HOM_VAR;
        } else if ( g.isHomVar() ) {
            b = HOM_REF;
        } else if ( g.isHet() ) {
            b = HET;
        } else {
            b = NO_CALL;
        }

        return (byte) (b << (2*offset));
    }

    private boolean checkGQIsGood(Genotype genotype) {
        if ( genotype.hasGQ() ) {
            return genotype.getGQ() >= minGenotypeQuality;
        } else if ( genotype.hasLikelihoods() ) {
            double log10gq = GenotypeLikelihoods.getGQLog10FromLikelihoods(genotype.getType().ordinal()-1,genotype.getLikelihoods().getAsVector());
            return MathUtils.log10ProbabilityToPhredScale(log10gq) >= minGenotypeQuality;
        }

        return minGenotypeQuality <= 0;
    }

    private static String getID(VariantContext v) {
        if ( v.hasID() ) {
            return v.getID();
        } else {
            return String.format("Var-%s-%d",v.getChr(),v.getStart());
        }
    }

    private double getAF(Object o) {
        if ( (o instanceof String) ) {
            return Double.parseDouble((String) o);
        } else if ( (o instanceof Double) ) {
            return (Double) o;
        } else {
            throw new UserException("Allele frequency appears to be neither String nor Double. Please check the header of your VCF.");
        }
    }

    private void writeBedHeader() {
        // write magic bits into the ped file
        try {
            outBed.write(new byte[] { (byte) 0x6c, (byte) 0x1b, (byte) (mode == OutputMode.INDIVIDUAL_MAJOR ? 0x0 : 0x1)});
            // ultimately, the bed will be in individual-major mode
        } catch (IOException e) {
            throw new ReviewedStingException("error writing to output file.");
        }
    }

    private Map<String,Map<String,String>> parseMetaData() {
        // write to the fam file, the first six columns of the standard ped file
        // first, load data from the input meta data file
        Map<String,Map<String,String>> metaValues = new HashMap<String,Map<String,String>>();
        logger.debug("Reading in metadata...");
        try {
            if ( metaDataFile.getAbsolutePath().endsWith(".fam") ) {
                for ( String line : new XReadLines(metaDataFile) ) {
                    String[] famSplit = line.split("\\s+");
                    if ( famSplit.length != 6 ) {
                        throw new UserException("Line of the fam file is malformatted. Expected 6 entries. Line is "+line);
                    }
                    String sid = famSplit[1];
                    String fid = famSplit[0];
                    String mom = famSplit[2];
                    String dad = famSplit[3];
                    String sex = famSplit[4];
                    String pheno = famSplit[5];
                    HashMap<String,String> values = new HashMap<String, String>();
                    values.put("mom",mom);
                    values.put("dad",dad);
                    values.put("fid",fid);
                    values.put("sex",sex);
                    values.put("phenotype",pheno);
                    metaValues.put(sid,values);
                    outFam.printf("%s%n",line);
                }
            } else {
                for ( String line : new XReadLines(metaDataFile) ) {
                    logger.debug(line);
                    String[] split = line.split("\\s+");
                    String sampleID = split[0];
                    String keyVals = split[1];
                    HashMap<String,String> values = new HashMap<String, String>();
                    for ( String kvp : keyVals.split(";") ) {
                        String[] kvp_split = kvp.split("=");
                        values.put(kvp_split[0],kvp_split[1]);
                    }
                    metaValues.put(sampleID,values);
                }
            }
        } catch (FileNotFoundException e) {
            throw new UserException("Meta data file not found: "+metaDataFile.getAbsolutePath(),e);
        }

        return metaValues;
    }

    private void validateVariantSite(VariantContext vc, ReferenceContext ref, AlignmentContext context) {
        final Allele reportedRefAllele = vc.getReference();
        final int refLength = reportedRefAllele.length();
        if ( refLength > 100 ) {
            logger.info(String.format("Reference allele is too long (%d) at position %s:%d; skipping that record.", refLength, vc.getChr(), vc.getStart()));
            return;
        }

        final byte[] observedRefBases = new byte[refLength];
        System.arraycopy(ref.getBases(), 0, observedRefBases, 0, refLength);
        final Allele observedRefAllele = Allele.create(observedRefBases);
        vc.validateReferenceBases(reportedRefAllele, observedRefAllele);
        vc.validateAlternateAlleles();
    }

    private String getReferenceAllele(VariantContext vc) {
        if ( vc.isSimpleInsertion() ) {
            // bi-allelic, so we just have "-" for ped output
            return PLINK_DELETION_MARKER;
        }
        if ( vc.isSymbolic() ) {
            // either symbolic or really long alleles. Plink alleles are allowed to be 1 or 2. Reference will just be 1.
            return "1";
        }
        if ( vc.isSimpleDeletion() ) {
            // bi-allelic. Want to take the standard representation and strip off the leading base.
            return vc.getReference().getBaseString().substring(1);
        }
        // snp or mnp
        return vc.getReference().getBaseString();
    }

    private String getAlternateAllele(VariantContext vc ) {
        if ( vc.isSimpleInsertion() ) {
            // bi-allelic. Want to take the standard representation and strip off the leading base.
            return vc.getAlternateAllele(0).getBaseString().substring(1);
        }
        if ( vc.isSymbolic() ) {
            // either symbolic or really long alleles. Plink alleles are allowed to be 1 or 2. Alt will just be 2.
            return "2";
        }
        if ( vc.isSimpleDeletion() ) {
            // bi-allelic, so we just have "-" for ped output
            return PLINK_DELETION_MARKER;
        }
        // snp or mnp
        return vc.getAlternateAllele(0).getBaseString();
    }
}
