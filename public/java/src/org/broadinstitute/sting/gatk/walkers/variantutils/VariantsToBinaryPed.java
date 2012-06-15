package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.R.RScriptExecutorException;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.*;
import java.util.*;

/**
 * Yet another VCF to Ped converter. The world actually does need one that will
 * work efficiently on large VCFs (or at least give a progress bar). This
 * produces a binary ped file in individual major mode.
 */
public class VariantsToBinaryPed extends RodWalker<Integer,Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    @Input(shortName="m",fullName = "metaData",required=true,doc="Sample metadata file. You may specify a .fam file " +
            "(in which case it will be copied to the file you provide as fam output).")
    File metaDataFile;

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

    private ValidateVariants vv = new ValidateVariants();

    private static double APPROX_CM_PER_BP = 1000000.0/750000.0;

    private static final byte HOM_REF = 0x0;
    private static final byte HOM_VAR = 0x3;
    private static final byte HET =     0x2;
    private static final byte NO_CALL = 0x1;

    private static final int BUFFER_SIZE = 1000; //4k genotypes per sample = Nmb for N*1000 samples

    // note that HET and NO_CALL are flipped from the documentation: that's because
    // plink actually reads these in backwards; and we want to use a shift operator
    // to put these in the appropriate location

    private Map<String,OutputStream> printMap = new HashMap<String,OutputStream>();
    private Map<String,File> tempFiles = new HashMap<String,File>();
    private Map<String,byte[]> genotypeBuffer = new HashMap<String,byte[]>();
    private int genotypeCount = 0;
    private int byteCount = 0;
    private List<String> famOrder = new ArrayList<String>();

    public void initialize() {
        vv.variantCollection = variantCollection;
        vv.dbsnp = dbsnp;
        vv.DO_NOT_VALIDATE_FILTERED = true;
        vv.type = ValidateVariants.ValidationType.REF;
        // create temporary output streams and buffers

        // write magic bits into the ped file
        try {
            outBed.write(new byte[] { (byte) 0x6c, (byte) 0x1b, 0x0});
            // ultimately, the bed will be in individual-major mode
        } catch (IOException e) {
            throw new ReviewedStingException("error writing to output file.");
        }
        // write to the fam file, the first six columns of the standard ped file
        // first, load data from the input meta data file
        Map<String,Map<String,String>> metaValues = new HashMap<String,Map<String,String>>();
        Set<String> samplesToUse = new HashSet<String>();
        logger.debug("Reading in metadata...");
        try {
            if ( metaDataFile.getAbsolutePath().endsWith(".fam") ) {
                for ( String line : new XReadLines(metaDataFile) ) {
                    String[] famSplit = line.split("\\t");
                    String sid = famSplit[1];
                    outFam.printf("%s%n",line);
                }
            } else {
                for ( String line : new XReadLines(metaDataFile) ) {
                    logger.debug(line);
                    String[] split = line.split("\\t");
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
        // family ID, individual ID, Paternal ID, Maternal ID, Sex, Phenotype
        int dummyID = 0; // increments for dummy parental and family IDs used
        // want to be especially careful to maintain order here
        Map<String,VCFHeader> headers = VCFUtils.getVCFHeadersFromRods(getToolkit());
        for ( Map.Entry<String,VCFHeader> header : headers.entrySet() ) {
            if ( ! header.getKey().equals(variantCollection.variants.getName()) && ! metaDataFile.getAbsolutePath().endsWith(".fam") ) {
                continue;
            }
            for ( String sample : header.getValue().getGenotypeSamples() ) {
                Map<String,String> mVals = metaValues.get(sample);
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
                try {
                    File temp = File.createTempFile(sample, ".tmp");
                    printMap.put(sample,new PrintStream(temp));
                    tempFiles.put(sample,temp);
                } catch (IOException e) {
                    throw new ReviewedStingException("Error creating temporary file",e);
                }
                genotypeBuffer.put(sample,new byte[BUFFER_SIZE]);
                famOrder.add(sample);
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ! tracker.hasValues(variantCollection.variants) ||
                tracker.getFirstValue(variantCollection.variants).isFiltered() ||
                ! tracker.getFirstValue(variantCollection.variants).isSNP() ||
                ! tracker.getFirstValue(variantCollection.variants).isBiallelic()) {
            return 0;
        }
        try {
            vv.map(tracker,ref,context);
        } catch (UserException e) {
            throw new UserException("Input VCF file is invalid; we cannot guarantee the resulting ped file. "+
            "Please run ValidateVariants for more detailed information.");
        }

        VariantContext vc = tracker.getFirstValue(variantCollection.variants);
        String refOut;
        String altOut;
        boolean altMajor;
        if ( majorAlleleFirst ) {
            // want to use the major allele as ref
            HashMap<String,Object> ats = new HashMap<String,Object>(vc.getAttributes());
            if ( ! vc.hasAttribute("AF") ) {
                VariantContextUtils.calculateChromosomeCounts(vc,ats,true);
            }
            if ( getAF(ats.get("AF")) > 0.5 ) {
                refOut = vc.getAlternateAllele(0).getBaseString();
                altOut = vc.getReference().getBaseString();
                altMajor = true;
            } else {
                refOut = vc.getReference().getBaseString();
                altOut = vc.getAlternateAllele(0).getBaseString();
                altMajor = false;
            }
        } else {
            refOut = vc.getReference().getBaseString();
            altOut = vc.getAlternateAllele(0).getBaseString();
            altMajor = false;
        }
        // write an entry into the map file
        outBim.printf("%s\t%s\t%.2f\t%d\t%s\t%s%n",vc.getChr(),getID(vc),APPROX_CM_PER_BP*vc.getStart(),vc.getStart(),
                refOut,altOut);
        // store genotypes per sample into the buffer
        for ( Genotype g : vc.getGenotypes() ) {
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
            }
            genotypeCount = 0;
        }

        return 1;
    }

    public Integer reduce(Integer m, Integer r) {
        return r + m;
    }

    public Integer reduceInit() {
        return 0;
    }

    public void onTraversalDone(Integer numSites) {
        logger.info(String.format("%d sites processed!",numSites));
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
                }
                if ( ttr > 0 ) {
                    byte[] readGenotypes = new byte[ttr];
                    inStream.read(readGenotypes);
                    outBed.write(readGenotypes);
                }
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
        if ( g.hasGQ() && g.getGQ() < minGenotypeQuality ) {
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
        if ( g.hasGQ() && g.getGQ() < minGenotypeQuality ) {
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
}
