package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Yet another VCF to Ped converter. The world actually does need one that will
 * work efficiently on large VCFs (or at least give a progress bar). This
 * produces a binary ped file in SNP-major mode.
 */
public class VariantsToPed extends RodWalker<Integer,Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    @Input(shortName="m",fullName = "metaData",required=true,doc="Sample metadata file. You may specify a .fam file (in which case it will be copied to the file you provide as fam output)")
    File metaDataFile;

    @Output(shortName="bed",fullName = "bed",required=true,doc="output ped file")
    PrintStream outBed;

    @Output(shortName="bim",fullName="bim",required=true,doc="output map file")
    PrintStream outBim;

    @Output(shortName="fam",fullName="fam",required=true,doc="output fam file")
    PrintStream outFam;

    private ValidateVariants vv = new ValidateVariants();

    private static double APPROX_CM_PER_BP = 1000000.0/750000.0;

    private static final byte HOM_REF = 0x0;
    private static final byte HOM_VAR = 0x3;
    private static final byte HET =     0x2;
    private static final byte NO_CALL = 0x1;

    // note that HET and NO_CALL are flippd from the documentation: that's because
    // plink actually reads these in backwards; and we want to use a shift operator
    // to put these in the appropriate location

    public void initialize() {
        vv.variantCollection = variantCollection;
        vv.dbsnp = dbsnp;
        vv.DO_NOT_VALIDATE_FILTERED = true;
        vv.type = ValidateVariants.ValidationType.REF;
        // write magic bits into the ped file
        try {
            outBed.write(new byte[] { (byte) 0x6c, (byte) 0x1b, 0x1 });
        } catch (IOException e) {
            throw new ReviewedStingException("error writing to output file.");
        }
        // write to the fam file, the first six columns of the standard ped file
        // first, load data from the input meta data file
        Map<String,Map<String,String>> metaValues = new HashMap<String,Map<String,String>>();
        try {
            if ( metaDataFile.getAbsolutePath().endsWith(".fam") ) {
                for ( String line : new XReadLines(metaDataFile) ) {
                    outFam.printf("%s%n",line);
                }
            } else {
                for ( String line : new XReadLines(metaDataFile) ) {
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
                outFam.printf("%s\t%s\t%s\t%s\t%s\t%s%n",fid,pid,sample,mid,sex,pheno);
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
        // write an entry into the map file
        outBim.printf("%s\t%s\t%.2f\t%d\t%s\t%s%n",vc.getChr(),getID(vc),APPROX_CM_PER_BP*vc.getStart(),vc.getStart(),
                vc.getReference().getBaseString(),vc.getAlternateAllele(0).getBaseString());
        // write an entry into the bed file
        int buf = 0;
        int idx = 0;
        byte out = 0x0;
        byte[] toWrite = new byte[1+(vc.getNSamples()/4)];
        for (Genotype g : vc.getGenotypes() ) {
            out |= getEncoding(g,buf);
            if ( buf == 3 ) {
                toWrite[idx] = out;
                buf = 0;
                out = 0x0;
                idx++;
            } else {
                buf++;
            }
        }
        if ( out != 0x0 ) {
            toWrite[idx]=out;
        }
        try {
            outBed.write(toWrite);
        } catch (IOException e) {
            throw new ReviewedStingException("Error writing to output file");
        }

        return 1;
    }

    public Integer reduce(Integer m, Integer r) {
        return r + m;
    }

    public Integer reduceInit() {
        return 0;
    }

    private static byte getEncoding(Genotype g, int offset) {
        byte b;
        if ( g.isHomRef() ) {
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

    private static String getID(VariantContext v) {
        if ( v.hasID() ) {
            return v.getID();
        } else {
            return String.format("SNP-%s-%d",v.getChr(),v.getStart());
        }
    }
}
