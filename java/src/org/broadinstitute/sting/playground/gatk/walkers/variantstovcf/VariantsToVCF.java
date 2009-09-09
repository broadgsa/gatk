package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;

import java.io.*;
import java.util.*;

public class VariantsToVCF extends RefWalker<Integer, Integer> {
    @Argument(fullName="vcfout", shortName="VO", doc="The output VCF file") public File VCF_OUT;
    @Argument(fullName="verbose", shortName="V", doc="Show extended output", required=false) public boolean VERBOSE = false;
    @Argument(fullName="suppress_multistate_alleles", shortName="SMA", doc="When multi-sample genotypes imply something other than a biallele state, suppress it.", required=false) public boolean SUPPRESS_MULTISTATE = false;

    private VCFWriter vcfwriter = null;
    private VCFHeader vcfheader = null;
    private TreeMap<String, String> sampleNames = null;

    public void initialize() {
        sampleNames = new TreeMap<String, String>();
        GATKArgumentCollection args = this.getToolkit().getArguments();

        for (String rodName : args.RODBindings) {
            //out.println(rodName);
            String[] rodPieces = rodName.split(",");
            String sampleName = rodPieces[0];

            if (sampleName.startsWith("NA"))
                sampleNames.put(sampleName.toUpperCase(), sampleName.toUpperCase());
        }

        vcfheader = getHeader(args, sampleNames.keySet());
        vcfwriter = new VCFWriter(vcfheader, VCF_OUT);
    }

    public static VCFHeader getHeader(GATKArgumentCollection args, Set<String> sampleNames) {
        Map<String, String> metaData = new HashMap<String, String>();
        List<String> additionalColumns = new ArrayList<String>();

        Calendar cal = Calendar.getInstance();

        metaData.put("format", "VCRv3.2");
        metaData.put("fileDate", String.format("%d%02d%02d", cal.get(Calendar.YEAR), cal.get(Calendar.MONTH), cal.get(Calendar.DAY_OF_MONTH)));
        metaData.put("source", "VariantsToVCF");
        metaData.put("reference", args.referenceFile.getAbsolutePath());

        additionalColumns.add("FORMAT");
        additionalColumns.addAll(sampleNames);

        return new VCFHeader(metaData, additionalColumns);
    }

    public boolean filter(RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        if (BaseUtils.simpleBaseToBaseIndex(ref) > -1) { return true; }

        for (ReferenceOrderedDatum rod : tracker.getAllRods()) {
            if (rod != null && sampleNames.keySet().contains(rod.getName().toUpperCase())) {
                return true;
            }
        }

        return false;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<VCFGenotypeRecord> gt = new ArrayList<VCFGenotypeRecord>();
        Map<VCFHeader.HEADER_FIELDS,String> map = new HashMap<VCFHeader.HEADER_FIELDS,String>();
        if ( generateVCFRecord(tracker, ref, context, vcfheader, gt, map, sampleNames, out, SUPPRESS_MULTISTATE, VERBOSE) ) {
            vcfwriter.addRecord(new VCFRecord(vcfheader, map, "GT:GQ:DP", gt));
            //vcfwriter.addRecord(new VCFRecord(vcfheader, map, "GT", gt));
            return  1;
        }
        return 0;
    }

    public static boolean generateVCFRecord(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context,
                                            VCFHeader vcfheader, List<VCFGenotypeRecord> gt, Map<VCFHeader.HEADER_FIELDS,String> map,
                                            Map<String, String> sampleNamesToRods, PrintStream out, boolean SUPPRESS_MULTISTATE, boolean VERBOSE) {
        int[] alleleCounts = new int[4];
        int numSNPs = 0;
        int numRefs = 0;
        int[] alleleNames = { 0, 1, 2, 3 };
        double snpQual = 0.0;
        int refbase = BaseUtils.simpleBaseToBaseIndex(ref.getBase());

        for (String name : vcfheader.getGenotypeSamples()) {
            ReferenceOrderedDatum rod = tracker.lookup(sampleNamesToRods.get(name), null);
            if (rod != null) {
                AllelicVariant av = (AllelicVariant) rod;
                String lod = String.format("%d", av.getVariationConfidence() > 99 ? 99 : (int) av.getVariationConfidence());
                int depth = 0;

                if (rod instanceof RodGeliText) {
                    RodGeliText rv = (RodGeliText) rod;
                    depth = rv.depth;
                }

                Map<String,String> str = new HashMap<String,String>();
                if (av.getGenotype().get(0).charAt(0) == av.getGenotype().get(0).charAt(1)) {
                    str.put("key","1/1:" + lod + (depth > 0 ? ":" + depth : ""));
                    //str.put("key","1/1");
                } else {
                    str.put("key","0/1:" + lod + (depth > 0 ? ":" + depth : ""));
                    //str.put("key","0/1");
                }
                
                List<String> alleles = av.getGenotype();

                int allele1 = BaseUtils.simpleBaseToBaseIndex(alleles.get(0).charAt(0));
                if (allele1 >= 0 && allele1 != refbase) { alleleCounts[allele1]++; }

                int allele2 = BaseUtils.simpleBaseToBaseIndex(alleles.get(0).charAt(1));
                if (allele2 >= 0 && allele2 != refbase) { alleleCounts[allele2]++; }

                gt.add(new VCFGenotypeRecord(name, str, alleles, VCFGenotypeRecord.PHASE.UNPHASED, ref.getBase()));

                numSNPs++;
                snpQual += av.getVariationConfidence();
            } else {
                Map<String,String> str = new HashMap<String,String>();
                str.put("key","0/0");

                List<String> alleles = new ArrayList<String>();
                alleles.add(ref.getBase() + "" + ref.getBase());

                gt.add(new VCFGenotypeRecord(name, str, alleles, VCFGenotypeRecord.PHASE.UNPHASED, ref.getBase()));
                
                numRefs++;
            }
        }

        if (numSNPs == 0)
            return false;


        Integer[] perm = Utils.SortPermutation(alleleCounts);
        int[] sortedCounts = Utils.PermuteArray(alleleCounts, perm);
        int[] sortedNames = Utils.PermuteArray(alleleNames, perm);

        rodDbSNP dbsnp = (rodDbSNP) tracker.lookup("dbsnp", null);

        String infoString = String.format("locus=%s ref=%c allele_count=( %c:%d %c:%d %c:%d %c:%d )",
            context.getLocation(),
            ref.getBase(),
            BaseUtils.baseIndexToSimpleBase(sortedNames[0]), sortedCounts[0],
            BaseUtils.baseIndexToSimpleBase(sortedNames[1]), sortedCounts[1],
            BaseUtils.baseIndexToSimpleBase(sortedNames[2]), sortedCounts[2],
            BaseUtils.baseIndexToSimpleBase(sortedNames[3]), sortedCounts[3]
        );

        if (SUPPRESS_MULTISTATE && sortedCounts[2] > 0) {
            out.println("[multistate] " + infoString);
            return false;
        } else {
            if (VERBOSE) {
                out.println("[locus_info] " + infoString);
            }
        }

        for (VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values()) {
            map.put(field,String.valueOf(1));

            if (field == VCFHeader.HEADER_FIELDS.CHROM) {
                map.put(field, context.getContig());
            } else if (field == VCFHeader.HEADER_FIELDS.POS) {
                map.put(field, String.valueOf(context.getPosition()));
            } else if (field == VCFHeader.HEADER_FIELDS.REF) {
                map.put(field, String.valueOf(ref.getBases()));
            } else if (field == VCFHeader.HEADER_FIELDS.ALT) {
                map.put(field, String.valueOf(BaseUtils.baseIndexToSimpleBase(sortedNames[3])));
            } else if (field == VCFHeader.HEADER_FIELDS.ID) {
                map.put(field, (dbsnp == null) ? "." : dbsnp.name);
            } else if (field == VCFHeader.HEADER_FIELDS.QUAL) {
                map.put(field, String.format("%d", snpQual > 99 ? 99 : (int) snpQual));
            } else if (field == VCFHeader.HEADER_FIELDS.FILTER) {
                map.put(field, String.valueOf("0"));
            } else if (field == VCFHeader.HEADER_FIELDS.INFO) {
                String infostr = ".";
                ArrayList<String> info = new ArrayList<String>();

                if (dbsnp != null) { info.add("DB=1"); }
                if (dbsnp != null && dbsnp.isHapmap()) { info.add("H2=1"); }

                for (int i = 0; i < info.size(); i++) {
                    if (i == 0) { infostr = ""; }

                    infostr += info.get(i);

                    if (i < info.size() - 1) {
                        infostr += ";";
                    }
                }

                map.put(field, infostr);
            }
        }
        return true;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer sum) {
        out.println("Processed " + sum + " variants.");

        vcfwriter.close();
    }
}
