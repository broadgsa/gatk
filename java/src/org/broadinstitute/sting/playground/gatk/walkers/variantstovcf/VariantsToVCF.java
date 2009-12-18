package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodGeliText;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Converts variants from other file formats (anything that implements the Variation interface) to VCF format.
 */
public class VariantsToVCF extends RefWalker<Integer, Integer> {
    @Argument(fullName="vcfout", shortName="VO", doc="The output VCF file") public File VCF_OUT;
    @Argument(fullName="verbose", shortName="V", doc="Show extended output", required=false) public boolean VERBOSE = false;
    @Argument(fullName="suppress_multistate_alleles", shortName="SMA", doc="When multi-sample genotypes imply something other than a biallele state, suppress it.", required=false) public boolean SUPPRESS_MULTISTATE = false;

    private VCFWriter vcfwriter = null;
    private VCFHeader vcfheader = null;
    private TreeMap<String, String> sampleNames = null;
    private static String format = "GT:GQ:DP";

    public void initialize() {
        sampleNames = new TreeMap<String, String>();
        GATKArgumentCollection args = this.getToolkit().getArguments();

        for (String rodName : args.RODBindings) {
            //out.println(rodName);
            String[] rodPieces = rodName.split(",");
            String sampleName = rodPieces[0];

            if (! ( sampleName.equals("dbSNP") || sampleName.equals("interval") ) ) {
                sampleNames.put(sampleName, sampleName);
            }

        }

        vcfheader = getHeader(args, sampleNames.keySet());
        vcfwriter = new VCFWriter(VCF_OUT);
        vcfwriter.writeHeader(vcfheader);
    }

    public static VCFHeader getHeader(GATKArgumentCollection args, Set<String> sampleNames) {

        // Don't output the data for now because it kills our unit test MD5s and is optional
        // TODO - figure out what to do here
        //Calendar cal = Calendar.getInstance();
        //metaData.put("fileDate", String.format("%d%02d%02d", cal.get(Calendar.YEAR), cal.get(Calendar.MONTH), cal.get(Calendar.DAY_OF_MONTH)));

        Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
        metaData.add(new VCFHeaderLine("source", "VariantsToVCF"));
        metaData.add(new VCFHeaderLine("reference", args.referenceFile.getAbsolutePath()));

        Set<String> additionalColumns = new HashSet<String>();
        additionalColumns.add("FORMAT");
        additionalColumns.addAll(sampleNames);

        return new VCFHeader(metaData, additionalColumns);
    }

    public boolean filter(RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        if (BaseUtils.simpleBaseToBaseIndex(ref) > -1) {
            return true;
        }

        for (ReferenceOrderedDatum rod : tracker.getAllRods()) {
            if (rod != null && sampleNames.keySet().contains(rod.getName().toUpperCase())) {
                return true;
            }
        }

        return false;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<VCFGenotypeRecord> gt = new ArrayList<VCFGenotypeRecord>();
        Map<VCFHeader.HEADER_FIELDS, String> map = new HashMap<VCFHeader.HEADER_FIELDS, String>();
        VCFRecord rec = generateVCFRecord(tracker, ref, context, vcfheader, gt, map, sampleNames, out, SUPPRESS_MULTISTATE, VERBOSE);
        if (rec != null) {
            vcfwriter.addRecord(rec);
            return 1;
        }
        return 0;
    }

    public static VCFRecord generateVCFRecord(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context,
                                              VCFHeader vcfheader, List<VCFGenotypeRecord> gt, Map<VCFHeader.HEADER_FIELDS, String> map,
                                              Map<String, String> sampleNamesToRods, PrintStream out, boolean SUPPRESS_MULTISTATE, boolean VERBOSE) {
        int[] alleleCounts = new int[4];
        int numSNPs = 0;
        int numRefs = 0;
        int[] alleleNames = {0, 1, 2, 3};
        double snpQual = 0.0;
        int refbase = BaseUtils.simpleBaseToBaseIndex(ref.getBase());
        List<VCFGenotypeEncoding> alts = new ArrayList<VCFGenotypeEncoding>();
        for (String name : vcfheader.getGenotypeSamples()) {
            ReferenceOrderedDatum rod = tracker.lookup(sampleNamesToRods.get(name), null);
            if (rod != null) {
                Variation av = (Variation) rod;
                String lod = String.format("%d", av.getNegLog10PError() > 99 ? 99 : (int) av.getNegLog10PError());
                int depth = 0;

                if (rod instanceof RodGeliText) {
                    RodGeliText rv = (RodGeliText) rod;
                    depth = rv.depth;
                }
                if (!(rod instanceof VariantBackedByGenotype))
                    throw new IllegalArgumentException("The passed in variant type must be backed by genotype data");
                Genotype genotype = ((VariantBackedByGenotype) rod).getCalledGenotype();
                List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
                for (char base : genotype.getBases().toCharArray()) {
                    alleles.add(new VCFGenotypeEncoding(String.valueOf(base)));
                    if (base != ref.getBase() && !alts.contains(String.valueOf(base))) alts.add(new VCFGenotypeEncoding(String.valueOf(base)));
                }
                int allele1 = BaseUtils.simpleBaseToBaseIndex(genotype.getBases().charAt(0));
                int allele2 = BaseUtils.simpleBaseToBaseIndex(genotype.getBases().charAt(1));
                if (allele1 >= 0 && allele1 != refbase) {
                    alleleCounts[allele1]++;
                }
                if (allele2 >= 0 && allele2 != refbase) {
                    alleleCounts[allele2]++;
                }
                Map<String, String> str = new HashMap<String, String>();
                str.put("GQ", lod);
                if (depth > 0) str.put("DP", String.valueOf(depth));

                gt.add(new VCFGenotypeRecord(name, alleles, VCFGenotypeRecord.PHASE.UNPHASED, str));

                numSNPs++;
                snpQual += av.getNegLog10PError();
            } else if (BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1) {
                Map<String, String> str = new HashMap<String, String>();
                List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
                alleles.add(new VCFGenotypeEncoding(String.valueOf(ref.getBase())));
                alleles.add(new VCFGenotypeEncoding(String.valueOf(ref.getBase())));
                gt.add(new VCFGenotypeRecord(name, alleles, VCFGenotypeRecord.PHASE.UNPHASED, str));

                numRefs++;
            }
        }

        if (numSNPs == 0)
            return null;


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
            return null;
        } else {
            if (VERBOSE) {
                out.println("[locus_info] " + infoString);
            }
        }

        Map<String,String> info = new HashMap<String,String>();
        if (dbsnp != null) info.put("DB","1");
        if (dbsnp != null && dbsnp.isHapmap()) info.put("H2","1");

        return new VCFRecord(ref.getBase(),
                             context.getContig(),
                             (int) context.getPosition(),
                             (dbsnp == null) ? "." : dbsnp.getRS_ID(),
                             alts,
                             snpQual > 99 ? 99 : (int) snpQual,
                             "0",
                             info,
                             format,
                             gt);

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
