package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import com.google.java.contract.Requires;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.*;
import java.util.*;

/**
 * Allows us to write out and read in information about exact calls (site, alleles, PLs, etc) in tabular format
 *
 * Once opened, calls can be writen to disk with printCallInfo
 */
public class ExactCallLogger implements Cloneable {
    private PrintStream callReport = null;

    /**
     * Create a new ExactCallLogger writing it's output to outputFile
     *
     * @param outputFile
     */
    public ExactCallLogger(final File outputFile) {
        try {
            callReport = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile), 10000000));
            callReport.println(Utils.join("\t", Arrays.asList("loc", "variable", "key", "value")));
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Summarizes information about an exact call that happened
     */
    public static class ExactCall {
        final VariantContext vc;
        final long runtime;
        final AFCalcResult originalCall;

        public ExactCall(VariantContext vc, final long runtime, final AFCalcResult originalCall) {
            this.vc = vc;
            this.runtime = runtime;
            this.originalCall = originalCall;
        }

        @Override
        public String toString() {
            return String.format("ExactCall %s:%d alleles=%s nSamples=%s orig.pNonRef=%.2f orig.runtime=%s",
                    vc.getChr(), vc.getStart(), vc.getAlleles(), vc.getNSamples(),
                    originalCall.getLog10PosteriorOfAFGT0(),
                    new AutoFormattingTime(runtime / 1e9).toString());
        }
    }

    protected final void printCallInfo(final VariantContext vc,
                                       final double[] log10AlleleFrequencyPriors,
                                       final long runtimeNano,
                                       final AFCalcResult result) {
        printCallElement(vc, "type", "ignore", vc.getType());

        int allelei = 0;
        for (final Allele a : vc.getAlleles())
            printCallElement(vc, "allele", allelei++, a.getDisplayString());

        for (final Genotype g : vc.getGenotypes())
            printCallElement(vc, "PL", g.getSampleName(), g.getLikelihoodsString());

        for (int priorI = 0; priorI < log10AlleleFrequencyPriors.length; priorI++)
            printCallElement(vc, "priorI", priorI, log10AlleleFrequencyPriors[priorI]);

        printCallElement(vc, "runtime.nano", "ignore", runtimeNano);
        printCallElement(vc, "log10PosteriorOfAFEq0", "ignore", result.getLog10PosteriorOfAFEq0());
        printCallElement(vc, "log10PosteriorOfAFGt0", "ignore", result.getLog10PosteriorOfAFGT0());

        for ( final Allele allele : result.getAllelesUsedInGenotyping() ) {
            if ( allele.isNonReference() ) {
                printCallElement(vc, "MLE", allele, result.getAlleleCountAtMLE(allele));
                printCallElement(vc, "pRefByAllele", allele, result.getLog10PosteriorOfAFEq0ForAllele(allele));
            }
        }

        callReport.flush();
    }

    @Requires({"vc != null", "variable != null", "key != null", "value != null", "callReport != null"})
    private void printCallElement(final VariantContext vc,
                                  final Object variable,
                                  final Object key,
                                  final Object value) {
        final String loc = String.format("%s:%d", vc.getChr(), vc.getStart());
        callReport.println(Utils.join("\t", Arrays.asList(loc, variable, key, value)));
    }

    /**
     * Read in a list of ExactCall objects from reader, keeping only those
     * with starts in startsToKeep or all sites (if this is empty)
     *
     * @param reader a just-opened reader sitting at the start of the file
     * @param startsToKeep a list of start position of the calls to keep, or empty if all calls should be kept
     * @param parser a genome loc parser to create genome locs
     * @return a list of ExactCall objects in reader
     * @throws IOException
     */
    public static List<ExactCall> readExactLog(final BufferedReader reader, final List<Integer> startsToKeep, GenomeLocParser parser) throws IOException {
        if ( reader == null ) throw new IllegalArgumentException("reader cannot be null");
        if ( startsToKeep == null ) throw new IllegalArgumentException("startsToKeep cannot be null");
        if ( parser == null ) throw new IllegalArgumentException("GenomeLocParser cannot be null");

        List<ExactCall> calls = new LinkedList<ExactCall>();

        // skip the header line
        reader.readLine();

        // skip the first "type" line
        reader.readLine();

        while (true) {
            final VariantContextBuilder builder = new VariantContextBuilder();
            final List<Allele> alleles = new ArrayList<Allele>();
            final List<Genotype> genotypes = new ArrayList<Genotype>();
            final double[] posteriors = new double[2];
            final double[] priors = MathUtils.normalizeFromLog10(new double[]{0.5, 0.5}, true);
            final List<Integer> mle = new ArrayList<Integer>();
            final Map<Allele, Double> log10pRefByAllele = new HashMap<Allele, Double>();
            long runtimeNano = -1;

            GenomeLoc currentLoc = null;
            while (true) {
                final String line = reader.readLine();
                if (line == null)
                    return calls;

                final String[] parts = line.split("\t");
                final GenomeLoc lineLoc = parser.parseGenomeLoc(parts[0]);
                final String variable = parts[1];
                final String key = parts[2];
                final String value = parts[3];

                if (currentLoc == null)
                    currentLoc = lineLoc;

                if (variable.equals("type")) {
                    if (startsToKeep.isEmpty() || startsToKeep.contains(currentLoc.getStart())) {
                        builder.alleles(alleles);
                        final int stop = currentLoc.getStart() + alleles.get(0).length() - 1;
                        builder.chr(currentLoc.getContig()).start(currentLoc.getStart()).stop(stop);
                        builder.genotypes(genotypes);
                        final int[] mleInts = ArrayUtils.toPrimitive(mle.toArray(new Integer[]{}));
                        final AFCalcResult result = new AFCalcResult(mleInts, 1, alleles, posteriors, priors, log10pRefByAllele);
                        calls.add(new ExactCall(builder.make(), runtimeNano, result));
                    }
                    break;
                } else if (variable.equals("allele")) {
                    final boolean isRef = key.equals("0");
                    alleles.add(Allele.create(value, isRef));
                } else if (variable.equals("PL")) {
                    final GenotypeBuilder gb = new GenotypeBuilder(key);
                    gb.PL(GenotypeLikelihoods.fromPLField(value).getAsPLs());
                    genotypes.add(gb.make());
                } else if (variable.equals("log10PosteriorOfAFEq0")) {
                    posteriors[0] = Double.valueOf(value);
                } else if (variable.equals("log10PosteriorOfAFGt0")) {
                    posteriors[1] = Double.valueOf(value);
                } else if (variable.equals("MLE")) {
                    mle.add(Integer.valueOf(value));
                } else if (variable.equals("pRefByAllele")) {
                    final Allele a = Allele.create(key);
                    log10pRefByAllele.put(a, Double.valueOf(value));
                } else if (variable.equals("runtime.nano")) {
                    runtimeNano = Long.valueOf(value);
                } else {
                    // nothing to do
                }
            }
        }
    }
}