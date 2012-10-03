package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 10/2/12
 * Time: 10:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class ExactAFCalculationPerformanceTest {
    final static Logger logger = Logger.getLogger(ExactAFCalculationPerformanceTest.class);

    private static abstract class Analysis {
        final GATKReport report;

        public Analysis(final String name, final List<String> columns) {
            report = GATKReport.newSimpleReport(name, columns);
        }

        public abstract void run(final ExactAFCalculationTestBuilder testBuilder,
                                 final List<Object> coreColumns);

        public String getName() {
            return getTable().getTableName();
        }

        public GATKReportTable getTable() {
            return report.getTables().iterator().next();
        }
    }

    private static class AnalyzeByACAndPL extends Analysis {
        public AnalyzeByACAndPL(final List<String> columns) {
            super("AnalyzeByACAndPL", Utils.append(columns, "non.type.pls", "ac", "n.alt.seg", "other.ac"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(10, 100, 1000) ) {
                final ExactAFCalculation calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                for ( int[] ACs : makeACs(testBuilder.numAltAlleles, testBuilder.nSamples*2) ) {
                    final VariantContext vc = testBuilder.makeACTest(ACs, nonTypePL);

                    timer.start();
                    final AlleleFrequencyCalculationResult result = calc.getLog10PNonRef(vc, priors);
                    final long runtime = timer.getElapsedTimeNano();

                    int otherAC = 0;
                    int nAltSeg = 0;
                    for ( int i = 0; i < ACs.length; i++ ) {
                        nAltSeg += ACs[i] > 0 ? 1 : 0;
                        if ( i > 0 ) otherAC += ACs[i];
                    }

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, result.getnEvaluations(), nonTypePL, ACs[0], nAltSeg, otherAC));
                    report.addRowList(columns);
                }
            }
        }

        private List<int[]> makeACs(final int nAltAlleles, final int nChrom) {
            if ( nAltAlleles > 2 ) throw new IllegalArgumentException("nAltAlleles must be < 3");

            final List<int[]> ACs = new LinkedList<int[]>();

            if ( nAltAlleles == 1 )
                for ( int i = 0; i < nChrom; i++ ) {
                    ACs.add(new int[]{i});
            } else if ( nAltAlleles == 2 ) {
                for ( int i = 0; i < nChrom; i++ ) {
                    for ( int j : Arrays.asList(0, 1, 5, 10, 50, 100, 1000, 10000, 100000) ) {
                        if ( j < nChrom - i )
                            ACs.add(new int[]{i, j});
                    }
                }
            } else {
                throw new IllegalStateException("cannot get here");
            }

            return ACs;
        }
    }

    private static class AnalyzeBySingletonPosition extends Analysis {
        public AnalyzeBySingletonPosition(final List<String> columns) {
            super("AnalyzeBySingletonPosition", Utils.append(columns, "non.type.pls", "position.of.singleton"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final ExactAFCalculation calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                final int[] ac = new int[testBuilder.numAltAlleles];
                ac[0] = 1;
                final VariantContext vc = testBuilder.makeACTest(ac, nonTypePL);

                for ( int position = 0; position < vc.getNSamples(); position++ ) {
                    final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                    final List<Genotype> genotypes = new ArrayList<Genotype>(vc.getGenotypes());
                    Collections.rotate(genotypes, position);
                    vcb.genotypes(genotypes);

                    timer.start();
                    final AlleleFrequencyCalculationResult result = calc.getLog10PNonRef(vcb.make(), priors);
                    final long runtime = timer.getElapsedTimeNano();

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, result.getnEvaluations(), nonTypePL, position));
                    report.addRowList(columns);
                }
            }
        }
    }

    private static class AnalyzeByNonInformative extends Analysis {
        public AnalyzeByNonInformative(final List<String> columns) {
            super("AnalyzeByNonInformative", Utils.append(columns, "non.type.pls", "n.non.informative"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final ExactAFCalculation calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                final int[] ac = new int[testBuilder.numAltAlleles];
                ac[0] = 1;
                final VariantContext vc = testBuilder.makeACTest(ac, nonTypePL);
                final Genotype nonInformative = testBuilder.makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0);

                for ( int nNonInformative = 0; nNonInformative < vc.getNSamples(); nNonInformative++ ) {
                    final VariantContextBuilder vcb = new VariantContextBuilder(vc);

                    final List<Genotype> genotypes = new ArrayList<Genotype>();
                    genotypes.addAll(vc.getGenotypes().subList(0, nNonInformative + 1));
                    genotypes.addAll(Collections.nCopies(vc.getNSamples() - nNonInformative, nonInformative));
                    vcb.genotypes(genotypes);

                    timer.start();
                    final AlleleFrequencyCalculationResult result = calc.getLog10PNonRef(vcb.make(), priors);
                    final long runtime = timer.getElapsedTimeNano();

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, result.getnEvaluations(), nonTypePL, nNonInformative));
                    report.addRowList(columns);
                }
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        logger.addAppender(new ConsoleAppender(new SimpleLayout()));

        final List<String> coreColumns = Arrays.asList("iteration", "n.alt.alleles", "n.samples",
                "exact.model", "prior.type", "runtime", "n.evaluations");

        final PrintStream out = new PrintStream(new FileOutputStream(args[0]));

        final boolean USE_GENERAL = false;
        final List<ExactAFCalculationTestBuilder.ModelType> modelTypes = USE_GENERAL
                ? Arrays.asList(ExactAFCalculationTestBuilder.ModelType.values())
                : Arrays.asList(ExactAFCalculationTestBuilder.ModelType.DiploidExact);

        final boolean ONLY_HUMAN_PRIORS = false;
        final List<ExactAFCalculationTestBuilder.PriorType> priorTypes = ONLY_HUMAN_PRIORS
                ? Arrays.asList(ExactAFCalculationTestBuilder.PriorType.values())
                : Arrays.asList(ExactAFCalculationTestBuilder.PriorType.human);

        final int MAX_N_SAMPLES_FOR_MULTI_ALLELIC = 100;

        final List<Analysis> analyzes = new ArrayList<Analysis>();
        analyzes.add(new AnalyzeByACAndPL(coreColumns));
        analyzes.add(new AnalyzeBySingletonPosition(coreColumns));
        analyzes.add(new AnalyzeByNonInformative(coreColumns));

        for ( int iteration = 0; iteration < 1; iteration++ ) {
            for ( final int nAltAlleles : Arrays.asList(1, 2) ) {
                for ( final int nSamples : Arrays.asList(1, 10, 100, 1000, 10000) ) {
                    if ( nSamples > MAX_N_SAMPLES_FOR_MULTI_ALLELIC && nAltAlleles > 1 )
                        continue; // skip things that will take forever!

                    for ( final ExactAFCalculationTestBuilder.ModelType modelType : modelTypes ) {
                        for ( final ExactAFCalculationTestBuilder.PriorType priorType : priorTypes ) {
                            final ExactAFCalculationTestBuilder testBuilder
                                    = new ExactAFCalculationTestBuilder(nSamples, nAltAlleles, modelType, priorType);

                            for ( final Analysis analysis : analyzes ) {
                                logger.info(Utils.join("\t", Arrays.asList(iteration, nAltAlleles, nSamples, modelType, priorType, analysis.getName())));
                                final List<?> values = Arrays.asList(iteration, nAltAlleles, nSamples, modelType, priorType);
                                analysis.run(testBuilder, (List<Object>)values);
                            }
                        }
                    }
                }
            }
        }

        final GATKReport report = new GATKReport();
        for ( final Analysis analysis : analyzes )
            report.addTable(analysis.getTable());
        report.print(out);
        out.close();
    }
}