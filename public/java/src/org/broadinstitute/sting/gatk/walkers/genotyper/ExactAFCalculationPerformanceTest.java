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
            super("AnalyzeByACAndPL", Utils.append(columns, "non.type.pls", "ac"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(10, 100, 1000) ) {
                final ExactAFCalculation calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                for ( int ac = 0; ac < testBuilder.getnSamples(); ac++ ) {
                    final VariantContext vc = testBuilder.makeACTest(ac, nonTypePL);

                    timer.start();
                    final AlleleFrequencyCalculationResult result = calc.getLog10PNonRef(vc, priors);
                    final long runtime = timer.getElapsedTimeNano();

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, result.getnEvaluations(), nonTypePL, ac));
                    report.addRowList(columns);
                }
            }
        }
    }

    private static class AnalyzeBySingletonPosition extends Analysis {
        public AnalyzeBySingletonPosition(final List<String> columns) {
            super("AnalyzeBySingletonPosition", Utils.append(columns, "non.type.pls", "position.of.singleton"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(10, 100, 1000) ) {
                final ExactAFCalculation calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                int ac = 1;
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

            for ( final int nonTypePL : Arrays.asList(10, 100, 1000) ) {
                final ExactAFCalculation calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                int ac = 1;
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

        final List<Analysis> analyzes = new ArrayList<Analysis>();
        analyzes.add(new AnalyzeByACAndPL(coreColumns));
        analyzes.add(new AnalyzeBySingletonPosition(coreColumns));
        analyzes.add(new AnalyzeByNonInformative(coreColumns));

        for ( int iteration = 0; iteration < 1; iteration++ ) {
            for ( final int nAltAlleles : Arrays.asList(1) ) {
                for ( final int nSamples : Arrays.asList(1, 10, 100) ) {
                    for ( final ExactAFCalculationTestBuilder.ModelType modelType : modelTypes ) {
                        for ( final ExactAFCalculationTestBuilder.PriorType priorType : priorTypes ) {
                            final ExactAFCalculationTestBuilder testBuilder
                                    = new ExactAFCalculationTestBuilder(nSamples, 1, modelType, priorType);

                            for ( final Analysis analysis : analyzes ) {
                                logger.info(Utils.join("\t", Arrays.asList(iteration, nSamples, modelType, priorType, analysis.getName())));
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