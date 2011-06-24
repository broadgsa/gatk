package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Given a VCF file and one or more scenarios for affected individuals, calculates the probability that a given site's genotypes
 * are consistent with the expected pattern for a given disease model.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
public class AffectedConsistencyWalker extends RodWalker<Integer, Integer> {
    public enum DiseaseModel { DOMINANT, RECESSIVE }

    @Output
    public PrintStream out;

    @Argument(fullName="affected", shortName="A", doc="A scenario file (or files) for affected individuals.  Scenarios are specified with an identifier and a comma-separated list of samples (e.g. Pedigree_1 sample1,sample2,sample3).  Each line is another scenario.", required=true)
    public String[] AFFECTED_SAMPLE_SCENARIOS;

    @Argument(fullName="diseaseModel", shortName="DM", doc="The disease model (DOMINANT or RECESSIVE)", required=true)
    public DiseaseModel DISEASE_MODEL;

    @Argument(fullName="verbose", shortName="V", doc="If specified, enable verbose mode with a lot of output useful for debugging", required=false)
    public PrintStream VERBOSE_WRITER = null;

    public Map<String, Set<String>> sampleScenarios;
    public GATKReport consistencyReport;
    public Set<String> availableSamples;

    private Map<String, Set<String>> loadAffectedSampleScenarios() {
        // Load all the specified sample scenarios specified in one or more files
        ArrayList<String> scenarioStrings = new ArrayList<String>();

        for (String affectedSampleScenario : AFFECTED_SAMPLE_SCENARIOS) {
            File affectedSampleScenarioFile = new File(affectedSampleScenario);

            try {
                XReadLines lineReader = new XReadLines(affectedSampleScenarioFile);

                for (String line : lineReader) {
                    // Ignore commented-out lines
                    if (!line.contains("#")) {
                        scenarioStrings.add(line);
                    }
                }
            } catch (FileNotFoundException e) {
                throw new UserException(String.format("The scenario file '%s' was not found", affectedSampleScenarioFile.getAbsolutePath()));
            }
        }

        // Parse all the sample scenario strings (comma- or white-space-separated sample lists)
        Map<String, Set<String>> scenarios = new HashMap<String, Set<String>>();

        for (String scenarioString : scenarioStrings) {
            String[] pieces = scenarioString.split("[\\s]+");

            if (pieces.length != 2) {
                throw new UserException(
                        String.format("The scenario line '%s' could not be understood.  Please make sure that your " +
                                      "scenario file has only two columns: the first being an arbitrary scenario id " +
                                      "(e.g. 'Pedigree_1') and the second being a comma-separated list of samples " +
                                      "(e.g. 'sample1,sample2,sample3')",
                                      scenarioString
                                      )
                        );
            }

            String scenarioId = pieces[0];

            String[] sampleNames = pieces[1].split(",");

            Set<String> samples = new HashSet<String>();
            for (String sample : sampleNames) {
                if (!availableSamples.contains(sample)) {
                    throw new UserException(
                            String.format("The sample '%s' was not found in the ROD bound as the 'variant' track " +
                                          "(i.e. the file that was supplied via '-B:variant,VCF /path/to/my.vcf'). " +
                                          "Please make sure all samples specified for processing are present in " +
                                          "your VCF.",
                                          sample)
                    );
                } else {
                    samples.add(sample);
                }
            }

            scenarios.put(scenarioId, samples);
        }

        if (scenarios.size() == 0) {
            throw new UserException("There were no scenarios specified.  Please specify at least one set of affected samples.");
        }

        return scenarios;
    }

    public void initialize() {
        // Figure out what samples I can possibly have (from the bound VCF file)
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add("variant");

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        availableSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        // Load the scenarios to consider
        sampleScenarios = loadAffectedSampleScenarios();

        // Prepare the output report
        consistencyReport = new GATKReport();
        consistencyReport.addTable("AffectedConsistency", "Table of results indicating if the observed genotypes matched the expected genotypes");

        GATKReportTable table = consistencyReport.getTable("AffectedConsistency");
        table.addPrimaryKey("locus_and_scenario", false);
        table.addColumn("chr", "unknown");
        table.addColumn("start", 0);
        table.addColumn("scenario", "unknown");
        table.addColumn("P_of_C_given_DM_is_true", 0.0);
        table.addColumn("P_of_C_given_DM_is_false", 0.0);
        table.addColumn("OR_DM_is_true_vs_DM_is_false", 0.0);

        for ( String sample : availableSamples ) {
            table.addColumn(sample, "unknown");
        }

        if (VERBOSE_WRITER != null) {
            VERBOSE_WRITER.println("This is a test of the verbose writer");
        }

        System.exit(0);
    }

    private VariantContext getExpectedGenotypeConfiguration(Set<String> affectedSamples, VariantContext obs) {
        List<Allele> homRefAlleles = new ArrayList<Allele>();
        homRefAlleles.add(obs.getReference());
        homRefAlleles.add(obs.getReference());

        List<Allele> hetAlleles = new ArrayList<Allele>();
        hetAlleles.add(obs.getReference());
        hetAlleles.add(obs.getAlternateAllele(0));

        List<Allele> homVarAlleles = new ArrayList<Allele>();
        homVarAlleles.add(obs.getAlternateAllele(0));
        homVarAlleles.add(obs.getAlternateAllele(0));

        Collection<Genotype> expectedGenotypes = new ArrayList<Genotype>();
        for ( String sample : obs.getSampleNames() ) {
            Genotype expectedGenotype = new Genotype(sample, homRefAlleles);
            if (affectedSamples.contains(sample)) {
                expectedGenotype = (DISEASE_MODEL == DiseaseModel.DOMINANT) ? new Genotype(sample, hetAlleles) : new Genotype(sample, homVarAlleles);
            }

            expectedGenotypes.add(expectedGenotype);
        }

        return new VariantContext("expected", obs.getChr(), obs.getStart(), obs.getEnd(), obs.getAlleles(), expectedGenotypes);
    }

    private double getLogLikelihoodOfDiseaseModelHypothesis(VariantContext obs, VariantContext exp, boolean diseaseModelIsSupported) {
        return getLogLikelihoodOfDiseaseModelHypothesis(obs, exp, diseaseModelIsSupported, 0, 0.0);
    }

    private double getLogLikelihoodOfDiseaseModelHypothesis(VariantContext obs, VariantContext exp, boolean diseaseModelIsSupported, int sampleIndex, double logLikelihoodSoFar) {
        if (sampleIndex < exp.getNSamples()) {
            Genotype expGenotype = exp.getGenotype(sampleIndex);
            Genotype obsGenotype = obs.getGenotype(sampleIndex);

            if (obsGenotype.hasLikelihoods()) {
                double[] normalizedLikelihoods = MathUtils.normalizeFromLog10(obsGenotype.getLikelihoods().getAsVector());
                boolean[] expectedGenotypes = { expGenotype.isHomRef(), expGenotype.isHet(), expGenotype.isHomVar() };

                for (int i = 0; i < 3; i++) {
                    if (expectedGenotypes[i] == diseaseModelIsSupported) {
                        return getLogLikelihoodOfDiseaseModelHypothesis(obs, exp, diseaseModelIsSupported, sampleIndex + 1, logLikelihoodSoFar + Math.log10(normalizedLikelihoods[i]));
                    }
                }
            } else {
                return getLogLikelihoodOfDiseaseModelHypothesis(obs, exp, diseaseModelIsSupported, sampleIndex + 1, logLikelihoodSoFar);
            }
        }

        return logLikelihoodSoFar;
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, null, ref.getLocus(), true, true);

            if (vcs.size() == 1) {
                VariantContext obs = vcs.iterator().next();

                for (String scenarioId : sampleScenarios.keySet()) {
                    Set<String> affectedSamples = sampleScenarios.get(scenarioId);

                    VariantContext exp = getExpectedGenotypeConfiguration(affectedSamples, obs);

                    /*
                    GATKReport report = new GATKReport();

                    String reportName = String.format("GenotypeTable_%s_%s_%s", scenarioId, ref.getLocus().getContig(), ref.getLocus().getStart());
                    String reportDesc = String.format("Info for scenario %s at locus %s", scenarioId, ref.getLocus());
                    report.addTable(reportName, reportDesc);

                    GATKReportTable table = report.getTable(reportName);

                    table.addPrimaryKey("sample_pk", false);
                    table.addColumn("table", "unknown");
                    table.addColumn("sample", "unknown");
                    table.addColumn("affected", false);
                    table.addColumn("homref_prob", "unknown");
                    table.addColumn("het_prob", "unknown");
                    table.addColumn("homvar_prob", "unknown");
                    table.addColumn("observed_genotype", "unknown");
                    table.addColumn("expected_genotype", "unknown");

                    for (String sample : obs.getSampleNames()) {
                        double[] normalizedLikelihoods = {0.0, 0.0, 0.0};
                        if (obs.getGenotype(sample).hasLikelihoods()) {
                            normalizedLikelihoods = MathUtils.normalizeFromLog10(obs.getGenotype(sample).getLikelihoods().getAsVector());
                        }

                        table.set(sample, "table", reportName);
                        table.set(sample, "sample", sample);
                        table.set(sample, "affected", affectedSamples.contains(sample));
                        table.set(sample, "homref_prob", normalizedLikelihoods[0]);
                        table.set(sample, "het_prob", normalizedLikelihoods[1]);
                        table.set(sample, "homvar_prob", normalizedLikelihoods[2]);
                        table.set(sample, "observed_genotype", obs.getGenotype(sample).getGenotypeString());
                        table.set(sample, "expected_genotype", exp.getGenotype(sample).getGenotypeString());
                    }

                    report.print(out);
                    */

                    double logLikelihoodThatDiseaseModelIsSupported = getLogLikelihoodOfDiseaseModelHypothesis(obs, exp, true);
                    double logLikelihoodThatDiseaseModelIsNotSupported = getLogLikelihoodOfDiseaseModelHypothesis(obs, exp, false);
                    double logOddsRatioThatDiseaseModelIsSupported = logLikelihoodThatDiseaseModelIsSupported / logLikelihoodThatDiseaseModelIsNotSupported;

                    String key = String.format("%s_%s_%s", ref.getLocus().getContig(), ref.getLocus().getStart(), scenarioId);

                    consistencyReport.getTable("AffectedConsistency").set(key, "scenario", scenarioId);
                    consistencyReport.getTable("AffectedConsistency").set(key, "chr", ref.getLocus().getContig());
                    consistencyReport.getTable("AffectedConsistency").set(key, "start", ref.getLocus().getStart());
                    consistencyReport.getTable("AffectedConsistency").set(key, "P_of_C_given_DM_is_true", logLikelihoodThatDiseaseModelIsSupported);
                    consistencyReport.getTable("AffectedConsistency").set(key, "P_of_C_given_DM_is_false", logLikelihoodThatDiseaseModelIsNotSupported);
                    consistencyReport.getTable("AffectedConsistency").set(key, "OR_DM_is_true_vs_DM_is_false", logOddsRatioThatDiseaseModelIsSupported);

                    for ( String sample : availableSamples ) {
                        String obsAndExpectedGenotypes = String.format("%s;%s", obs.getGenotype(sample).getGenotypeString(), exp.getGenotype(sample).getGenotypeString());
                        consistencyReport.getTable("AffectedConsistency").set(key, sample, obsAndExpectedGenotypes);
                    }
                }
            }
        }

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer result) {
        consistencyReport.print(out);
    }
}
