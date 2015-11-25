/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.varianteval.util;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.StandardEval;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.RequiredStratification;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.StandardStratification;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.util.*;

public class VariantEvalUtils {
    private final VariantEval variantEvalWalker;
    Logger logger;

    public VariantEvalUtils(VariantEval variantEvalWalker) {
        this.variantEvalWalker = variantEvalWalker;
        this.logger = variantEvalWalker.getLogger();
    }

    /**
     * List all of the available evaluation modules, then exit successfully
     */
    public void listModulesAndExit() {
        List<Class<? extends VariantStratifier>> vsClasses = new PluginManager<VariantStratifier>(VariantStratifier.class).getPlugins();
        List<Class<? extends VariantEvaluator>> veClasses = new PluginManager<VariantEvaluator>(VariantEvaluator.class).getPlugins();

        logger.info("Available stratification modules:");
        logger.info("(Standard modules are starred)");
        for (Class<? extends VariantStratifier> vsClass : vsClasses) {
            logger.info("\t" + vsClass.getSimpleName() + (RequiredStratification.class.isAssignableFrom(vsClass) || StandardStratification.class.isAssignableFrom(vsClass) ? "*" : ""));
        }
        logger.info("");

        logger.info("Available evaluation modules:");
        logger.info("(Standard modules are starred)");
        for (Class<? extends VariantEvaluator> veClass : veClasses) {
            logger.info("\t" + veClass.getSimpleName() + (StandardEval.class.isAssignableFrom(veClass) ? "*" : ""));
        }
        logger.info("");

        System.exit(0);
    }

    /**
     * Initialize required, standard and user-specified stratification objects
     *
     * @param noStandardStrats  don't use the standard stratifications
     * @param modulesToUse      the list of stratification modules to use
     * @return set of stratifications to use
     */
    public List<VariantStratifier> initializeStratificationObjects(boolean noStandardStrats, String[] modulesToUse) {
        TreeSet<VariantStratifier> strats = new TreeSet<VariantStratifier>();
        Set<String> stratsToUse = new HashSet<String>();

        // Create a map for all stratification modules for easy lookup.
        HashMap<String, Class<? extends VariantStratifier>> classMap = new HashMap<String, Class<? extends VariantStratifier>>();
        for (Class<? extends VariantStratifier> c : new PluginManager<VariantStratifier>(VariantStratifier.class).getPlugins()) {
            classMap.put(c.getSimpleName(), c);
        }

        // We must use all required stratification modules.
        for (Class<? extends RequiredStratification> reqClass : new PluginManager<RequiredStratification>(RequiredStratification.class).getPlugins()) {
            if (classMap.containsKey(reqClass.getSimpleName())) {
                stratsToUse.add(reqClass.getSimpleName());
            }
        }

        // By default, use standard stratification modules.
        if (!noStandardStrats) {
            for (Class<? extends StandardStratification> stdClass : new PluginManager<StandardStratification>(StandardStratification.class).getPlugins()) {
                if (classMap.containsKey(stdClass.getSimpleName())) {
                    stratsToUse.add(stdClass.getSimpleName());
                }
            }
        }

        // Now add the user-selected modules
        stratsToUse.addAll(Arrays.asList(modulesToUse));

        // Instantiate the stratifications
        for (String module : stratsToUse) {
            if (!classMap.containsKey(module)) {
                throw new UserException.CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if (classMap.containsKey(module)) {
                Class<? extends VariantStratifier> c = classMap.get(module);

                try {
                    VariantStratifier vs = c.newInstance();
                    vs.setVariantEvalWalker(variantEvalWalker);
                    vs.initialize();

                    strats.add(vs);
                } catch (InstantiationException e) {
                    throw new GATKException("Unable to instantiate stratification module '" + c.getSimpleName() + "'");
                } catch (IllegalAccessException e) {
                    throw new GATKException("Illegal access error when trying to instantiate stratification module '" + c.getSimpleName() + "'");
                }
            }
        }

        return new ArrayList<VariantStratifier>(strats);
    }

    /**
     * Initialize required, standard and user-specified evaluation objects
     *
     * @param noStandardEvals don't use the standard evaluations
     * @param modulesToUse    the list of evaluation modules to use
     * @return set of evaluations to use
     */
    public Set<Class<? extends VariantEvaluator>> initializeEvaluationObjects(boolean noStandardEvals, String[] modulesToUse) {
        Set<Class<? extends VariantEvaluator>> evals = new HashSet<Class<? extends VariantEvaluator>>();

        // Create a map for all eval modules for easy lookup.
        HashMap<String, Class<? extends VariantEvaluator>> classMap = new HashMap<String, Class<? extends VariantEvaluator>>();
        for (Class<? extends VariantEvaluator> c : new PluginManager<VariantEvaluator>(VariantEvaluator.class).getPlugins()) {
            classMap.put(c.getSimpleName(), c);
        }

        // By default, use standard eval modules.
        if (!noStandardEvals) {
            for (Class<? extends StandardEval> stdClass : new PluginManager<StandardEval>(StandardEval.class).getPlugins()) {
                if (classMap.containsKey(stdClass.getSimpleName())) {
                    evals.add(classMap.get(stdClass.getSimpleName()));
                }
            }
        }

        // Get the specific classes provided.
        for (String module : modulesToUse) {
            if (!classMap.containsKey(module)) {
                throw new UserException.CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if (classMap.containsKey(module)) {
                evals.add(classMap.get(module));
            }
        }

        //add MetricsCollection if required modules are included

        if(evals.contains(classMap.get("CompOverlap")) && evals.contains(classMap.get("IndelSummary")) && evals.contains(classMap.get("TiTvVariantEvaluator")) && evals.contains(classMap.get("CountVariants")) && evals.contains(classMap.get("MultiallelicSummary")) )
            evals.add(classMap.get("MetricsCollection"));

        return evals;
    }

    /**
     * Subset a VariantContext to a single sample
     *
     * @param vc         the VariantContext object containing multiple samples
     * @param sampleName the sample to pull out of the VariantContext
     * @return a new VariantContext with just the requested sample
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, String sampleName) {
        return getSubsetOfVariantContext(vc, Collections.singleton(sampleName));
    }

    /**
     * Subset a VariantContext to a set of samples
     *
     * @param vc          the VariantContext object containing multiple samples
     * @param sampleNames the samples to pull out of the VariantContext
     * @return a new VariantContext with just the requested samples
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, Set<String> sampleNames) {
        // if we want to preserve AC0 sites as polymorphic we need to not rederive alleles
        final boolean deriveAlleles = variantEvalWalker.ignoreAC0Sites();
        return ensureAnnotations(vc, vc.subContextFromSamples(sampleNames, deriveAlleles));
    }

    public VariantContext ensureAnnotations(final VariantContext vc, final VariantContext vcsub) {
        final int originalAlleleCount = vc.getHetCount() + 2 * vc.getHomVarCount();
        final int newAlleleCount = vcsub.getHetCount() + 2 * vcsub.getHomVarCount();
        final boolean isSingleton = originalAlleleCount == newAlleleCount && newAlleleCount == 1;
        final boolean hasChrCountAnnotations = vcsub.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) &&
                vcsub.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) &&
                vcsub.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY);

        if ( ! isSingleton && hasChrCountAnnotations ) {
            // nothing to update
            return vcsub;
        } else {
            // have to do the work
            VariantContextBuilder builder = new VariantContextBuilder(vcsub);

            if ( isSingleton )
                builder.attribute(VariantEval.IS_SINGLETON_KEY, true);

            if ( ! hasChrCountAnnotations )
                VariantContextUtils.calculateChromosomeCounts(builder, true);

            return builder.make();
        }
    }

    /**
     * For a list of track names, bind the variant contexts to a trackName->sampleName->VariantContext mapping.
     * Additional variant contexts per sample are automatically generated and added to the map unless the sample name
     * matches the ALL_SAMPLE_NAME constant.
     *
     * @param tracker        the metadata tracker
     * @param ref            the reference context
     * @param tracks         the list of tracks to process
     * @param byFilter       if false, only accept PASSing VariantContexts.  Otherwise, accept both PASSing and filtered
     *                       sites
     * @param subsetBySample if false, do not separate the track into per-sample VCs
     * @param trackPerSample if false, don't stratify per sample (and don't cut up the VariantContext like we would need
     *                       to do this)
     * @return the mapping of track to VC list that should be populated
     */
    public HashMap<RodBinding<VariantContext>, HashMap<String, Collection<VariantContext>>>
    bindVariantContexts(RefMetaDataTracker tracker,
                        ReferenceContext ref,
                        List<RodBinding<VariantContext>> tracks,
                        boolean byFilter,
                        boolean subsetBySample,
                        boolean trackPerSample,
                        boolean trackPerFamily,
                        boolean mergeTracks) {
        if (tracker == null)
            return null;

        HashMap<RodBinding<VariantContext>, HashMap<String, Collection<VariantContext>>> bindings = new HashMap<RodBinding<VariantContext>, HashMap<String, Collection<VariantContext>>>();

        RodBinding<VariantContext> firstTrack = tracks.isEmpty() ? null : tracks.get(0);
        for (RodBinding<VariantContext> track : tracks) {
            HashMap<String, Collection<VariantContext>> mapping = new HashMap<String, Collection<VariantContext>>();

            for (VariantContext vc : tracker.getValues(track, ref.getLocus())) {

                // First, filter the VariantContext to represent only the samples for evaluation
                VariantContext vcsub = vc;

                if ((subsetBySample) && vc.hasGenotypes())
                    vcsub = getSubsetOfVariantContext(vc, variantEvalWalker.getSampleNamesForEvaluation());

                //always add a mapping for all samples together
                if ((byFilter || !vcsub.isFiltered())) {
                    addMapping(mapping, VariantEval.getAllSampleName(), vcsub);
                }

                // Now, if stratifying, split the subsetted vc per sample and add each as a new context
                if (vc.hasGenotypes() && trackPerSample) {
                    for (String sampleName : variantEvalWalker.getSampleNamesForEvaluation()) {
                        VariantContext samplevc = getSubsetOfVariantContext(vc, sampleName);

                        if (byFilter || !samplevc.isFiltered()) {
                            addMapping(mapping, sampleName, samplevc);
                        }
                    }
                }
                else if (vc.hasGenotypes() && trackPerFamily) {
                    for (final String familyName : variantEvalWalker.getFamilyNamesForEvaluation()) {
                        Set<String> familyMemberNames = new HashSet<>();
                        //if the current stratification family name is "all", then add all the families to the VC for evaluation here
                        if (familyName.equals(VariantEval.getAllFamilyName())) {
                            familyMemberNames = variantEvalWalker.getSampleNamesForEvaluation();
                        }
                        else {
                            Set<Sample> familyMembers = variantEvalWalker.getToolkit().getSampleDB().getFamily(familyName);
                            for (final Sample s : familyMembers) {
                                familyMemberNames.add(s.getID());
                            }
                        }
                        VariantContext samplevc = getSubsetOfVariantContext(vc, familyMemberNames);

                        if (byFilter || !samplevc.isFiltered()) {
                            addMapping(mapping, familyName, samplevc);
                        }
                    }
                }
            }

            if (mergeTracks && bindings.containsKey(firstTrack)) {
                // go through each binding of sample -> value and add all of the bindings from this entry
                HashMap<String, Collection<VariantContext>> firstMapping = bindings.get(firstTrack);
                for (Map.Entry<String, Collection<VariantContext>> elt : mapping.entrySet()) {
                    Collection<VariantContext> firstMappingSet = firstMapping.get(elt.getKey());
                    if (firstMappingSet != null) {
                        firstMappingSet.addAll(elt.getValue());
                    } else {
                        firstMapping.put(elt.getKey(), elt.getValue());
                    }
                }
            } else {
                bindings.put(track, mapping);
            }
        }

        return bindings;
    }

    private void addMapping(HashMap<String, Collection<VariantContext>> mappings, String sample, VariantContext vc) {
        if (!mappings.containsKey(sample))
            mappings.put(sample, new ArrayList<VariantContext>(1));
        mappings.get(sample).add(vc);
    }
}