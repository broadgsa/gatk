/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEval;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.StandardEval;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.RequiredStratification;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.StandardStratification;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.VariantContextUtils;

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
                    throw new StingException("Unable to instantiate stratification module '" + c.getSimpleName() + "'");
                } catch (IllegalAccessException e) {
                    throw new StingException("Illegal access error when trying to instantiate stratification module '" + c.getSimpleName() + "'");
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

                if (subsetBySample && vc.hasGenotypes())
                    vcsub = getSubsetOfVariantContext(vc, variantEvalWalker.getSampleNamesForEvaluation());

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