package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotationType;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.PackageUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;

import java.util.*;


public class VariantAnnotatorEngine {

    private ArrayList<InfoFieldAnnotation> requestedInfoAnnotations;
    private ArrayList<GenotypeAnnotation> requestedGenotypeAnnotations;

    // should we annotate dbsnp?
    private boolean annotateDbsnp = false;
    // how about hapmap2?
    private boolean annotateHapmap2 = false;
    // how about hapmap3?
    private boolean annotateHapmap3 = false;


    // use this constructor if you want all possible annotations
    public VariantAnnotatorEngine(GenomeAnalysisEngine engine) {
        List<Class<? extends InfoFieldAnnotation>> infoAnnotationClasses = PackageUtils.getClassesImplementingInterface(InfoFieldAnnotation.class);
        requestedInfoAnnotations = getInstances(infoAnnotationClasses);
        List<Class<? extends GenotypeAnnotation>> genotypeAnnotationClasses = PackageUtils.getClassesImplementingInterface(GenotypeAnnotation.class);
        requestedGenotypeAnnotations = getInstances(genotypeAnnotationClasses);

        initialize(engine);
    }

    // use this constructor if you want to select specific annotations (and/or interfaces)
    public VariantAnnotatorEngine(GenomeAnalysisEngine engine, String[] annotationClassesToUse, String[] annotationsToUse) {
        // create a map for all annotation classes which implement our top-level interfaces
        HashMap<String, Class> classMap = new HashMap<String, Class>();
        for ( Class c : PackageUtils.getClassesImplementingInterface(InfoFieldAnnotation.class) )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : PackageUtils.getClassesImplementingInterface(GenotypeAnnotation.class) )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : PackageUtils.getInterfacesExtendingInterface(AnnotationType.class) )
            classMap.put(c.getSimpleName(), c);

        HashSet<Class> classes = new HashSet<Class>();
        // get the classes from the provided groups (interfaces)
        for ( String group : annotationClassesToUse ) {
            Class interfaceClass = classMap.get(group);
            if ( interfaceClass == null )
                interfaceClass = classMap.get(group + "Annotation");
            if ( interfaceClass == null )
                throw new StingException("Class " + group + " is not found; please check that you have specified the class name correctly");
            classes.addAll(PackageUtils.getClassesImplementingInterface(interfaceClass));
        }
        // get the specific classes provided
        for ( String annotation : annotationsToUse ) {
            Class annotationClass = classMap.get(annotation);
            if ( annotationClass == null )
                annotationClass = classMap.get(annotation + "Annotation");
            if ( annotationClass == null )
                throw new StingException("Class " + annotation + " is not found; please check that you have specified the class name correctly");
            classes.add(annotationClass);
        }

        // get the instances
        requestedInfoAnnotations = new ArrayList<InfoFieldAnnotation>();
        requestedGenotypeAnnotations = new ArrayList<GenotypeAnnotation>();

        for ( Class c : classes ) {
            // note that technically an annotation can work on both the INFO and FORMAT fields
            if ( InfoFieldAnnotation.class.isAssignableFrom(c) )
                requestedInfoAnnotations.add((InfoFieldAnnotation)getInstance(c));
            if ( GenotypeAnnotation.class.isAssignableFrom(c) )
                requestedGenotypeAnnotations.add((GenotypeAnnotation)getInstance(c));
        }

        initialize(engine);
    }

    private static <T> ArrayList<T> getInstances(List<Class<? extends T>> classes) {
        ArrayList<T> objects = new ArrayList<T>();
        for ( Class c : classes )
            objects.add((T)getInstance(c));
        return objects;
    }

    private static <T> T getInstance(Class<T> c) {
        try {
            return c.newInstance();
        } catch (InstantiationException e) {
            throw new StingException(String.format("Cannot instantiate annotation class '%s': must be concrete class", c.getSimpleName()));
        } catch (IllegalAccessException e) {
            throw new StingException(String.format("Cannot instantiate annotation class '%s': must have no-arg constructor", c.getSimpleName()));
        }
    }

    private void initialize(GenomeAnalysisEngine engine) {

        // check to see whether a dbsnp rod was included
        List<ReferenceOrderedDataSource> dataSources = engine.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(rodDbSNP.class) ) {
                annotateDbsnp = true;
            }
            if ( rod.getName().equals("hapmap2") ) {
                annotateHapmap2 = true;
            }
            if ( rod.getName().equals("hapmap3") ) {
                annotateHapmap3 = true;
            }
        }
    }

    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {

        Set<VCFHeaderLine> descriptions = new HashSet<VCFHeaderLine>();

        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations )
            descriptions.add(annotation.getDescription());
        for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations )
            descriptions.add(annotation.getDescription());
        if ( annotateDbsnp )
            descriptions.add(new VCFInfoHeaderLine(VCFRecord.DBSNP_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "dbSNP membership"));
        if ( annotateHapmap2 )
            descriptions.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP2_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Hapmap 2 membership"));
        if ( annotateHapmap3 )
            descriptions.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP3_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Hapmap 3 membership"));

        return descriptions;
    }

    public VariantContext annotateContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        Map<String, Object> infoAnnotations = new HashMap<String, Object>(vc.getAttributes());

        // annotate dbsnp occurrence
        if ( annotateDbsnp ) {
            rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getReferenceMetaData("dbsnp"));
            infoAnnotations.put(VCFRecord.DBSNP_KEY, dbsnp == null ? "0" : "1");
            // annotate dbsnp id if available and not already there
            if ( dbsnp != null && !vc.hasAttribute("ID") )
                infoAnnotations.put("ID", dbsnp.getRS_ID());
        }

        if ( annotateHapmap2 ) {
            List<Object> hapmap2 = tracker.getReferenceMetaData("hapmap2");
            infoAnnotations.put(VCFRecord.HAPMAP2_KEY, hapmap2.size() == 0 ? "0" : "1");
        }

        if ( annotateHapmap3 ) {
            List<Object> hapmap3 = tracker.getReferenceMetaData("hapmap3");
            infoAnnotations.put(VCFRecord.HAPMAP3_KEY, hapmap3.size() == 0 ? "0" : "1");
        }

        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations ) {
            Map<String, Object> result = annotation.annotate(tracker, ref, stratifiedContexts, vc);
            if ( result != null )
                infoAnnotations.putAll(result);
        }

        Map<String, Genotype> genotypes;
        if ( requestedGenotypeAnnotations.size() == 0 ) {
            genotypes = vc.getGenotypes();
        } else {
            genotypes = new HashMap<String, Genotype>(vc.getNSamples());
            for ( Map.Entry<String, Genotype> g : vc.getGenotypes().entrySet() ) {
                Genotype genotype = g.getValue();
                StratifiedAlignmentContext context = stratifiedContexts.get(g.getKey());
                if ( context == null ) {
                    genotypes.put(g.getKey(), genotype);
                    continue;
                }

                Map<String, Object> genotypeAnnotations = new HashMap<String, Object>(genotype.getAttributes());
                for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
                    Map<String, Object> result = annotation.annotate(tracker, ref, context, vc, genotype);
                    if ( result != null )
                        genotypeAnnotations.putAll(result);
                }
                genotypes.put(g.getKey(), new Genotype(g.getKey(), genotype.getAlleles(), genotype.getNegLog10PError(), genotype.getFilters(), genotypeAnnotations, genotype.genotypesArePhased()));
            }
        }

        return new VariantContext(vc.getName(), vc.getLocation(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.getFilters(), infoAnnotations);
    }
}