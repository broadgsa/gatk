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
import org.broadinstitute.sting.playground.gatk.walkers.annotator.GenomicAnnotation;
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

    // command-line option used for GenomicAnnotation.
    private Map<String, Set<String>> requestedColumnsMap;

    // command-line option used for GenomicAnnotation.
    private boolean oneToMany;


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
            descriptions.add(new VCFInfoHeaderLine(VCFRecord.DBSNP_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "dbSNP Membership"));
        if ( annotateHapmap2 )
            descriptions.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP2_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Hapmap 2 Membership"));
        if ( annotateHapmap3 )
            descriptions.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP3_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Hapmap 3 Membership"));

        return descriptions;
    }

    public Collection<VariantContext> annotateContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

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


        //Process the info field
        List<Map<String, Object>> infoAnnotationOutputsList = new LinkedList<Map<String, Object>>(); //each element in infoAnnotationOutputs corresponds to a single line in the output VCF file
        infoAnnotationOutputsList.add(new HashMap<String, Object>(vc.getAttributes())); //keep the existing info-field annotations. After this infoAnnotationOutputsList.size() == 1, which means the output VCF file gains 1 line.

        //go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations )
        {
            Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, ref, stratifiedContexts, vc);
            if ( annotationsFromCurrentType == null ) {
                continue;
            }

            if(annotationType instanceof GenomicAnnotation)
            {
                //go through the annotations returned by GenericAnnotation for each -B input file.
                for( Map.Entry<String, Object> annotationsFromInputFile : annotationsFromCurrentType.entrySet() )
                {
                    final String inputFileBindingName = annotationsFromInputFile.getKey();
                    final List<Map<String, String>> matchingRecords = (List<Map<String, String>>) annotationsFromInputFile.getValue();

                    if( matchingRecords.size() > 1 && oneToMany)
                    {
                        //More than one record matched in this file. After this, infoAnnotationOutputsList.size() will be infoAnnotationOutputsList.size()*matchingRecords.size().
                        infoAnnotationOutputsList = explodeInfoAnnotationOutputsList( infoAnnotationOutputsList, matchingRecords, inputFileBindingName);
                    }
                    else
                    {
                        //This doesn't change infoAnnotationOutputsList.size(). If more than one record matched, their annotations will
                        //all be added to the same output line, with keys disambiguated by appending _i .
                        addToExistingAnnotationOutputs( infoAnnotationOutputsList, matchingRecords, inputFileBindingName);
                    }
                }
            }
            else
            {
                //add the annotations to each output line.
                for(Map<String, Object> infoAnnotationOutput : infoAnnotationOutputsList) {
                    infoAnnotationOutput.putAll(annotationsFromCurrentType);
                }
            }
        }

        //Process genotypes
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

      //Create a separate VariantContext (aka. output line) for each element in infoAnnotationOutputsList
        Collection<VariantContext> returnValue = new LinkedList<VariantContext>();
        for(Map<String, Object> infoAnnotationOutput : infoAnnotationOutputsList) {
            returnValue.add( new VariantContext(vc.getName(), vc.getLocation(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.getFilters(), infoAnnotationOutput) );
        }

        return returnValue;
    }


    /**
     * Implements non-explode mode, where the output lines have a one-to-one relationship
     * with the input variants, and all multiple-match records are collapsed into the single info field.
     * The collapsing is done by appending an _i to each key name (where 'i' is a record counter).
     *
     * @param infoAnnotationOutputsList
     * @param matchingRecords
     * @param bindingName
     */
    private void addToExistingAnnotationOutputs(
            final List<Map<String, Object>> infoAnnotationOutputsList,
            final List<Map<String, String>> matchingRecords,
            final String bindingName) {
        //For each matching record, just add its annotations to all existing output lines.
        final boolean renameKeys = matchingRecords.size() > 1;
        for(int i = 0; i < matchingRecords.size(); i++) {
            Map<String,String> annotationsForRecord = matchingRecords.get(i);
            annotationsForRecord = selectColumnsFromRecord(bindingName, annotationsForRecord); //use only those columns that the user specifically requested.

            if(renameKeys) {
                //Rename keys to avoid naming conflicts (eg. if you have multiple dbsnp matches,
                // dbSNP.avHet=value1 from record 1 and dbSNP.avHet=value2 from record 2 will become dbSNP.avHet_1=value1 and dbSNP.avHet_2=value2 )
                Map<String,String> annotationsForRecordWithRenamedKeys = new HashMap<String, String>();
                for(Map.Entry<String, String> annotation : annotationsForRecord.entrySet()) {
                    annotationsForRecordWithRenamedKeys.put(annotation.getKey() + "_" + i, annotation.getValue());
                }

                annotationsForRecord = annotationsForRecordWithRenamedKeys;
            }

            //Add the annotations from this record to each output line.
            for(Map<String, Object> infoAnnotationOutput : infoAnnotationOutputsList) {
                infoAnnotationOutput.putAll(annotationsForRecord);
            }
        }
    }

    /**
     * Implements "explode" mode. Takes the current list of
     * infoAnnotationOutputs (each element of will end up in a different line
     * of the output VCF file), and generates/returns a new list of infoAnnotationOutputs
     * which contain one copy of the current infoAnnotationOutputs for each record
     * in matching records. The returned list will have size:
     *
     * infoAnnotationOutputsList.size() * matchingRecords.size()
     *
     * See class-level comments for more details.
     *
     * @param infoAnnotationOutputsList
     * @param matchingRecords
     * @param bindingName
     * @return
     */
    private List<Map<String, Object>> explodeInfoAnnotationOutputsList(
            final List<Map<String, Object>> infoAnnotationOutputsList,
            final List<Map<String, String>> matchingRecords,
            final String bindingName) {


        //This is the return value. It represents the new list of lines in the output VCF file.
        final List<Map<String, Object>> newInfoAnnotationOutputsList = new LinkedList<Map<String, Object>>();

        //For each matching record, generate a new output line
        for(int i = 0; i < matchingRecords.size(); i++) {
            Map<String,String> annotationsForRecord = matchingRecords.get(i);
            annotationsForRecord = selectColumnsFromRecord(bindingName, annotationsForRecord); //use only those columns that the user specifically requested.

            //Add the annotations from this record to each output line.
            for(Map<String, Object> infoAnnotationOutput : infoAnnotationOutputsList) {
                Map<String, Object> infoAnnotationOutputCopy = new HashMap<String, Object>(infoAnnotationOutput); //create a new copy of this line.
                infoAnnotationOutputCopy.putAll(annotationsForRecord); //Adds the column-value pairs from this record to this line.

                newInfoAnnotationOutputsList.add(infoAnnotationOutputCopy); //Add the line to the new list of lines.
            }
        }

        return newInfoAnnotationOutputsList;
    }



    /**
     * Takes a list of key-value pairs and returns a new Map containing only the columns which were requested by the user
     * via the -s arg. If there was no -s arg that referenced the given bindingName, all annotationsForRecord returned untouched.
     *
     * @param bindingName The binding name for a particular ROD input file.
     * @param annotationsForRecord The list of column_name -> value pairs for a particular record from the given input file.
     *
     * @return Map - see above.
     */
    private Map<String, String> selectColumnsFromRecord( String bindingName, Map<String, String> annotationsForRecord) {
        if(requestedColumnsMap == null || !requestedColumnsMap.containsKey(bindingName)) {
            return annotationsForRecord;
        }

        Set<String> requestedColumns = requestedColumnsMap.get(bindingName);
        Map<String, String> subsettedAnnotations = new HashMap<String, String>();
        for(Map.Entry<String, String> e : annotationsForRecord.entrySet() ) {
            if(requestedColumns.contains(e.getKey())) {
                subsettedAnnotations.put(e.getKey(), e.getValue());
            }
        }

        if(subsettedAnnotations.isEmpty()) {
            throw new StingException("Invalid -s argument for the '" + bindingName + "' input file. " +
                    "It caused all columns in the file to be rejected. Please check to make sure the -s column " +
                    "names match the column names in the '" + bindingName + "' file's HEADER line.");
        }

        return subsettedAnnotations;
    }



    /**
     * Determines how the engine will handle the case where multiple records in a ROD file
     * overlap a particular single locus. If explode is set to true, the output will be
     * one-to-many, so that each locus in the input VCF file could result in multiple
     * entries in the output VCF file. Otherwise, the output will be one-to-one, and
     * all multiple-match records will be collapsed into the single info field.
     * The collapsing is done by appending an _i to each key name (where 'i' is a
     * record counter).
     *
     * See class-level comments for more details.
     *
     * @param oneToMany
     */
    public void setOneToMany(boolean oneToMany) {
        this.oneToMany = oneToMany;
    }

    /**
     * Sets the columns that will be used for the info annotation field.
     * Column names should be of the form bindingName.columnName (eg. dbsnp.avHet).
     *
     * @param columns An array of strings where each string is a comma-separated list
     * of columnNames (eg ["dbsnp.avHet,dbsnp.valid", "file2.col1,file3.col1"] ).
     */
    public void setRequestedColumns(String[] columns) {
        if(columns == null) {
            throw new IllegalArgumentException("columns arg is null. Please check the -s command-line arg.");
        }

        //System.err.println("COLUMNS:  "+Arrays.asList(columns).toString());

        this.requestedColumnsMap = GenomicAnnotation.parseColumnsArg(columns);
    }
}
