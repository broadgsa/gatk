/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.annotator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotationType;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.playground.gatk.walkers.annotator.GenomicAnnotation;
import org.broadinstitute.sting.playground.gatk.walkers.annotator.JoinTable;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.classloader.PackageUtils;


public class VariantAnnotatorEngine {

    public static final String dbPrefix = "comp";

    private List<InfoFieldAnnotation> requestedInfoAnnotations;
    private List<GenotypeAnnotation> requestedGenotypeAnnotations;

    private HashMap<String, String> dbAnnotations = new HashMap<String, String>();

        // command-line option from GenomicAnnotator.
    private Map<String, Set<String>> requestedColumnsMap;

    // command-line option from GenomicAnnotator.
    private boolean oneToMany;

    // command-line option from GenomicAnnotator.
    private List<JoinTable> joinTables;

    // used by GenomicAnnotator. Maps binding name to number of output VCF records
    // annotated with records from the input table with this binding name. Only used for
    // printing out stats at the end.
    private Map<String, Integer> inputTableHitCounter = new HashMap<String, Integer>();


    // use this constructor if you want all possible annotations
    public VariantAnnotatorEngine(GenomeAnalysisEngine engine) {
        requestedInfoAnnotations = PackageUtils.getInstancesOfClassesImplementingInterface(InfoFieldAnnotation.class);
        requestedGenotypeAnnotations = PackageUtils.getInstancesOfClassesImplementingInterface(GenotypeAnnotation.class);
        initialize(engine);
    }

    // use this constructor if you want to select specific annotations (and/or interfaces)
    public VariantAnnotatorEngine(GenomeAnalysisEngine engine, String[] annotationGroupsToUse, String[] annotationsToUse) {

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
        for ( String group : annotationGroupsToUse ) {
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
                requestedInfoAnnotations.add((InfoFieldAnnotation)PackageUtils.getSimpleInstance(c));
            if ( GenotypeAnnotation.class.isAssignableFrom(c) )
                requestedGenotypeAnnotations.add((GenotypeAnnotation)PackageUtils.getSimpleInstance(c));
        }

        initialize(engine);
    }

    private void initialize(GenomeAnalysisEngine engine) {

        // check to see whether comp rods were included
        List<ReferenceOrderedDataSource> dataSources = engine.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            if ( source.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                dbAnnotations.put(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, VCFConstants.DBSNP_KEY);
            }
            else if ( source.getName().startsWith(dbPrefix) ) {
                dbAnnotations.put(source.getName(), source.getName().substring(dbPrefix.length()));
            }
        }
    }

    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {

        Set<VCFHeaderLine> descriptions = new HashSet<VCFHeaderLine>();

        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( Map.Entry<String, String> dbSet : dbAnnotations.entrySet() )
            descriptions.add(new VCFInfoHeaderLine(dbSet.getValue(), 0, VCFHeaderLineType.Flag, (dbSet.getKey().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ? "dbSNP" : dbSet.getValue()) + " Membership"));

        return descriptions;
    }

    // A slightly simplified interface for when you don't have any reads, so the stratifiedContexts aren't necessary, and
    // you only permit a single return value
    public VariantContext annotateContext(RefMetaDataTracker tracker, ReferenceContext ref, VariantContext vc) {
        Collection<VariantContext> results = this.annotateContext(tracker, ref, EMPTY_STRATIFIED_ALIGNMENT_CONTEXT, vc);

        if ( results.size() != 1 )
            throw new StingException("BUG: annotateContext call requires 1 resulting annotated VC, but got " + results);

        return results.iterator().next();

    }
    private static final Map<String, StratifiedAlignmentContext> EMPTY_STRATIFIED_ALIGNMENT_CONTEXT = (Map<String, StratifiedAlignmentContext>)Collections.EMPTY_MAP;

    public Collection<VariantContext> annotateContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        Map<String, Object> infoAnnotations = new LinkedHashMap<String, Object>(vc.getAttributes());

        // annotate db occurrences
        for ( Map.Entry<String, String> dbSet : dbAnnotations.entrySet() ) {
            if ( dbSet.getKey().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                DbSNPFeature dbsnp = DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
                infoAnnotations.put(VCFConstants.DBSNP_KEY, dbsnp == null ? false : true);
                // annotate dbsnp id if available and not already there
                if ( dbsnp != null && (!vc.hasAttribute(VariantContext.ID_KEY) || vc.getAttribute(VariantContext.ID_KEY).equals(VCFConstants.EMPTY_ID_FIELD)) )
                    infoAnnotations.put(VariantContext.ID_KEY, dbsnp.getRsID());
            } else {
                boolean overlapsComp = false;
                for ( VariantContext comp : tracker.getVariantContexts(ref, dbSet.getKey(), null, ref.getLocus(), false, false) ) {
                    if ( !comp.isFiltered() ) {
                        overlapsComp = true;
                        break;
                    }
                }
                infoAnnotations.put(dbSet.getValue(), overlapsComp ? true : false);
            }
        }

        //Process the info field
        List<Map<String, Object>> infoAnnotationOutputsList = new LinkedList<Map<String, Object>>(); //each element in infoAnnotationOutputs corresponds to a single line in the output VCF file
        infoAnnotationOutputsList.add(new LinkedHashMap<String, Object>(vc.getAttributes())); //keep the existing info-field annotations. After this infoAnnotationOutputsList.size() == 1, which means the output VCF file has 1 additional line.
        infoAnnotationOutputsList.get(0).putAll(infoAnnotations);  // put the DB membership info in

        //go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations )
        {
            Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, ref, stratifiedContexts, vc);
            if ( annotationsFromCurrentType == null ) {
                continue;
            }

            if(annotationType instanceof GenomicAnnotation)
            {
                infoAnnotationOutputsList = processGenomicAnnotation( infoAnnotationOutputsList, annotationsFromCurrentType );
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
            returnValue.add( new VariantContext(vc.getName(), vc.getLocation(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, infoAnnotationOutput) );
        }

        return returnValue;
    }

    // Finish processing data from GenomicAnnotation.
    private List<Map<String, Object>> processGenomicAnnotation( List<Map<String, Object>> infoAnnotationOutputsList, Map<String, Object> annotationsForCurrentLocusFromAllAnnotatorInputTables)
    {

        //process the map returned by GenomicAnnotation. This completes processing of the -B args.
        for( Map.Entry<String, Object> annotationsFromOneInputTable : annotationsForCurrentLocusFromAllAnnotatorInputTables.entrySet() )
        {
            final String inputTableBindingName = annotationsFromOneInputTable.getKey();
            final List<Map<String, String>> matchingRecords = (List<Map<String, String>>) annotationsFromOneInputTable.getValue();

            if( matchingRecords.size() > 1 && oneToMany)
            {
                //More than one record matched in this file. After this, infoAnnotationOutputsList.size() will be infoAnnotationOutputsList.size()*matchingRecords.size().
                infoAnnotationOutputsList = explodeInfoAnnotationOutputsList( infoAnnotationOutputsList, matchingRecords, inputTableBindingName );
            }
            else
            {
                //This doesn't change infoAnnotationOutputsList.size(). If more than one record matched, their annotations will
                //all be added to the same output line, with keys disambiguated by appending _i .
                addToExistingAnnotationOutputs( infoAnnotationOutputsList, matchingRecords, inputTableBindingName );
            }
        }

        //process -J args
        if(joinTables != null)
        {
            //for each joinTable, join it with the data in the info-field of each output line.
            for(JoinTable joinTable : joinTables)
            {
                //for each info field, join it to the current join table
                final List<Map<String, Object>> previousInfoAnnotationOutputsList = new LinkedList<Map<String, Object>>(infoAnnotationOutputsList);  //create a shallow copy because infoAnnotationOutputsList will change during the iteration.
                for(Map<String, Object> outputRecordInfoField : previousInfoAnnotationOutputsList)
                {
                    infoAnnotationOutputsList = performJoin( infoAnnotationOutputsList, outputRecordInfoField, joinTable );
                }
            }
        }

        //apply -S args last to select the columns requested by the user
        if(requestedColumnsMap != null) {
            infoAnnotationOutputsList = applySelectArg(infoAnnotationOutputsList);
        }

        return infoAnnotationOutputsList;
    }

   // Performs a join between the an info field record represented by outputRecordInfoField and the infoAnnotationOutputsList.
   private List<Map<String, Object>> performJoin( List<Map<String, Object>> infoAnnotationOutputsList,  Map<String, Object> outputRecordInfoField, JoinTable joinTable)
    {
        //System.err.println("Looking at: " + joinTable.getLocalBindingName()+ "- join to " + joinTable.getExternalBindingName() + "." + joinTable.getExternalColumnName() );
        //for the current joinTable, for each output line, find the externalJoinColumnValue and see if it matches the joinColumnValue of any record(s) in this joinTable.
        final String externalBindingName = joinTable.getExternalBindingName();
        final String externalColumnName = joinTable.getExternalColumnName();
        final String fullyQualifiedExternalColumnName = GenomicAnnotation.generateInfoFieldKey(externalBindingName, externalColumnName);

        //find the externalJoinColumnValue in the current info field, and then look up any joinTable records that have this value for the localJoinColumnValue
        List<ArrayList<String>> matchingJoinTableRecords = null; //record(s) in the join table whose joinColumnValue(s) matches the joinColumnValue inside the current outputRecordInfoField. Since the join keys don't have to be unique, there may be more than one record in the join table thtat matches.
        final Object numInfoFieldKeysToCheckObj = outputRecordInfoField.get(GenomicAnnotation.generateInfoFieldKey(externalBindingName, GenomicAnnotation.NUM_MATCHES_SPECIAL_INFO_FIELD));
        if(numInfoFieldKeysToCheckObj == null) {
            //only 1 record in the externalBindingName -B AnnotatoInfoTable overlapped the current position
            Object externalColumnValue = outputRecordInfoField.get(fullyQualifiedExternalColumnName);
            if(externalColumnValue != null) {
                matchingJoinTableRecords = joinTable.get(externalColumnValue.toString());
                //System.err.println("Found matching record(s) in join table for record: " + outputRecordInfoField + " where " + fullyQualifiedExternalColumnName  + "==" + externalColumnValue +  ":  " + matchingJoinTableRecords);
            }
        } else {
            //multiple records in the externalBindingName -B AnnotatoInfoTable overlapped the current position
            final int numInfoFieldKeysToCheck = Integer.parseInt(numInfoFieldKeysToCheckObj.toString());
            for(int i = 0; i < numInfoFieldKeysToCheck; i++) {
                final Object externalColumnValue = outputRecordInfoField.get(fullyQualifiedExternalColumnName + "_" + i);
                if(externalColumnValue != null) {
                    matchingJoinTableRecords = joinTable.get(externalColumnValue.toString());
                    if(matchingJoinTableRecords != null) {
                        //System.err.println("Found matching record(s) in join table for record: " + outputRecordInfoField + " where " + fullyQualifiedExternalColumnName  + "==" + externalColumnValue +  ":  " + matchingJoinTableRecords);
                        break;
                    }
                }
            }
        }

        //if match(s) for the externalJoinColumnValue in the current outputRecordInfoField have been found in the join table, perform the join.
        if(matchingJoinTableRecords != null)
        {
            final String joinTableBindingName = joinTable.getLocalBindingName();

            //convert the List<ArrayList<String>> to List<Map<String, String>> by hashing the values from the ArrayList<String> by their column names.
            final List<Map<String, String>> matchingJoinTableRecordsConverted = new LinkedList<Map<String,String>>();
            for(ArrayList<String> columnValues : matchingJoinTableRecords) {
                final List<String> columnNames = joinTable.getColumnNames();

                final Map<String, String> matchingRecord = new LinkedHashMap<String, String>();
                for(int i = 0; i < columnNames.size(); i++) {
                    matchingRecord.put(columnNames.get(i), columnValues.get(i));
                }

                matchingJoinTableRecordsConverted.add(GenomicAnnotation.convertRecordToAnnotations(joinTableBindingName, matchingRecord));
            }



            // do the join between the outputRecordInfoField and the matchingJoinTableRecords, then add the results to to infoAnnotationOutputsList
            List<Map<String, Object>> tempList = new LinkedList<Map<String, Object>>();
            tempList.add(outputRecordInfoField);
            if( matchingJoinTableRecordsConverted.size() > 1 && oneToMany)
            {
                //More than one record in the joinTable matched the current info field. After this, infoAnnotationOutputsList.size() will be infoAnnotationOutputsList.size()*matchingRecords.size().
                tempList = explodeInfoAnnotationOutputsList( tempList, matchingJoinTableRecordsConverted, joinTableBindingName );
            }
            else
            {
                //This doesn't change infoAnnotationOutputsList.size(). If more than one record matched, their annotations will
                //all be added to the same output line, with keys disambiguated by appending _i .
                addToExistingAnnotationOutputs( tempList, matchingJoinTableRecordsConverted, joinTableBindingName );
            }

            infoAnnotationOutputsList.remove(outputRecordInfoField); //remove the old info field
            infoAnnotationOutputsList.addAll(tempList); //add the new info field(s) that have been joined with the matchingJoinTableRecords
        }
        return infoAnnotationOutputsList;
    }


    // Implements not-oneToMany mode, where the output lines have a one-to-one relationship
    // with the input variants, and all multiple-match records are collapsed into the single info field.
    // The collapsing is done by appending an _i to each key name (where 'i' is a record counter), as well
    // as a special bindingName.numMatchingRecords=n key-value pair which specifies the upper limit of the counter.
    private void addToExistingAnnotationOutputs(
            final List<Map<String, Object>> infoAnnotationOutputsList,
            final List<Map<String, String>> matchingRecords,
            final String bindingName) {
        //For each matching record, just add its annotations to all existing output lines.
        final boolean renameKeys = matchingRecords.size() > 1;
        for(int i = 0; i < matchingRecords.size(); i++) {
            Map<String,String> currentRecord = matchingRecords.get(i);

            if(renameKeys) {
                //Rename keys to avoid naming conflicts. After this all keys from the i'th matching record will have _i appended to them.
                // (This solves the following problem: if you have multiple dbsnp matches - such as dbSNP.avHet=value1 from record 1 and
                //  dbSNP.avHet=value2 from record 2, the keys will be renamed to dbSNP.avHet_1=value1 and dbSNP.avHet_2=value2 )
                Map<String,String> currentRecordWithRenamedKeys = new LinkedHashMap<String, String>();
                for(final Map.Entry<String, String> annotation : currentRecord.entrySet()) {
                    currentRecordWithRenamedKeys.put(annotation.getKey() + "_" + (i + 1), annotation.getValue());
                }
                currentRecordWithRenamedKeys.put(GenomicAnnotation.generateInfoFieldKey(bindingName, GenomicAnnotation.NUM_MATCHES_SPECIAL_INFO_FIELD),
                        Integer.toString(matchingRecords.size())); //add the special field that specifies how many matchingRecords there were.
                currentRecord = currentRecordWithRenamedKeys;
            }

            //Add the annotations from this record to each output line.
            for(Map<String, Object> outputRecordInfoField : infoAnnotationOutputsList) {
                outputRecordInfoField.putAll(currentRecord);
            }
        }

        incrementStatsCounter(bindingName, infoAnnotationOutputsList.size());
    }

    /**
     * Records statistics that will be printed when GenomicAnnotator finishes.
     *
     * @param bindingName The table from which annotations were gotten
     * @param numNewRecords The number of new output VCF records created with annotations from this table
     */
    private void incrementStatsCounter( final String bindingName, int numNewRecords) {
        //record some stats - there were infoAnnotationOutputsList.size() output VCF records annotated with data from the 'bindingName' input table.
        Integer counter = inputTableHitCounter.get(bindingName);
        if( counter == null ) {
            inputTableHitCounter.put(bindingName, numNewRecords); //init the counter
        } else {
            inputTableHitCounter.put(bindingName, counter + numNewRecords); //increment the counter
        }
    }

    // Implements oneToMany mode. Takes the current infoAnnotationOutputsList
    // (where each element represents a line in the output VCF file), and
    // generates a new infoAnnotationOutputsList which contains one copy of the current
    // infoAnnotationOutputs for each record matchingRecords.
    // The returned list will have size:
    // infoAnnotationOutputsList.size() * matchingRecords.size()
    private List<Map<String, Object>> explodeInfoAnnotationOutputsList(
            final List<Map<String, Object>> infoAnnotationOutputsList,
            final List<Map<String, String>> matchingRecords,
            final String bindingName) {


        //This is the return value. It represents the new list of lines in the output VCF file.
        final List<Map<String, Object>> newInfoAnnotationOutputsList = new LinkedList<Map<String, Object>>();

        //For each matching record, generate a new output line
        for(int i = 0; i < matchingRecords.size(); i++) {
            Map<String,String> annotationsForRecord = matchingRecords.get(i);

            //Add the annotations from this record to each output line.
            for(Map<String, Object> outputRecordInfoField : infoAnnotationOutputsList) {
                Map<String, Object> outputRecordInfoFieldCopy = new LinkedHashMap<String, Object>(outputRecordInfoField); //create a new copy of this line.
                outputRecordInfoFieldCopy.putAll(annotationsForRecord); //Adds the column-value pairs from this record to this line.

                newInfoAnnotationOutputsList.add(outputRecordInfoFieldCopy); //Add the line to the new list of lines.
            }
        }

        recordStats(bindingName, newInfoAnnotationOutputsList.size(), infoAnnotationOutputsList, matchingRecords.size());

        return newInfoAnnotationOutputsList;
    }


    /**
     * Records statistics for the explodeInfoAnnotationOutputsList(..) calculation.
     * @param bindingName The table from which annotations were gotten
     * @param numNewVCFRecordsAnnotatedWithBindingNameData The number of new output VCF records created with annotations from this table
     * @param infoAnnotationOutputsList output list
     * @param matchingRecordsSize  matching records size
     */
    private void recordStats( final String bindingName, int numNewVCFRecordsAnnotatedWithBindingNameData, final List<Map<String, Object>> infoAnnotationOutputsList, int matchingRecordsSize ) {

        //update stats for the 'bindingName' table
        incrementStatsCounter(bindingName, numNewVCFRecordsAnnotatedWithBindingNameData); //All records in newInfoAnnotationOutputsList were annotated with data from bindingName.

        //update stats for all other tables besides 'bindingName'
        for(String otherBindingName : inputTableHitCounter.keySet()) {
            if(otherBindingName.equals(bindingName)) {
                continue;
            }

            //count how many records in the initial infoAnnotationOutputsList were annotated with data from otherBindingName
            int numAnnotatedWithOtherBindingNameData = 0;
            for(Map<String, Object> outputRecordInfoField : infoAnnotationOutputsList) {
                for(String outputRecordInfoFieldKey : outputRecordInfoField.keySet()) {
                    if(outputRecordInfoFieldKey.contains(otherBindingName)) {
                        //this record has some annotations from the otherBindingName table
                        numAnnotatedWithOtherBindingNameData++;
                        break;
                    }
                }
            }

            if(numAnnotatedWithOtherBindingNameData > 0) {
                //numAnnotatedWithOtherBindingNameData * (matchingRecordsSize - 1) is how many additional output VCF records were created with annotations from otherBindingName
                incrementStatsCounter(otherBindingName, numAnnotatedWithOtherBindingNameData * (matchingRecordsSize - 1));
            }
        }
    }


    // Applies the -S arg to the results
    private List<Map<String, Object>> applySelectArg( final List<Map<String, Object>> infoAnnotationOutputsList )
    {
        final List<Map<String, Object>> newInfoAnnotationOutputList = new LinkedList<Map<String, Object>>();
        for(final Map<String, Object> outputRecordInfoField : infoAnnotationOutputsList) {
            final Map<String, Object> newOutputRecordInfoField = new LinkedHashMap<String, Object>();
            for(final Entry<String, Object>  keyValue : outputRecordInfoField.entrySet()) {
                if(!isKeyFilteredOutBySelectArg(keyValue.getKey())) {
                    newOutputRecordInfoField.put(keyValue.getKey(), keyValue.getValue());
                }
            }
            newInfoAnnotationOutputList.add(newOutputRecordInfoField);
        }

        return newInfoAnnotationOutputList;
    }


    /**
     * Determines whether to exclude the given column from the annotations.
     * @param key The fully qualified columnName
     * @return Whether the -S arg specifies that this column should be included in the annotations.
     *
     * TODO this function can be optimized through memoization
     */
    private boolean isKeyFilteredOutBySelectArg(String key)
    {
        for(final String bindingName : requestedColumnsMap.keySet()) {

            if(key.contains(bindingName)) {
                final Set<String> selectArgsWithThisBindingName = requestedColumnsMap.get(bindingName);
                for(final String selectArgWithThisBindingName : selectArgsWithThisBindingName) {
                    if(key.contains(selectArgWithThisBindingName)) {
                        return false; //this key matches one of the -s args, so the user explicitly requested this key
                    }
                }
                if(!selectArgsWithThisBindingName.isEmpty()) {
                    return true; //the -S arg contains some keys with this binding name, but doesn't include this key
                }
            }
        }

        return false; //the -S arg doesn't have anything with the same binding name as this key, so the user implicitly requested this key
    }




    /**
     * Determines how the engine will handle the case where multiple records in a ROD file
     * overlap a particular single locus. If oneToMany is set to true, the output will be
     * one-to-many, so that each locus in the input VCF file could result in multiple
     * entries in the output VCF file. Otherwise, the output will be one-to-one, and
     * all multiple-match records will be collapsed into the single info field.
     * The collapsing is done by appending an _i to each key name (where 'i' is a
     * record counter).
     *
     * See class-level comments for more details.
     *
     * @param oneToMany true if we should break out from one to many
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

        this.requestedColumnsMap = parseColumnsArg(columns);
    }


    /**
     * Passes in a pointer to the JoinTables.
     *
     * @param joinTables The list of JoinTables. There should be one JoinTable object for each -J arg.
     */
    public void setJoinTables(List<JoinTable> joinTables) {
        this.joinTables = joinTables;
    }


    /**
     * Parses the columns arg and returns a Map of columns hashed by their binding name.
     * For example:
     *   The command line:
     *      -s dbSnp.valid,dbsnp.avHet -s refGene.txStart,refGene.txEnd
     *
     *   will be passed to this method as:
     *       ["dbSnp.valid,dbsnp.avHet", "refGene.txStart,refGene.txEnd"]
     *
     *   resulting in a return value of:
     *      {
     *       "dbSnp" -> "dbSnp.valid" ,
     *       "dbSnp" -> "dbsnp.avHet" ,
     *       "refGene" -> "refGene.txStart",
     *       "refGene" -> "refGene.txEnd"
     *      }
     *
     * @param columnsArg The -s command line arg value.
     *
     * @return Map representing a parsed version of this arg - see above.
     */
    private static Map<String, Set<String>> parseColumnsArg(String[] columnsArg) {
        Map<String, Set<String>> result = new HashMap<String, Set<String>>();

        for(String s : columnsArg) {
            for(String columnSpecifier : s.split(",") ) {
                String[] rodNameColumnName = columnSpecifier.split("\\.");
                if(rodNameColumnName.length != 2) {
                    throw new IllegalArgumentException("The following column specifier in the -s arg is invalid: [" + columnSpecifier + "]. It must be of the form 'bindingName.columnName'.");
                }
                String rodName = rodNameColumnName[0];
                //String columnName = rodNameColumnName[1];

                Set<String> requestedColumns = result.get(rodName);
                if(requestedColumns == null) {
                    requestedColumns = new HashSet<String>();
                    result.put(rodName, requestedColumns);
                }
                requestedColumns.add(columnSpecifier);
            }
        }

        return result;
    }


    //Returns a map containing stats on how many output vcf records were annotated from each database
    public Map<String, Integer> getInputTableHitCounter() {
        return Collections.unmodifiableMap(inputTableHitCounter);
    }


}
