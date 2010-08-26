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


package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableCodec;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

/**
 * Annotates variant calls with information from user-specified tabular files.
 *
 * For details, see:  http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=VariantContext.class))
@By(DataSource.REFERENCE)
public class GenomicAnnotator extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="vcfOutput", shortName="vcf", doc="Please use --out instead", required=false)
    @Deprecated
    protected String oldOutArg;

    @Argument(fullName="sampleName", shortName="sample", doc="The sample (NA-ID) corresponding to the variant input (for non-VCF input only)", required=false)
    protected String sampleName = null;

    @Argument(fullName="select", shortName="s", doc="Optionally specifies which subset of columns from which -B inputs should be used for annotations. For example, -B:mydbsnp,AnnotatorInputTable /path/to/mydbsnp.txt -B:mytable,AnnotatorInputTable /path/mytable.txt -s mydbsnp.avHet,mydbsnp.name,mytable.column3 will cause annotations to only be generated from the 3 columns specified using -s.", required=false)
    protected String[] SELECT_COLUMNS = {};

    @Argument(fullName="join", shortName="J", doc="Optionally specifies a file and column within that file that should be LEFT-JOIN'ed to a column in a previously-specified file. The file provided to -J must be tab-delimited, with the first non-comment/non-empty line containing column names. (example: -B name,AnnotatorInputTable,/path/to/file1   -J name2,/path/to/file2,name.columnName=name2.columnName2  - this will join the table in file2 to the table in file1) ", required=false)
    protected String[] JOIN_ARGS = {};

    @Argument(fullName="oneToMany", shortName="m", doc="If more than one record from the same file matches a particular locus (for example, multiple dbSNP records with the same position), create multiple entries in the ouptut VCF file - one for each match. If a particular tabular file has J matches, and another tabular file has K matches for a given locus, then J*K output VCF records will be generated - one for each pair of K, J.   If this flag is not provided, the multiple records are still generated, but they are stored in the INFO field of a single output VCF record, with their annotation keys differentiated by appending '_i' with i varying from 1 to K*J. ", required=false)
    protected Boolean ONE_TO_MANY = false;

    private VariantAnnotatorEngine engine;

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        //read all ROD file headers and construct a set of all column names to be used for validation of command-line args
        final Set<String> allFullyQualifiedColumnNames = new LinkedHashSet<String>();
        final Set<String> allBindingNames = new LinkedHashSet<String>();
        try {
            for(ReferenceOrderedDataSource ds : getToolkit().getRodDataSources()) {
                if(! ds.getReferenceOrderedData().getType().equals(AnnotatorInputTableCodec.class)) {
                    continue; //skip all non-AnnotatorInputTable files.
                }
                final String bindingName = ds.getName();
                allBindingNames.add(bindingName);
                final ArrayList<String> header = AnnotatorInputTableCodec.readHeader(ds.getReferenceOrderedData().getFile());
                for(String columnName : header) {
                    allFullyQualifiedColumnNames.add(bindingName + "." + columnName);
                }
            }
        } catch(IOException e) {
            throw new StingException("Failed when attempting to read file header. ", e);
        }

        //parse the JOIN_COLUMNS args, read in the specified files, and validate column names in the = relation. This end result of this loop is to populate the List of joinTables with one entry per -J arg.
        final List<JoinTable> joinTables = new LinkedList<JoinTable>();
        for(String joinArg : JOIN_ARGS) {

            //parse the tokens
            final String[] arg = joinArg.split(",");
            if(arg.length != 3) {
                throw new StingException("The following -J arg: \"" + joinArg + "\" must contain 3 comma-separated values. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }
            final String bindingName = arg[0];
            final String filename = arg[1];
            final String columnsToJoin = arg[2];

            if(allBindingNames.contains(bindingName)) {
                throw new StingException("The name \"" + bindingName + "\" in the -J arg: \"" + joinArg + "\" has already been used in another binding.");
            }

            String[] splitOnEquals = columnsToJoin.split("=+");
            if(splitOnEquals.length != 2) {
                throw new StingException("The -J arg: \"" + joinArg + "\" must specify the columns to join on. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }

            String[] splitOnDot1 = splitOnEquals[0].split("\\.");
            String[] splitOnDot2 = splitOnEquals[1].split("\\.");
            if(splitOnDot1.length != 2 || splitOnDot2.length != 2) {
                throw new StingException("The -J arg: \"" + joinArg + "\" must fully specify the columns to join on. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }

            final String bindingName1 = splitOnDot1[0];
            final String columnName1 = splitOnDot1[1];
            final String bindingName2 = splitOnDot2[0];
            final String columnName2 = splitOnDot2[1];

            //figure out which of the 2 binding names within the = relation matches the -J bindingName
            final String localBindingName = bindingName; //alias
            final String localColumnName;
            final String externalBindingName;
            final String externalColumnName;
            if(bindingName1.equals(bindingName)) {
                localColumnName = columnName1;
                externalBindingName = bindingName2;
                externalColumnName = columnName2;
            } else if(bindingName2.equals(bindingName)) {
                localColumnName = columnName2;
                externalBindingName = bindingName1;
                externalColumnName = columnName1;
            } else {
                throw new StingException("The name \"" + bindingName + "\" in the -J arg: \"" + joinArg + "\" must be specified in one the columns to join on. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }

            //validate externalColumnName
            final String fullyQualifiedExternalColumnName = externalBindingName + '.' + externalColumnName;
            if( !allFullyQualifiedColumnNames.contains(fullyQualifiedExternalColumnName) ) {
                throw new StingException("The -J arg: \"" + joinArg + "\" specifies an unknown column name: \"" + fullyQualifiedExternalColumnName + "\"");
            }

            //read in the file contents into a JoinTable object
            final JoinTable joinTable = new JoinTable();
            joinTable.parseFromFile(filename, localBindingName, localColumnName, externalBindingName, externalColumnName);
            joinTables.add(joinTable);

            //validate localColumnName, and add all column names in this file to the list of allFullyQualifiedColumnNames so that they can be referenced from subsequent -J args.
            final List<String> columnNames = joinTable.getColumnNames();
            final List<String> fullyQualifiedColumnNames = new LinkedList<String>();
            boolean found = false;
            for ( String columnName : columnNames ) {
                 if ( columnName.equals(localColumnName) )
                     found = true;
                 fullyQualifiedColumnNames.add(localBindingName + '.' + columnName);
            }
            if ( !found )
                throw new StingException("The -J arg: \"" + joinArg + "\" specifies an unknown column name: \"" + localColumnName + "\". It's not one of the column names in the header " + columnNames + " of the file: " + filename);

            allFullyQualifiedColumnNames.addAll(fullyQualifiedColumnNames);
        }

        //parse the SELECT_COLUMNS arg and validate the column names
        List<String> parsedSelectColumns = new LinkedList<String>();
        for ( String token : SELECT_COLUMNS )
            parsedSelectColumns.addAll(Arrays.asList(token.split(",")));
        SELECT_COLUMNS = parsedSelectColumns.toArray(SELECT_COLUMNS);

        for ( String columnName : SELECT_COLUMNS ) {
            if ( !allFullyQualifiedColumnNames.contains(columnName) )
                throw new StingException("The column name '" + columnName + "' provided to -s doesn't match any of the column names in any of the -B files. Here is the list of available column names: " + allFullyQualifiedColumnNames);
        }

        //instantiate the VariantAnnotatorEngine
        ArrayList<String> annotationsToUse = new ArrayList<String>();
        annotationsToUse.add("GenomicAnnotation");
        engine = new VariantAnnotatorEngine(getToolkit(), new ArrayList<String>(), annotationsToUse);
        engine.setOneToMany(ONE_TO_MANY);
        engine.setRequestedColumns(SELECT_COLUMNS);
        engine.setJoinTables(joinTables);

        // set up the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), Arrays.asList("variant")));
        hInfo.addAll(engine.getVCFAnnotationDescriptions());

        Set<String> rodName = new HashSet<String>();
        rodName.add("variant");
        TreeSet<String> samples = new TreeSet<String>(SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName));
        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     * We want reads that span deletions
     *
     * @return true
     */
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * For each site of interest, annotate based on the requested annotation types
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Set<VariantContext> results = new LinkedHashSet<VariantContext>();
        for (VariantContext vc : tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, false)) {
            if ( vc.isFiltered() ||
                    (vc.isVariant() && (!vc.isSNP() || !vc.isBiallelic())) ) {
                results.add(vc);
            } else {
                Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
                if ( stratifiedContexts != null )
                    results.addAll(engine.annotateContext(tracker, ref, stratifiedContexts, vc));
                else
                    results.add(vc);
            }
        }

        for ( VariantContext vc : results )
            vcfWriter.add(vc ,ref.getBase());

        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer sum) {

        //out.printf("Generated %d annotated VCF records.\n", totalOutputVCFRecords);
        Map<String, Integer> inputTableHitCounter = engine.getInputTableHitCounter();
        for ( Entry<String, Integer> e : inputTableHitCounter.entrySet() ) {
            final String bindingName = e.getKey();
            final int counter = e.getValue();
            //final float percent = 100 * counter /(float) totalOutputVCFRecords;
            //out.printf(" %-6.1f%%   (%d) annotated with %s.\n", percent, counter, bindingName );
            System.out.printf(" %d annotated with %s.\n", counter, bindingName );
        }
    }
}

