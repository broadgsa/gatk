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

package org.broadinstitute.gatk.engine.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.engine.filters.BadCigarFilter;
import org.broadinstitute.gatk.engine.filters.MalformedReadFilter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.baq.BAQ;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.engine.recalibration.BQSRMode;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 1:53:31 PM
 * To change this template use File | Settings | File Templates.
 */
@ReadFilters({MalformedReadFilter.class,BadCigarFilter.class})
@PartitionBy(PartitionType.NONE)
@Downsample(by = DownsampleType.NONE)
@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = ReadTransformer.ApplicationTime.ON_INPUT)
@BQSRMode(ApplicationTime = ReadTransformer.ApplicationTime.ON_INPUT)
@DocumentedGATKFeature(groupName = "Uncategorized", extraDocs = {CommandLineGATK.class})
public abstract class Walker<MapType, ReduceType> {
    final protected static Logger logger = Logger.getLogger(Walker.class);
    private GenomeAnalysisEngine toolkit;

    protected Walker() {
    }

    /**
     * Set the toolkit, for peering into internal structures that can't
     * otherwise be read.
     * @param toolkit The genome analysis toolkit.
     */
    public void setToolkit(GenomeAnalysisEngine toolkit) {
        this.toolkit = toolkit;
    }

    /**
     * Retrieve the toolkit, for peering into internal structures that can't
     * otherwise be read.  Use sparingly, and discuss uses with software engineering
     * team.
     * @return The genome analysis toolkit.
     */
    protected GenomeAnalysisEngine getToolkit() {
        return toolkit;
    }

    /**
     * Gets the master sequence dictionary for this walker
     * @link GenomeAnalysisEngine.getMasterSequenceDictionary
     * @return the master sequence dictionary or null if no genome analysis toolkit.
     */
    protected SAMSequenceDictionary getMasterSequenceDictionary() {
        if ( toolkit == null )
            return null;
        else
            return toolkit.getMasterSequenceDictionary();
    }

    /**
     * Gets the GATK argument collection
     * @link GenomeAnalysisEngine.getArguments
     * @return the GATK argument collection or null if no genome analysis toolkit.
     */
    public GATKArgumentCollection getArguments(){
        if ( toolkit == null )
            return null;
        else
            return toolkit.getArguments();
    }

    /**
     * Gets the GATK samples database
     * @link GenomeAnalysisEngine.getSampleDB
     * @return the GATK samples database or null if no genome analysis toolkit.
     */
    public SampleDB getSampleDB() {
        if ( toolkit == null )
            return null;
        else
            return toolkit.getSampleDB();
    }

    /**
     * Gets a sample from the GATK samples database
     * @param id the sample ID
     * @return the sample from the GATK samples database or null if no genome analysis toolkit or samples database.
     */
    protected Sample getSample(final String id) {
        if ( getSampleDB() == null )
            return null;
        else
            return getSampleDB().getSample(id);
    }

    /**
     * (conceptual static) method that states whether you want to see reads piling up at a locus
     * that contain a deletion at the locus.
     *
     * ref:   ATCTGA
     * read1: ATCTGA
     * read2: AT--GA
     *
     * Normally, the locus iterator only returns a list of read1 at this locus at position 3, but
     * if this function returns true, then the system will return (read1, read2) with offsets
     * of (3, -1).  The -1 offset indicates a deletion in the read.
     *
     * @return false if you don't want to see deletions, or true if you do
     */
    public boolean includeReadsWithDeletionAtLoci() { 
        return false;
    }

    public void initialize() { }

    /**
     * A function for overloading in subclasses providing a mechanism to abort early from a walker.
     *
     * If this ever returns true, then the Traversal engine will stop executing map calls
     * and start the process of shutting down the walker in an orderly fashion.
     * @return
     */
    public boolean isDone() {
        return false;
    }

    /**
     * Provide an initial value for reduce computations.
     * @return Initial value of reduce.
     */
    public abstract ReduceType reduceInit();

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     * @param value result of the map.
     * @param sum accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public abstract ReduceType reduce(MapType value, ReduceType sum);    

    public void onTraversalDone(ReduceType result) {
        logger.info("[REDUCE RESULT] Traversal result is: " + result);
    }

    /**
     * General interval reduce routine called after all of the traversals are done
     * @param results interval reduce results
     */
    public void onTraversalDone(List<Pair<GenomeLoc, ReduceType>> results) {
        for ( Pair<GenomeLoc, ReduceType> result : results ) {
            logger.info(String.format("[INTERVAL REDUCE RESULT] at %s ", result.getFirst()));
            this.onTraversalDone(result.getSecond());
        }
    }

    /**
     * Return true if your walker wants to reduce each interval separately.  Default is false.
     *
     * If you set this flag, several things will happen.
     *
     * The system will invoke reduceInit() once for each interval being processed, starting a fresh reduce
     * Reduce will accumulate normally at each map unit in the interval
     * However, onTraversalDone(reduce) will be called after each interval is processed.
     * The system will call onTraversalDone( GenomeLoc -> reduce ), after all reductions are done,
     *   which is overloaded here to call onTraversalDone(reduce) for each location
     *
     * @return true if your walker wants to reduce each interval separately.
     */
    public boolean isReduceByInterval() {
        return false;
    }
}
