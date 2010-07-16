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

package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * A set of static utility methods for common operations on VCF files/records.
 */
public class VCFUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private VCFUtils() { }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit) {

        // keep a map of sample name to occurrences encountered
        TreeSet<VCFHeaderLine> fields = new TreeSet<VCFHeaderLine>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getRecordType().equals(VCFRecord.class) ) {
                fields.addAll(rod.getHeader(VCFHeader.class).getMetaData());                
            }
        }

        return fields;
    }

    /**
     * Merges various vcf records into a single one using the mapping from rodNamesToSampleNames to get unique sample names
     *
     * @param rods                   the vcf rods
     * @param rodNamesToSampleNames  mapping of rod/sample pairs to new uniquified sample names
     * @return the new merged vcf record
     */
    public static VCFRecord mergeRecords(Map<VCFRecord,String> rods, Map<Pair<String, String>, String> rodNamesToSampleNames) {

        VCFParameters params = new VCFParameters();
        params.addFormatItem(VCFConstants.GENOTYPE_KEY);

        // keep track of the data so we can merge them intelligently
        double maxConfidence = 0.0;
        String id = null;
        Map<String, String> infoFields = new HashMap<String, String>();
        List<String> filters = new ArrayList<String>();

        for ( VCFRecord rod : rods.keySet() ) {
            List<VCFGenotypeRecord> myGenotypes = rod.getVCFGenotypeRecords();
            for ( VCFGenotypeRecord call : myGenotypes ) {
                // set the name to be the new uniquified name and add it to the list of genotypes
                call.setSampleName(rodNamesToSampleNames.get(new Pair<String, String>(rods.get(rod), call.getSampleName())));
                if ( params.getPosition() < 1 )
                    params.setLocations(GenomeLocParser.createGenomeLoc(rod.getChr(), rod.getStart()), call.getReference());
                params.addGenotypeRecord(createVCFGenotypeRecord(params, call, rod));
            }

            // set the overall confidence to be the max entry we see
            double confidence = 10.0 * rod.getNegLog10PError();
            if ( confidence > maxConfidence )
                maxConfidence = confidence;

            if ( rod.getID() != null )
                id = rod.getID();

            if ( rod.isFiltered() )
                filters.add(rod.getFilterString());

            // just take the last value we see for a given key
            infoFields.putAll(rod.getInfoValues());
        }

        return new VCFRecord(params.getReferenceBases(),
                params.getContig(),
                params.getPosition(),
                (id != null ? id : "."),
                params.getAlternateBases(),
                maxConfidence,
                filters.size() == 0 ? "0" : Utils.join(";", filters),
                infoFields,
                params.getFormatString(),
                params.getGenotypeRecords());
    }

    public static Set<VCFHeaderLine> smartMergeHeaders(Collection<VCFHeader> headers, Logger logger) throws IllegalStateException {
        HashMap<String, VCFHeaderLine> map = new HashMap<String, VCFHeaderLine>(); // from KEY.NAME -> line
        HashSet<VCFHeaderLine> lines = new HashSet<VCFHeaderLine>();

        // todo -- needs to remove all version headers from sources and add its own VCF version line
        for ( VCFHeader source : headers ) {
            //System.out.printf("Merging in header %s%n", source);
            for ( VCFHeaderLine line : source.getMetaData()) {
                String key = line.getKey();

                if ( line instanceof VCFNamedHeaderLine )
                    key = key + "." + ((VCFNamedHeaderLine) line).getName();

                if ( map.containsKey(key) ) {
                    VCFHeaderLine other = map.get(key);
                    if ( line.equals(other) || line.getVersion() != VCFHeaderVersion.VCF4_0 ) // todo -- remove me when everything is 4
                        continue;
                    else if ( ! line.getClass().equals(other.getClass()) )
                        throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    else if ( line instanceof VCFFilterHeaderLine ) {
                        String lineName = ((VCFFilterHeaderLine) line).getName();
                        String otherName = ((VCFFilterHeaderLine) other).getName();
                        if ( ! lineName.equals(otherName) )
                            throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    } else if ( line instanceof VCFCompoundHeaderLine ) {
                        VCFCompoundHeaderLine compLine = (VCFCompoundHeaderLine)line;
                        VCFCompoundHeaderLine compOther = (VCFCompoundHeaderLine)other;

                        // if the names are the same, but the values are different, we need to quit
                        if (! (compLine).equalsExcludingDescription(compOther) ) {
                            if ( compLine.getType().equals(compOther.getType()) ) {
                                // The Number entry is an Integer that describes the number of values that can be
                                // included with the INFO field. For example, if the INFO field contains a single
                                // number, then this value should be 1. However, if the INFO field describes a pair
                                // of numbers, then this value should be 2 and so on. If the number of possible
                                // values varies, is unknown, or is unbounded, then this value should be '.'.
                                if ( logger != null ) logger.warn("Promoting header field Number to . due to number differences in header lines: " + line + " " + other);
                                compOther.setNumberToUnbounded();
                            } else {
                                throw new IllegalStateException("Incompatible header types, collision between these two types: " + line + " " + other );
                            }
                        }
                        if ( ! compLine.getDescription().equals(compOther) )
                            if ( logger != null ) logger.warn(String.format("Allowing unequal description fields through: keeping " + compOther + " excluding " + compLine));
                    } else {
                        // we are not equal, but we're not anything special either
                        if ( logger != null ) logger.warn(String.format("Ignoring header line already in map: this header line = " + line + " already present header = " + other));
                    }
                } else {
                    line.setVersion(VCFHeaderVersion.VCF4_0); // todo -- remove this when we finally have vcf3/4 unified headers
                    map.put(key, line);
                    //System.out.printf("Adding header line %s%n", line);
                }
            }
        }

        return new HashSet<VCFHeaderLine>(map.values());
    }

    /**
     * create the VCF genotype record
     *
     * @param params the VCF parameters object
     * @param gtype  the genotype
     * @param vcfrecord  the VCF record
     *
     * @return a VCFGenotypeRecord
     */
    public static VCFGenotypeRecord createVCFGenotypeRecord(VCFParameters params, VCFGenotypeRecord gtype, VCFRecord vcfrecord) {

        List<VCFGenotypeEncoding> alleles = createAlleleArray(gtype);
        for (VCFGenotypeEncoding allele : alleles) {
            params.addAlternateBase(allele);
        }

        VCFGenotypeRecord record = new VCFGenotypeRecord(gtype.getSampleName(), alleles, VCFGenotypeRecord.PHASE.UNPHASED);
        for ( Map.Entry<String, String> entry : gtype.getFields().entrySet() ) {
            record.setField(entry.getKey(), entry.getValue());
            params.addFormatItem(entry.getKey());
        }

        record.setVCFRecord(vcfrecord);
        return record;
    }

    /**
     * create the allele array?
     *
     * @param gtype the gentoype object
     *
     * @return a list of string representing the string array of alleles
     */
    private static List<VCFGenotypeEncoding> createAlleleArray(VCFGenotypeRecord gtype) {
        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
        for (char allele : gtype.getBases().toCharArray()) {
            alleles.add(new VCFGenotypeEncoding(String.valueOf(allele)));
        }
        return alleles;
    }
}