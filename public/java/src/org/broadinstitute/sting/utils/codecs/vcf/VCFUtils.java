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

package org.broadinstitute.sting.utils.codecs.vcf;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * A set of static utility methods for common operations on VCF files/records.
 */
public class VCFUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private VCFUtils() { }

    public static Map<String, VCFHeader> getVCFHeadersFromRods(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        Map<String, VCFHeader> data = new HashMap<String, VCFHeader>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if it's not in our list
            if ( rodNames != null && !rodNames.contains(source.getName()) )
                continue;

            if ( source.getHeader() != null && source.getHeader() instanceof VCFHeader )
                data.put(source.getName(), (VCFHeader)source.getHeader());
        }

        return data;
    }

    public static Map<String,VCFHeader> getVCFHeadersFromRodPrefix(GenomeAnalysisEngine toolkit,String prefix) {
        Map<String, VCFHeader> data = new HashMap<String, VCFHeader>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if lacks the prefix
            if ( ! source.getName().startsWith(prefix) )
                continue;

            if ( source.getHeader() != null && source.getHeader() instanceof VCFHeader )
                data.put(source.getName(), (VCFHeader)source.getHeader());
        }

        return data;
    }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit) {
        return getHeaderFields(toolkit, null);
    }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     * @param rodNames   names of rods to use, or null if we should use all possible ones
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {

        // keep a map of sample name to occurrences encountered
        TreeSet<VCFHeaderLine> fields = new TreeSet<VCFHeaderLine>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if it's not in our list
            if ( rodNames != null && !rodNames.contains(source.getName()) )
                continue;

            if ( source.getRecordType().equals(VariantContext.class)) {
                VCFHeader header = (VCFHeader)source.getHeader();
                if ( header != null )
                    fields.addAll(header.getMetaData());
            }
        }

        return fields;
    }

    /** Only displays a warning if a logger is provided and an identical warning hasn't been already issued */
    private static final class HeaderConflictWarner {
        Logger logger;
        Set<String> alreadyIssued = new HashSet<String>();

        private HeaderConflictWarner(final Logger logger) {
            this.logger = logger;
        }

        public void warn(final VCFHeaderLine line, final String msg) {
            if ( logger != null && ! alreadyIssued.contains(line.getKey()) ) {
                alreadyIssued.add(line.getKey());
                logger.warn(msg);
            }
        }
    }

    public static Set<VCFHeaderLine> smartMergeHeaders(Collection<VCFHeader> headers, Logger logger) throws IllegalStateException {
        HashMap<String, VCFHeaderLine> map = new HashMap<String, VCFHeaderLine>(); // from KEY.NAME -> line
        HeaderConflictWarner conflictWarner = new HeaderConflictWarner(logger);

        // todo -- needs to remove all version headers from sources and add its own VCF version line
        for ( VCFHeader source : headers ) {
            //System.out.printf("Merging in header %s%n", source);
            for ( VCFHeaderLine line : source.getMetaData()) {
                String key = line.getKey();

                if ( line instanceof VCFNamedHeaderLine)
                    key = key + "" + ((VCFNamedHeaderLine) line).getName();

                if ( map.containsKey(key) ) {
                    VCFHeaderLine other = map.get(key);
                    if ( line.equals(other) )
                        continue;
                    else if ( ! line.getClass().equals(other.getClass()) )
                        throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    else if ( line instanceof VCFFilterHeaderLine) {
                        String lineName = ((VCFFilterHeaderLine) line).getName();                                                                                                         String otherName = ((VCFFilterHeaderLine) other).getName();
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
                                conflictWarner.warn(line, "Promoting header field Number to . due to number differences in header lines: " + line + " " + other);
                                compOther.setNumberToUnbounded();
                            } else if ( compLine.getType() == VCFHeaderLineType.Integer && compOther.getType() == VCFHeaderLineType.Float ) {
                                // promote key to Float
                                conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                                map.put(key, compOther);
                            } else if ( compLine.getType() == VCFHeaderLineType.Float && compOther.getType() == VCFHeaderLineType.Integer ) {
                                // promote key to Float
                                conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                            } else {
                                throw new IllegalStateException("Incompatible header types, collision between these two types: " + line + " " + other );
                            }
                        }
                        if ( ! compLine.getDescription().equals(compOther) )
                            conflictWarner.warn(line, "Allowing unequal description fields through: keeping " + compOther + " excluding " + compLine);
                    } else {
                        // we are not equal, but we're not anything special either
                        conflictWarner.warn(line, "Ignoring header line already in map: this header line = " + line + " already present header = " + other);
                    }
                } else {
                    map.put(key, line);
                    //System.out.printf("Adding header line %s%n", line);
                }
            }
        }

        return new HashSet<VCFHeaderLine>(map.values());
    }

    public static String rsIDOfFirstRealVariant(List<VariantContext> VCs, VariantContext.Type type) {
        if ( VCs == null )
            return null;

        String rsID = null;
        for ( VariantContext vc : VCs ) {
            if ( vc.getType() == type ) {
                rsID = vc.getID();
                break;
            }
        }

        return rsID;
    }
}