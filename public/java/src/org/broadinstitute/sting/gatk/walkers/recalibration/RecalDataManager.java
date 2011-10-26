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

package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 6, 2009
 *
 * This helper class holds the data HashMap as well as submaps that represent the marginal distributions collapsed over all needed dimensions.
 * It also has static methods that are used to perform the various solid recalibration modes that attempt to correct the reference bias.
 * This class holds the parsing methods that are shared between CountCovariates and TableRecalibration.
 */

public class RecalDataManager {

    public final NestedHashMap data; // The full dataset
    private final NestedHashMap dataCollapsedReadGroup; // Table where everything except read group has been collapsed
    private final NestedHashMap dataCollapsedQualityScore; // Table where everything except read group and quality score has been collapsed
    private final ArrayList<NestedHashMap> dataCollapsedByCovariate; // Tables where everything except read group, quality score, and given covariate has been collapsed

    public final static String ORIGINAL_QUAL_ATTRIBUTE_TAG = "OQ"; // The tag that holds the original quality scores
    public final static String COLOR_SPACE_QUAL_ATTRIBUTE_TAG = "CQ"; // The tag that holds the color space quality scores for SOLID bams
    public final static String COLOR_SPACE_ATTRIBUTE_TAG = "CS"; // The tag that holds the color space for SOLID bams
    public final static String COLOR_SPACE_INCONSISTENCY_TAG = "ZC"; // A new tag made up for the recalibrator which will hold an array of ints which say if this base is inconsistent with its color
    private static boolean warnUserNullReadGroup = false;
    private static boolean warnUserNullPlatform = false;

    public enum SOLID_RECAL_MODE {
        /** Treat reference inserted bases as reference matching bases. Very unsafe! */
        DO_NOTHING,
        /** Set reference inserted bases and the previous base (because of color space alignment details) to Q0. This is the default option. */
        SET_Q_ZERO,
        /** In addition to setting the quality scores to zero, also set the base itself to 'N'. This is useful to visualize in IGV. */
        SET_Q_ZERO_BASE_N,
        /** Look at the color quality scores and probabilistically decide to change the reference inserted base to be the base which is implied by the original color space instead of the reference. */
        REMOVE_REF_BIAS
    }

    public enum SOLID_NOCALL_STRATEGY {
        /** When a no call is detected throw an exception to alert the user that recalibrating this SOLiD data is unsafe. This is the default option. */
        THROW_EXCEPTION,
        /** Leave the read in the output bam completely untouched. This mode is only okay if the no calls are very rare. */
        LEAVE_READ_UNRECALIBRATED,
        /** Mark these reads as failing vendor quality checks so they can be filtered out by downstream analyses. */
        PURGE_READ
    }

    RecalDataManager() {
        data = new NestedHashMap();
        dataCollapsedReadGroup = null;
        dataCollapsedQualityScore = null;
        dataCollapsedByCovariate = null;
    }

    RecalDataManager( final boolean createCollapsedTables, final int numCovariates ) {
        if( createCollapsedTables ) { // Initialize all the collapsed tables, only used by TableRecalibrationWalker
            data = null;
            dataCollapsedReadGroup = new NestedHashMap();
            dataCollapsedQualityScore = new NestedHashMap();
            dataCollapsedByCovariate = new ArrayList<NestedHashMap>();
            for( int iii = 0; iii < numCovariates - 2; iii++ ) { // readGroup and QualityScore aren't counted here, their tables are separate
                dataCollapsedByCovariate.add( new NestedHashMap() );
            }
        } else {
            data = new NestedHashMap();
            dataCollapsedReadGroup = null;
            dataCollapsedQualityScore = null;
            dataCollapsedByCovariate = null;
        }
    }

    /**
     * Add the given mapping to all of the collapsed hash tables
     * @param key The list of comparables that is the key for this mapping
     * @param fullDatum The RecalDatum which is the data for this mapping
     * @param PRESERVE_QSCORES_LESS_THAN The threshold in report quality for adding to the aggregate collapsed table
     */
    public final void addToAllTables( final Object[] key, final RecalDatum fullDatum, final int PRESERVE_QSCORES_LESS_THAN ) {

        // The full dataset isn't actually ever used for anything because of the sequential calculation so no need to keep the full data HashMap around
        //data.put(key, thisDatum); // add the mapping to the main table

        final int qualityScore = Integer.parseInt( key[1].toString() );
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] qualityScoreCollapsedKey = new Object[2];
        final Object[] covariateCollapsedKey = new Object[3];
        RecalDatum collapsedDatum;

        // Create dataCollapsedReadGroup, the table where everything except read group has been collapsed
        if( qualityScore >= PRESERVE_QSCORES_LESS_THAN ) {
            readGroupCollapsedKey[0] = key[0]; // Make a new key with just the read group
            collapsedDatum = (RecalDatum) dataCollapsedReadGroup.get( readGroupCollapsedKey );
            if( collapsedDatum == null ) {
                dataCollapsedReadGroup.put( new RecalDatum(fullDatum), readGroupCollapsedKey );
            } else {
                collapsedDatum.combine( fullDatum ); // using combine instead of increment in order to calculate overall aggregateQReported
            }
        }

        // Create dataCollapsedQuality, the table where everything except read group and quality score has been collapsed
        qualityScoreCollapsedKey[0] = key[0]; // Make a new key with the read group ...
        qualityScoreCollapsedKey[1] = key[1]; //                                    and quality score
        collapsedDatum = (RecalDatum) dataCollapsedQualityScore.get( qualityScoreCollapsedKey );
        if( collapsedDatum == null ) {
            dataCollapsedQualityScore.put( new RecalDatum(fullDatum), qualityScoreCollapsedKey );
        } else {
            collapsedDatum.increment( fullDatum );
        }

        // Create dataCollapsedByCovariate's, the tables where everything except read group, quality score, and given covariate has been collapsed
        for( int iii = 0; iii < dataCollapsedByCovariate.size(); iii++ ) {
            covariateCollapsedKey[0] = key[0]; // Make a new key with the read group ...
            covariateCollapsedKey[1] = key[1]; //                                    and quality score ...
            final Object theCovariateElement = key[iii + 2]; //                                        and the given covariate
            if( theCovariateElement != null ) {
                covariateCollapsedKey[2] = theCovariateElement;
                collapsedDatum = (RecalDatum) dataCollapsedByCovariate.get(iii).get( covariateCollapsedKey );
                if( collapsedDatum == null ) {
                    dataCollapsedByCovariate.get(iii).put( new RecalDatum(fullDatum), covariateCollapsedKey );
                } else {
                    collapsedDatum.increment( fullDatum );
                }
            }
        }
    }

    /**
     * Loop over all the collapsed tables and turn the recalDatums found there into an empirical quality score
     *   that will be used in the sequential calculation in TableRecalibrationWalker
     * @param smoothing The smoothing parameter that goes into empirical quality score calculation
     * @param maxQual At which value to cap the quality scores
     */
    public final void generateEmpiricalQualities( final int smoothing, final int maxQual ) {

        recursivelyGenerateEmpiricalQualities(dataCollapsedReadGroup.data, smoothing, maxQual);
        recursivelyGenerateEmpiricalQualities(dataCollapsedQualityScore.data, smoothing, maxQual);
        for( NestedHashMap map : dataCollapsedByCovariate ) {
            recursivelyGenerateEmpiricalQualities(map.data, smoothing, maxQual);
            checkForSingletons(map.data);
        }
    }

    private void recursivelyGenerateEmpiricalQualities( final Map data, final int smoothing, final int maxQual ) {

        for( Object comp : data.keySet() ) {
            final Object val = data.get(comp);
            if( val instanceof RecalDatum ) { // We are at the end of the nested hash maps
                ((RecalDatum)val).calcCombinedEmpiricalQuality(smoothing, maxQual);
            } else { // Another layer in the nested hash map
                recursivelyGenerateEmpiricalQualities( (Map) val, smoothing, maxQual);
            }
        }
    }

    private void checkForSingletons( final Map data ) {
        // todo -- this looks like it's better just as a data.valueSet() call?
        for( Object comp : data.keySet() ) {
            final Object val = data.get(comp);
            if( val instanceof RecalDatum ) { // We are at the end of the nested hash maps
                if( data.keySet().size() == 1) {
                    data.clear(); // don't TableRecalibrate a non-required covariate if it only has one element because that correction has already been done ...
                                                                                                    // in a previous step of the sequential calculation model
                }
            } else { // Another layer in the nested hash map
                checkForSingletons( (Map) val );
            }
        }
    }

    /**
     * Get the appropriate collapsed table out of the set of all the tables held by this Object
     * @param covariate Which covariate indexes the desired collapsed HashMap
     * @return The desired collapsed HashMap
     */
    public final NestedHashMap getCollapsedTable( final int covariate ) {
        if( covariate == 0) {
            return dataCollapsedReadGroup; // Table where everything except read group has been collapsed
        } else if( covariate == 1 ) {
            return dataCollapsedQualityScore; // Table where everything except read group and quality score has been collapsed
        } else {
            return dataCollapsedByCovariate.get( covariate - 2 ); // Table where everything except read group, quality score, and given covariate has been collapsed
        }
    }

    /**
     * Section of code shared between the two recalibration walkers which uses the command line arguments to adjust attributes of the read such as quals or platform string
     * @param read The read to adjust
     * @param RAC The list of shared command line arguments
     */
    public static void parseSAMRecord( final SAMRecord read, final RecalibrationArgumentCollection RAC ) {
        GATKSAMReadGroupRecord readGroup = ((GATKSAMRecord)read).getReadGroup();

        // If there are no read groups we have to default to something, and that something could be specified by the user using command line arguments
        if( readGroup == null ) {
            if( RAC.DEFAULT_READ_GROUP != null && RAC.DEFAULT_PLATFORM != null) {
                if( !warnUserNullReadGroup && RAC.FORCE_READ_GROUP == null ) {
                    Utils.warnUser("The input .bam file contains reads with no read group. " +
                                    "Defaulting to read group ID = " + RAC.DEFAULT_READ_GROUP + " and platform = " + RAC.DEFAULT_PLATFORM + ". " +
                                    "First observed at read with name = " + read.getReadName() );
                    warnUserNullReadGroup = true;
                }
                // There is no readGroup so defaulting to these values
                readGroup = new GATKSAMReadGroupRecord( RAC.DEFAULT_READ_GROUP );
                readGroup.setPlatform( RAC.DEFAULT_PLATFORM );
                ((GATKSAMRecord)read).setReadGroup( readGroup );
            } else {
                throw new UserException.MalformedBAM(read, "The input .bam file contains reads with no read group. First observed at read with name = " + read.getReadName() );
            }
        }

        if( RAC.FORCE_READ_GROUP != null && !readGroup.getReadGroupId().equals(RAC.FORCE_READ_GROUP) ) { // Collapse all the read groups into a single common String provided by the user
            final String oldPlatform = readGroup.getPlatform();
            readGroup = new GATKSAMReadGroupRecord( RAC.FORCE_READ_GROUP );
            readGroup.setPlatform( oldPlatform );
            ((GATKSAMRecord)read).setReadGroup( readGroup );
        }

        if( RAC.FORCE_PLATFORM != null && (readGroup.getPlatform() == null || !readGroup.getPlatform().equals(RAC.FORCE_PLATFORM))) {
            readGroup.setPlatform( RAC.FORCE_PLATFORM );
        }

        if ( readGroup.getPlatform() == null ) {
            if( RAC.DEFAULT_PLATFORM != null ) {
                if( !warnUserNullPlatform ) {
                    Utils.warnUser("The input .bam file contains reads with no platform information. " +
                                        "Defaulting to platform = " + RAC.DEFAULT_PLATFORM + ". " +
                                        "First observed at read with name = " + read.getReadName() );
                    warnUserNullPlatform = true;
                }
                readGroup.setPlatform( RAC.DEFAULT_PLATFORM );
            } else {
                throw new UserException.MalformedBAM(read, "The input .bam file contains reads with no platform information. First observed at read with name = " + read.getReadName() );
            }
        }
    }

    /**
     * Parse through the color space of the read and add a new tag to the SAMRecord that says which bases are inconsistent with the color space
     * @param read The SAMRecord to parse
     */
    public static void parseColorSpace( final SAMRecord read ) {

        // If this is a SOLID read then we have to check if the color space is inconsistent. This is our only sign that SOLID has inserted the reference base
        if( read.getReadGroup().getPlatform().toUpperCase().contains("SOLID") ) {
            if( read.getAttribute(RecalDataManager.COLOR_SPACE_INCONSISTENCY_TAG) == null ) { // Haven't calculated the inconsistency array yet for this read
                final Object attr = read.getAttribute(RecalDataManager.COLOR_SPACE_ATTRIBUTE_TAG);
                if( attr != null ) {
                    byte[] colorSpace;
                    if( attr instanceof String ) {
                        colorSpace = ((String)attr).getBytes();
                    } else {
                        throw new UserException.MalformedBAM(read, String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.COLOR_SPACE_ATTRIBUTE_TAG, read.getReadName()));
                    }

                    // Loop over the read and calculate first the inferred bases from the color and then check if it is consistent with the read
                    byte[] readBases = read.getReadBases();
                    if( read.getReadNegativeStrandFlag() ) {
                        readBases = BaseUtils.simpleReverseComplement( read.getReadBases() );
                    }
                    final byte[] inconsistency = new byte[readBases.length];
                    int iii;
                    byte prevBase = colorSpace[0]; // The sentinel
                    for( iii = 0; iii < readBases.length; iii++ ) {
                        final byte thisBase = getNextBaseFromColor( read, prevBase, colorSpace[iii + 1] );
                        inconsistency[iii] = (byte)( thisBase == readBases[iii] ? 0 : 1 );
                        prevBase = readBases[iii];
                    }
                    read.setAttribute( RecalDataManager.COLOR_SPACE_INCONSISTENCY_TAG, inconsistency );

                } else {
                    throw new UserException.MalformedBAM(read, "Unable to find color space information in SOLiD read. First observed at read with name = " + read.getReadName() +
                                            " Unfortunately this .bam file can not be recalibrated without color space information because of potential reference bias.");
                }
            }
        }
    }

    /**
     * Parse through the color space of the read and apply the desired --solid_recal_mode correction to the bases
     * This method doesn't add the inconsistent tag to the read like parseColorSpace does
     * @param read The SAMRecord to parse
     * @param originalQualScores The array of original quality scores to modify during the correction
     * @param solidRecalMode Which mode of solid recalibration to apply
     * @param refBases The reference for this read
     * @return A new array of quality scores that have been ref bias corrected
     */
    public static byte[] calcColorSpace( final SAMRecord read, byte[] originalQualScores, final SOLID_RECAL_MODE solidRecalMode, final byte[] refBases ) {

        final Object attr = read.getAttribute(RecalDataManager.COLOR_SPACE_ATTRIBUTE_TAG);
        if( attr != null ) {
            byte[] colorSpace;
            if( attr instanceof String ) {
                colorSpace = ((String)attr).getBytes();
            } else {
                throw new ReviewedStingException(String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.COLOR_SPACE_ATTRIBUTE_TAG, read.getReadName()));
            }

            // Loop over the read and calculate first the inferred bases from the color and then check if it is consistent with the read
            byte[] readBases = read.getReadBases();
            final byte[] colorImpliedBases = readBases.clone();
            byte[] refBasesDirRead = AlignmentUtils.alignmentToByteArray( read.getCigar(), read.getReadBases(), refBases ); //BUGBUG: This needs to change when read walkers are changed to give the aligned refBases
            if( read.getReadNegativeStrandFlag() ) {
                readBases = BaseUtils.simpleReverseComplement( read.getReadBases() );
                refBasesDirRead = BaseUtils.simpleReverseComplement( refBasesDirRead.clone() );
            }
            final int[] inconsistency = new int[readBases.length];
            byte prevBase = colorSpace[0]; // The sentinel
            for( int iii = 0; iii < readBases.length; iii++ ) {
                final byte thisBase = getNextBaseFromColor( read, prevBase, colorSpace[iii + 1] );
                colorImpliedBases[iii] = thisBase;
                inconsistency[iii] = ( thisBase == readBases[iii] ? 0 : 1 );
                prevBase = readBases[iii];
            }

            // Now that we have the inconsistency array apply the desired correction to the inconsistent bases
            if( solidRecalMode == SOLID_RECAL_MODE.SET_Q_ZERO ) { // Set inconsistent bases and the one before it to Q0
                final boolean setBaseN = false;
                originalQualScores = solidRecalSetToQZero(read, readBases, inconsistency, originalQualScores, refBasesDirRead, setBaseN);
            } else if( solidRecalMode == SOLID_RECAL_MODE.SET_Q_ZERO_BASE_N ) {
                final boolean setBaseN = true;
                originalQualScores = solidRecalSetToQZero(read, readBases, inconsistency, originalQualScores, refBasesDirRead, setBaseN);
            } else if( solidRecalMode == SOLID_RECAL_MODE.REMOVE_REF_BIAS ) { // Use the color space quality to probabilistically remove ref bases at inconsistent color space bases
                solidRecalRemoveRefBias(read, readBases, inconsistency, colorImpliedBases, refBasesDirRead);
            }

        } else {
            throw new UserException.MalformedBAM(read, "Unable to find color space information in SOLiD read. First observed at read with name = " + read.getReadName() +
                    " Unfortunately this .bam file can not be recalibrated without color space information because of potential reference bias.");
        }

        return originalQualScores;
    }

    public static boolean checkNoCallColorSpace( final SAMRecord read ) {
        if( read.getReadGroup().getPlatform().toUpperCase().contains("SOLID") ) {
            final Object attr = read.getAttribute(RecalDataManager.COLOR_SPACE_ATTRIBUTE_TAG);
            if( attr != null ) {
                byte[] colorSpace;
                if( attr instanceof String ) {
                    colorSpace = ((String)attr).substring(1).getBytes(); // trim off the Sentinel
                } else {
                    throw new ReviewedStingException(String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.COLOR_SPACE_ATTRIBUTE_TAG, read.getReadName()));
                }

                for( byte color : colorSpace ) {
                    if( color != (byte)'0' && color != (byte)'1' && color != (byte)'2' && color != (byte)'3' ) {
                        return true; // There is a bad color in this SOLiD read and the user wants to skip over it
                    }
                }

            } else {
                throw new UserException.MalformedBAM(read, "Unable to find color space information in SOLiD read. First observed at read with name = " + read.getReadName() +
                                        " Unfortunately this .bam file can not be recalibrated without color space information because of potential reference bias.");
            }
        }

        return false; // There aren't any color no calls in this SOLiD read
    }

    /**
     * Perform the SET_Q_ZERO solid recalibration. Inconsistent color space bases and their previous base are set to quality zero
     * @param read The SAMRecord to recalibrate
     * @param readBases The bases in the read which have been RC'd if necessary
     * @param inconsistency The array of 1/0 that says if this base is inconsistent with its color
     * @param originalQualScores The array of original quality scores to set to zero if needed
     * @param refBases The reference which has been RC'd if necessary
     * @param setBaseN Should we also set the base to N as well as quality zero in order to visualize in IGV or something similar
     * @return The byte array of original quality scores some of which might have been set to zero
     */
    private static byte[] solidRecalSetToQZero( final SAMRecord read, byte[] readBases, final int[] inconsistency, final byte[] originalQualScores,
                                                final byte[] refBases, final boolean setBaseN ) {

        final boolean negStrand = read.getReadNegativeStrandFlag();
        for( int iii = 1; iii < originalQualScores.length; iii++ ) {
            if( inconsistency[iii] == 1 ) {
                if( readBases[iii] == refBases[iii] ) {
                    if( negStrand ) { originalQualScores[originalQualScores.length-(iii+1)] = (byte)0; }
                    else { originalQualScores[iii] = (byte)0; }
                    if( setBaseN ) { readBases[iii] = (byte)'N'; }
                }
                // Set the prev base to Q0 as well
                if( readBases[iii-1] == refBases[iii-1] ) {
                    if( negStrand ) { originalQualScores[originalQualScores.length-iii] = (byte)0; }
                    else { originalQualScores[iii-1] = (byte)0; }
                    if( setBaseN ) { readBases[iii-1] = (byte)'N'; }
                }
            }
        }
        if( negStrand ) {
            readBases = BaseUtils.simpleReverseComplement( readBases.clone() ); // Put the bases back in reverse order to stuff them back in the read
        }
        read.setReadBases( readBases );

        return originalQualScores;
    }

    /**
     * Peform the REMOVE_REF_BIAS solid recalibration. Look at the color space qualities and probabilistically decide if the base should be change to match the color or left as reference
     * @param read The SAMRecord to recalibrate
     * @param readBases The bases in the read which have been RC'd if necessary
     * @param inconsistency The array of 1/0 that says if this base is inconsistent with its color
     * @param colorImpliedBases The bases implied by the color space, RC'd if necessary
     * @param refBases The reference which has been RC'd if necessary
     */
    private static void solidRecalRemoveRefBias( final SAMRecord read, byte[] readBases, final int[] inconsistency, final byte[] colorImpliedBases,
                                                 final byte[] refBases) {

        final Object attr = read.getAttribute(RecalDataManager.COLOR_SPACE_QUAL_ATTRIBUTE_TAG);
        if( attr != null ) {
            byte[] colorSpaceQuals;
            if( attr instanceof String ) {
                String x = (String)attr;
                colorSpaceQuals = x.getBytes();
                SAMUtils.fastqToPhred(colorSpaceQuals);
            } else {
                throw new ReviewedStingException(String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.COLOR_SPACE_QUAL_ATTRIBUTE_TAG, read.getReadName()));
            }

            for( int iii = 1; iii < inconsistency.length - 1; iii++ ) {
                if( inconsistency[iii] == 1 ) {
                    for( int jjj = iii - 1; jjj <= iii; jjj++ ) { // Correct this base and the one before it along the direction of the read
                        if( jjj == iii || inconsistency[jjj] == 0 ) { // Don't want to correct the previous base a second time if it was already corrected in the previous step
                            if( readBases[jjj] == refBases[jjj] ) {
                                if( colorSpaceQuals[jjj] == colorSpaceQuals[jjj+1] ) { // Equal evidence for the color implied base and the reference base, so flip a coin
                                    final int rand = GenomeAnalysisEngine.getRandomGenerator().nextInt( 2 );
                                    if( rand == 0 ) { // The color implied base won the coin flip
                                        readBases[jjj] = colorImpliedBases[jjj];
                                    }
                                } else {
                                    final int maxQuality = Math.max((int)colorSpaceQuals[jjj], (int)colorSpaceQuals[jjj+1]);
                                    final int minQuality = Math.min((int)colorSpaceQuals[jjj], (int)colorSpaceQuals[jjj+1]);
                                    int diffInQuality = maxQuality - minQuality;
                                    int numLow = minQuality;
                                    if( numLow == 0 ) {
                                        numLow++;
                                        diffInQuality++;
                                    }
                                    final int numHigh = Math.round( numLow * (float)Math.pow(10.0f, (float) diffInQuality / 10.0f) ); // The color with higher quality is exponentially more likely
                                    final int rand = GenomeAnalysisEngine.getRandomGenerator().nextInt( numLow + numHigh );
                                    if( rand >= numLow ) { // higher q score won
                                        if( maxQuality == (int)colorSpaceQuals[jjj] ) {
                                            readBases[jjj] = colorImpliedBases[jjj];
                                        } // else ref color had higher q score, and won out, so nothing to do here
                                    } else { // lower q score won
                                        if( minQuality == (int)colorSpaceQuals[jjj] ) {
                                            readBases[jjj] = colorImpliedBases[jjj];
                                        } // else ref color had lower q score, and won out, so nothing to do here
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if( read.getReadNegativeStrandFlag() ) {
                readBases = BaseUtils.simpleReverseComplement( readBases.clone() ); // Put the bases back in reverse order to stuff them back in the read
            }
            read.setReadBases( readBases );
        } else { // No color space quality tag in file
            throw new UserException.MalformedBAM(read, "REMOVE_REF_BIAS recal mode requires color space qualities but they can't be found for read: " + read.getReadName());
        }
    }

    /**
     * Given the base and the color calculate the next base in the sequence
     * @param prevBase The base
     * @param color The color
     * @return The next base in the sequence
     */
    private static byte getNextBaseFromColor( SAMRecord read, final byte prevBase, final byte color ) {
        switch(color) {
            case '0':
                return prevBase;
            case '1':
                return performColorOne( prevBase );
            case '2':
                return performColorTwo( prevBase );
            case '3':
                return performColorThree( prevBase );
            default:
                throw new UserException.MalformedBAM(read, "Unrecognized color space in SOLID read, color = " + (char)color +
                                          " Unfortunately this bam file can not be recalibrated without full color space information because of potential reference bias.");
        }
    }

    /**
     * Check if this base is inconsistent with its color space. If it is then SOLID inserted the reference here and we should reduce the quality
     * @param read The read which contains the color space to check against
     * @param offset The offset in the read at which to check
     * @return Returns true if the base was inconsistent with the color space
     */
    public static boolean isInconsistentColorSpace( final SAMRecord read, final int offset ) {
        final Object attr = read.getAttribute(RecalDataManager.COLOR_SPACE_INCONSISTENCY_TAG);
        if( attr != null ) {
            final byte[] inconsistency = (byte[])attr;
            // NOTE: The inconsistency array is in the direction of the read, not aligned to the reference!
            if( read.getReadNegativeStrandFlag() ) { // Negative direction
                return inconsistency[inconsistency.length - offset - 1] != (byte)0;
            } else { // Forward direction
                return inconsistency[offset] != (byte)0;
            }

            // This block of code is for if you want to check both the offset and the next base for color space inconsistency
            //if( read.getReadNegativeStrandFlag() ) { // Negative direction
            //    if( offset == 0 ) {
            //        return inconsistency[0] != 0;
            //    } else {
            //        return (inconsistency[inconsistency.length - offset - 1] != 0) || (inconsistency[inconsistency.length - offset] != 0);
            //    }
            //} else { // Forward direction
            //    if( offset == inconsistency.length - 1 ) {
            //        return inconsistency[inconsistency.length - 1] != 0;
            //    } else {
            //        return (inconsistency[offset] != 0) || (inconsistency[offset + 1] != 0);
            //    }
            //}

        } else { // No inconsistency array, so nothing is inconsistent
            return false;
        }
    }

    /**
     * Computes all requested covariates for every offset in the given read
     * by calling covariate.getValues(..).
     *
     * @param gatkRead The read for which to compute covariate values.
     * @param requestedCovariates The list of requested covariates.
     * @return An array of covariate values where result[i][j] is the covariate
     * value for the ith position in the read and the jth covariate in
     * reqeustedCovariates list.
     */
     public static Comparable[][] computeCovariates(final GATKSAMRecord gatkRead, final List<Covariate> requestedCovariates) {
         //compute all covariates for this read
         final List<Covariate> requestedCovariatesRef = requestedCovariates;
         final int numRequestedCovariates = requestedCovariatesRef.size();
         final int readLength = gatkRead.getReadLength();

         final Comparable[][] covariateValues_offset_x_covar = new Comparable[readLength][numRequestedCovariates];
         final Comparable[] tempCovariateValuesHolder = new Comparable[readLength];

         // Loop through the list of requested covariates and compute the values of each covariate for all positions in this read
         for( int i = 0; i < numRequestedCovariates; i++ ) {
             requestedCovariatesRef.get(i).getValues( gatkRead, tempCovariateValuesHolder );
             for(int j = 0; j < readLength; j++) {
                 //copy values into a 2D array that allows all covar types to be extracted at once for
                 //an offset j by doing covariateValues_offset_x_covar[j]. This avoids the need to later iterate over covar types.
                 covariateValues_offset_x_covar[j][i] = tempCovariateValuesHolder[j];
             }
         }

         return covariateValues_offset_x_covar;
     }

    /**
     * Perform a ceratin transversion (A <-> C or G <-> T) on the base.
     *
     * @param base the base [AaCcGgTt]
     * @return the transversion of the base, or the input base if it's not one of the understood ones
     */
    private static byte performColorOne(byte base) {
        switch (base) {
            case 'A':
            case 'a': return 'C';
            case 'C':
            case 'c': return 'A';
            case 'G':
            case 'g': return 'T';
            case 'T':
            case 't': return 'G';
            default: return base;
        }
    }

    /**
     * Perform a transition (A <-> G or C <-> T) on the base.
     *
     * @param base the base [AaCcGgTt]
     * @return the transition of the base, or the input base if it's not one of the understood ones
     */
    private static byte performColorTwo(byte base) {
        switch (base) {
            case 'A':
            case 'a': return 'G';
            case 'C':
            case 'c': return 'T';
            case 'G':
            case 'g': return 'A';
            case 'T':
            case 't': return 'C';
            default: return base;
        }
    }

    /**
     * Return the complement (A <-> T or C <-> G) of a base.
     *
     * @param base the base [AaCcGgTt]
     * @return the complementary base, or the input base if it's not one of the understood ones
     */
    private static byte performColorThree(byte base) {
        switch (base) {
            case 'A':
            case 'a': return 'T';
            case 'C':
            case 'c': return 'G';
            case 'G':
            case 'g': return 'C';
            case 'T':
            case 't': return 'A';
            default: return base;
        }
    }
}
