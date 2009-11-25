package org.broadinstitute.sting.utils;

import net.sf.samtools.*;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
@Deprecated
abstract public class BasicPileup {
    public static final char DELETION_CHAR = 'D';

    abstract GenomeLoc getLocation();
    abstract char getRef();
    abstract int size();


    /**
     * This is the right way to get bases
     *
     * @return
     */
    byte[] getBases() { return null; }

    /**
     * This is the right way to get quals
     *
     * @return
     */
    byte[] getQuals()  { return null; }

    /**
     * This is a terrible way to get bases.  Use getBases() or getBasesAsArrayList()
     *
     * @return
     */
    @Deprecated
    String getBasesAsString()  { return null; }

    /**
     * This is a terrible way to get quals.  Use getQuals() or getQualsAsArrayList()
     *
     * @return
     */
    @Deprecated
    String getQualsAsString()  { return null; }

    public String getPileupString() {
        return String.format("%s: %s %s %s", getLocation(), getRef(), getBasesAsString(), getQualsAsString());
    }

    public static String basePileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        StringBuilder bases = new StringBuilder();
        for ( byte base : getBasesAsArrayList(reads, offsets)) {
            bases.append((char)base);
        }
        return bases.toString();
    }

    public static String baseWithStrandPileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        StringBuilder bases = new StringBuilder();

        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            char base;
            if ( offset == -1 ) {
                base = DELETION_CHAR;
            } else {
                base = (char) read.getReadBases()[offset];
            }

            base = Character.toUpperCase(base);
            if (read.getReadNegativeStrandFlag()) {
                base = Character.toLowerCase(base);
            }

            bases.append(base);
        }

        return bases.toString();
    }

    //
    // byte[] methods
    //
    public static byte[] getBases( List<SAMRecord> reads, List<Integer> offsets ) {
        return getBasesAsArray(reads,offsets);
    }

    public static byte[] getQuals( List<SAMRecord> reads, List<Integer> offsets ) {
        return getQualsAsArray( reads, offsets );
    }

    //
    // ArrayList<Byte> methods
    //
    public static byte[] getBasesAsArray( List<SAMRecord> reads, List<Integer> offsets ) {
        byte array[] = new byte[reads.size()];
        int index = 0;
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            if ( offset == -1 ) {
               array[index++] = ((byte)DELETION_CHAR);
            } else {
                array[index++] = read.getReadBases()[offset];
            }
        }
        return array;
     }


    public static ArrayList<Byte> getBasesAsArrayList( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> bases = new ArrayList<Byte>(reads.size());
        for (byte value : getBasesAsArray(reads, offsets))
            bases.add(value);
        return bases;
     }

    public static ArrayList<Byte> getQualsAsArrayList( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> quals = new ArrayList<Byte>(reads.size());
        for (byte value : getQualsAsArray(reads, offsets))
            quals.add(value);
        return quals;
    }

    public static byte[] getQualsAsArray( List<SAMRecord> reads, List<Integer> offsets ) {
        byte array[] = new byte[reads.size()];
        int index = 0;
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            // skip deletion sites
            if ( offset == -1 ) {
                array[index++] = ((byte)0);
            } else {
                array[index++] = read.getBaseQualities()[offset];
            }
        }
        return array;
    }

    public static ArrayList<Byte> mappingQualPileup( List<SAMRecord> reads) {
        ArrayList<Byte> quals = new ArrayList<Byte>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            byte qual = (byte)read.getMappingQuality();
            quals.add(qual);
        }
        return quals;
    }

    public static String mappingQualPileupAsString( List<SAMRecord> reads) {
        return quals2String(mappingQualPileup(reads));
    }

    public static String quals2String( List<Byte> quals ) {
        StringBuilder qualStr = new StringBuilder();
        for ( int qual : quals ) {
            qual = Math.min(qual, 63);              // todo: fixme, this isn't a good idea
            char qualChar = (char) (33 + qual);     // todo: warning, this is illegal for qual > 63
            qualStr.append(qualChar);
        }

        return qualStr.toString();
    }

    public static String qualPileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        return quals2String(getQualsAsArrayList(reads, offsets));
    }


    public static ArrayList<Byte> getSecondaryBasesAsArrayList( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> bases2 = new ArrayList<Byte>(reads.size());
        boolean hasAtLeastOneSQorE2Field = false;

        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            byte base2 = BaseUtils.getSecondBase(read, offset);
            hasAtLeastOneSQorE2Field = hasAtLeastOneSQorE2Field || BaseUtils.simpleBaseToBaseIndex((char)base2) != -1;
            bases2.add(base2);
        }
        
        return (hasAtLeastOneSQorE2Field ? bases2 : null);
    }

    public static String getSecondaryBasePileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        StringBuilder bases2 = new StringBuilder();
        ArrayList<Byte> sbases = getSecondaryBasesAsArrayList(reads, offsets);

        if (sbases == null) { return null; }

        ArrayList<Byte> pbases = getBasesAsArrayList(reads, offsets);

        //Random generator = new Random();

        if ( sbases.size() != pbases.size() ) {
            throw new StingException("BUG in conversion of secondary bases: primary and secondary base vectors are different sizes!");
        }

        for (int pileupIndex = 0; pileupIndex < sbases.size(); pileupIndex++) {
            byte pbase = pbases.get(pileupIndex);
            byte sbase = sbases.get(pileupIndex);

            if ( sbase == pbase ) {
                throw new StingException("BUG in conversion of secondary bases!");
            }

//            while (sbase == pbase) { // TODO why is here?
//                sbase = (byte) BaseUtils.baseIndexToSimpleBase(generator.nextInt(4));
//            }

            bases2.append((char) sbase);
        }

        return bases2.toString();
    }

    @Deprecated // todo -- delete me
    public static String[] indelPileup( List<SAMRecord> reads, List<Integer> offsets ) 
    {
        String[] indels = new String[reads.size()];

        for (int i = 0; i < reads.size(); i++)
        {
            SAMRecord read = reads.get(i);
            Cigar cigar    = read.getCigar();
	    int offset     = offsets.get(i);

            String cigar_string = read.getCigarString();	 
			if (! (cigar_string.contains("I") || cigar_string.contains("D"))) { indels[i] = "null"; continue; }

			//System.out.printf("%s:%d %s %s %s ", read.getReferenceName(), read.getAlignmentStart(), read.getReadName(), read.getReadString(), cigar_string);
            int k = 0;
            for (int j = 0; j < cigar.numCigarElements(); j++)
            {
                CigarOperator operator = cigar.getCigarElement(j).getOperator();
                int           length   = cigar.getCigarElement(j).getLength();
                if (operator == CigarOperator.M) 
                { 
                    k += length; 
                }
                else if ((k == offset+1) && (operator == CigarOperator.I))
                {
                    // this insertion is associated with this offset (kinda ;) ).
  					indels[i] = read.getReadString().substring(k, k+length);	 
					//System.out.printf("(I,%d,%d)", k, offset);
					break;
                }
                else if ((k != offset+1) && (operator == CigarOperator.I)) 
                { 
					//System.out.printf("(i,%d,%d)", k, offset);
                    k += length; 
                }
                else if ((k == offset) && (operator == CigarOperator.D))
                {
					// this deletion is associated with this offset.
  					indels[i] = length + "D";
					//System.out.printf("(D,%d,%d)", k, offset);
   					break;
                }
                else if (k >= offset) 
                {
                    // no indel here.
                    indels[i] = "null";
					//System.out.printf("(N,%d,%d)", k, offset);
                    break;
                }
            }
			if (indels[i] == null) { indels[i] = "null"; }
			//System.out.printf("\n");
        }

        return indels;
    }
}


