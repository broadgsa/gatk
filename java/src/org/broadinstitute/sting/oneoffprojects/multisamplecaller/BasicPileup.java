package org.broadinstitute.sting.oneoffprojects.multisamplecaller;

import org.broadinstitute.sting.utils.*;

import net.sf.samtools.*;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class BasicPileup {
    public static final char DELETION_CHAR = 'D';

    @Deprecated
    abstract GenomeLoc getLocation();
    @Deprecated
    abstract char getRef();
    @Deprecated
    abstract int size();


    /**
     * This is the right way to get bases
     *
     * @return
     */
    @Deprecated
    byte[] getBases() { return null; }

    /**
     * This is the right way to get quals
     *
     * @return
     */
    @Deprecated
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

    @Deprecated
    public static ArrayList<Byte> getBasesAsArrayList( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> bases = new ArrayList<Byte>(reads.size());
        for (byte value : getBasesAsArray(reads, offsets))
            bases.add(value);
        return bases;
     }

    @Deprecated
    public static ArrayList<Byte> getQualsAsArrayList( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> quals = new ArrayList<Byte>(reads.size());
        for (byte value : getQualsAsArray(reads, offsets))
            quals.add(value);
        return quals;
    }

    @Deprecated
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

    @Deprecated
    public static ArrayList<Byte> mappingQualPileup( List<SAMRecord> reads) {
        ArrayList<Byte> quals = new ArrayList<Byte>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            byte qual = (byte)read.getMappingQuality();
            quals.add(qual);
        }
        return quals;
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


