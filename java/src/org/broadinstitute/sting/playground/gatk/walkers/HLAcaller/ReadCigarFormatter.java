/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

/**
 * Returns formatted read given read string and cigar string
 * Essentially removes header bases, soft clipped bases, and currently removes insertions
 * Deletions coded as "D"
 *
 * @author shermanjia
 */
public class ReadCigarFormatter {
    public String FormatRead(String cigar, String read){
        // returns a cigar-formatted sequence (removes insertions, inserts 'D' to where deletions occur
        String formattedRead = ""; char c; String count;
        int cigarPlaceholder = 0; int subcigarLength = 0;
        int readPlaceholder = 0; int subreadLength = 0;

        //reads cigar string
        for (int i = 0; i < cigar.length(); i++){
            c = cigar.charAt(i);
            if (c == 'M'){
                //If reach M for match/mismatch, get number immediately preceeding 'M' and tack on that many characters to sequence
                subcigarLength = i-cigarPlaceholder;
                count = cigar.substring(cigarPlaceholder, i);

                subreadLength = Integer.parseInt(count);
                formattedRead = formattedRead + read.substring(readPlaceholder, readPlaceholder+subreadLength);

                //increment placeholders
                cigarPlaceholder = i+1;
                readPlaceholder = readPlaceholder + subreadLength;
            } else if (c == 'I'){
                //***NOTE: To be modified later if needed (insertions removed here)***

                //If reaches I for insertion, get number before 'I' and skip that many characters in sequence
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);

                //increment placeholders without adding inserted bases to sequence (effectively removes insertion).
                cigarPlaceholder = i+1;
                readPlaceholder = readPlaceholder + subreadLength;
            } else if (c == 'H' || c == 'S'){
                //(H = Headers or S = Soft clipped removed here)***

                //If reaches H for insertion, get number before 'H' and skip that many characters in sequence
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);

                //increment cigar placeholder without adding inserted bases to sequence (effectively removes insertion).
                cigarPlaceholder = i+1;
            } else if (c == 'D'){
                //If reaches D for deletion, insert 'D' into sequence as placeholder
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);

                //Add one 'D' for each deleted base
                String deletion = "";
                for (int j = 1; j <= subreadLength; j++){
                    deletion = deletion + "D";
                }

                //update placeholders
                formattedRead = formattedRead + deletion;
                cigarPlaceholder = i+1;
            }

        }
        return formattedRead;
    }

}
