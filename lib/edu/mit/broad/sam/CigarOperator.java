package edu.mit.broad.sam;/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

/**
 * The operators that can appear in a cigar string.
 */
public enum CigarOperator {
    M,
    I,
    D,
    N,
    S,
    H,
    P,
    C; // I don't know what C means, but it is in the BAM spec

    // Readable synonyms of the above enums
    public static final CigarOperator MATCH_OR_MISMATCH = M;
    public static final CigarOperator INSERTION = I;
    public static final CigarOperator DELETION = D;
    public static final CigarOperator SKIPPED_REGION = N;
    public static final CigarOperator SOFT_CLIP = S;
    public static final CigarOperator HARD_CLIP = H;
    public static final CigarOperator PADDING = P;

    // Representation of CigarOperator in BAM file
    private static final byte OP_M = 0;
    private static final byte OP_I = 1;
    private static final byte OP_D = 2;
    private static final byte OP_N = 3;
    private static final byte OP_S = 4;
    private static final byte OP_H = 5;
    private static final byte OP_P = 6;
    private static final byte OP_C = 7;



    public static CigarOperator characterToEnum(final int b) {
        switch (b) {
        case 'M':
            return M;
        case 'I':
            return I;
        case 'D':
            return D;
        case 'N':
            return N;
        case 'S':
            return S;
        case 'H':
            return H;
        case 'P':
            return P;
        case 'C':
            return C;
        default:
            throw new IllegalArgumentException("Unrecognized CigarOperator: " + b);
        }
    }

    public static CigarOperator binaryToEnum(final int i) {
        switch(i) {
            case OP_M:
                return M;
            case OP_I:
                return I;
            case OP_D:
                return D;
            case OP_N:
                return N;
            case OP_S:
                return S;
            case OP_H:
                return H;
            case OP_P:
                return P;
            case OP_C:
                return C;
            default:
                throw new IllegalArgumentException("Unrecognized CigarOperator: " + i);
        }
    }

    public static int enumToBinary(final CigarOperator e) {
        switch(e) {
            case M:
                return OP_M;
            case I:
                return OP_I;
            case D:
                return OP_D;
            case N:
                return OP_N;
            case S:
                return OP_S;
            case H:
                return OP_H;
            case P:
                return OP_P;
            case C:
                return OP_C;
            default:
                throw new IllegalArgumentException("Unrecognized CigarOperator: " + e);
        }
    }
}
