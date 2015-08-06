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

package org.broadinstitute.gatk.engine.crypt;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.io.IOUtils;

import java.io.*;
import java.security.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Class to represent a GATK user key.
 *
 * A GATK user key contains an email address and a cryptographic signature.
 * The signature is the SHA-1 hash of the email address encrypted using
 * the GATK master private key. The GATK master public key (distributed
 * with the GATK) is used to decrypt the signature and validate the key
 * at the start of each GATK run that requires a key.
 *
 * Keys are cryptographically secure in that valid keys definitely come
 * from us and cannot be fabricated, however nothing prevents keys from
 * being shared between users.
 *
 * GATK user keys have the following on-disk format:
 *
 *     GZIP Container:
 *         Email address
 *         NUL byte (delimiter)
 *         Cryptographic Signature (encrypted SHA-1 hash of email address)
 *
 * The key data is wrapped within a GZIP container to placate over-zealous
 * email filters (since keys must often be emailed) and also to provide an
 * additional integrity check via the built-in GZIP CRC.
 *
 * @author David Roazen
 */
public class GATKKey {

    /**
     * Private key used to sign the GATK key. Required only when creating a new
     * key from scratch, not when loading an existing key from disk.
     */
    private PrivateKey privateKey;

    /**
     * Public key used to validate the GATK key.
     */
    private PublicKey publicKey;

    /**
     * The user's email address, stored within the key and signed.
     */
    private String emailAddress;

    /**
     * The cryptographic signature of the email address. By default, this is
     * the SHA-1 hash of the email address encrypted using the RSA algorithm.
     */
    private byte[] signature;

    /**
     * The combination of hash/encryption algorithms to use to generate the signature.
     * By default this is "SHA1withRSA"
     */
    private String signingAlgorithm;

    /**
     * Default hash/encryption algorithms to use to sign the key.
     */
    public static final String DEFAULT_SIGNING_ALGORITHM = "SHA1withRSA";

    /**
     * Byte value used to separate the email address from its signature in the key file.
     */
    public static final byte GATK_KEY_SECTIONAL_DELIMITER = 0;


    // -----------------------
    // Constructors:
    // -----------------------

    /**
     * Constructor to create a new GATK key from scratch using an email address
     * and public/private key pair. The private key is used for signing, and the
     * public key is used to validate the newly-created key.
     *
     * @param privateKey Private key used to sign the new GATK key
     * @param publicKey Public key used to validate the new GATK key
     * @param emailAddress The user's email address, which we will store in the key and sign
     */
    public GATKKey ( PrivateKey privateKey, PublicKey publicKey, String emailAddress ) {
        this(privateKey, publicKey, emailAddress, DEFAULT_SIGNING_ALGORITHM);
    }

    /**
     * Constructor to create a new GATK key from scratch using an email address
     * and public/private key pair, and additionally specify the signing algorithm
     * to use. The private key is used for signing, and the public key is used to
     * validate the newly-created key.
     *
     * @param privateKey Private key used to sign the new GATK key
     * @param publicKey Public key used to validate the new GATK key
     * @param emailAddress The user's email address, which we will store in the key and sign
     * @param signingAlgorithm The combination of hash and encryption algorithms to use to sign the key
     */
    public GATKKey ( PrivateKey privateKey, PublicKey publicKey, String emailAddress, String signingAlgorithm ) {
        if ( privateKey == null || publicKey == null || emailAddress == null || emailAddress.length() == 0 || signingAlgorithm == null ) {
            throw new ReviewedGATKException("Cannot construct GATKKey using null/empty arguments");
        }

        this.privateKey = privateKey;
        this.publicKey = publicKey;
        this.emailAddress = emailAddress;
        this.signingAlgorithm = signingAlgorithm;

        validateEmailAddress();
        generateSignature();

        if ( ! isValid() ) {
            throw new ReviewedGATKException("Newly-generated GATK key fails validation -- this should never happen!");
        }
    }

    /**
     * Constructor to load an existing GATK key from a file.
     *
     * During loading, the key file is checked for integrity, but not cryptographic
     * validity (which must be done through a subsequent call to isValid()).
     *
     * @param publicKey Public key that will be used to validate the loaded GATK key
     *                  in subsequent calls to isValid()
     * @param keyFile File containing the GATK key to load
     */
    public GATKKey ( PublicKey publicKey, File keyFile ) {
        this(publicKey, keyFile, DEFAULT_SIGNING_ALGORITHM);
    }

    /**
     * Constructor to load an existing GATK key from a file, and additionally specify
     * the signing algorithm used to sign the key being loaded.
     *
     * During loading, the key file is checked for integrity, but not cryptographic
     * validity (which must be done through a subsequent call to isValid()).
     *
     * @param publicKey Public key that will be used to validate the loaded GATK key
     *                  in subsequent calls to isValid()
     * @param keyFile File containing the GATK key to load
     * @param signingAlgorithm The combination of hash and encryption algorithms used to sign the key
     */
    public GATKKey ( PublicKey publicKey, File keyFile, String signingAlgorithm ) {
        if ( publicKey == null || keyFile == null || signingAlgorithm == null ) {
            throw new ReviewedGATKException("Cannot construct GATKKey using null arguments");
        }

        this.publicKey = publicKey;
        this.signingAlgorithm = signingAlgorithm;

        readKey(keyFile);
    }

    // -----------------------
    // Public API Methods:
    // -----------------------

    /**
     * Writes out this key to a file in the format described at the top of this class,
     * encapsulating the key within a GZIP container.
     *
     * @param destination File to write the key to
     */
    public void writeKey ( File destination ) {
        try {
            byte[] keyBytes = marshalKeyData();
            IOUtils.writeByteArrayToStream(keyBytes, new GZIPOutputStream(new FileOutputStream(destination)));
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(destination, e);
        }
    }

    /**
     * Checks whether the signature of this key is cryptographically valid (ie., can be
     * decrypted by the public key to produce a valid SHA-1 hash of the email address
     * in the key).
     *
     * @return True if the key's signature passes validation, otherwise false
     */
    public boolean isValid() {
        try {
            Signature sig = Signature.getInstance(signingAlgorithm);
            sig.initVerify(publicKey);
            sig.update(emailAddress.getBytes());
            return sig.verify(signature);
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Signing algorithm %s not found", signingAlgorithm), e);
        }
        catch ( InvalidKeyException e ) {
            // If the GATK public key is invalid, it's likely our problem, not the user's:
            throw new ReviewedGATKException(String.format("Public key %s is invalid", publicKey), e);
        }
        catch ( SignatureException e ) {
            throw new UserException.UnreadableKeyException("Signature is invalid or signing algorithm was unable to process the input data", e);
        }
    }

    // -----------------------
    // Private Helper Methods:
    // -----------------------

    /**
     * Helper method that creates a signature for this key using the combination of
     * hash/encryption algorithms specified at construction time.
     */
    private void generateSignature() {
        try {
            Signature sig = Signature.getInstance(signingAlgorithm);
            sig.initSign(privateKey, CryptUtils.createRandomnessSource());
            sig.update(emailAddress.getBytes());
            signature = sig.sign();
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Signing algorithm %s not found", signingAlgorithm), e);
        }
        catch ( InvalidKeyException e ) {
            throw new ReviewedGATKException(String.format("Private key %s is invalid", privateKey), e);
        }
        catch ( SignatureException e ) {
            throw new ReviewedGATKException(String.format("Error creating signature for email address %s", emailAddress), e);
        }
    }

    /**
     * Helper method that reads in a GATK key from a file. Should not be called directly --
     * use the appropriate constructor above.
     *
     * @param source File to read the key from
     */
    private void readKey ( File source ) {
        try {
            byte[] keyBytes = IOUtils.readStreamIntoByteArray(new GZIPInputStream(new FileInputStream(source)));

            // As a sanity check, compare the number of bytes read to the uncompressed file size
            // stored in the GZIP ISIZE field. If they don't match, the key must be corrupt:
            if ( keyBytes.length != IOUtils.getGZIPFileUncompressedSize(source) ) {
                throw new UserException.UnreadableKeyException("Number of bytes read does not match the uncompressed size specified in the GZIP ISIZE field");
            }

            unmarshalKeyData(keyBytes);
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(source, e);
        }
        catch ( IOException e ) {
            throw new UserException.UnreadableKeyException(source, e);
        }
        catch ( UserException.CouldNotReadInputFile e ) {
            throw new UserException.UnreadableKeyException(source, e);
        }
    }

    /**
     * Helper method that assembles the email address and signature into a format
     * suitable for writing to disk.
     *
     * @return The aggregated key data, ready to be written to disk
     */
    private byte[] marshalKeyData() {
        byte[] emailAddressBytes = emailAddress.getBytes();
        byte[] assembledKey = new byte[emailAddressBytes.length + 1 + signature.length];

        System.arraycopy(emailAddressBytes, 0, assembledKey, 0, emailAddressBytes.length);
        assembledKey[emailAddressBytes.length] = GATK_KEY_SECTIONAL_DELIMITER;
        System.arraycopy(signature, 0, assembledKey, emailAddressBytes.length + 1, signature.length);

        return assembledKey;
    }

    /**
     * Helper method that parses the raw key data from disk into its component
     * email address and signature. Performs some basic validation in the process.
     *
     * @param keyBytes The raw, uncompressed key data read from disk
     */
    private void unmarshalKeyData ( byte[] keyBytes ) {
        int delimiterPosition = -1;

        for ( int i = 0; i < keyBytes.length; i++ ) {
            if ( keyBytes[i] == GATK_KEY_SECTIONAL_DELIMITER ) {
                delimiterPosition = i;
                break;
            }
        }

        if ( delimiterPosition == -1 ) {
            throw new UserException.UnreadableKeyException("Malformed GATK key contains no sectional delimiter");
        }
        else if ( delimiterPosition == 0 ) {
            throw new UserException.UnreadableKeyException("Malformed GATK key contains no email address");
        }
        else if ( delimiterPosition == keyBytes.length - 1 ) {
            throw new UserException.UnreadableKeyException("Malformed GATK key contains no signature");
        }

        byte[] emailAddressBytes = new byte[delimiterPosition];
        System.arraycopy(keyBytes, 0, emailAddressBytes, 0, delimiterPosition);
        emailAddress = new String(emailAddressBytes);

        signature = new byte[keyBytes.length - delimiterPosition - 1];
        System.arraycopy(keyBytes, delimiterPosition + 1, signature, 0, keyBytes.length - delimiterPosition - 1);
    }

    /**
     * Helper method that ensures that the user's email address does not contain the NUL byte, which we
     * reserve as a delimiter within each key file.
     */
    private void validateEmailAddress() {
        for ( byte b : emailAddress.getBytes() ) {
            if ( b == GATK_KEY_SECTIONAL_DELIMITER ) {
                throw new UserException(String.format("Email address must not contain a byte with value %d", GATK_KEY_SECTIONAL_DELIMITER));
            }
        }
    }
}
