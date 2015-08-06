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
import org.broadinstitute.gatk.utils.io.IOUtils;

import javax.crypto.Cipher;
import java.io.File;
import java.io.InputStream;
import java.security.*;
import java.security.spec.InvalidKeySpecException;
import java.security.spec.KeySpec;
import java.security.spec.PKCS8EncodedKeySpec;
import java.security.spec.X509EncodedKeySpec;
import java.util.Arrays;

/**
 * A set of cryptographic utility methods and constants.
 *
 * Contains methods to:
 *
 * -Create a public/private key pair
 * -Read and write public/private keys to/from files/streams
 * -Load the GATK master private/public keys
 * -Encrypt/decrypt data
 *
 * Also contains constants that control the cryptographic defaults
 * throughout the GATK.
 *
 * @author David Roazen
 */
public class CryptUtils {

    // ---------------------------------------------------------------------------------
    // Constants (these control the default cryptographic settings throughout the GATK):
    // ---------------------------------------------------------------------------------

    /**
     * Default key length in bits of newly-created keys. 2048 bits provides a good balance between
     * security and speed.
     */
    public static final int DEFAULT_KEY_LENGTH = 2048;

    /**
     * Default encryption algorithm to use, when none is specified.
     */
    public static final String DEFAULT_ENCRYPTION_ALGORITHM = "RSA";

    /**
     * Default random-number generation algorithm to use, when none is specified.
     */
    public static final String DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM = "SHA1PRNG";

    /**
     * Name of the public key file distributed with the GATK. This file is packaged
     * into the GATK jar, and we use the system ClassLoader to find it.
     */
    public static final String GATK_DISTRIBUTED_PUBLIC_KEY_FILE_NAME = "GATK_public.key";

    /**
     * Location of the master copy of the GATK private key.
     */
    public static final String GATK_MASTER_PRIVATE_KEY_FILE = "/humgen/gsa-hpprojects/GATK/data/gatk_master_keys/GATK_private.key";

    /**
     * Location of the master copy of the GATK public key. This file should always be the same as
     * the public key file distributed with the GATK (and there are automated tests to ensure that it is).
     */
    public static final String GATK_MASTER_PUBLIC_KEY_FILE =  "/humgen/gsa-hpprojects/GATK/data/gatk_master_keys/GATK_public.key";

    /**
     * Directory where generated GATK user keys are stored. See the GATKKey class for more information.
     */
    public static final String GATK_USER_KEY_DIRECTORY =      "/humgen/gsa-hpprojects/GATK/data/gatk_user_keys/";


    // -----------------------
    // Utility Methods:
    // -----------------------

    /**
     * Generate a new public/private key pair using the default encryption settings defined above.
     *
     * @return A new public/private key pair created using the default settings
     */
    public static KeyPair generateKeyPair() {
        return generateKeyPair(DEFAULT_KEY_LENGTH, DEFAULT_ENCRYPTION_ALGORITHM, DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM);
    }

    /**
     * Generate a new public/private key pair using custom encryption settings.
     *
     * @param keyLength Length of the key in bits
     * @param encryptionAlgorithm Encryption algorithm to use
     * @param randNumberAlgorithm Random-number generation algorithm to use
     * @return A new public/private key pair, created according to the specified parameters
     */
    public static KeyPair generateKeyPair( int keyLength, String encryptionAlgorithm, String randNumberAlgorithm ) {
        try {
            KeyPairGenerator keyGen = KeyPairGenerator.getInstance(encryptionAlgorithm);
            SecureRandom randomnessSource = createRandomnessSource(randNumberAlgorithm);

            keyGen.initialize(keyLength, randomnessSource);
            return keyGen.generateKeyPair();
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Could not find an implementation of the requested encryption algorithm %s", encryptionAlgorithm), e);
        }
        catch ( Exception e ) {
            throw new ReviewedGATKException("Error while generating key pair", e);
        }
    }

    /**
     * Create a source of randomness using the default random-number generation algorithm.
     *
     * @return A randomness source that uses the default algorithm
     */
    public static SecureRandom createRandomnessSource() {
        return createRandomnessSource(DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM);
    }

    /**
     * Create a source of randomness using a custom random-number generation algorithm.
     *
     * @param randAlgorithm The random-number generation algorithm to use
     * @return A randomness sources that uses the specified algorithm
     */
    public static SecureRandom createRandomnessSource ( String randAlgorithm ) {
        try {
            return SecureRandom.getInstance(randAlgorithm);
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Could not find an implementation of the requested random-number generation algorithm %s", randAlgorithm), e);
        }
    }

    /**
     * Writes a public/private key pair to disk
     *
     * @param keyPair The key pair we're writing to disk
     * @param privateKeyFile Location to write the private key
     * @param publicKeyFile Location to write the public key
     */
    public static void writeKeyPair ( KeyPair keyPair, File privateKeyFile, File publicKeyFile ) {
        writeKey(keyPair.getPrivate(), privateKeyFile);
        writeKey(keyPair.getPublic(), publicKeyFile);
    }

    /**
     * Writes an arbitrary key to disk
     *
     * @param key The key to write
     * @param destination Location to write the key to
     */
    public static void writeKey ( Key key, File destination ) {
        IOUtils.writeByteArrayToFile(key.getEncoded(), destination);
    }

    /**
     * Reads in a public key created using the default encryption algorithm from a file.
     *
     * @param source File containing the public key
     * @return The public key read
     */
    public static PublicKey readPublicKey ( File source ) {
        return decodePublicKey(IOUtils.readFileIntoByteArray(source), DEFAULT_ENCRYPTION_ALGORITHM);
    }

    /**
     * Reads in a public key created using the default encryption algorithm from a stream.
     *
     * @param source Stream attached to the public key
     * @return The public key read
     */
    public static PublicKey readPublicKey ( InputStream source ) {
        return decodePublicKey(IOUtils.readStreamIntoByteArray(source), DEFAULT_ENCRYPTION_ALGORITHM);
    }

    /**
     * Decodes the raw bytes of a public key into a usable object.
     *
     * @param rawKey The encoded bytes of a public key as read from, eg., a file. The
     *               key must be in the standard X.509 format for a public key.
     * @param encryptionAlgorithm The encryption algorithm used to create the public key
     * @return The public key as a usable object
     */
    public static PublicKey decodePublicKey ( byte[] rawKey, String encryptionAlgorithm ) {
        try {
            KeySpec keySpec = new X509EncodedKeySpec(rawKey);
            KeyFactory keyFactory = KeyFactory.getInstance(encryptionAlgorithm);
            return keyFactory.generatePublic(keySpec);
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Could not find an implementation of the requested encryption algorithm %s", encryptionAlgorithm), e);
        }
        catch ( InvalidKeySpecException e ) {
            throw new ReviewedGATKException("Unable to use X.509 key specification to decode the given key", e);
        }
    }

    /**
     * Reads in a private key created using the default encryption algorithm from a file.
     *
     * @param source File containing the private key
     * @return The private key read
     */
    public static PrivateKey readPrivateKey ( File source ) {
        return decodePrivateKey(IOUtils.readFileIntoByteArray(source), DEFAULT_ENCRYPTION_ALGORITHM);
    }

    /**
     * Reads in a private key created using the default encryption algorithm from a stream.
     *
     * @param source Stream attached to the private key
     * @return The private key read
     */
    public static PrivateKey readPrivateKey ( InputStream source ) {
        return decodePrivateKey(IOUtils.readStreamIntoByteArray(source), DEFAULT_ENCRYPTION_ALGORITHM);
    }

    /**
     * Decodes the raw bytes of a private key into a usable object.
     *
     * @param rawKey The encoded bytes of a private key as read from, eg., a file. The
     *               key must be in the standard PKCS #8 format for a private key.
     * @param encryptionAlgorithm The encryption algorithm used to create the private key
     * @return The private key as a usable object
     */
    public static PrivateKey decodePrivateKey ( byte[] rawKey, String encryptionAlgorithm ) {
        try {
            KeySpec keySpec = new PKCS8EncodedKeySpec(rawKey);
            KeyFactory keyFactory = KeyFactory.getInstance(encryptionAlgorithm);
            return keyFactory.generatePrivate(keySpec);
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Could not find an implementation of the requested encryption algorithm %s", encryptionAlgorithm), e);
        }
        catch ( InvalidKeySpecException e ) {
            throw new ReviewedGATKException("Unable to use the PKCS #8 key specification to decode the given key", e);
        }
    }

    /**
     * Loads the copy of the GATK public key that is distributed with the GATK. Uses the system
     * ClassLoader to locate the public key file, which should be stored at the root of the GATK
     * jar file.
     *
     * @return The GATK public key as a usable object
     */
    public static PublicKey loadGATKDistributedPublicKey() {
        InputStream publicKeyInputStream = ClassLoader.getSystemResourceAsStream(GATK_DISTRIBUTED_PUBLIC_KEY_FILE_NAME);

        if ( publicKeyInputStream == null ) {
            throw new ReviewedGATKException(String.format("Could not locate the GATK public key %s in the classpath",
                                                           GATK_DISTRIBUTED_PUBLIC_KEY_FILE_NAME));
        }

        return readPublicKey(publicKeyInputStream);
    }

    /**
     * Loads the master copy of the GATK private key. You must have the appropriate UNIX permissions
     * to do this!
     *
     * @return The GATK master private key as a usable object
     */
    public static PrivateKey loadGATKMasterPrivateKey() {
        return readPrivateKey(new File(GATK_MASTER_PRIVATE_KEY_FILE));
    }

    /**
     * Loads the master copy of the GATK public key. This should always be the same as the
     * public key distributed with the GATK returned by loadGATKDistributedPublicKey().
     *
     * @return The GATK master public key as a usable object
     */
    public static PublicKey loadGATKMasterPublicKey() {
        return readPublicKey(new File(GATK_MASTER_PUBLIC_KEY_FILE));
    }

    /**
     * Encrypts the given data using the key provided.
     *
     * @param data The data to encrypt, as a byte array
     * @param encryptKey The key with which to encrypt the data
     * @return The encrypted version of the provided data
     */
    public static byte[] encryptData ( byte[] data, Key encryptKey ) {
        return transformDataUsingCipher(data, encryptKey, Cipher.ENCRYPT_MODE);
    }

    /**
     * Decrypts the given data using the key provided.
     *
     * @param encryptedData Data to decrypt, as a byte array
     * @param decryptKey The key with which to decrypt the data
     * @return The decrypted version of the provided data
     */
    public static byte[] decryptData ( byte[] encryptedData, Key decryptKey ) {
        return transformDataUsingCipher(encryptedData, decryptKey, Cipher.DECRYPT_MODE);
    }

    /**
     * Helper method for encryption/decryption that takes data and processes it using
     * the given key
     *
     * @param data Data to encrypt/decrypt
     * @param key Key to use to encrypt/decrypt the data
     * @param cipherMode Specifies whether we are encrypting or decrypting
     * @return The encrypted/decrypted data
     */
    private static byte[] transformDataUsingCipher ( byte[] data, Key key, int cipherMode ) {
        try {
            Cipher cipher = Cipher.getInstance(key.getAlgorithm());
            cipher.init(cipherMode, key);
            return cipher.doFinal(data);
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new ReviewedGATKException(String.format("Could not find an implementation of the requested algorithm %s",
                                             key.getAlgorithm()), e);
        }
        catch ( InvalidKeyException e ) {
            throw new ReviewedGATKException("Key is invalid", e);
        }
        catch ( GeneralSecurityException e ) {
            throw new ReviewedGATKException("Error during encryption", e);
        }
    }

    /**
     * Tests whether the public/private keys provided can each decrypt data encrypted by
     * the other key -- ie., tests whether these two keys are part of the same public/private
     * key pair.
     *
     * @param privateKey The private key to test
     * @param publicKey The public key to test
     * @return True if the keys are part of the same key pair and can decrypt each other's
     *         encrypted data, otherwise false.
     */
    public static boolean keysDecryptEachOther ( PrivateKey privateKey, PublicKey publicKey ) {
        byte[] plainText = "Test PlainText".getBytes();

        byte[] dataEncryptedUsingPrivateKey = CryptUtils.encryptData(plainText, privateKey);
        byte[] dataEncryptedUsingPublicKey = CryptUtils.encryptData(plainText, publicKey);

        byte[] privateKeyDataDecryptedWithPublicKey = CryptUtils.decryptData(dataEncryptedUsingPrivateKey, publicKey);
        byte[] publicKeyDataDecryptedWithPrivateKey = CryptUtils.decryptData(dataEncryptedUsingPublicKey, privateKey);

        // Make sure we actually transformed the data during encryption:
        if ( Arrays.equals(plainText, dataEncryptedUsingPrivateKey) ||
             Arrays.equals(plainText, dataEncryptedUsingPublicKey) ||
             Arrays.equals(dataEncryptedUsingPrivateKey, dataEncryptedUsingPublicKey) ) {
            return false;
        }

        // Make sure that we were able to recreate the original plaintext using
        // both the public key on the private-key-encrypted data and the private
        // key on the public-key-encrypted data:
        if ( ! Arrays.equals(plainText, privateKeyDataDecryptedWithPublicKey) ||
             ! Arrays.equals(plainText, publicKeyDataDecryptedWithPrivateKey) ) {
            return false;
        }

        return true;
    }
}
