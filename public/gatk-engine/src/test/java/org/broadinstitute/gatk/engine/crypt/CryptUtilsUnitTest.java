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

import org.broadinstitute.gatk.engine.crypt.CryptUtils;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.security.Key;
import java.security.KeyPair;
import java.security.PrivateKey;
import java.security.PublicKey;
import java.util.Arrays;

public class CryptUtilsUnitTest extends BaseTest {

    @Test
    public void testGenerateValidKeyPairWithDefaultSettings() {
        KeyPair keyPair = CryptUtils.generateKeyPair();
        Assert.assertTrue(CryptUtils.keysDecryptEachOther(keyPair.getPrivate(), keyPair.getPublic()));
    }

    @DataProvider( name = "InvalidKeyPairSettings" )
    public Object[][] invalidKeyPairSettingsDataProvider() {
        return new Object[][] {
            { -1, CryptUtils.DEFAULT_ENCRYPTION_ALGORITHM, CryptUtils.DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM},
            { CryptUtils.DEFAULT_KEY_LENGTH, "Made-up algorithm", CryptUtils.DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM},
            { CryptUtils.DEFAULT_KEY_LENGTH, CryptUtils.DEFAULT_ENCRYPTION_ALGORITHM, "Made-up algorithm"}
        };
    }

    @Test( dataProvider = "InvalidKeyPairSettings", expectedExceptions = ReviewedGATKException.class )
    public void testGenerateKeyPairWithInvalidSettings( int keyLength, String encryptionAlgorithm, String randomNumberGenerationAlgorithm ) {
        KeyPair keyPair = CryptUtils.generateKeyPair(keyLength, encryptionAlgorithm, randomNumberGenerationAlgorithm);
    }

    @Test
    public void testGATKMasterKeyPairMutualDecryption() {
        if ( gatkPrivateKeyExistsButReadPermissionDenied() ) {
            throw new SkipException(String.format("Skipping test %s because we do not have permission to read the GATK private key",
                                    "testGATKMasterKeyPairMutualDecryption"));
        }

        Assert.assertTrue(CryptUtils.keysDecryptEachOther(CryptUtils.loadGATKMasterPrivateKey(), CryptUtils.loadGATKMasterPublicKey()));
    }

    @Test
    public void testGATKMasterPrivateKeyWithDistributedPublicKeyMutualDecryption() {
        if ( gatkPrivateKeyExistsButReadPermissionDenied() ) {
            throw new SkipException(String.format("Skipping test %s because we do not have permission to read the GATK private key",
                                    "testGATKMasterPrivateKeyWithDistributedPublicKeyMutualDecryption"));
        }

        Assert.assertTrue(CryptUtils.keysDecryptEachOther(CryptUtils.loadGATKMasterPrivateKey(), CryptUtils.loadGATKDistributedPublicKey()));
    }

    @Test
    public void testKeyPairWriteThenRead() {
        KeyPair keyPair = CryptUtils.generateKeyPair();
        File privateKeyFile = createTempFile("testKeyPairWriteThenRead_private", "key");
        File publicKeyFile = createTempFile("testKeyPairWriteThenRead_public", "key");

        CryptUtils.writeKeyPair(keyPair, privateKeyFile, publicKeyFile);

        assertKeysAreEqual(keyPair.getPrivate(), CryptUtils.readPrivateKey(privateKeyFile));
        assertKeysAreEqual(keyPair.getPublic(), CryptUtils.readPublicKey(publicKeyFile));
    }

    @Test
    public void testPublicKeyWriteThenReadFromFile() {
        File keyFile = createTempFile("testPublicKeyWriteThenReadFromFile", "key");
        PublicKey publicKey = CryptUtils.generateKeyPair().getPublic();

        CryptUtils.writeKey(publicKey, keyFile);

        assertKeysAreEqual(publicKey, CryptUtils.readPublicKey(keyFile));
    }

    @Test
    public void testPublicKeyWriteThenReadFromStream() throws IOException {
        File keyFile = createTempFile("testPublicKeyWriteThenReadFromStream", "key");
        PublicKey publicKey = CryptUtils.generateKeyPair().getPublic();

        CryptUtils.writeKey(publicKey, keyFile);

        assertKeysAreEqual(publicKey, CryptUtils.readPublicKey(new FileInputStream(keyFile)));
    }

    @Test
    public void testPrivateKeyWriteThenReadFromFile() {
        File keyFile = createTempFile("testPrivateKeyWriteThenReadFromFile", "key");
        PrivateKey privateKey = CryptUtils.generateKeyPair().getPrivate();

        CryptUtils.writeKey(privateKey, keyFile);

        assertKeysAreEqual(privateKey, CryptUtils.readPrivateKey(keyFile));
    }

    @Test
    public void testPrivateKeyWriteThenReadFromStream() throws IOException {
        File keyFile = createTempFile("testPrivateKeyWriteThenReadFromStream", "key");
        PrivateKey privateKey = CryptUtils.generateKeyPair().getPrivate();

        CryptUtils.writeKey(privateKey, keyFile);

        assertKeysAreEqual(privateKey, CryptUtils.readPrivateKey(new FileInputStream(keyFile)));
    }

    @Test( expectedExceptions = UserException.CouldNotReadInputFile.class )
    public void testReadNonExistentPublicKey() {
        File nonExistentFile = new File("jdshgkdfhg.key");
        Assert.assertFalse(nonExistentFile.exists());

        CryptUtils.readPublicKey(nonExistentFile);
    }

    @Test( expectedExceptions = UserException.CouldNotReadInputFile.class )
    public void testReadNonExistentPrivateKey() {
        File nonExistentFile = new File("jdshgkdfhg.key");
        Assert.assertFalse(nonExistentFile.exists());

        CryptUtils.readPrivateKey(nonExistentFile);
    }

    @Test
    public void testDecodePublicKey() {
        PublicKey originalKey = CryptUtils.generateKeyPair().getPublic();
        PublicKey decodedKey = CryptUtils.decodePublicKey(originalKey.getEncoded(), CryptUtils.DEFAULT_ENCRYPTION_ALGORITHM);
        assertKeysAreEqual(originalKey, decodedKey);
    }

    @Test
    public void testDecodePrivateKey() {
        PrivateKey originalKey = CryptUtils.generateKeyPair().getPrivate();
        PrivateKey decodedKey = CryptUtils.decodePrivateKey(originalKey.getEncoded(), CryptUtils.DEFAULT_ENCRYPTION_ALGORITHM);
        assertKeysAreEqual(originalKey, decodedKey);
    }

    @Test
    public void testLoadGATKMasterPrivateKey() {
        if ( gatkPrivateKeyExistsButReadPermissionDenied() ) {
            throw new SkipException(String.format("Skipping test %s because we do not have permission to read the GATK private key",
                                    "testLoadGATKMasterPrivateKey"));
        }

        PrivateKey gatkMasterPrivateKey = CryptUtils.loadGATKMasterPrivateKey();
    }

    @Test
    public void testLoadGATKMasterPublicKey() {
        PublicKey gatkMasterPublicKey = CryptUtils.loadGATKMasterPublicKey();
    }

    @Test
    public void testLoadGATKDistributedPublicKey() {
        PublicKey gatkDistributedPublicKey = CryptUtils.loadGATKDistributedPublicKey();
    }

    private void assertKeysAreEqual( Key originalKey, Key keyFromDisk ) {
        Assert.assertTrue(Arrays.equals(originalKey.getEncoded(), keyFromDisk.getEncoded()));
        Assert.assertEquals(originalKey.getAlgorithm(), keyFromDisk.getAlgorithm());
        Assert.assertEquals(originalKey.getFormat(), keyFromDisk.getFormat());
    }

    private boolean gatkPrivateKeyExistsButReadPermissionDenied() {
        File gatkPrivateKey = new File(CryptUtils.GATK_MASTER_PRIVATE_KEY_FILE);
        return gatkPrivateKey.exists() && ! gatkPrivateKey.canRead();
    }
}
