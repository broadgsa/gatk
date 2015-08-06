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

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.SkipException;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.security.KeyPair;
import java.security.PrivateKey;
import java.security.PublicKey;

public class GATKKeyUnitTest extends BaseTest {

    @Test
    public void testCreateGATKKeyUsingMasterKeyPair() {
        if ( gatkPrivateKeyExistsButReadPermissionDenied() ) {
            throw new SkipException(String.format("Skipping test %s because we do not have permission to read the GATK private key",
                                    "testCreateGATKKeyUsingMasterKeyPair"));
        }

        PrivateKey masterPrivateKey = CryptUtils.loadGATKMasterPrivateKey();
        PublicKey masterPublicKey = CryptUtils.loadGATKMasterPublicKey();

        // We should be able to create a valid GATKKey using our master key pair:
        GATKKey key = new GATKKey(masterPrivateKey, masterPublicKey, "foo@bar.com");
        Assert.assertTrue(key.isValid());
    }

    @Test
    public void testCreateGATKKeyUsingMasterPrivateKeyAndDistributedPublicKey() {
        if ( gatkPrivateKeyExistsButReadPermissionDenied() ) {
            throw new SkipException(String.format("Skipping test %s because we do not have permission to read the GATK private key",
                                    "testCreateGATKKeyUsingMasterPrivateKeyAndDistributedPublicKey"));
        }

        PrivateKey masterPrivateKey = CryptUtils.loadGATKMasterPrivateKey();
        PublicKey distributedPublicKey = CryptUtils.loadGATKDistributedPublicKey();

        // We should also be able to create a valid GATKKey using our master private
        // key and the public key we distribute with the GATK:
        GATKKey key = new GATKKey(masterPrivateKey, distributedPublicKey, "foo@bar.com");
        Assert.assertTrue(key.isValid());
    }

    @Test( expectedExceptions = ReviewedGATKException.class )
    public void testKeyPairMismatch() {
        KeyPair firstKeyPair = CryptUtils.generateKeyPair();
        KeyPair secondKeyPair = CryptUtils.generateKeyPair();

        // Attempting to create a GATK Key with private and public keys that aren't part of the
        // same key pair should immediately trigger a validation failure:
        GATKKey key = new GATKKey(firstKeyPair.getPrivate(), secondKeyPair.getPublic(), "foo@bar.com");
    }

    @Test( expectedExceptions = ReviewedGATKException.class )
    public void testEncryptionAlgorithmMismatch() {
        KeyPair keyPair = CryptUtils.generateKeyPair(CryptUtils.DEFAULT_KEY_LENGTH, "DSA", CryptUtils.DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM);

        // Attempting to use a DSA private key to create an RSA signature should throw an error:
        GATKKey key = new GATKKey(keyPair.getPrivate(), keyPair.getPublic(), "foo@bar.com", "SHA1withRSA");
    }

    @Test( expectedExceptions = UserException.class )
    public void testInvalidEmailAddress() {
        String emailAddressWithNulByte = new String(new byte[] { 0 });
        KeyPair keyPair = CryptUtils.generateKeyPair();

        // Email addresses cannot contain the NUL byte, since it's used as a sectional delimiter in the key file:
        GATKKey key = new GATKKey(keyPair.getPrivate(), keyPair.getPublic(), emailAddressWithNulByte);
    }

    @Test
    public void testCreateGATKKeyFromValidKeyFile() {
        GATKKey key = new GATKKey(CryptUtils.loadGATKDistributedPublicKey(), new File(keysDataLocation + "valid.key"));
        Assert.assertTrue(key.isValid());
    }

    @Test( expectedExceptions = UserException.UnreadableKeyException.class )
    public void testCreateGATKKeyFromCorruptKeyFile() {
        GATKKey key = new GATKKey(CryptUtils.loadGATKDistributedPublicKey(), new File(keysDataLocation + "corrupt_random_contents.key"));
    }

    @Test
    public void testCreateGATKKeyFromRevokedKeyFile() {
        GATKKey key = new GATKKey(CryptUtils.loadGATKDistributedPublicKey(), new File(keysDataLocation + "revoked.key"));
        Assert.assertFalse(key.isValid());
    }

    @Test( expectedExceptions = UserException.CouldNotReadInputFile.class )
    public void testCreateGATKKeyFromNonExistentFile() {
        File nonExistentFile = new File("ghfdkgsdhg.key");
        Assert.assertFalse(nonExistentFile.exists());

        GATKKey key = new GATKKey(CryptUtils.loadGATKDistributedPublicKey(), nonExistentFile);
    }

    private boolean gatkPrivateKeyExistsButReadPermissionDenied() {
        File gatkPrivateKey = new File(CryptUtils.GATK_MASTER_PRIVATE_KEY_FILE);
        return gatkPrivateKey.exists() && ! gatkPrivateKey.canRead();
    }
}
