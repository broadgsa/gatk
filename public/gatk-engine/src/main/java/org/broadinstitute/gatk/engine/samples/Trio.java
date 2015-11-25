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

package org.broadinstitute.gatk.engine.samples;

/**
 * A class for imposing a trio structure on three samples; a common paradigm
 *
 * todo -- there should probably be an interface or abstract class "Pedigree" that generalizes the notion of
 *      -- imposing structure on samples. But given how complex pedigrees can quickly become, it's not
 *      -- clear the best way to do this.
 */
public class Trio {
    private Sample mother;
    private Sample father;
    private Sample child;

    public Trio(Sample mom, Sample dad, Sample spawn) {
        assert mom.getID().equals(spawn.getMaternalID()) && dad.getID().equals(spawn.getPaternalID()) : "Samples passed to trio constructor do not form a trio";
        mother = mom;
        father = dad;
        child = spawn;
    }

    public Sample getMother() {
        return mother;
    }

    public String getMaternalID() {
        return mother.getID();
    }

    public Sample getFather() {
        return father;
    }

    public String getPaternalID() {
        return father.getID();
    }

    public Sample getChild() {
        return child;
    }

    public String getChildID() {
        return child.getID();
    }
}
