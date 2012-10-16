package org.broadinstitute.sting.gatk.samples;

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
