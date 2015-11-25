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

package org.broadinstitute.gatk.tools.walkers.annotator;

import org.broadinstitute.gatk.tools.walkers.annotator.SnpEff.EffectType;

import java.util.*;
/**
 * Created with IntelliJ IDEA.
 * User: farjoun
 * Date: 6/5/13
 * Time: 12:06 PM
 * To change this template use File | Settings | File Templates.
 */

/* This class holds a tree representation of the annotations used in snpEff, and provides a mechanism for telling if a
given annotation is a descendant of another.
The idea is to be able to stratify effects by large branches and not only the specific
snpEff annotation that a variant might have. For example if we want to know whether a variant is in CDS
but if it's marked SYNONYMOUS_CODING or NON_SYNONYMOUS_CODING (or many other options) still imply that its in the CDS.

The hierarchy was determined by Yossi Farjoun with input from Pablo (SNPEFF) and Tim Fennel.
*/


public class SnpEffUtil {

    // A map holding for every child, it's parent.
    // A node that isn't a key node is a root node.
    static private final Map<EffectType,EffectType> snpEffectGraph = new HashMap<>();

    //A map from each value of EffectType to a set of it's ancestors
    static private final Map<EffectType,Set<EffectType>> snpEffectAncestorSet = new HashMap<>();

    static {


        //INTERGENIC
        snpEffectGraph.put(EffectType.UPSTREAM,EffectType.INTERGENIC);
        snpEffectGraph.put(EffectType.DOWNSTREAM,EffectType.INTERGENIC);
        snpEffectGraph.put(EffectType.INTERGENIC_CONSERVED,EffectType.INTERGENIC);

        //INTRON
        snpEffectGraph.put(EffectType.INTRON_CONSERVED,EffectType.INTRON);
        snpEffectGraph.put(EffectType.SPLICE_SITE_ACCEPTOR,EffectType.INTRON);
        snpEffectGraph.put(EffectType.SPLICE_SITE_DONOR,EffectType.INTRON);

        //CDS
        snpEffectGraph.put(EffectType.EXON_DELETED,EffectType.CDS);
        snpEffectGraph.put(EffectType.SYNONYMOUS_CODING,EffectType.CDS);
        snpEffectGraph.put(EffectType.NON_SYNONYMOUS_CODING,EffectType.CDS);

        //SYNONYMOUS_CODING
        snpEffectGraph.put(EffectType.SYNONYMOUS_STOP,EffectType.SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.SYNONYMOUS_START,EffectType.SYNONYMOUS_CODING);

        //NON_SYNONYMOUS_CODING
        snpEffectGraph.put(EffectType.START_LOST,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.STOP_GAINED,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.STOP_LOST,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.CODON_CHANGE,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.CODON_INSERTION,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.CODON_DELETION,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.CODON_CHANGE_PLUS_CODON_DELETION,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.CODON_CHANGE_PLUS_CODON_INSERTION,EffectType.NON_SYNONYMOUS_CODING);
        snpEffectGraph.put(EffectType.FRAME_SHIFT,EffectType.NON_SYNONYMOUS_CODING);

        //UTRs
        snpEffectGraph.put(EffectType.UTR_5_DELETED,EffectType.UTR_5_PRIME);
        snpEffectGraph.put(EffectType.UTR_3_DELETED,EffectType.UTR_3_PRIME);
        snpEffectGraph.put(EffectType.START_GAINED,EffectType.UTR_5_PRIME);

        //EXON
        snpEffectGraph.put(EffectType.UTR_5_PRIME,EffectType.EXON);
        snpEffectGraph.put(EffectType.UTR_3_PRIME,EffectType.EXON);
        snpEffectGraph.put(EffectType.CDS,EffectType.EXON);


        //TRANSCRIPT
        snpEffectGraph.put(EffectType.INTRON,EffectType.TRANSCRIPT);
        snpEffectGraph.put(EffectType.EXON,EffectType.TRANSCRIPT);

        //GENE
        snpEffectGraph.put(EffectType.TRANSCRIPT,EffectType.GENE);
        snpEffectGraph.put(EffectType.REGULATION,EffectType.GENE);

        //CHROMOSOME
        snpEffectGraph.put(EffectType.GENE,EffectType.CHROMOSOME);
        snpEffectGraph.put(EffectType.INTERGENIC,EffectType.CHROMOSOME);
    }

    //A helper function that gets the parent set of the set of children
    private static Set<EffectType> getParentSet(final Set<EffectType> children){
        final Set<EffectType> parents=new HashSet<>();
        for(EffectType child:children){
            final EffectType parent = snpEffectGraph.get(child);
            if(parent!=null) parents.add(parent);
        }
        return parents;
    }

    //builds the total set of ancestors of a given node
    private static Set<EffectType> getAncestorSet(final EffectType child, final boolean isSelfIncluded){

        final Set<EffectType> ancestors=new HashSet<>();
        if(isSelfIncluded) ancestors.add(child);

        Set<EffectType> untraversedNodes=Collections.singleton(child);

        while(!untraversedNodes.isEmpty()){
            final Set<EffectType> putativeParents = getParentSet(untraversedNodes); //get immediate parents of unexamined set
            putativeParents.removeAll(ancestors); //remove all known parents, remaining with previously unknown parents
            ancestors.addAll(putativeParents); // add these parents to growing list of ancestors
            untraversedNodes=putativeParents; //still need to traverse parents of these nodes
        }
        return ancestors;
    }

    //returns true if the child effect is a subType of the parentEffect (including itself)
    public static boolean isSubTypeOf(final SnpEff.EffectType childEffect, final SnpEff.EffectType parentEffect){

        Set<EffectType> ancestorSet=snpEffectAncestorSet.get(childEffect);

        if(ancestorSet==null) {  //lazy population of map.
            ancestorSet = new HashSet<>();
            ancestorSet.addAll(getAncestorSet(childEffect, true)); //"true" so that a type is considered a subtype of itself
            snpEffectAncestorSet.put(childEffect, ancestorSet);
        }
        return ancestorSet.contains(parentEffect);
    }
}
