package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.wiggle.WiggleHeader;
import org.broadinstitute.sting.utils.wiggle.WiggleWriter;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashSet;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Jul 21, 2010
 */
public class CreateTiTvTrack extends RodWalker<VariantContext,TiTvWindow> {
    @Argument(shortName="size",doc="Size of the window",required = true)
    int size = -1;

    private WiggleWriter writer;

    public TiTvWindow reduceInit() {
        writer = new WiggleWriter(out);
        writer.writeHeader(new WiggleHeader("TiTv",String.format("The Transition Transversion rate across the genome using variants windows of size %d",size)));
        return new TiTvWindow(size);
    }

   // public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
    //    return tracker == null ||  tracker.getVariantContext(ref, "variant", null, context.getLocation(), true).isFiltered();
    //}

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc;
        if ( tracker != null ) {
            vc = tracker.getVariantContext(ref, "variant", null, context.getLocation(), true);
            return vc;
        } else {
            return null;
        }
    }

    public TiTvWindow reduce(VariantContext vc, TiTvWindow window) {
        if ( vc == null || ! vc.isSNP() || vc.getAlternateAlleles().size() > 1) {
            return window;
        }
        
        window.update(vc.isTransition());
        if ( window.getTiTv() != null ) {
            writer.writeData(vc.getLocation(),window.getTiTv());
        }

        return window;
    }

    public void onTraversalDone() {
        
    }

}

class TiTvWindow {
    long nTi;
    long nTv;
    ArrayList<Boolean> variants;
    int maxSize;

    public TiTvWindow(int size) {
        maxSize = size;
        variants = new ArrayList<Boolean>(size);
    }

    public void update(Boolean isTi) {
        if ( variants.size() == maxSize ) {
            Boolean first = variants.remove(0);
            if ( first ) {
                nTi--;
            } else {
                nTv--;
            }
        }

        variants.add(isTi);
        if ( isTi ) {
            nTi++;
        } else {
            nTv++;
        }

        //System.out.println(variants.size());
    }

    public Double getTiTv() {
        if ( variants.size() == maxSize ) {
            return ( nTi + 1.0 )/(nTv + 1.0);
        } else {
            return null; // window not full
        }
    }
}
