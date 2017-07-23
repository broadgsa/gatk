## Seeing deletion spanning reads in LocusWalkers

http://gatkforums.broadinstitute.org/gatk/discussion/1348/seeing-deletion-spanning-reads-in-locuswalkers

<h2>1. Introduction</h2>
<p>The <code>LocusTraversal</code> now supports passing walkers reads that have deletions spanning the current locus.  This is useful in many situation where you want to calculate coverage, call variants and need to avoid calling variants where there are a lot of deletions, etc.  </p>
<p>Currently, the system by default will not pass you deletion-spanning reads.  In order to see them, you need to overload the function:</p>
<pre><code class="pre_md">/**
 * (conceptual static) method that states whether you want to see reads piling up at a locus
 * that contain a deletion at the locus.
 *
 * ref:   ATCTGA
 * read1: ATCTGA
 * read2: AT--GA
 *
 * Normally, the locus iterator only returns a list of read1 at this locus at position 3, but
 * if this function returns true, then the system will return (read1, read2) with offsets
 * of (3, -1).  The -1 offset indicates a deletion in the read.
 *
 * @return false if you don't want to see deletions, or true if you do
 */
public boolean includeReadsWithDeletionAtLoci() { return true; }</code class="pre_md"></pre>
<p>in your walker.  Now you will start seeing deletion-spanning reads in your walker.  These reads are flagged with offsets of -1, so that you can:</p>
<pre><code class="pre_md">    for ( int i = 0; i &lt; context.getReads().size(); i++ ) {
        SAMRecord read = context.getReads().get(i);
        int offset = context.getOffsets().get(i);

       if ( offset == -1 ) 
               nDeletionReads++;
        else 
               nCleanReads++;
    }</code class="pre_md"></pre>
<p>There are also two convenience functions in <code>AlignmentContext</code> to extract subsets of the reads with and without spanning deletions:</p>
<pre><code class="pre_md">/**
 * Returns only the reads in ac that do not contain spanning deletions of this locus
 * 
 * @param ac
 * @return
 */
public static AlignmentContext withoutSpanningDeletions( AlignmentContext ac );

/**
 * Returns only the reads in ac that do contain spanning deletions of this locus
 * 
 * @param ac
 * @return
 */
public static AlignmentContext withSpanningDeletions( AlignmentContext ac );</code class="pre_md"></pre>