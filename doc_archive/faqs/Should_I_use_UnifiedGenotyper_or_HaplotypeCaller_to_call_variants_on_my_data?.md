## Should I use UnifiedGenotyper or HaplotypeCaller to call variants on my data?

http://gatkforums.broadinstitute.org/gatk/discussion/3151/should-i-use-unifiedgenotyper-or-haplotypecaller-to-call-variants-on-my-data

<p><strong>Use HaplotypeCaller!</strong></p>
<p>The HaplotypeCaller is a more recent and sophisticated tool than the UnifiedGenotyper. Its ability to call SNPs is equivalent to that of the UnifiedGenotyper, its ability to call indels is far superior, and it is now capable of calling non-diploid samples. It also comprises several unique functionalities such as the reference confidence model (which enables efficient and incremental variant discovery on ridiculously large cohorts) and special settings for RNAseq data. </p>
<p><strong>As of GATK version 3.3, we recommend using HaplotypeCaller in all cases, with no exceptions.</strong></p>
<p><em>Caveats for older versions</em></p>
<p>If you are limited to older versions for project continuity, you may opt to use UnifiedGenotyper in the following cases:</p>
<ul>
<li>If you are working with non-diploid organisms (UG can handle different levels of ploidy while older versions of HC cannot)  </li>
<li>If you are working with pooled samples (also due to the HC’s limitation regarding ploidy)  </li>
<li>If you want to analyze more than 100 samples at a time (for performance reasons) (versions 2.x) </li>
</ul>