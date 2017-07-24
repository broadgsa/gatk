## What is the difference between QUAL and GQ annotations?

http://gatkforums.broadinstitute.org/gatk/discussion/4860/what-is-the-difference-between-qual-and-gq-annotations

<p>There has been a lot of confusion about the difference between QUAL and GQ, and we hope this FAQ will clarify the difference.</p>
<p>The basic difference is that QUAL refers to the variant site whereas GQ refers to a specific sample's GT. </p>
<ul>
<li>
<p>QUAL tells you how confident we are that there is some kind of variation at a given site. The variation may be present in one or more samples. </p>
</li>
<li>GQ tells you how confident we are that the genotype we assigned to a particular sample is correct. It is simply the second lowest PL, because it is the difference between the second lowest PL and the lowest PL (always 0).</li>
</ul>
<p>QUAL (or more importantly, its normalized form, QD) is mostly useful in multisample context. When you are recalibrating a cohort callset, you're going to be looking exclusively at site-level annotations like QD, because at that point what you're looking for is evidence of variation overall. That way you don't rely too much on individual sample calls, which are less robust.</p>
<p>In fact, many cohort studies don't even really care about individual genotype assignments, so they only use site annotations for their entire analysis.</p>
<p>Conversely, QUAL may seem redundant if you have only one sample. Especially if it has a good GQ (and more importantly, well separated PLs) then admittedly you don't really need to look at the QUAL -- you know what you have. If the GQ is not good, you can typically rely on the PLs to tell you whether you do probably have a variant, but we're just not sure if it's het or hom-var. If hom-ref is also a possibility, the call may be a potential false positive.</p>
<p>That said, it is more effective to filter on site-level annotations first, then refine and filter genotypes as appropriate. That's the workflow we recommend, based on years of experience doing this at fairly large scales...</p>