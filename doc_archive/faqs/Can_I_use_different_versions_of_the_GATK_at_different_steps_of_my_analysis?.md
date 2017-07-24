## Can I use different versions of the GATK at different steps of my analysis?

http://gatkforums.broadinstitute.org/gatk/discussion/3536/can-i-use-different-versions-of-the-gatk-at-different-steps-of-my-analysis

<p>Short answer: NO. </p>
<p>Medium answer: no, at least not if you want to run a low-risk pipeline.</p>
<p>Long answer: see below for details.</p>
<hr />
<p><strong>The rationale</strong></p>
<p>There are several reasons why you might want to do this: you're using the latest version of GATK and one of the tools has a show-stopping bug, so you'd like to use an older, pre-bug version of that tool, but still use the latest version of all the other tools; or maybe you've been using an older version of GATK and you'd like to use a new tool, but keep using the rest in the version that you've been using to process hundreds of samples already.</p>
<p><strong>The problem: compatibility is not guaranteed</strong></p>
<p>In many cases, when we modify one tool in the GATK, we need to make adjustments to other tools that interact either directly or indirectly with the data consumed or produced by the upgraded tool. If you mix and match tools from different versions of GATK, you risk running into compatibility issues. For example, HaplotypeCaller expects a BAM compressed by Reduce Reads to have its data annotated in a certain way. If the information is formatted differently than what the HC expects (because that's how the corresponding RR from the same version does it), it can blow up -- or worse, do the wrong thing but not tell you there's a problem.</p>
<p><strong>But what if the tools/tasks are in unrelated workflows?</strong></p>
<p>Would it really be so bad to use CountReads from GATK version 2.7 for a quick QC check that's not actually part of my pipeline, which uses version 2.5? Well, maaaaybe not, but we still think it's a source of error, and we do our damnedest to eliminate those.</p>
<p><strong>The conclusion</strong></p>
<p>You shouldn't use tools from different versions within the same workflow, that's for sure. We don't think it's worth the risks. If there's a show-stopping bug, let us know and we promise to fix it as soon as (humanly) possible. For the rest, either accept that you're stuck with the version you started your study with (we may be able to help with workarounds for known issues), or upgrade your entire workflow and start your analysis from scratch. Depending on how far along you are one of those options will be less painful to you; go with that. </p>
<p><strong>The plea bargain, and a warning</strong></p>
<p>If despite our dire warnings you're still going to mix and match tool versions, fine, we can't stop you. But be really careful, and check every version release notes document ever. And keep in mind that when things go wrong, we will deny you support if we think you've been reckless. </p>