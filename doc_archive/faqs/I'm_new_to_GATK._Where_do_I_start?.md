## I'm new to GATK. Where do I start?

http://gatkforums.broadinstitute.org/gatk/discussion/4863/im-new-to-gatk-where-do-i-start

<p>If this is your first rodeo, you're probably asking yourself:</p>
<ul>
<li>
<p><strong>What can GATK do for me?</strong>
Identify variants in a bunch of sample sequences, with great sensitivity and specificity.</p>
</li>
<li>
<p><strong>How do I get GATK to do that?</strong>
You run the recommended <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices</a> steps, one by one, from start to finish, as described in the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices documentation</a>.</p>
</li>
<li>
<p><strong>No but really, how do I know what to do?</strong>
For each step in the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices</a>, there is a tutorial that details how to run the tools involved, with example commands. The idea is to daisy-chain all thosee tutorials in the order that they're referenced in the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices</a> doc into a pipeline.</p>
</li>
<li>
<p><strong>Oh, you mean I can just copy/paste all the tutorial commands as they are?</strong>
Not quite, because there are a few things that need to be tweaked. For example, the tutorials use the <code>-L/--intervals</code> argument to restrict analysis for demo purposes, but depending on your data and experimental design, you may need to remove it (e.g. for WGS) or adapt it (for WEx). Hopefully it's explained clearly enough in the tutorials.</p>
</li>
<li>
<p><strong>Why don't you just provide one script that runs all the tools?</strong>
It's really hard to build and maintain a one-size-fits-all pipeline solution. Really really hard. And not nearly as much fun as developing new analysis methods. We do provide a pipelining program called Queue that has the advantage of understanding GATK argument syntax natively, but you still have to actually write scripts yourself in Scala to use it. Sorry. Maybe one day we will be able to offer GATK analysis on the Cloud. But not today. </p>
</li>
<li>
<p><strong>What if I want to know what a command line argument does or change a parameter?</strong>
First, check out the <a href="https://www.broadinstitute.org/gatk/guide/article?id=4669">basic GATK command syntax FAQ</a> if it's your first time using GATK, then consult the relevant <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/index">Tool Documentation</a> page. Keep in mind that some arguments are &quot;engine parameters&quot; that are shared by many tools, and are listed in a <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php">separate document</a>. Also, you can always use the search box to find an argument description really quickly. </p>
</li>
<li>
<p><strong>The documentation seems chaotic. Is there any logic to how it's organized?</strong>
Sort of. (And, ouch. Tough crowd.) The main category names should be obvious enough (if not, see the &quot;Documentation Categories&quot; tab). Within categories, everything is just in alphabetical order. In future, we're going to try to provide more use-case based structure, but for now this is what we have. The best way to find practical information is to either go from the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices</a> doc (which provide links to all FAQs, method articles and tutorials directly related to a given step), or use the search box and search-by-tag functions (see the &quot;Search tab&quot;). Be sure to also check out the <a href="https://www.broadinstitute.org/gatk/guide/presentations">Presentations section</a>, which provides workshop materials and videos that explain a lot of the motivation and methods behind the Best Practices. </p>
</li>
<li>
<p><strong>Does GATK include other tools beside the ones in the Best Practices?</strong>
Oh sure, there's a whole bunch of them, all listed in the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/index">Tool Documentation</a> section, categorized by type of analysis. But be aware that anything that's not part of the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices</a> is most likely either a tool that was written for a one-off analysis years ago, an experimental feature that we're still not sure is actually useful, or an accessory utility that can be used in many different ways and takes expert inside knowledge to use properly. All these may be buggy, insufficiently documented, or both. We provide support for them as well as humanly possible but ultimately, you use them at your own risk. </p>
</li>
<li>
<p><strong>Why do the answers to these questions keep getting longer and longer?</strong>
I don't know what you're talking about. </p>
</li>
<li><strong>What else should I know before I start?</strong>
You should probably browse the titles of the <a href="https://www.broadinstitute.org/gatk/guide/topic?name=faqs">Frequently Asked Questions</a> -- there will be at least a handful you'll want to read, but it's hard for us to predict which ones.</li>
</ul>