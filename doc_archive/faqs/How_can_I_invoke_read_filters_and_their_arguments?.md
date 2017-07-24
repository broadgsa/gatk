## How can I invoke read filters and their arguments?

http://gatkforums.broadinstitute.org/gatk/discussion/2338/how-can-i-invoke-read-filters-and-their-arguments

<p>Most GATK tools apply several read filters by default. You can look up exactly what are the defaults for each tool in their respective <a href="http://www.broadinstitute.org/gatk/gatkdocs/">Technical Documentation</a> pages. </p>
<p>But sometimes you want to specify additional filters yourself (and before you ask, no, you cannot disable the default read filters used by a given tool). This is how you do it:</p>
<p>The <code>--read-filter</code> argument (or <code>-rf</code> for short) allows you to apply whatever read filters you'd like.  For example, to add the <code>MaxReadLengthFilter</code> filter above to <code>PrintReads</code>, you just add this to your command line:</p>
<pre><code class="pre_md">--read_filter MaxReadLength </code class="pre_md"></pre>
<h4>Notice that when you specify a read filter, you need to strip the Filter part of its name off!</h4>
<p>The read filter will be applied with its default value (which you can also look up in the Tech Docs for that filter). Now, if you want to specify a different value from the default, you pass the relevant argument by adding this right after the read filter:</p>
<pre><code class="pre_md">--read_filter MaxReadLength -maxReadLength 76</code class="pre_md"></pre>
<p>It's important that you pass the argument right after the filter itself, otherwise the command line parser won't know that they're supposed to go together.</p>
<p>And of course, you can add as many filters as you like by using multiple copies of the <code>--read_filter</code> parameter:</p>
<pre><code class="pre_md">--read_filter MaxReadLength --maxReadLength 76 --read_filter ZeroMappingQualityRead</code class="pre_md"></pre>