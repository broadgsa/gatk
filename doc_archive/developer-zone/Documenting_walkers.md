## Documenting walkers

http://gatkforums.broadinstitute.org/gatk/discussion/1346/documenting-walkers

<p>The GATK discovers walker documentation by reading it out of the Javadoc, Sun's design pattern for providing documentation for packages and classes.  This page will provide an extremely brief explanation of how to write Javadoc; more information on writing javadoc comments can be found in <a href="http://www.oracle.com/technetwork/java/javase/documentation/index-137868.html">Sun's documentation</a>.</p>
<h2>1. Adding walker and package descriptions to the help text</h2>
<p>The GATK's build system uses the javadoc parser to extract the javadoc for classes and packages and embed the contents of that javadoc in the help system.  If you add Javadoc to your package or walker, it will automatically appear in the help.  The javadoc parser will pick up on 'standard' javadoc comments, such as the following, taken from PrintReadsWalker:</p>
<pre><code class="pre_md">/**
 * This walker prints out the input reads in SAM format.  Alternatively, the walker can write reads into a specified BAM file.
 */</code class="pre_md"></pre>
<p>You can add javadoc to your package by creating a special file, <code>package-info.java</code>, in the package directory.  This file should consist of the javadoc for your package plus a package descriptor line.  One such example follows:</p>
<pre><code class="pre_md">/**
 * @help.display.name Miscellaneous walkers (experimental)
 */
package org.broadinstitute.sting.playground.gatk.walkers;</code class="pre_md"></pre>
<p>Additionally, the GATK provides two extra custom tags for overriding the information that ultimately makes it into the help.</p>
<ul>
<li>
<p><code>@help.display.name</code> Changes the name of the package as it appears in help.  Note that the name of the walker cannot be changed as it is required to be passed verbatim to the <code>-T</code> argument.</p>
</li>
<li>
<p><code>@help.summary</code> Changes the description which appears on the right-hand column of the help text.  This is useful if you'd like to provide a more concise description of the walker that should appear in the help.</p>
</li>
<li><code>@help.description</code> Changes the description which appears at the bottom of the help text with <code>-T &lt;your walker&gt; --help</code> is specified.  This is useful if you'd like to present a more complete description of your walker.</li>
</ul>
<h2>2. Hiding experimental walkers (use sparingly, please!)</h2>
<p>Walkers can be hidden from the documentation system by adding the <code>@Hidden</code> annotation to the top of each walker.  <code>@Hidden</code> walkers can still be run from the command-line, but their documentation will not be visible to end users.  Please use this functionality sparingly to avoid walkers with hidden command-line options that are required for production use.</p>
<h2>3. Disabling building of help</h2>
<p>Because the building of our help text is actually heavyweight and can dramatically increase compile time on some systems, we have a mechanism to disable help generation.</p>
<p>Compile with the following command:</p>
<pre><code class="pre_md">ant -Ddisable.help=true</code class="pre_md"></pre>
<p>to disable generation of help.</p>