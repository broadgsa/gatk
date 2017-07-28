## Sting to GATK renaming

http://gatkforums.broadinstitute.org/gatk/discussion/4173/sting-to-gatk-renaming

<h1>Overview</h1>
<p>The GATK 3.2 source code uses new java package names, directory paths, and executable jars. Post GATK 3.2, any patches submitted via pull requests should also include classes moved to the appropriate artifact.</p>
<p>Note that the document includes references to the <code>private</code> module, which is part of our internal development codebase but is not available to the general public. </p>
<h1>Summary</h1>
<p>A long term ideal of the GATK is to separate out reusable parts and eventually make them available as compiled libraries via centralized binary repositories. Ahead of publishing a number of steps must be completed. One of the larger steps has been completed for GATK 3.2, where the code base rebranded all references of Sting to GATK.</p>
<p>Currently implemented changes include:</p>
<ul>
<li>Java/Scala package names changed from org.broadinstitute.sting to org.broadinstitute.gatk</li>
<li>Renamed Maven artifacts including new directories</li>
</ul>
<p>As of May 16, 2014, remaining TODOs ahead of publishing to central include:</p>
<ul>
<li>Uploading all transitive GATK dependencies to central repositories</li>
<li>Separating a bit more of the intertwined utility, engine, and tool classes</li>
</ul>
<p>Now that the new package names and Maven artifacts are available, any pull request should include ensuring that updated classes are also moved into the correct GATK Maven artifact. While there are a significant number of classes, cleaning up as we go along will allow the larger task to be completed in a distributed fashion.</p>
<p>The full lists of new Maven artifacts and renamed packages are below under [Renamed Artifact Directories]. For those developers in the middle of a <code>git rebase</code> around commits before and after 3.2, here is an abridged mapping of renamed directories for those trying to locate files:</p>
<table class="table table-striped">
<thead>
<tr>
<th>Old Maven Artifact</th>
<th>New Maven Artifact</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>public/sting-root</code></td>
<td><code>public/gatk-root</code></td>
</tr>
<tr>
<td><code>public/sting-utils</code></td>
<td><code>public/gatk-utils</code></td>
</tr>
<tr>
<td><code>public/gatk-framework</code></td>
<td><code>public/gatk-tools-public</code></td>
</tr>
<tr>
<td><code>public/queue-framework</code></td>
<td><code>public/gatk-queue</code></td>
</tr>
<tr>
<td><code>protected/gatk-protected</code></td>
<td><code>protected/gatk-tools-protected</code></td>
</tr>
<tr>
<td><code>private/gatk-private</code></td>
<td><code>private/gatk-tools-private</code></td>
</tr>
<tr>
<td><code>private/queue-private</code></td>
<td><code>private/gatk-queue-private</code></td>
</tr>
</tbody>
</table>
<p>QScripts are no longer located with the Queue engine, and instead are now located with the GATK wrappers implemented as Queue extensions. See [Separated Queue Extensions] for more info.</p>
<h1>Changes</h1>
<h2>Separating the GATK Engine and Tools</h2>
<p>Starting with GATK 3.2, separate Maven utility artifacts exist to separate reusable portions of the GATK engine apart from tool specific implementations. The biggest impact this will have on developers is the separation of the walkers packages.</p>
<p>In GATK versions &lt;= 3.1 there was one package for both the base classes and the implementations of walkers:</p>
<ul>
<li>org.broadinstitute.sting.gatk.walkers</li>
</ul>
<p>In GATK versions &gt;= 3.2 threre are two packages. The first contains the base interfaces, annotations, etc. The latter package is for the concrete tools implemented as walkers:</p>
<ul>
<li>
<p>org.broadinstitute.<strong>gatk.engine</strong>.walkers</p>
<ul>
<li>Ex: ReadWalker, LocusWalker, @PartitionBy, @Requires, etc.</li>
</ul>
</li>
<li>org.broadinstitute.<strong>gatk.tools</strong>.walkers
<ul>
<li>Ex: PrintReads, VariantEval, IndelRealigner, HaplotypeCaller, etc.</li>
</ul></li>
</ul>
<h2>Renamed Binary Packages</h2>
<p>Previously, depending on how the source code was compiled, the executable gatk-package-3.1.jar and queue-package-3.1.jar (aka GenomeAnalysisTK.jar and Queue.jar) contained various mixes of public/protected/private code. For example, if the private directory was present when the source code was compiled, the same artifact named gatk-package-3.1.jar might, or might not contain private code.</p>
<p>Starting with 3.2, there are two versions of the jar created, each with specific file contents.</p>
<table class="table table-striped">
<thead>
<tr>
<th>New Maven Artifact</th>
<th>Alias in the /target folder</th>
<th>Packaged contents</th>
</tr>
</thead>
<tbody>
<tr>
<td>gatk-package-distribution-3.2.jar</td>
<td>GenomeAnalysisTK.jar</td>
<td>public,protected</td>
</tr>
<tr>
<td>gatk-package-internal-3.2.jar</td>
<td>GenomeAnalysisTK-internal.jar</td>
<td>public,protected,private</td>
</tr>
<tr>
<td>gatk-queue-package-distribution-3.2.jar</td>
<td>Queue.jar</td>
<td>public,protected</td>
</tr>
<tr>
<td>gatk-queue-package-internal-3.2.jar</td>
<td>Queue-internal.jar</td>
<td>public,protected,private</td>
</tr>
</tbody>
</table>
<h2>Separated Queue Extensions</h2>
<p>When creating a packaged version of Queue, the GATKExtensionsGenerator builds Queue engine compatible command line wrappers around each GATK walker. Previously, the wrappers were generated during the compilation of the Queue framework. Similar to the binary packages, depending on who built the source code, queue-framework-3.1.jar would contain various mixes of public/protected/private wrappers.</p>
<p>Starting with GATK 3.2, the gatk-queue-3.2.jar only contains code for the Queue engine. Generated and manually created extensions for wrapping any other command line programs are all included in separate artifacts. Due to a current limitation regarding how the generator uses reflection, the generator cannot build wrappers for just private classes without also generating protected and public classes. Thus, there are three different Maven artifacts generated, that contain different mixes of public, protected and private wrappers.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Extensions Artifact</th>
<th>Generated wrappers for GATK tools</th>
</tr>
</thead>
<tbody>
<tr>
<td>gatk-queue-extensions-public-3.2.jar</td>
<td>public <em>only</em></td>
</tr>
<tr>
<td>gatk-queue-extensions-distribution-3.2.jar</td>
<td>public,protected</td>
</tr>
<tr>
<td>gatk-queue-extensions-internal-3.2.jar</td>
<td>public,protected,private</td>
</tr>
</tbody>
</table>
<p>As for QScripts that used to be located with the framework, they are now located with the generated wrappers.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Old QScripts Artifact Directory</th>
<th>New QScripts Artifact Directory</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>public/queue-framework/src/main/qscripts</code></td>
<td><code>public/gatk-queue-extensions-public/src/main/qscripts</code></td>
</tr>
<tr>
<td><code>private/queue-private/src/main/qscripts</code></td>
<td><code>private/gatk-queue-extensions-internal/src/main/qscripts</code></td>
</tr>
</tbody>
</table>
<h2>Renamed Artifact Directories</h2>
<p>The following list shows the mapping of artifact names pre and post GATK 3.2. In addition to the engine changes, the packaging updates and extensions changes above also affected Maven artifact refactoring. The packaging artifacts have split from a single public to protected and private versions, and new queue extensions artifacts have been added as well.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Maven Artifact &lt;= GATK 3.1</th>
<th>Maven Artifact &gt;= GATK 3.2</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>/pom.xml</code> <em>(sting-aggregator)</em></td>
<td><code>/pom.xml</code> _(gatk<em>aggregator)</em></td>
</tr>
<tr>
<td><code>public/sting-root</code></td>
<td><code>public/gatk-root</code></td>
</tr>
<tr>
<td><code>public/sting-utils</code></td>
<td><code>public/gatk-utils</code></td>
</tr>
<tr>
<td><em>none</em></td>
<td><code>public/gatk-engine</code></td>
</tr>
<tr>
<td><code>public/gatk-framework</code></td>
<td><code>public/gatk-tools-public</code></td>
</tr>
<tr>
<td><code>public/queue-framework</code></td>
<td><code>public/gatk-queue</code></td>
</tr>
<tr>
<td><code>public/gatk-queue-extgen</code></td>
<td><code>public/gatk-queue-extensions-generator</code></td>
</tr>
<tr>
<td><code>protected/gatk-protected</code></td>
<td><code>protected/gatk-tools-protected</code></td>
</tr>
<tr>
<td><code>private/gatk-private</code></td>
<td><code>private/gatk-tools-private</code></td>
</tr>
<tr>
<td><code>private/queue-private</code></td>
<td><code>private/gatk-queue-private</code></td>
</tr>
<tr>
<td><code>public/gatk-package</code></td>
<td><code>protected/gatk-package-distribution</code></td>
</tr>
<tr>
<td><code>public/queue-package</code></td>
<td><code>protected/gatk-queue-package-distribution</code></td>
</tr>
<tr>
<td><em>none</em></td>
<td><code>private/gatk-package-internal</code></td>
</tr>
<tr>
<td><em>none</em></td>
<td><code>private/gatk-queue-package-internal</code></td>
</tr>
<tr>
<td><em>none</em></td>
<td><code>public/gatk-queue-extensions-public</code></td>
</tr>
<tr>
<td><em>none</em></td>
<td><code>protected/gatk-queue-extensions-distribution</code></td>
</tr>
<tr>
<td><em>none</em></td>
<td><code>private/gatk-queue-extensions-internal</code></td>
</tr>
</tbody>
</table>
<p><em>A note regarding the aggregator:</em></p>
<p>The aggregator is the pom.xml in the top directory level of the GATK source code. When someone clones the GATK source code and runs <code>mvn</code> in the top level directory, the aggregator the pom.xml executed.</p>
<p>The root is a pom.xml that contains all common Maven configuration. There are a couple dependent pom.xml files that inherit configuration from the root, but are <em>NOT</em> aggregated during normal source compilation.</p>
<p>As of GATK 3.2, these un-aggregated child artifacts are VectorPairHMM and picard-maven. They should not run by default with each instance of <code>mvn</code> run on the GATK source code.</p>
<p>For more clarification on Maven Inheritance vs. Aggregation, see the Maven <a href="http://maven.apache.org/guides/introduction/introduction-to-the-pom.html#Project_Inheritance_vs_Project_Aggregation">introduction to the pom</a>.</p>
<h2>Renamed Java/Scala Package Names</h2>
<p>In GATK 3.2, except for classes with Sting in the name, all file names are still the same. To locate migrated files under new java package names, developers should either use <a href="http://www.jetbrains.com/idea/webhelp/navigating-to-class-file-or-symbol-by-name.html">Intellij IDEA Navigation</a> or <code>/bin/find</code> to locate the same file they used previously.</p>
<p>The biggest change most developers will face is the new package names for GATK classes. Code entanglement does not permit simply moving the classes into the correct Maven artifacts, as a few number of lines of code must be edited inside a large number of files. So post renaming only a very small number of classes were moved out of the incorrect Maven artifacts as examples.</p>
<p>As of the May 16, 2014, the migrated GATK package distribution is as follows. This list includes only main classes. The table excludes all tests, renamed files such as StingException, certain private Queue wrappers, and qscripts renamed to end in *.scala.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Scope</th>
<th>Type</th>
<th>&lt;= 3.1 Artifact</th>
<th>&lt;= 3.1 Package</th>
<th>&gt;= GATK 3.2 Artifact</th>
<th>&gt;= 3.2 GATK Package</th>
<th style="text-align: right;">Files</th>
</tr>
</thead>
<tbody>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s</td>
<td>gatk-utils</td>
<td>o.b.g</td>
<td style="text-align: right;">4</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s.gatk</td>
<td>gatk-engine</td>
<td>o.b.g.engine</td>
<td style="text-align: right;">2</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s</td>
<td>gatk-tools-public</td>
<td>o.b.g</td>
<td style="text-align: right;">202</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s</td>
<td>gatk-tools-public</td>
<td>o.b.g.utils</td>
<td style="text-align: right;">49</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s</td>
<td>gatk-tools-public</td>
<td>o.b.g.engine</td>
<td style="text-align: right;">34</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-public</td>
<td>o.b.g.engine</td>
<td style="text-align: right;">244</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-public</td>
<td>o.b.g.tools</td>
<td style="text-align: right;">134</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-framework</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-public</td>
<td>o.b.g.tools.walkers</td>
<td style="text-align: right;">2</td>
</tr>
<tr>
<td>protected</td>
<td>java</td>
<td>gatk-protected</td>
<td>o.b.s</td>
<td>gatk-tools-protected</td>
<td>o.b.g</td>
<td style="text-align: right;">44</td>
</tr>
<tr>
<td>protected</td>
<td>java</td>
<td>gatk-protected</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-protected</td>
<td>o.b.g.engine</td>
<td style="text-align: right;">1</td>
</tr>
<tr>
<td>protected</td>
<td>java</td>
<td>gatk-protected</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-protected</td>
<td>o.b.g.tools</td>
<td style="text-align: right;">209</td>
</tr>
<tr>
<td>private</td>
<td>java</td>
<td>gatk-private</td>
<td>o.b.s</td>
<td>gatk-tools-private</td>
<td>o.b.g</td>
<td style="text-align: right;">23</td>
</tr>
<tr>
<td>private</td>
<td>java</td>
<td>gatk-private</td>
<td>o.b.s</td>
<td>gatk-tools-private</td>
<td>o.b.g.utils</td>
<td style="text-align: right;">7</td>
</tr>
<tr>
<td>private</td>
<td>java</td>
<td>gatk-private</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-private</td>
<td>o.b.g.engine</td>
<td style="text-align: right;">5</td>
</tr>
<tr>
<td>private</td>
<td>java</td>
<td>gatk-private</td>
<td>o.b.s.gatk</td>
<td>gatk-tools-private</td>
<td>o.b.g.tools</td>
<td style="text-align: right;">133</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>queue-framework</td>
<td>o.b.s</td>
<td>gatk-queue</td>
<td>o.b.g</td>
<td style="text-align: right;">2</td>
</tr>
<tr>
<td>public</td>
<td>scala</td>
<td>queue-framework</td>
<td>o.b.s</td>
<td>gatk-queue</td>
<td>o.b.g</td>
<td style="text-align: right;">72</td>
</tr>
<tr>
<td>public</td>
<td>scala</td>
<td>queue-framework</td>
<td>o.b.s</td>
<td>gatk-queue-extensions-public</td>
<td>o.b.g</td>
<td style="text-align: right;">31</td>
</tr>
<tr>
<td>public</td>
<td>qscripts</td>
<td>queue-framework</td>
<td>o.b.s</td>
<td>gatk-queue-extensions-public</td>
<td>o.b.g</td>
<td style="text-align: right;">12</td>
</tr>
<tr>
<td>private</td>
<td>scala</td>
<td>queue-private</td>
<td>o.b.s</td>
<td>gatk-queue-private</td>
<td>o.b.g</td>
<td style="text-align: right;">2</td>
</tr>
<tr>
<td>private</td>
<td>qscripts</td>
<td>queue-private</td>
<td>o.b.s</td>
<td>gatk-queue-extensions-internal</td>
<td>o.b.g</td>
<td style="text-align: right;">118</td>
</tr>
</tbody>
</table>
<p><strong>During all future code modifications and pull requests, classes should be refactored to correct artifacts and package as follows.</strong></p>
<p>All non-engine tools should be in the tools artifacts, with appropriate sub-package names.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Scope</th>
<th>Type</th>
<th>Artifact</th>
<th>Package(s)</th>
</tr>
</thead>
<tbody>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-utils</td>
<td>o.b.g.utils</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-engine</td>
<td>o.b.g.engine</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-tools-public</td>
<td>o.b.g.tools.walkers</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-tools-public</td>
<td>o.b.g.tools.*</td>
</tr>
<tr>
<td>protected</td>
<td>java</td>
<td>gatk-tools-protected</td>
<td>o.b.g.tools.walkers</td>
</tr>
<tr>
<td>protected</td>
<td>java</td>
<td>gatk-tools-protected</td>
<td>o.b.g.tools.*</td>
</tr>
<tr>
<td>private</td>
<td>java</td>
<td>gatk-tools-private</td>
<td>o.b.g.tools.walkers</td>
</tr>
<tr>
<td>private</td>
<td>java</td>
<td>gatk-tools-private</td>
<td>o.b.g.tools.*</td>
</tr>
<tr>
<td>public</td>
<td>java</td>
<td>gatk-queue</td>
<td>o.b.g.queue</td>
</tr>
<tr>
<td>public</td>
<td>scala</td>
<td>gatk-queue</td>
<td>o.b.g.queue</td>
</tr>
<tr>
<td>public</td>
<td>scala</td>
<td>gatk-queue-extensions-public</td>
<td>o.b.g.queue.extensions</td>
</tr>
<tr>
<td>public</td>
<td>qscripts</td>
<td>gatk-queue-extensions-public</td>
<td>o.b.g.queue.qscripts</td>
</tr>
<tr>
<td>private</td>
<td>scala</td>
<td>gatk-queue-private</td>
<td>o.b.g.queue</td>
</tr>
<tr>
<td>private</td>
<td>qscripts</td>
<td>gatk-queue-extensions-internal</td>
<td>o.b.g.queue.qscripts</td>
</tr>
</tbody>
</table>
<h2>Renamed Classes</h2>
<p>The following class names were updated to replace Sting with GATK.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Old Sting class</th>
<th>New GATK class</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>ArtificialStingSAMFileWriter</code></td>
<td><code>ArtificialGATKSAMFileWriter</code></td>
</tr>
<tr>
<td><code>ReviewedStingException</code></td>
<td><code>ReviewedGATKException</code></td>
</tr>
<tr>
<td><code>StingException</code></td>
<td><code>GATKException</code></td>
</tr>
<tr>
<td><code>StingSAMFileWriter</code></td>
<td><code>GATKSAMFileWriter</code></td>
</tr>
<tr>
<td><code>StingSAMIterator</code></td>
<td><code>GATKSAMIterator</code></td>
</tr>
<tr>
<td><code>StingSAMIteratorAdapter</code></td>
<td><code>GATKSAMIteratorAdapter</code></td>
</tr>
<tr>
<td><code>StingSAMRecordIterator</code></td>
<td><code>GATKSAMRecordIterator</code></td>
</tr>
<tr>
<td><code>StingTextReporter</code></td>
<td><code>GATKTextReporter</code></td>
</tr>
</tbody>
</table>
<h1>Common Git/Maven Issues</h1>
<h2>Renamed files</h2>
<p>The 3.2 renaming patch is actually split into two commits. The first commit renames the files without making any content changes, while the second changes the contents of the files without changing any file paths.</p>
<p>When dealing with renamed files, it is best to work with a clean directory during rebasing. It will be easier for you track files that you may not have added to git.</p>
<p>After running a git rebase or merge, you may first run into problems with files that you renamed and were moved during the GATK 3.2 package renaming. As a general rule, the renaming only changes directory names. The exception to this rule are classes such as StingException that are renamed to GATKException, and are listed under [Renamed Classes]. The workflow for resolving these merge issues is to find the list of your renamed files, put your content in the correct location, then register the changes with git.</p>
<p>To obtain the list of renamed directories and files:</p>
<ol>
<li>Use <code>git status</code> to get a list of affected files</li>
<li>Find the common old directory and file name under &quot;both deleted&quot;</li>
<li>Find your new file name under &quot;added by them&quot; (yes, you are &quot;them&quot;)</li>
<li>Find the new directory under &quot;added by us&quot;</li>
</ol>
<p>Then, to resolve the issue for each file:</p>
<ol>
<li>Move your copy of your renamed file to the new directory</li>
<li><code>git rm</code> the old paths as appropriate</li>
<li><code>git add</code> the new path</li>
<li>Repeat for other files until git status shows &quot;all conflicts fixed&quot;</li>
</ol>
<p>Upon first rebasing you will see a lot of text. <strong>At this moment, you can ignore most of it, and use git status instead.</strong></p>
<p>For the purposes of illustration, while running <code>git rebase</code> it is perfectly normal to see something similar to:</p>
<pre><code class="pre_md">$ git rebase master
First, rewinding head to replay your work on top of it...
Applying: &lt;&lt;&lt; Your first commit message here &gt;&gt;&gt;
Using index info to reconstruct a base tree...
A   protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngine.java
A   protected/gatk-protected/src/test/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngineUnitTest.java
&lt;&lt;&lt;Other files that you renamed.&gt;&gt;&gt;
warning: squelched 12 whitespace errors
warning: 34 lines add whitespace errors.
Falling back to patching base and 3-way merge...
CONFLICT (rename/rename): Rename "protected/gatk-protected/src/test/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngineUnitTest.java"-&gt;"protected/gatk-tools-protected/src/test/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngineUnitTest.java" in branch "HEAD" rename "protected/gatk-protected/src/test/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngineUnitTest.java"-&gt;"protected/gatk-protected/src/test/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/HaplotypeCallerGenotypingEngineUnitTest.java" in "&lt;&lt;&lt; Your first commit message here &gt;&gt;&gt;"
CONFLICT (rename/rename): Rename "protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngine.java"-&gt;"protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngine.java" in branch "HEAD" rename "protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngine.java"-&gt;"protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java" in "&lt;&lt;&lt; Your first commit message here &gt;&gt;&gt;"
Failed to merge in the changes.
Patch failed at 0001 Example conflict.
The copy of the patch that failed is found in:
   /Users/zzuser/src/gsa-unstable/.git/rebase-apply/patch

When you have resolved this problem, run "git rebase --continue".
If you prefer to skip this patch, run "git rebase --skip" instead.
To check out the original branch and stop rebasing, run "git rebase --abort".

$</code class="pre_md"></pre>
<p>While everything you need to resolve the issue is technically in the message above, it may be much easier to track what's going on using <code>git status</code>.</p>
<pre><code class="pre_md">$ git status
rebase in progress; onto cba4321
You are currently rebasing branch 'zz_renaming_haplotypecallergenotypingengine' on 'cba4321'.
  (fix conflicts and then run "git rebase --continue")
  (use "git rebase --skip" to skip this patch)
  (use "git rebase --abort" to check out the original branch)

Unmerged paths:
  (use "git reset HEAD &lt;file&gt;..." to unstage)
  (use "git add/rm &lt;file&gt;..." as appropriate to mark resolution)

    added by them:      protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java
    both deleted:       protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngine.java
    added by them:      protected/gatk-protected/src/test/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/HaplotypeCallerGenotypingEngineUnitTest.java
    both deleted:       protected/gatk-protected/src/test/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngineUnitTest.java
    added by us:        protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngine.java
    added by us:        protected/gatk-tools-protected/src/test/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngineUnitTest.java

Untracked files:
  (use "git add &lt;file&gt;..." to include in what will be committed)

&lt;&lt;&lt; possible untracked files if your working directory is not clean&gt;&gt;&gt;

no changes added to commit (use "git add" and/or "git commit -a")
$ </code class="pre_md"></pre>
<p>Let's look at the main java file as an example. If you are having issues figuring out the new directory and new file name, they are all listed in the output.</p>
<pre><code class="pre_md">Path in the common ancestor branch:
 |      old source directory       |                     old package name                     |   old file name     |
  protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/GenotypingEngine.java

Path in the new master branch before merge:
 |           new source directory             |                 new package name                    |   old file name     |
  protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngine.java

Path in your branch before merge:
 |      old source directory       |                     old package name                     |           new file name            |
  protected/gatk-protected/src/main/java/org/broadinstitute/sting/gatk/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java

Path in your branch post merge:
 |           new source directory             |                 new package name                    |           new file name            |
  protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java    </code class="pre_md"></pre>
<p>After identifying the new paths for use post merge, use the following workflow for each file:</p>
<ol>
<li>Move or copy your version of the renamed file to the new directory</li>
<li><code>git rm</code> the three old file paths: common ancestor, old directory with new file name, and new directory with old file name</li>
<li><code>git add</code> the new file name in the new directory</li>
</ol>
<p>After you process all files correctly, in the output of <code>git status</code> you should see the &quot;all conflicts fixed&quot; and all your files renamed.</p>
<pre><code class="pre_md">$ git status
rebase in progress; onto cba4321
You are currently rebasing branch 'zz_renaming_haplotypecallergenotypingengine' on 'cba4321'.
  (all conflicts fixed: run "git rebase --continue")

Changes to be committed:
  (use "git reset HEAD &lt;file&gt;..." to unstage)

    renamed:    protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngine.java -&gt; protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java
    renamed:    protected/gatk-tools-protected/src/test/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/GenotypingEngineUnitTest.java -&gt; protected/gatk-tools-protected/src/test/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/HaplotypeCallerGenotypingEngineUnitTest.java

Untracked files:
  (use "git add &lt;file&gt;..." to include in what will be committed)

&lt;&lt;&lt; possible untracked files if your working directory is not clean&gt;&gt;&gt;

$</code class="pre_md"></pre>
<p>Continue your rebase, handling other merges as normal.</p>
<pre><code class="pre_md">$ git rebase --continue</code class="pre_md"></pre>
<h2>Fixing imports</h2>
<p>Because all the packages names are different in 3.2, while rebasing you may run into conflicts due to imports you also changed. Use your favorite editor to fix the imports within the files. Then try recompiling, and repeat as necessary until your code works.</p>
<p>While editing the files with conflicts with a basic text editor may work, IntelliJ IDEA also offers a special merge tool that may help via the menu:</p>
<pre><code class="pre_md">VCS &gt; Git &gt; Resolve Conflicts...</code class="pre_md"></pre>
<p>For each file, click on the &quot;Merge&quot; button in the first dialog. Use the various buttons in the <a href="https://www.jetbrains.com/idea/webhelp/resolving-conflicts.html">Conflict Resolution Tool</a> to automatically accept any changes that are not in conflict. Then find any edit any remaining conflicts that require further manual intervention.</p>
<p>Once you begin editing the import statements in the three way merge tool, another IntelliJ IDEA 13.1 feature that may speed up repairing blocks of import statements is <a href="http://blog.jetbrains.com/idea/2014/03/intellij-idea-13-1-rc-introduces-sublime-text-style-multiple-selections/">Multiple Selections</a>. Find a block of import lines that need the same changes. Hold down the option key as you drag your cursor vertically down the edit point on each import line. Then begin typing or deleting text from the multiple lines.</p>
<h2>Switching branches</h2>
<p>Even after a successful merge, you may still run into stale GATK code or links from modifications before and after the 3.2 package renaming. To significantly reduce these chances, run <code>mvn clean</code> <em>before</em> and then again <em>after</em> switching branches.</p>
<p>If this doesn't work, run <code>mvn clean &amp;&amp; git status</code>, looking for any directories you don't that shouldn't be in the current branch. It is possible that some files were not correctly moved, including classes or test resources. Find the file still in the old directories via a command such as <code>find public/gatk-framework -type f</code>. Then move them to the correct new directories and commit them into git.</p>
<h2>Slow Builds with Queue and Private</h2>
<p>Due to the [Renamed Binary Packages], the separate artifacts including and excluding private code are now packaged during the Maven package build lifecycle.</p>
<p>When building packages, to significantly speed up the default packaging time, if you only require the GATK tools run <code>mvn verify -P\!queue</code>.</p>
<p>Alternatively, if you do not require building private source, then disable private compiling via <code>mvn verify -P\!private</code>. </p>
<p>The two may be combined as well via: <code>mvn verify -P\!queue,\!private</code>. </p>
<p>The exclamation mark is a shell command that must be escaped, in the above case with a backslash. Shell quotes may also be used: <code>mvn verify -P'!queue,!private'</code>.</p>
<p>Alternatively, developers with access to private may often want to disable packaging the protected distributions. In this case, use the <code>gsadev</code> profile. This may be done via <code>mvn verify -Pgsadev</code> or, excluding Queue, <code>mvn verify -Pgsadev,\!queue</code>.</p>
<h2>Stale symlinks</h2>
<p>Users see errors from maven when an unclean repo in git is updated.
Because BaseTest.java currently hardcodes relative paths to
&quot;public/testdata&quot;, maven creates these symbolic links all over the
file system to help the various tests in different modules find the
relative path &quot;<current module>/public/testdata&quot;.</p>
<p>However, our Maven support has evolved from 2.8, to 3.0, to now the
3.2 renaming, each time has changed the symbolic link's target
directory. Whenever a stale symbolic link to an old testdata directory
remains in the users folder, maven is saying it will not remove the
link, because maven basically doesn't know why the link is pointing to
the wrong folder (answer, the link is from an old git checkout) and
thinks it's a bug in the build.</p>
<p>If one doesn't have an stale / unclean maven repo when updating git
via merge/rebase/checkout, you will never see this issue.</p>
<p>The script that can remove the stale symlinks, <code>public/src/main/scripts/shell/delete_maven_links.sh</code>, should run automatically during a <code>mvn test-compile</code> or <code>mvn verify</code>.</p>