## Setting up your dev environment: Maven and IntelliJ for GATK 3+

http://gatkforums.broadinstitute.org/gatk/discussion/4023/setting-up-your-dev-environment-maven-and-intellij-for-gatk-3

<h3>Overview</h3>
<p>Since GATK 3.0, we use Apache Maven (instead of Ant) as our build system, and IntelliJ as our IDE (Integrated Development Environment). This document describes how to get set up to use Maven as well as how to create an IntelliJ project around our Maven project structure.</p>
<h3>Before you start</h3>
<ul>
<li>Ensure that you have git clones of our repositories on your machine. See <a href="http://www.broadinstitute.org/gatk/guide/article?id=4022">this document</a> for details on obtaining the GATK source code from our Git repos.</li>
</ul>
<h3>Setting up Maven</h3>
<ol>
<li>
<p>Check whether you can run <code>mvn --version</code> on your machine. If you can't, install Maven from <a href="http://maven.apache.org/">here</a>.</p>
</li>
<li>
<p>Ensure that the JAVA_HOME environment variable is properly set. If it's not, add the appropriate line to your shell's startup file:</p>
<p>for tcsh: </p>
<pre><code class="pre_md">setenv JAVA_HOME  \`/usr/libexec/java_home\`</code class="pre_md"></pre>
<p>for bash: </p>
<pre><code class="pre_md">export JAVA_HOME=\`/usr/libexec/java_home\`</code class="pre_md"></pre>
</li>
</ol>
<p>Note that the commands above use backticks, not single quotes.</p>
<h3>Basic Maven usage</h3>
<ol>
<li>
<p>To compile everything, type:</p>
<pre><code class="pre_md">mvn verify</code class="pre_md"></pre>
</li>
<li>
<p>To compile the GATK but not Queue (much faster!), the command is:</p>
<pre><code class="pre_md">mvn verify -P\!queue</code class="pre_md"></pre>
<p>Note that the <code>!</code> needs to be escaped with a backslash to avoid interpretation by the shell.</p>
</li>
<li>
<p>To obtain a clean working directory, type:</p>
<pre><code class="pre_md">mvn clean</code class="pre_md"></pre>
</li>
<li>
<p>If you're used to using ant to compile the GATK, you should be able to feed your old ant commands to the <code>ant-bridge.sh</code> script in the root directory. For example:</p>
<pre><code class="pre_md">./ant-bridge.sh test -Dsingle=MyTestClass</code class="pre_md"></pre>
</li>
</ol>
<h3>Setting up IntelliJ</h3>
<ol>
<li>
<p>Run <code>mvn test-compile</code> in your git clone's root directory.</p>
</li>
<li>
<p>Open IntelliJ</p>
</li>
<li>
<p>File -&gt; import project, select your git clone directory, then click &quot;ok&quot;</p>
</li>
<li>
<p>On the next screen, select &quot;import project from external model&quot;, then &quot;maven&quot;, then click &quot;next&quot;</p>
</li>
<li>
<p>Click &quot;next&quot; on the next screen without changing any defaults -- in particular:</p>
<ul>
<li>DON'T check &quot;Import maven projects automatically&quot;       </li>
<li>DON'T check &quot;Create module groups for multi-module maven projects&quot;</li>
</ul>
</li>
<li>
<p>On the &quot;Select Profiles&quot; screen, make sure private and protected ARE checked, then click &quot;next&quot;.</p>
</li>
<li>
<p>On the next screen, the &quot;gatk-aggregator&quot; project should already be checked for you -- if not, then check it. Click &quot;next&quot;.</p>
</li>
<li>
<p>Select the 1.7 SDK, then click &quot;next&quot;.</p>
</li>
<li>
<p>Select an appropriate project name (can be anything), then click &quot;next&quot; (or &quot;finish&quot;, depending on your version of IntelliJ).</p>
</li>
<li>
<p>Click &quot;Finish&quot; to create the new IntelliJ project.</p>
</li>
<li>
<p>That's it! Due to Maven magic, everything else will be set up for you automatically, including modules, libraries, Scala facets, etc.</p>
</li>
<li>You will see a popup &quot;Maven projects need to be imported&quot; on every IntelliJ startup. You should click import unless you're working on the actual pom files that make up the build system.</li>
</ol>