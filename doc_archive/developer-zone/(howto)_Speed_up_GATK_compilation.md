## (howto) Speed up GATK compilation

http://gatkforums.broadinstitute.org/gatk/discussion/5784/howto-speed-up-gatk-compilation

<hr />
<p>TL;DR: <code>mvn -Ddisable.shadepackage verify</code></p>
<hr />
<h3>Background</h3>
<p>In addition to Queue's GATK-wrapper codegen, relatively slow scala compilation, etc. there's still a lot of legacy compatibility from our <code>ant</code> days in the Maven scripts. Our <code>mvn verify</code> behaves more like when one runs <code>ant</code>, and builds <em>everything</em> needed to bundle the GATK.</p>
<p>As of GATK 3.4, by default the build for the &quot;protected&quot; code generates jar files that contains every class needed for running, one for the GATK and one for Queue. This is done by the <a href="https://maven.apache.org/plugins/maven-shade-plugin/">Maven shade plugin</a>, and are each called the &quot;package jar&quot;. But, there's a way to generate a jar file that only contains <code>META-INF/MANIFEST.MF</code> pointers to the dependency jar files, instead of zipping/shading them up. These are each the &quot;executable jar&quot;, and FYI are always generated as it takes seconds, not minutes.</p>
<hr />
<h3>Instructions for fast compilation</h3>
<p>While developing and recompiling Queue, disable the shaded jar with <code>-Ddisable.shadepackage</code>. Then run <code>java -jar target/executable/Queue.jar ...</code> If you need to transfer this jar to another machine / directory, you can't copy (or rsync) just the jar, you'll need the entire executable directory.</p>
<pre><code class="pre_md"># Total expected time, on a local disk, with Queue:
#   ~5.0 min from clean
#   ~1.5 min per recompile
mvn -Ddisable.shadepackage verify

# always available
java -jar target/executable/Queue.jar --help

# not found when shade disabled
java -jar target/package/Queue.jar --help</code class="pre_md"></pre>
<p>If one is only developing for the GATK, skip Queue by adding  <code>-P\!queue</code> also.</p>
<pre><code class="pre_md">mvn -Ddisable.shadepackage -P\!queue verify

# always available
java -jar target/executable/GenomeAnalysisTK.jar --help

# not found when queue profile disabled
java -jar target/executable/Queue.jar --help</code class="pre_md"></pre>