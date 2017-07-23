## Migration from Apache Ant to Apache Maven

http://gatkforums.broadinstitute.org/gatk/discussion/3437/migration-from-apache-ant-to-apache-maven

<h1>Overview</h1>
<hr />
<p><strong>We're replacing Ant with Maven. To build, run <code>mvn verify</code>.</strong></p>
<h2>Background</h2>
<p>In the early days of the Genome Analysis Toolkit (GATK), the code base separated the GATK genomics engine from the core java utilities, encompassed in a wider project called Sting. During this time, the build tool of choice was the relatively flexible Java build tool <a href="http://ant.apache.org">Apache Ant</a>, run via the command <code>ant</code>.</p>
<p>As our code base expanded to more and more packages, groups internal and external to GSA, and the Broad, have expressed interest in using portions of Sting/GATK as modules in larger projects. Unfortunately over time, many parts of the GATK and Sting intermingled, producing the current situation where developers finds it easier to copy the monolithic GATK instead, or individual java files, instead of using the tools as libraries.</p>
<p>The goal of this first stage is to split the parts of the monolithic Sting/GATK into easily recognizable sub artifacts. The tool used to accomplish this task is <a href="http://maven.apache.org">Apache Maven</a>, also known as <em>Maven</em>, and run via the command <code>mvn</code>. Maven convention encourages developers to separate code, and accompanying resources, into a hierarchical structure of reusable artifacts. Maven attempts to avoid build configuration, preferring source repositories to lay out code in a conventional structure. When needed, a Maven configuration file called <em>pom.xml</em> specifies each artifact's build configuration, that one may think of as similar to an Ant <em>build.xml</em>.</p>
<p>The actual migration consisted of zero changes to the contents of existing Java source files, easing git merges and rebasing. The Java files from public, protected, and private have all moved into Maven conventional child artifacts, with each artifact containing a separate <em>pom.xml</em>.</p>
<h1>Examples</h1>
<h2>Obtaining the GATK with Maven support</h2>
<p>Clone the repository:</p>
<p><code>git clone ssh://git@github.com/broadinstitute/gsa-unstable.git cd gsa-unstable</code></p>
<h2>Building GATK and Queue</h2>
<p>Clone the repository:</p>
<p><code>git clone ssh://git@github.com/broadinstitute/gsa-unstable.git cd gsa-unstable</code></p>
<p>If running on a Broad server, add maven to your environment via the dotkit:</p>
<p><code>reuse Maven-3.0.3</code></p>
<p>Build all of Sting, including packaged versions of the GATK and Queue:</p>
<p><code>mvn verify</code></p>
<p>The packaged, executable jar files will be output to:</p>
<p><code>public/gatk-package/target/gatk-package-2.8-SNAPSHOT.jar public/queue-package/target/queue-package-2.8-SNAPSHOT.jar</code></p>
<p>Find equivalent maven commands for existing ant targets:</p>
<p><code>./ant-bridge.sh &lt;target&gt; &lt;properties&gt;</code></p>
<p>Example output:</p>
<p><code>$ ./ant-bridge.sh fasttest -Dsingle=GATKKeyUnitTest Equivalent maven command mvn verify -Dsting.committests.skipped=false -pl private/gatk-private -am -Dresource.bundle.skip=true -Dit.test=disabled -Dtest=GATKKeyUnitTest $</code></p>
<h2>Running the GATK and Queue</h2>
<p>To run the GATK, or copy the compiled jar, find the packaged jar under public/gatk-package/target</p>
<p><code>public/gatk-package/target/gatk-package-2.8-SNAPSHOT.jar</code></p>
<p>To run Queue, the jar is under the similarly named public/queue-package/target</p>
<p><code>public/queue-package/target/queue-package-2.8-SNAPSHOT.jar</code></p>
<p><strong>NOTE:</strong> Unlike builds with Ant, you <em>cannot</em> execute the jar file built by the gatk-framework module. This is because maven does not include dependent artifacts in the target folder with assembled framework jar. Instead, use the packaged jars, listed above, that contain all the classes and resources needed to run the GATK, or Queue.</p>
<h2>Excluding Queue</h2>
<p><em>NOTE:</em> If you make changes to sting-utils, gatk-framework, or any other dependencies <em>and</em> disable queue, you may accidentally end up breaking the full repository build without knowing.</p>
<p>The Queue build contributes a majority portion of the Sting project build time. To exclude Queue from your build, run maven with either (the already shell escaped) <code>-P\!queue</code> or <code>-Ddisable.queue</code>. Currently the latter property also disables the maven queue profile. This allows one other semi-permanent option to disable building Queue as part of the Sting repository. Configure your local Maven settings to always pass the property <code>-Ddisable.queue</code> by adding and activating a custom profile in your local ~/.m2/settings.xml</p>
<p>```$ cat ~/.m2/settings.xml</p>
<settings>
  <!--
  Other settings.xml changes...
  -->
  <!--
  Define a new profile to set disable.queue
  -->
  <profiles>
    <profile>
      <id>disable.queue</id>
      <properties>
        <disable.queue>true</disable.queue>
      </properties>
    </profile>
  </profiles>
  <!--
  Activate the profile defined above
  -->
  <activeProfiles>
    <activeProfile>disable.queue</activeProfile>
  </activeProfiles>
</settings>
<p>$```</p>
<h2>Using the GATK framework as a module</h2>
<p>Currently the GATK artifacts are not available via any centralized repository. To build code using the GATK you must still have a checkout of the GATK source code, and install the artifacts to your local mvn repository (by default ~/.m2/repository). The installation copies the artifacts to your local repo such that it may be used by your external project. The checkout of the local repo provides several artifacts under <code>public/repo</code> that will be required for your project.</p>
<p>After updating to the latest version of the Sting source code, install the Sting artifacts via:</p>
<p><code>mvn install</code></p>
<p>After the GATK has been installed locally, in your own source repository, include the artifact gatk-framework as a library.</p>
<p>In Apache Maven add this dependency:</p>
<p>```<dependency></p>
<groupId>org.broadinstitute.sting</groupId>
<pre><code class="pre_md">&lt;artifactId&gt;gatk-framework&lt;/artifactId&gt;
&lt;version&gt;2.8-SNAPSHOT&lt;/version&gt;</code class="pre_md"></pre>
<p></dependency>```</p>
<p>For Apache Ivy, you may need to specify <code>~/.m2/repository</code> as a local repo. Once the local repository has been configured, ivy may find the dependency via:</p>
<p><code>&lt;dependency org="org.broadinstitute.sting" name="gatk-framework" rev="2.8-SNAPSHOT" /&gt;</code></p>
<p>If you decide to also use Maven to build your project, your source code should go under the conventional directory <code>src/main/java</code>. The <code>pom.xml</code> contains any special configuration for your project. To see an example pom.xml and maven conventional project structure in:</p>
<p><code>public/external-example</code></p>
<h2>Moved directories</h2>
<p>If you have an old git branch that needs to be merged, you may need to know where to move files in order for your classes to now build with Maven. In general, most directories were moved with minimal or no changes.</p>
<table class="table table-striped">
<thead>
<tr>
<th><strong>Old directory</strong></th>
<th><strong>New maven directory</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td>private/java/src/</td>
<td>private/gatk-private/src/main/java/</td>
</tr>
<tr>
<td>private/R/scripts/</td>
<td>private/gatk-private/src/main/resources/</td>
</tr>
<tr>
<td>private/java/test/</td>
<td>private/gatk-private/src/test/java/</td>
</tr>
<tr>
<td>private/testdata/</td>
<td>private/gatk-private/src/test/resources/</td>
</tr>
<tr>
<td>private/scala/qscript/</td>
<td>private/queue-private/src/main/qscripts/</td>
</tr>
<tr>
<td>private/scala/src/</td>
<td>private/queue-private/src/main/scala/</td>
</tr>
<tr>
<td>private/scala/test/</td>
<td>private/queue-private/src/test/scala/</td>
</tr>
<tr>
<td>protected/java/src/</td>
<td>protected/gatk-protected/src/main/java/</td>
</tr>
<tr>
<td>protected/java/test/</td>
<td>protected/gatk-protected/src/test/java/</td>
</tr>
<tr>
<td>public/java/src/</td>
<td>public/gatk-framework/src/main/java/</td>
</tr>
<tr>
<td>public/java/test/</td>
<td>public/gatk-framework/src/test/java/</td>
</tr>
<tr>
<td>public/testdata/</td>
<td>public/gatk-framework/src/test/resources/</td>
</tr>
<tr>
<td>public/scala/qscript/</td>
<td>public/queue-framework/src/main/qscripts/</td>
</tr>
<tr>
<td>public/scala/src/</td>
<td>public/queue-framework/src/main/scala/</td>
</tr>
<tr>
<td>public/scala/test/</td>
<td>public/queue-framework/src/test/scala/</td>
</tr>
</tbody>
</table>
<h1>Future Directions</h1>
<h2>Further segregate source code</h2>
<p>Currently, the artifacts sting-utils and the gatk-framework contain intertwined code bases. This leads to the current setup where all sting-utils code is actually found in the gatk-framework artifact, including generic utilities that could be used by other software modules. In the future, all elements under <code>org.broadinstitute.sting.gatk</code> will be located the gatk-framework, while all other packages under <code>org.broadinstitut.sting</code> will be evaluated and then separated under the gatk-framework or sting-utils artifacts.</p>
<h2>Publishing artifacts</h2>
<p>Tangentially related to segregating sting-utils and the gatk-framework, the current Sting and GATK artifacts are ineligible to be pushed to the <a href="http://search.maven.org">Maven Central Repository</a>, due to several other issues:</p>
<ul>
<li>Need to provide trivial workflow for Picard, and possibly SnpEff, to submit to central</li>
<li>Missing <a href="https://docs.sonatype.org/display/Repository/Sonatype+OSS+Maven+Repository+Usage+Guide#SonatypeOSSMavenRepositoryUsageGuide-6.CentralSyncRequirement">meta files</a> for the jars:
<ul>
<li>*-sources.jar</li>
<li>*-javadoc.jar</li>
<li>*.md5</li>
<li>*.sha1</li>
</ul></li>
</ul>
<p><em>NOTE:</em> Artifact jars do NOT need to actually be in Central, and may be available as pom reference only, for example <a href="http://central.maven.org/maven2/com/oracle/ojdbc14/">Oracle ojdbc</a>.</p>
<p>In the near term, we could use a private repos based on <a href="http://www.jfrog.com/home/v_artifactorycloud_overview">Artifactory</a> or <a href="http://www.sonatype.org/nexus">Nexus</a> (<a href="http://docs.codehaus.org/display/MAVENUSER/Maven+Repository+Manager+Feature+Matrix">comparison</a>). After more work of adding, cleaning up, or centrally publishing all the dependencies for Sting, we may then publish into the basic Central repo. Or, we could move to a social service like <a href="https://bintray.com">BinTray</a> (think GitHub vs. Git).</p>
<h1>Status Updates</h1>
<h2>February 13, 2014</h2>
<p>Maven is now the default in gsa-unstable's master branch. For GATK developers, the git migration is effectively complete. Software engineers are resolving a few remaining issues related to the automated build and testing infrastructure, but the basic workflow for developers should now be up to date.</p>
<h2>January 30, 2014</h2>
<p>The migration to to maven has begun in the <a href="https://github.com/broadinstitute/gsa-unstable">gsa-unstable repository</a> on the ks_new_maven_build_system branch.</p>
<h2>November 5, 2013</h2>
<p>The maven port of the existing ant build resides in the <a href="https://github.com/broadinstitute/gsa-qc">gsa-qc repository</a>.</p>
<p>This is an old branch of Sting/GATK, with the existing files relocated to Maven appropriate locations, pom.xml files added, along with basic resources to assist in artifact generation.</p>