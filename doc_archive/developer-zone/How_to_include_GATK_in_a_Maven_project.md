## How to include GATK in a Maven project

http://gatkforums.broadinstitute.org/gatk/discussion/6214/how-to-include-gatk-in-a-maven-project

<p>GATK 3.x releases are not currently published to Central. But it is possible to install the GATK into your local repository, where Maven can then pick up the GATK as a dependency.</p>
<hr />
<p><strong>TL;DR</strong> Clone GATK 3.4, <code>mvn install</code>, then use the GATK as any other artifact.</p>
<hr />
<p>The repository you should use depends on what is your goal.</p>
<p>If you want to build your own analysis tools on top of the GATK engine (not including the GATK analysis tools), with the option of distributing your project to others, you should clone the <a href="https://github.com/broadgsa/gatk"><code>gatk</code></a> repo.</p>
<p>If you want to integrate the full GATK into a project for in-house purposes (redistribution is not allowed under the licensing terms), in which your tools can call GATK tools directly internally, you should clone <a href="https://github.com/broadgsa/gatk-protected"><code>gatk-protected</code></a>. This can be done by running the following code:</p>
<pre><code>: 'GATK 3.4 code has known issues with the Java 8 compiler. Make sure you are using Java 7.'
java -version

: 'The entire GATK repo is relatively large. This only downloads 3.4.'
git clone -b 3.4 --depth 1 git@github.com:broadgsa/gatk-protected.git gatk-protected-3.4
cd gatk-protected-3.4

: 'Install the gatk into a the local ~/.m2/repository, where your project can then refer to the GATK.'
mvn install

: 'Build the "external example" as a demo of using the GATK as a library.'
cd public/external-example
mvn verify
java -jar target/external-example-1.0-SNAPSHOT.jar -T MyExampleWalker --help</code></pre>
<p>After the GATK is installed, add this dependency to your Maven artifact, and all other GATK dependencies will be included as well.</p>
<pre><code>&lt;dependency&gt;
    &lt;groupId&gt;org.broadinstitute.gatk&lt;/groupId&gt;
    &lt;artifactId&gt;gatk-tools-protected&lt;/artifactId&gt;
    &lt;version&gt;3.4&lt;/version&gt;
&lt;/dependency&gt;</code></pre>
<p>One thing you might run into is that the GATK artifacts, and hence the external-example, transitively depend on artifacts that are also not in Central. They are instead committed under the path <code>public/repo</code>.  Like in the <code>public/external-example/pom.xml</code>, your Maven project may need to include this directory as an additional repository. That being said <code>mvn install</code> <em>should</em> copy the artifacts to <code>~/.m2/repository</code> for you. For example, after the install, you should have a directory <code>~/.m2/repository/com/google/code/cofoja/cofoja</code>.</p>
<p>If you somehow need to add the GATK's public repo as a repository, use a repository element like the one below:</p>
<pre><code>&lt;repositories&gt;
    &lt;repository&gt;
        &lt;id&gt;gatk.public.repo.local&lt;/id&gt;
        &lt;name&gt;GATK Public Local Repository&lt;/name&gt;
        &lt;url&gt;file:/Users/someuser/src/gatk-protected-3.4/public/repo&lt;/url&gt;
    &lt;/repository&gt;
&lt;/repositories&gt;</code></pre>
<p>Since the GATK is not in Central, each developer will need to install the GATK 3.4 once. Or, as an advanced step, your may also want to explore publishing the GATK on one of your shared local systems. If you have a shared filesystem you'd like to use as a repository, publish the GATK 3.4 to the directory using <code>mvn install -Dmaven.repo.local=/mount/path/to/shared/repo</code>, and then add a repository element to your Maven project. If your team is using a Maven repository such as Artifactory or Nexus, we can't provide guidance for publishing &quot;third party&quot; artifacts. But it should theoretically be possible, with instructions hopefully available through either Maven or the repository manager's help forums.</p>