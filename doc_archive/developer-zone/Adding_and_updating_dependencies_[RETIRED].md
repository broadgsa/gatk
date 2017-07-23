## Adding and updating dependencies [RETIRED]

http://gatkforums.broadinstitute.org/gatk/discussion/1352/adding-and-updating-dependencies-retired

<h2>Adding Third-party Dependencies</h2>
<p>The GATK build system uses the <a href="http://ant.apache.org/ivy/">Ivy dependency manager</a> to make it easy for our users to add additional dependencies.  Ivy can pull the latest jars and their dependencies from the <a href="http://mvnrepository.com">Maven repository</a>, making adding or updating a dependency as simple as adding a new line to the <code>ivy.xml</code> file.</p>
<p>If your tool is available in the maven repository, add a line to the <code>ivy.xml</code> file similar to the following:</p>
<pre><code class="pre_md">&lt;dependency org="junit" name="junit" rev="4.4" /&gt;</code class="pre_md"></pre>
<p>If you would like to add a dependency to a tool not available in the maven repository, please email <a href="mailto:gsahelp@broadinstitute.org">gsahelp@broadinstitute.org</a></p>
<h2>Updating SAM-JDK and Picard</h2>
<p>Because we work so closely with the SAM-JDK/Picard team and are critically dependent on the code they produce, we have a special procedure for updating the SAM/Picard jars.  Please use the following procedure to when updating <code>sam-*.jar</code> or <code>picard-*.jar</code>.</p>
<ul>
<li>
<p>Download and build the latest versions of <a href="http://picard.svn.sourceforge.net/svnroot/picard/trunk/">Picard public</a> and <a href="https://svnrepos.broad.mit.edu/picard/trunk">Picard private</a> (restricted to Broad Institute users) from their respective svns.  </p>
</li>
<li>
<p>Get the latest svn versions for picard public and picard private by running the following commands:</p>
<p>svn info $PICARD_PUBLIC_HOME | grep &quot;Revision&quot;
svn info $PICARD_PRIVATE_HOME | grep &quot;Revision&quot;</p>
</li>
</ul>
<h3>Updating the Picard public jars</h3>
<ul>
<li>
<p>Rename the jars and xmls in <code>$STING_HOME/settings/repository/net.sf</code> to <code>{picard|sam}-$PICARD_PUBLIC_MAJOR_VERSION.$PICARD_PUBLIC_MINOR_VERSION.PICARD_PUBLIC_SVN_REV.{jar|xml}</code></p>
</li>
<li>
<p>Update the jars in <code>$STING_HOME/settings/repository/net.sf</code> with their newer equivalents in <code>$PICARD_PUBLIC_HOME/dist/picard_lib</code>.</p>
</li>
<li>Update the xmls in <code>$STING_HOME/settings/repository/net.sf</code> with the appropriate version number (<code>$PICARD_PUBLIC_MAJOR_VERSION.$PICARD_PUBLIC_MINOR_VERSION.$PICARD_PUBLIC_SVN_REV</code>).</li>
</ul>
<h3>Updating the Picard private jar</h3>
<ul>
<li>
<p>Create the picard private jar with the following command:</p>
<p>ant clean package -Dexecutable=PicardPrivate -Dpicard.dist.dir=${PICARD_PRIVATE_HOME}/dist</p>
</li>
<li>
<p>Rename <code>picard-private-parts-*.jar</code> in <code>$STING_HOME/settings/repository/edu.mit.broad</code> to <code>picard-private-parts-$PICARD_PRIVATE_SVN_REV.jar</code>.</p>
</li>
<li>
<p>Update <code>picard-private-parts-*.jar</code> in <code>$STING_HOME/settings/repository/edu.mit.broad</code> with the <code>picard-private-parts.jar</code> in <code>$STING_HOME/dist/packages/picard-private-parts</code>.</p>
</li>
<li>Update the xml in <code>$STING_HOME/settings/repository/edu.mit.broad</code> to reflect the new revision and publication date.</li>
</ul>