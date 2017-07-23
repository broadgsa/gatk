## How to get and install Firepony

http://gatkforums.broadinstitute.org/gatk/discussion/6020/how-to-get-and-install-firepony

<p>Binary packages for various versions of Linux are available at <a href="http://packages.shadau.com/">http://packages.shadau.com/</a></p>
<p>Below are installation instructions for Debian, Ubunto, CentOS and Fedora. For other Linux distributions, the Firepony source code is available at <a href="https://github.com/broadinstitute/firepony">https://github.com/broadinstitute/firepony</a> along with compilation instructions.</p>
<hr />
<h3>On Debian or Ubuntu systems</h3>
<p>The following commands can be used to install Firepony:</p>
<pre><code class="pre_md">sudo apt-get install software-properties-common
sudo add-apt-repository http://packages.shadau.com/
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-key 285514D704F4CDB7
sudo apt-get update
sudo apt-get install firepony</code class="pre_md"></pre>
<p>Once this initial install is done, updates will be automatically installed as part of the standard Ubuntu/Debian update procedure.</p>
<hr />
<h3>On CentOS 7 and Fedora 21 systems</h3>
<p>On CentOS 7, the following commands can be used to install Firepony:</p>
<pre><code class="pre_md">sudo curl -o /etc/yum.repos.d/packages.shadau.com.repo \
    http://packages.shadau.com/rpm/centos-7/packages.shadau.com.repo
sudo yum install firepony</code class="pre_md"></pre>
<p>For Fedora 21, use the following sequence of commands:</p>
<pre><code class="pre_md">sudo curl -o /etc/yum.repos.d/packages.shadau.com.repo \
    http://packages.shadau.com/rpm/fedora-21/packages.shadau.com.repo
sudo yum install firepony</code class="pre_md"></pre>
<p>Any subsequent updates will automatically be installed when running ‘yum update’.</p>