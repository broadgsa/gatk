## How can I prepare a FASTA file to use as reference?

http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

<p>This article describes the steps necessary to prepare your reference file (if it's not one that you got from us). As a complement to this article, see the relevant <a href="http://www.broadinstitute.org/gatk/guide/article?id=2798">tutorial</a>.</p>
<h3>Why these steps are necessary</h3>
<p>The GATK uses two files to access and safety check access to the reference files: a <code>.dict</code> dictionary of the contig names and sizes and a <code>.fai</code> fasta index file to allow efficient random access to the reference bases. You have to generate these files in order to be able to use a Fasta file as reference.</p>
<p><strong>NOTE: Picard and samtools treat spaces in contig names differently. We recommend that you avoid using spaces in contig names.</strong></p>
<h3>Creating the fasta sequence dictionary file</h3>
<p>We use CreateSequenceDictionary.jar from Picard to create a .dict file from a fasta file. </p>
<pre><code class="pre_md">&gt; java -jar CreateSequenceDictionary.jar R= Homo_sapiens_assembly18.fasta O= Homo_sapiens_assembly18.dict
[Fri Jun 19 14:09:11 EDT 2009] net.sf.picard.sam.CreateSequenceDictionary R= Homo_sapiens_assembly18.fasta O= Homo_sapiens_assembly18.dict
[Fri Jun 19 14:09:58 EDT 2009] net.sf.picard.sam.CreateSequenceDictionary done.
Runtime.totalMemory()=2112487424
44.922u 2.308s 0:47.09 100.2%   0+0k 0+0io 2pf+0w</code class="pre_md"></pre>
<p>This produces a SAM-style header file describing the contents of our fasta file.</p>
<pre><code class="pre_md">&gt; cat Homo_sapiens_assembly18.dict 
@HD     VN:1.0  SO:unsorted
@SQ     SN:chrM LN:16571        UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:d2ed829b8a1628d16cbeee88e88e39eb
@SQ     SN:chr1 LN:247249719    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:9ebc6df9496613f373e73396d5b3b6b6
@SQ     SN:chr2 LN:242951149    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:b12c7373e3882120332983be99aeb18d
@SQ     SN:chr3 LN:199501827    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:0e48ed7f305877f66e6fd4addbae2b9a
@SQ     SN:chr4 LN:191273063    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:cf37020337904229dca8401907b626c2
@SQ     SN:chr5 LN:180857866    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:031c851664e31b2c17337fd6f9004858
@SQ     SN:chr6 LN:170899992    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:bfe8005c536131276d448ead33f1b583
@SQ     SN:chr7 LN:158821424    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:74239c5ceee3b28f0038123d958114cb
@SQ     SN:chr8 LN:146274826    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:1eb00fe1ce26ce6701d2cd75c35b5ccb
@SQ     SN:chr9 LN:140273252    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:ea244473e525dde0393d353ef94f974b
@SQ     SN:chr10        LN:135374737    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:4ca41bf2d7d33578d2cd7ee9411e1533
@SQ     SN:chr11        LN:134452384    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:425ba5eb6c95b60bafbf2874493a56c3
@SQ     SN:chr12        LN:132349534    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:d17d70060c56b4578fa570117bf19716
@SQ     SN:chr13        LN:114142980    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:c4f3084a20380a373bbbdb9ae30da587
@SQ     SN:chr14        LN:106368585    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:c1ff5d44683831e9c7c1db23f93fbb45
@SQ     SN:chr15        LN:100338915    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:5cd9622c459fe0a276b27f6ac06116d8
@SQ     SN:chr16        LN:88827254     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:3e81884229e8dc6b7f258169ec8da246
@SQ     SN:chr17        LN:78774742     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:2a5c95ed99c5298bb107f313c7044588
@SQ     SN:chr18        LN:76117153     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:3d11df432bcdc1407835d5ef2ce62634
@SQ     SN:chr19        LN:63811651     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:2f1a59077cfad51df907ac25723bff28
@SQ     SN:chr20        LN:62435964     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:f126cdf8a6e0c7f379d618ff66beb2da
@SQ     SN:chr21        LN:46944323     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:f1b74b7f9f4cdbaeb6832ee86cb426c6
@SQ     SN:chr22        LN:49691432     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:2041e6a0c914b48dd537922cca63acb8
@SQ     SN:chrX LN:154913754    UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:d7e626c80ad172a4d7c95aadb94d9040
@SQ     SN:chrY LN:57772954     UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:62f69d0e82a12af74bad85e2e4a8bd91
@SQ     SN:chr1_random  LN:1663265      UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:cc05cb1554258add2eb62e88c0746394
@SQ     SN:chr2_random  LN:185571       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:18ceab9e4667a25c8a1f67869a4356ea
@SQ     SN:chr3_random  LN:749256       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:9cc571e918ac18afa0b2053262cadab6
@SQ     SN:chr4_random  LN:842648       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:9cab2949ccf26ee0f69a875412c93740
@SQ     SN:chr5_random  LN:143687       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:05926bdbff978d4a0906862eb3f773d0
@SQ     SN:chr6_random  LN:1875562      UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:d62eb2919ba7b9c1d382c011c5218094
@SQ     SN:chr7_random  LN:549659       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:28ebfb89c858edbc4d71ff3f83d52231
@SQ     SN:chr8_random  LN:943810       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:0ed5b088d843d6f6e6b181465b9e82ed
@SQ     SN:chr9_random  LN:1146434      UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:1e3d2d2f141f0550fa28a8d0ed3fd1cf
@SQ     SN:chr10_random LN:113275       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:50be2d2c6720dabeff497ffb53189daa
@SQ     SN:chr11_random LN:215294       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:bfc93adc30c621d5c83eee3f0d841624
@SQ     SN:chr13_random LN:186858       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:563531689f3dbd691331fd6c5730a88b
@SQ     SN:chr15_random LN:784346       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:bf885e99940d2d439d83eba791804a48
@SQ     SN:chr16_random LN:105485       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:dd06ea813a80b59d9c626b31faf6ae7f
@SQ     SN:chr17_random LN:2617613      UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:34d5e2005dffdfaaced1d34f60ed8fc2
@SQ     SN:chr18_random LN:4262 UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:f3814841f1939d3ca19072d9e89f3fd7
@SQ     SN:chr19_random LN:301858       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:420ce95da035386cc8c63094288c49e2
@SQ     SN:chr21_random LN:1679693      UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:a7252115bfe5bb5525f34d039eecd096
@SQ     SN:chr22_random LN:257318       UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:4f2d259b82f7647d3b668063cf18378b
@SQ     SN:chrX_random  LN:1719168      UR:file:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/Homo_sapiens_assembly18.fasta      M5:f4d71e0758986c15e5455bf3e14e5d6f</code class="pre_md"></pre>
<h3>Creating the fasta index file</h3>
<p>We use the faidx command in samtools to prepare the fasta index file. This file describes byte offsets in the fasta file for each contig, allowing us to compute exactly where a particular reference base at contig:pos is in the fasta file.</p>
<pre><code class="pre_md">&gt; samtools faidx Homo_sapiens_assembly18.fasta 
108.446u 3.384s 2:44.61 67.9%   0+0k 0+0io 0pf+0w</code class="pre_md"></pre>
<p>This produces a text file with one record per line for each of the fasta contigs. Each record is of the: contig, size, location, basesPerLine, bytesPerLine. The index file produced above looks like:</p>
<pre><code class="pre_md">&gt; cat Homo_sapiens_assembly18.fasta.fai 
chrM    16571   6       50      51
chr1    247249719       16915   50      51
chr2    242951149       252211635       50      51
chr3    199501827       500021813       50      51
chr4    191273063       703513683       50      51
chr5    180857866       898612214       50      51
chr6    170899992       1083087244      50      51
chr7    158821424       1257405242      50      51
chr8    146274826       1419403101      50      51
chr9    140273252       1568603430      50      51
chr10   135374737       1711682155      50      51
chr11   134452384       1849764394      50      51
chr12   132349534       1986905833      50      51
chr13   114142980       2121902365      50      51
chr14   106368585       2238328212      50      51
chr15   100338915       2346824176      50      51
chr16   88827254        2449169877      50      51
chr17   78774742        2539773684      50      51
chr18   76117153        2620123928      50      51
chr19   63811651        2697763432      50      51
chr20   62435964        2762851324      50      51
chr21   46944323        2826536015      50      51
chr22   49691432        2874419232      50      51
chrX    154913754       2925104499      50      51
chrY    57772954        3083116535      50      51
chr1_random     1663265 3142044962      50      51
chr2_random     185571  3143741506      50      51
chr3_random     749256  3143930802      50      51
chr4_random     842648  3144695057      50      51
chr5_random     143687  3145554571      50      51
chr6_random     1875562 3145701145      50      51
chr7_random     549659  3147614232      50      51
chr8_random     943810  3148174898      50      51
chr9_random     1146434 3149137598      50      51
chr10_random    113275  3150306975      50      51
chr11_random    215294  3150422530      50      51
chr13_random    186858  3150642144      50      51
chr15_random    784346  3150832754      50      51
chr16_random    105485  3151632801      50      51
chr17_random    2617613 3151740410      50      51
chr18_random    4262    3154410390      50      51
chr19_random    301858  3154414752      50      51
chr21_random    1679693 3154722662      50      51
chr22_random    257318  3156435963      50      51
chrX_random     1719168 3156698441      50      51</code class="pre_md"></pre>