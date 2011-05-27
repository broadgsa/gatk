-- This script parses the table of genes from  http://geneticassociationdb.nih.gov/
-- and generates an annotated interval list for use with the GATK
--
-- Author: carneiro
-- Date: 5/24/2011


ref_table = {
  "ID",
  "AS",
  "PH",
  "DI",
  "DC",
  "DT",
  "CH",
  "CB",
  "GE",
  "ST",
  "SP",
  "PV",
  "RE",
  "PI",
  "AA",
  "AF",
  "PC",
  "GN",
  "RS",
  "PO",
  "GO",
  "SU",
  "LN",
  "UN",
  "NP",
  "MP",
  "JO",
  "TI",
  "RN",
  "OM",
  "YR",
  "CN",
  "SI",
  "EF",
  "GIGA",
  "GIAA",
  "GIGB",
  "GIAB",
  "GIGC",
  "GIAC",
  "GIAS",
  "GIEF"
}

chr_limits = {
  249250621,
  243199373,
  198022430,
  191154276,
  180915260,
  171115067,
  159138663,
  146364022,
  141213431,
  135534747,
  135006516,
  133851895,
  115169878,
  107349540,
  102531392,
  90354753,
  81195210,
  78077248,
  59128983,
  63025520,
  48129895,
  51304566,
  155270560,
  59373566,
  16569
}

local header = [[@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1b22b98cdeb4a9304cb5d48026a85128	SP:Homo Sapiens
@SQ	SN:2	LN:243199373	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:a0d9851da00400dec1098a9255ac712e	SP:Homo Sapiens
@SQ	SN:3	LN:198022430	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:fdfd811849cc2fadebc929bb925902e5	SP:Homo Sapiens
@SQ	SN:4	LN:191154276	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:23dccd106897542ad87d2765d28a19a1	SP:Homo Sapiens
@SQ	SN:5	LN:180915260	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:0740173db9ffd264d728f32784845cd7	SP:Homo Sapiens
@SQ	SN:6	LN:171115067	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1d3a93a248d92a729ee764823acbbc6b	SP:Homo Sapiens
@SQ	SN:7	LN:159138663	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:618366e953d6aaad97dbe4777c29375e	SP:Homo Sapiens
@SQ	SN:8	LN:146364022	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:96f514a9929e410c6651697bded59aec	SP:Homo Sapiens
@SQ	SN:9	LN:141213431	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:3e273117f15e0a400f01055d9f393768	SP:Homo Sapiens
@SQ	SN:10	LN:135534747	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:988c28e000e84c26d552359af1ea2e1d	SP:Homo Sapiens
@SQ	SN:11	LN:135006516	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:98c59049a2df285c76ffb1c6db8f8b96	SP:Homo Sapiens
@SQ	SN:12	LN:133851895	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:51851ac0e1a115847ad36449b0015864	SP:Homo Sapiens
@SQ	SN:13	LN:115169878	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:283f8d7892baa81b510a015719ca7b0b	SP:Homo Sapiens
@SQ	SN:14	LN:107349540	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:98f3cae32b2a2e9524bc19813927542e	SP:Homo Sapiens
@SQ	SN:15	LN:102531392	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:e5645a794a8238215b2cd77acb95a078	SP:Homo Sapiens
@SQ	SN:16	LN:90354753	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:fc9b1a7b42b97a864f56b348b06095e6	SP:Homo Sapiens
@SQ	SN:17	LN:81195210	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:351f64d4f4f9ddd45b35336ad97aa6de	SP:Homo Sapiens
@SQ	SN:18	LN:78077248	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:b15d4b2d29dde9d3e4f93d1d0f2cbc9c	SP:Homo Sapiens
@SQ	SN:19	LN:59128983	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1aacd71f30db8e561810913e0b72636d	SP:Homo Sapiens
@SQ	SN:20	LN:63025520	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:0dec9660ec1efaaf33281c0d5ea2560f	SP:Homo Sapiens
@SQ	SN:21	LN:48129895	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:2979a6085bfe28e3ad6f552f361ed74d	SP:Homo Sapiens
@SQ	SN:22	LN:51304566	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:a718acaa6135fdca8357d5bfe94211dd	SP:Homo Sapiens
@SQ	SN:X	LN:155270560	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:7e0e2e580297b7764e31dbc80c2540dd	SP:Homo Sapiens
@SQ	SN:Y	LN:59373566	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1fa3474750af0948bdf97d5a0ee52e51	SP:Homo Sapiens
@SQ	SN:MT	LN:16569	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:c68f52674c9fb33aef52dcf399755519	SP:Homo Sapiens
@SQ	SN:GL000207.1	LN:4262	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:f3814841f1939d3ca19072d9e89f3fd7	SP:Homo Sapiens
@SQ	SN:GL000226.1	LN:15008	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1c1b2cd1fccbc0a99b6a447fa24d1504	SP:Homo Sapiens
@SQ	SN:GL000229.1	LN:19913	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:d0f40ec87de311d8e715b52e4c7062e1	SP:Homo Sapiens
@SQ	SN:GL000231.1	LN:27386	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:ba8882ce3a1efa2080e5d29b956568a4	SP:Homo Sapiens
@SQ	SN:GL000210.1	LN:27682	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:851106a74238044126131ce2a8e5847c	SP:Homo Sapiens
@SQ	SN:GL000239.1	LN:33824	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:99795f15702caec4fa1c4e15f8a29c07	SP:Homo Sapiens
@SQ	SN:GL000235.1	LN:34474	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:118a25ca210cfbcdfb6c2ebb249f9680	SP:Homo Sapiens
@SQ	SN:GL000201.1	LN:36148	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:dfb7e7ec60ffdcb85cb359ea28454ee9	SP:Homo Sapiens
@SQ	SN:GL000247.1	LN:36422	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:7de00226bb7df1c57276ca6baabafd15	SP:Homo Sapiens
@SQ	SN:GL000245.1	LN:36651	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:89bc61960f37d94abf0df2d481ada0ec	SP:Homo Sapiens
@SQ	SN:GL000197.1	LN:37175	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:6f5efdd36643a9b8c8ccad6f2f1edc7b	SP:Homo Sapiens
@SQ	SN:GL000203.1	LN:37498	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:96358c325fe0e70bee73436e8bb14dbd	SP:Homo Sapiens
@SQ	SN:GL000246.1	LN:38154	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:e4afcd31912af9d9c2546acf1cb23af2	SP:Homo Sapiens
@SQ	SN:GL000249.1	LN:38502	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1d78abec37c15fe29a275eb08d5af236	SP:Homo Sapiens
@SQ	SN:GL000196.1	LN:38914	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:d92206d1bb4c3b4019c43c0875c06dc0	SP:Homo Sapiens
@SQ	SN:GL000248.1	LN:39786	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:5a8e43bec9be36c7b49c84d585107776	SP:Homo Sapiens
@SQ	SN:GL000244.1	LN:39929	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:0996b4475f353ca98bacb756ac479140	SP:Homo Sapiens
@SQ	SN:GL000238.1	LN:39939	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:131b1efc3270cc838686b54e7c34b17b	SP:Homo Sapiens
@SQ	SN:GL000202.1	LN:40103	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:06cbf126247d89664a4faebad130fe9c	SP:Homo Sapiens
@SQ	SN:GL000234.1	LN:40531	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:93f998536b61a56fd0ff47322a911d4b	SP:Homo Sapiens
@SQ	SN:GL000232.1	LN:40652	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:3e06b6741061ad93a8587531307057d8	SP:Homo Sapiens
@SQ	SN:GL000206.1	LN:41001	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:43f69e423533e948bfae5ce1d45bd3f1	SP:Homo Sapiens
@SQ	SN:GL000240.1	LN:41933	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:445a86173da9f237d7bcf41c6cb8cc62	SP:Homo Sapiens
@SQ	SN:GL000236.1	LN:41934	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:fdcd739913efa1fdc64b6c0cd7016779	SP:Homo Sapiens
@SQ	SN:GL000241.1	LN:42152	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:ef4258cdc5a45c206cea8fc3e1d858cf	SP:Homo Sapiens
@SQ	SN:GL000243.1	LN:43341	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:cc34279a7e353136741c9fce79bc4396	SP:Homo Sapiens
@SQ	SN:GL000242.1	LN:43523	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:2f8694fc47576bc81b5fe9e7de0ba49e	SP:Homo Sapiens
@SQ	SN:GL000230.1	LN:43691	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:b4eb71ee878d3706246b7c1dbef69299	SP:Homo Sapiens
@SQ	SN:GL000237.1	LN:45867	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:e0c82e7751df73f4f6d0ed30cdc853c0	SP:Homo Sapiens
@SQ	SN:GL000233.1	LN:45941	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:7fed60298a8d62ff808b74b6ce820001	SP:Homo Sapiens
@SQ	SN:GL000204.1	LN:81310	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:efc49c871536fa8d79cb0a06fa739722	SP:Homo Sapiens
@SQ	SN:GL000198.1	LN:90085	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:868e7784040da90d900d2d1b667a1383	SP:Homo Sapiens
@SQ	SN:GL000208.1	LN:92689	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:aa81be49bf3fe63a79bdc6a6f279abf6	SP:Homo Sapiens
@SQ	SN:GL000191.1	LN:106433	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:d75b436f50a8214ee9c2a51d30b2c2cc	SP:Homo Sapiens
@SQ	SN:GL000227.1	LN:128374	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:a4aead23f8053f2655e468bcc6ecdceb	SP:Homo Sapiens
@SQ	SN:GL000228.1	LN:129120	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:c5a17c97e2c1a0b6a9cc5a6b064b714f	SP:Homo Sapiens
@SQ	SN:GL000214.1	LN:137718	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:46c2032c37f2ed899eb41c0473319a69	SP:Homo Sapiens
@SQ	SN:GL000221.1	LN:155397	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:3238fb74ea87ae857f9c7508d315babb	SP:Homo Sapiens
@SQ	SN:GL000209.1	LN:159169	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:f40598e2a5a6b26e84a3775e0d1e2c81	SP:Homo Sapiens
@SQ	SN:GL000218.1	LN:161147	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:1d708b54644c26c7e01c2dad5426d38c	SP:Homo Sapiens
@SQ	SN:GL000220.1	LN:161802	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:fc35de963c57bf7648429e6454f1c9db	SP:Homo Sapiens
@SQ	SN:GL000213.1	LN:164239	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:9d424fdcc98866650b58f004080a992a	SP:Homo Sapiens
@SQ	SN:GL000211.1	LN:166566	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:7daaa45c66b288847b9b32b964e623d3	SP:Homo Sapiens
@SQ	SN:GL000199.1	LN:169874	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:569af3b73522fab4b40995ae4944e78e	SP:Homo Sapiens
@SQ	SN:GL000217.1	LN:172149	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:6d243e18dea1945fb7f2517615b8f52e	SP:Homo Sapiens
@SQ	SN:GL000216.1	LN:172294	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:642a232d91c486ac339263820aef7fe0	SP:Homo Sapiens
@SQ	SN:GL000215.1	LN:172545	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:5eb3b418480ae67a997957c909375a73	SP:Homo Sapiens
@SQ	SN:GL000205.1	LN:174588	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:d22441398d99caf673e9afb9a1908ec5	SP:Homo Sapiens
@SQ	SN:GL000219.1	LN:179198	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:f977edd13bac459cb2ed4a5457dba1b3	SP:Homo Sapiens
@SQ	SN:GL000224.1	LN:179693	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:d5b2fc04f6b41b212a4198a07f450e20	SP:Homo Sapiens
@SQ	SN:GL000223.1	LN:180455	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:399dfa03bf32022ab52a846f7ca35b30	SP:Homo Sapiens
@SQ	SN:GL000195.1	LN:182896	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:5d9ec007868d517e73543b005ba48535	SP:Homo Sapiens
@SQ	SN:GL000212.1	LN:186858	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:563531689f3dbd691331fd6c5730a88b	SP:Homo Sapiens
@SQ	SN:GL000222.1	LN:186861	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:6fe9abac455169f50470f5a6b01d0f59	SP:Homo Sapiens
@SQ	SN:GL000200.1	LN:187035	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:75e4c8d17cd4addf3917d1703cacaf25	SP:Homo Sapiens
@SQ	SN:GL000193.1	LN:189789	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:dbb6e8ece0b5de29da56601613007c2a	SP:Homo Sapiens
@SQ	SN:GL000194.1	LN:191469	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:6ac8f815bf8e845bb3031b73f812c012	SP:Homo Sapiens
@SQ	SN:GL000225.1	LN:211173	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:63945c3e6962f28ffd469719a747e73c	SP:Homo Sapiens
@SQ	SN:GL000192.1	LN:547496	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:325ba9e808f669dfeee210fdd7b470ac	SP:Homo Sapiens
@SQ	SN:NC_007605	LN:171823	AS:NC_007605.1	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:6743bd63b3ff2b5b8985d8933c53290a	SP:Epstein-Barr virus]]

local function comp(a, b)
  if a.chr == b.chr then return a.start < b.start end
  local x = tonumber(a.chr)
  local y = tonumber(b.chr)
  if x and y then
    return x < y end
  return a.chr < b.chr
end

local function cleanNumber (x)
  local t = x:gsub(",","")
  return tonumber(t)
end

local function cleanDisease (x)
  local d = x:gsub("[^%w%s]", "")
  return d:gsub(" ", "_")
end

local function isSameGene (g1, g2)
  return g1.chr == g2.chr and
         g1.start == g2.start and
         g1.finish == g2.finish
end

local function isValidGene(gene)
  local x
  if gene.chr == "X" then x = 23
  elseif gene.chr == "Y" then x = 24
  elseif gene.chr == "MT" then x = 25
  else x = tonumber(gene.chr) end

  return x >= 1 and x <= 25 and gene.start < gene.finish and chr_limits[x] > gene.finish
end


local function addGene(geneTable, gene)
  if isValidGene(gene) then  -- only adds valid genes
    local geneAdded = false
    for _, g in ipairs(geneTable) do
      if isSameGene(gene, g) then
         geneAdded = true
        if g.disease ~= gene.disease then
          g.disease = g.disease .. ","..gene.disease
        end
      end
    end
    if (not geneAdded) then
      table.insert(geneTable, gene)
    end
  end
end

local function cleanIntervals(f, t)
  for l in io.lines(f) do
    local counter = 1
    local id, c, b, e, disease
    for match in l:gmatch("(.-);;") do
  --    print("DEBUG: ", counter, match)
      if     counter == 1 then id = match
      elseif counter == 3 then disease = cleanDisease(match)
      elseif counter ==  7 then c = match
      elseif counter == 10 then b = cleanNumber(match)
      elseif counter == 11 then e = cleanNumber(match)
      end
      counter = counter + 1
    end
    if id and c and b and e and c~= "" then
      addGene(t, {id=id, chr=c, start=b, finish=e, disease=disease})
    end
  end
  return table.sort(t, comp)
end

local function printIntervals(t)
  print(header)
  for _,interval in ipairs(t) do
    print(interval.chr,interval.start,interval.finish, "+", interval.id..":::"..interval.disease)
  end
end

t = {}
cleanIntervals(arg[1], t)
printIntervals(t)

--[[
-- REFERENCE TABLE

1	ID
2	Association(Y/N)
3	Broad Phenotype
4	Disease Class
5	Disease Class Code
6	MeSH Disease Terms
7	Chromosom
8	Chr-Band
9	Gene
10	DNA Start
11	DNA End
12	P Value
13	Reference
14	Pubmed ID
15	Allele Author Description
16	Allele Functional Effects
17	Polymophism Class
18	Gene Name
19	RefSeq
20	Population
21	MeSH Geolocation
22	Submitter
23	Locus Number
24	Unigene
25	Narrow Phenotype
26	Mole. Phenotype
27	Journal
28	Title
29	rs Number
30	OMIM ID
31	Year
32	Conclusion
33	Study Info
34	Env. Factor
35	GI Gene A
36	GI Allele of Gene A
37	GI Gene B
38	GI Allele of Gene B
39	GI Gene C
40	GI Allele of Gene C
41	GI Association?
42	GI combine Env. Factor
--]]