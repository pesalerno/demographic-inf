# demographic inference using momi2
The following is the workflow for setting up the input files and running the demographic inference analyses in momi2 for the Chimantá ddRAD and the channel Islands RADseq datasets. We will prep input files and run demographic inference models in the following order: 

1. *Stefania ginesi* from Chimantá - 3 populations

2. *Tepuihyla edelcae* from Chimantá and Auyan - 4 populations
3. *Pseudacris regilla* from Channel Islands and mainland - 4 populations
4. *Xantusia riversiana* from Channel Islands - 3 populations

## Preparing the input SNP matrix (`vcf`)


In order to run a proper demographic analysis, we should have a filtered SNP matrix which excludes: 

- highly linked loci
- loci potentially under selection



1. Filtering loci potentially under selection.
----- 


>Add here details on how to run PCAdapt and code used.

First, we **generated the whitelists of the adaptive loci** directly from the PCAdapt outputs for Stefania and Tepuihyla. In Tepuihyla, 308 out of 4751 loci were found to be putatively adaptive. In Stefania, 680 out of 8722 loci were found to be putatively adaptive. 


PCAdapt outputs the "order" of the loci rather than the IDs of the loci themselves, so to exclude the outliers we generate a blacklist like this:

	##Fetch identities of Tepuihyla adaptive loci:
	
	awk '{print $16,$17,$34,$54,$60,$76,$129,$147,$179,$200,$218,$245,$305,$306,$308,$318,$319,$332,$333,$348,$356,$361,$369,$370,$371,$383,$393,$416,$423,$435,$447,$488,$491,$552,$558,$561,$566,$568,$585,$591,$604,$621,$638,$647,$660,$679,$680,$685,$686,$689,$693,$696,$702,$723,$733,$749,$764,$776,$787,$791,$816,$847,$851,$855,$901,$915,$918,$923,$963,$970,$990,$996,$1005,$1037,$1047,$1093,$1098,$1118,$1139,$1155,$1157,$1161,$1164,$1172,$1197,$1230,$1248,$1257,$1273,$1276,$1285,$1337,$1361,$1385,$1408,$1416,$1444,$1446,$1463,$1470,$1477,$1481,$1485,$1506,$1524,$1545,$1574,$1596,$1610,$1681,$1692,$1702,$1719,$1726,$1730,$1755,$1763,$1768,$1795,$1798,$1809,$1815,$1816,$1823,$1867,$1880,$1881,$1882,$1903,$1906,$1933,$1949,$1961,$1973,$2005,$2017,$2020,$2057,$2074,$2123,$2176,$2215,$2274,$2291,$2320,$2333,$2356,$2378,$2382,$2383,$2391,$2421,$2422,$2440,$2483,$2486,$2492,$2496,$2498,$2504,$2525,$2532,$2533,$2556,$2582,$2588,$2589,$2621,$2630,$2634,$2638,$2655,$2670,$2675,$2701,$2724,$2755,$2769,$2781,$2810,$2813,$2849,$2872,$2900,$2907,$2911,$2956,$2970,$3000,$3007,$3014,$3032,$3033,$3041,$3072,$3075,$3111,$3130,$3132,$3226,$3236,$3238,$3243,$3244,$3257,$3263,$3277,$3290,$3294,$3297,$3311,$3328,$3337,$3346,$3380,$3381,$3393,$3398,$3404,$3416,$3430,$3434,$3439,$3442,$3445,$3481,$3484,$3490,$3527,$3537,$3574,$3583,$3592,$3594,$3607,$3617,$3618,$3627,$3635,$3642,$3660,$3698,$3716,$3724,$3727,$3745,$3777,$3778,$3787,$3801,$3827,$3835,$3837,$3882,$3890,$3896,$3903,$3914,$3948,$3952,$3961,$3984,$3996,$3998,$4005,$4080,$4105,$4125,$4144,$4168,$4193,$4198,$4212,$4216,$4220,$4226,$4240,$4243,$4245,$4257,$4267,$4286,$4337,$4377,$4388,$4412,$4423,$4431,$4437,$4450,$4452,$4459,$4464,$4481,$4509,$4512,$4514,$4577,$4586,$4588,$4618,$4653,$4656,$4718,$4728,$4731,$4733,$4749}' Tepuihyla.stru > Tep_white-adaptive-loci.txt

	##Fetch identities of Stefania adaptive loci: 
	awk '{print $37,$48,$61,$62,$67,$74,$99,$109,$117,$126,$149,$157,$180,$184,$188,$199,$200,$211,$215,$225,$230,$257,$263,$265,$286,$288,$299,$308,$418,$456,$484,$517,$521,$536,$539,$544,$581,$623,$635,$644,$649,$675,$706,$711,$720,$740,$759,$767,$772,$773,$792,$821,$822,$830,$875,$885,$893,$901,$911,$921,$923,$927,$932,$941,$944,$967,$983,$993,$1017,$1026,$1034,$1040,$1042,$1050,$1052,$1054,$1065,$1110,$1112,$1118,$1122,$1130,$1131,$1134,$1149,$1161,$1181,$1190,$1198,$1230,$1236,$1255,$1271,$1281,$1304,$1313,$1338,$1349,$1352,$1359,$1375,$1386,$1390,$1408,$1410,$1434,$1437,$1440,$1441,$1446,$1449,$1468,$1473,$1497,$1499,$1502,$1516,$1530,$1533,$1543,$1565,$1623,$1626,$1629,$1635,$1639,$1649,$1654,$1656,$1662,$1669,$1690,$1698,$1709,$1725,$1731,$1739,$1740,$1792,$1800,$1801,$1806,$1816,$1828,$1832,$1846,$1850,$1857,$1860,$1902,$1918,$1934,$1935,$1943,$1960,$2003,$2027,$2032,$2047,$2061,$2063,$2078,$2106,$2111,$2113,$2124,$2127,$2144,$2161,$2163,$2202,$2231,$2243,$2257,$2263,$2266,$2284,$2289,$2292,$2336,$2376,$2379,$2389,$2392,$2405,$2408,$2419,$2426,$2446,$2457,$2469,$2476,$2489,$2511,$2517,$2518,$2542,$2545,$2546,$2558,$2576,$2580,$2605,$2613,$2628,$2646,$2653,$2663,$2665,$2675,$2687,$2693,$2718,$2722,$2730,$2739,$2744,$2746,$2747,$2752,$2762,$2768,$2783,$2791,$2795,$2801,$2804,$2807,$2815,$2826,$2827,$2828,$2855,$2864,$2880,$2887,$2890,$2904,$2924,$2925,$2960,$2965,$2979,$2990,$2999,$3005,$3007,$3025,$3037,$3039,$3044,$3046,$3065,$3090,$3095,$3127,$3157,$3166,$3186,$3194,$3202,$3247,$3261,$3275,$3283,$3332,$3333,$3340,$3356,$3407,$3418,$3426,$3456,$3465,$3476,$3499,$3517,$3545,$3565,$3609,$3655,$3674,$3696,$3715,$3719,$3724,$3751,$3803,$3807,$3845,$3861,$3885,$3903,$3908,$3916,$3918,$3925,$3976,$3981,$4006,$4026,$4028,$4043,$4044,$4052,$4060,$4068,$4078,$4081,$4084,$4087,$4118,$4119,$4141,$4174,$4190,$4192,$4198,$4206,$4211,$4233,$4237,$4239,$4250,$4278,$4340,$4346,$4348,$4365,$4369,$4459,$4474,$4482,$4490,$4491,$4534,$4550,$4555,$4558,$4581,$4582,$4587,$4588,$4597,$4606,$4616,$4619,$4632,$4635,$4638,$4646,$4650,$4655,$4656,$4671,$4676,$4681,$4682,$4686,$4696,$4749,$4756,$4767,$4770,$4791,$4819,$4820,$4829,$4834,$4838,$4841,$4850,$4881,$4892,$4896,$4898,$4901,$4904,$4922,$4937,$4945,$4989,$5001,$5022,$5042,$5064,$5070,$5138,$5141,$5174,$5175,$5186,$5188,$5199,$5251,$5267,$5273,$5274,$5285,$5308,$5316,$5349,$5350,$5352,$5357,$5366,$5448,$5460,$5466,$5472,$5476,$5477,$5481,$5500,$5517,$5523,$5527,$5550,$5570,$5599,$5603,$5606,$5612,$5628,$5634,$5658,$5660,$5662,$5677,$5679,$5697,$5705,$5711,$5721,$5728,$5732,$5763,$5764,$5828,$5834,$5861,$5888,$5895,$5904,$5909,$5930,$5954,$5971,$5973,$5995,$5998,$6000,$6009,$6018,$6037,$6050,$6052,$6090,$6097,$6104,$6152,$6154,$6180,$6194,$6204,$6211,$6219,$6230,$6234,$6240,$6245,$6250,$6272,$6282,$6291,$6292,$6293,$6308,$6318,$6334,$6344,$6345,$6352,$6363,$6387,$6398,$6400,$6405,$6414,$6423,$6435,$6441,$6446,$6453,$6464,$6467,$6498,$6499,$6503,$6510,$6533,$6538,$6564,$6566,$6579,$6592,$6614,$6620,$6627,$6630,$6639,$6650,$6655,$6666,$6671,$6697,$6722,$6727,$6740,$6741,$6759,$6766,$6772,$6785,$6788,$6797,$6802,$6804,$6816,$6819,$6844,$6846,$6848,$6852,$6856,$6871,$6877,$6891,$6892,$6895,$6897,$6903,$6910,$6940,$6953,$6976,$6980,$6993,$7004,$7006,$7019,$7058,$7090,$7109,$7112,$7119,$7128,$7129,$7152,$7169,$7175,$7177,$7206,$7231,$7235,$7238,$7254,$7256,$7260,$7271,$7275,$7283,$7308,$7315,$7318,$7324,$7335,$7336,$7370,$7383,$7398,$7457,$7464,$7511,$7512,$7548,$7550,$7551,$7569,$7577,$7584,$7599,$7611,$7612,$7627,$7645,$7646,$7649,$7655,$7656,$7657,$7668,$7675,$7683,$7698,$7702,$7732,$7755,$7783,$7785,$7795,$7808,$7820,$7857,$7861,$7908,$7913,$7918,$7923,$7936,$7945,$7954,$7958,$7972,$7974,$7990,$7995,$7998,$8010,$8020,$8047,$8048,$8059,$8062,$8076,$8082,$8171,$8179,$8190,$8200,$8227,$8228,$8251,$8266,$8289,$8290,$8292,$8297,$8298,$8308,$8333,$8341,$8386,$8394,$8396,$8398,$8411,$8436,$8438,$8442,$8455,$8460,$8463,$8465,$8476,$8496,$8550,$8552,$8586,$8607,$8624,$8629,$8640,$8666,$8673,$8674,$8676,$8699,$8714,$8721}' Stefania.stru > Stef_white-adaptive-loci.txt


We use the PCAdapt outlier loci (keeping both the locus and the SNP ID, such as: 257_57) to generate a **blacklist of loci to exclude in plink**: 

	./plink --file Tepuihyla --exclude Tep_adaptive_blacklist.txt --recode --out Tep_neutral --noweb


We then used the ***.map*** output from ***plink*** in Text Wrangler, and generated the populations whitelist using find and replace arguments with **grep**:

	search for \d\t(\d*)_\d*\t\d\t\d*$
	replace with \1

Based the **.irem** file from the second iteration of *plink* we also removed poorly genotyped individuals (that did not pass plink filter) from the popmap to use in populations input (i.e. individuals with >50% missing data). 

	populations -b 1 -P ./denovo-03-2017 -M ./popmap-Tep.txt  -t 36 -p 1 -r 0.5 -W Tep-whitelist --write_random_snp --phylip --structure --plink --vcf --genepop --fstats

Here are the final *Stefania*  [neutral](https://github.com/pesalerno/Chimanta-genomics/blob/master/Stef_neutral.stru) matrix, and the *Tepuihyla*  [neutral](https://github.com/pesalerno/Chimanta-genomics/blob/master/Tep_neutral.stru) matrix. 

---------------------

2. **Filtering loci under Linkage Disequilibrium** 
---

Using the `.ped` and `.map` input files (which were transformed using `PGDSpider`, we can check for LD using plink with the following code: 

	plink --noweb --file mydata --r2 --out outputfilename 
	
	
The pairs of loci are then exported into excel and we filter out one of the two - or more - loci that are linked to one another.

To filter out excess spaces in the file, we run the short command:

	cat Stef-LD.ld | tr -s ' ' > LD-stef-b
 
Then we can decide which SNPs to filter visualizing the matrix in excel. I filtered all SNPs that were more than 80% associated in each matrix, which resulted in a much higher numner of loci excluded in *Tepuihyla* (1173 SNPs) than in *Stefania* (89 SNPs). 

In order to create a blacklist of loci in LD to exclude from the SNP matrix, we used find/replace arguments in TextWrangler to create a list of SNPs like this:

	59185_31
	539339_84
	52224_62
	49662_45
	348531_41
	##       ...etc

And then using this file (named **LD-loci-list.txt**) we generate our LE SNP matrix using *plink*: 

	./plink --file filename --exclude LD-loci-list.txt --recode --out LD-filtered
	

--------------------

3. **Transforming and exporting the vcf file**
----

> Add details here on transforming `.ped` to `.vcf` using `vcftools/bcftools`

**Properly formatted VCF.** momi2 expects very large input files, so it needs them to be preprocessed to speed things up. The details of this preprocessing step are not very interesting, but we are basically compressing and indexing the VCF so it’s faster to search.


First, we perform a blockwise compression using `bgzip`. The -c flag directs bgzip to leave the original vcf file untouched and create a new file for the vcf.gz
	
	cp /rigel/edu/radcamp/users/work1/ipyrad-workshop/anoles_outfiles/	anoles.vcf anolis.vcf
	bgzip -c anolis.vcf > anolis.vcf.gz

Second, we index the file for searching using `tabix`: 

	tabix anolis.vcf.gz
	ls
	anolis.alleles.loci  anolis.loci      anolis.snps.phy	anolis.u.snps.phy
	anolis.geno	     anolis.nex       anolis_stats.txt	anolis.ustr
	anolis.gphocs	     anolis.phy       anolis.str	anolis.vcf.gz
	anolis.hdf5	     anolis.snps.map  anolis.u.geno	anolis.vcf.gz.tbi

------------------

## Preparing the population assignment file


Based on the original population localities for *Stefania* and *Tepuihyla*, we assigned individual specimen samples to one of the three Chimantá populations - Abakapá (`AB`), Churí (`CH`), and Eruoda (`ER`), or to Auyantepui (`AU`; only for *Tepuihyla*): the “North” population, and 8 samples to the “South” population:

	punc_ICST764    AB
	punc_MUFAL9635  AB
	punc_IBSPCRIB0361       ER
	punc_JFT773     ER
	punc_MTR05978   ER
	punc_MTR17744   CH
	punc_MTR21545   CH
	punc_MTR34414   CH
	punc_MTRX1468   CH
	punc_MTRX1478   CH
	.
	#etc...

----
## Preparing other input files for analyses


>esto es para que Alex termine de rellenar/hacer de acuerdo al [tutorial de RADcamp](https://radcamp.github.io/Yale2019/07_momi2_API.html). 

**BED file**: This file specifies genomic regions to include in when calculating the SFS. It is composed of 3 columns which specify ‘chrom’, ‘chromStart’, and ‘chromEnd’.

**The allele counts file**: The allele counts file is an intermediate file that we must generate on the way to constructing the SFS. momi2 provides a function for this.

**Genereate the SFS**: The culmination of all this housekeeping is the SFS file which we will use for demographic inference.

----
## First set of demographic inference models


We will begin by running basic demographic models, which are based both on sampling/distributions of the species and on the results from previous analyses. Also, the setup of these models is based a lot on [this publication](https://github.com/pesalerno/demographic-inf/blob/master/files/passerines-demography.pdf).

>**Basic *Stefania* Demographic Models setup.** For *Stefania*, we have three sampled populations atop Chimantá, and we observe two incongruent inferred histories between `mtDNA` and `nuDNA`, as observed in the below graphs: 
>
>![](https://github.com/pesalerno/demographic-inf/blob/master/files/stefania-trees.png)

Thus, the demographic models to test for ***Stefania*** are: 


1. **no-migration models:** here we will test various divergence scenarios for the three populations and assess AIC supports. The parameters to vary are time to divergence and population size (constant).


>model |	T1|	T2 |	N1 | N2 | N3 |
>-------------|------|------|-----|------|------|
>stef1_a |	T1|	T2|	N1 | N2 | N3 |
stef1_b |	T1|	T2|	N1 | N2 | N3 |
stef1_c |	T1|	T2|	N1 | N2 | N3 |
stef1_d |	T1|	T2|	N1 | N2 | N3 |
stef1_e |	T1|	T2|	N1 | N2 | N3 |
stef1_f |	T1|	T2|	N1 | N2 | N3 |
stef1_g |	T1|	T2|	N1 | N2 | N3 |

2. **migration models:** based on the previous results of best models, we will test various migration scenarios for the three populations and assess AIC supports. The parameters to vary are time to divergence, population size (constant), and migrations among populations.

>model |	T1|	T2 |	N1 | N2 | N3 | m1 | m2 | m3 |
>-------------|------|------|-----|------|------|------|------|------- |
stef2_a |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |
stef2_b |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |
stef2_c |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |
stef2_d |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |
stef2_e |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |
stef2_f |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |
stef2_g |	T1|	T2|	N1 | N2 | N3 | m1 | m2 | m3 |

--------------------------
>**Basic *Tepuihyla* Demographic Models setup.** For *Tepuihyla*, we have the same three sampled populations atop Chimantá, plus an additional sampled population on the neighboring formation of Auyantepui, and we mostly observe very little information/signal from `mtDNA` when compared to `nuDNA`, as observed in the below graphs: 
>
>![](https://github.com/pesalerno/demographic-inf/blob/master/files/tepuihyla-trees.png)