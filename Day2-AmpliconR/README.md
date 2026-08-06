# Day 2: Symbiosis and Amplicon analysis using Qiime2

## Morning Lecture (9am - 10am):
Studying Symbiosis through Sequencing (Liana Burghardt) - [pdf](MicrobiomeKickStart_Burghardt_Aug2026.pdf)

## Amplicon Analysis using Qiime2 (10am - 4pm): 
### Resources

- [Qiime2](https://amplicon-docs.qiime2.org/en/stable/)

- [Qiime2 view](https://view.qiime2.org/)

- [Guide for amplicon analysis](day2_qiime2-workflow.md)

- [mini lecture](Day2_amplicon-VRR.pdf)

### Analysis

1. Log on to the portal
http://portal.hpc.psu.edu/ 

2. Open the terminal; to do so, click on ">_Open in Terminal"

🖥️ in terminal

3. Request an interactive session in the terminal (here we are requesting 3hours):  

```
salloc -N 1 -n 2 -t 4:00:00 --account=one_sc_default
```
Wait for the terminal prompts that indicate you have the requested resources.

4. Load the qiime2 environment

```
module load anaconda/2023.09
conda activate /storage/group/one/default/workshop/2026/envs/rachis-qiime2-2026.7
```
5. Create a variable (DAY2) to make the commands shorter and list the files in the variable DAY2

```
DAY2=/storage/group/one/default/workshop/2026/day2_16S/qiime2
ls $DAY2
```

6. Follow the commands in the Qiime2 guide (day2_qiime2_workflow.md). Please run the scripts interactively, meaning copy and paste from the guide into the terminal.

Running the script interactively helps you spot errors and follow along. After you’ve tested everything, you can execute all commands from a single script (this file is on Roar and has the extension `*.slurm*`).

### Visualization

For today’s session, we’ll use the Qiime2 visualization files (`.qzv`). However, for greater reproducibility and flexibility in plotting, it is recommended to download the underlying feature tables and perform the visualization in R using packages such as **phyloseq** and **vegan**.

🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹🔹
   
## 🔹Extra Material🔹

### **Suggested Reading**

Weinroth MD, Belk AD, Dean C, et al. Considerations and best practices in animal science 16S ribosomal RNA gene sequencing microbiome studies. J Anim Sci. 2022;100(2):skab346. doi:10.1093/jas/skab346 [Link](https://academic.oup.com/jas/article/100/2/skab346/6519592)

Mirzayi C, Renson A; Genomic Standards Consortium; Reporting guidelines for human microbiome research: the STORMS checklist Nat Med. 2021 Nov;27(11):1885-1892. doi: 10.1038/s41591-021-01552-x. Epub 2021 Nov 17. PMID: 34789871; PMCID: PMC9105086.[Link](https://www.nature.com/articles/s41591-021-01552-x)

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. 2016;13(7):581-583. doi:10.1038/nmeth.3869 [Link](https://www.nature.com/articles/nmeth.3869)

Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome Datasets Are Compositional: And This Is Not Optional. Front Microbiol. 2017;8:2224. Published 2017 Nov 15. doi:10.3389/fmicb.2017.02224 [Link](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full)

### **Best practices for differential abundance analysis**

Article: [Beyond Normalization: Incorporating Scale Uncertainty in Microbiome and Gene Expression Analysis](https://www.biorxiv.org/content/10.1101/2024.04.01.587602v1.full)

### **Online "textbook" on fundamental concepts of microbiome**
[Meet the Metaorganism](https://metaorganism.app/en)

### **Online educational game on fundamental concepts of microbiome**
[Tiny Biome Tales](https://microbiome.gamelabgraz.at/)
Article: [‘Tiny biome tales’: playing a game to understand the human microbiome](https://www.the-microbiologist.com/news/tiny-biome-tales-playing-a-game-to-understand-the-human-microbiome/3858.article)

### **Relevant TED talks**
Laura Boykin: [How we're using DNA tech to help farmers fight crop diseases](https://www.ted.com/talks/laura_boykin_how_we_re_using_dna_tech_to_help_farmers_fight_crop_diseases?utm_campaign=tedspread&utm_medium=referral&utm_source=tedcomshare)
Steffanie Strathdee: [How Sewage Saved My Husband's Life from a Superbug](https://youtu.be/AbAZU8FqzX4?si=JDqz6OFD1Vb9Y0Wu) and [How does a virus eat bacteria?](https://www.cnn.com/videos/tv/2022/07/08/steffanie-strathdee-life-itself-wellness.cnn)
Rob Knight: [How our microbes make us who we are](https://www.ted.com/talks/rob_knight_how_our_microbes_make_us_who_we_are?utm_campaign=tedspread&utm_medium=referral&utm_source=tedcomshare)


### **Suggested Resources**

#### Background: Coursera - gut check
It is great to start grappling with concepts - videos are short<br>
https://www.coursera.org/lecture/microbiome/welcome-to-the-course-AuMZk

### Coding:

#### 🧑‍💻 The best way to get started with bioinformatics
For getting used to best practices in coding in the command line<br>
https://www.biostarhandbook.com/index.html

##### 🖥️ Getting started with R
1.	https://www.codecademy.com/learn/learn-r
2.	https://www.codecademy.com/search?query=R
   
##### 💻 Getting started with BASH
Before you do biostars handbook do this quick shell one
https://swcarpentry.github.io/shell-novice/

##### ✅ Scientific rigor, reproducibility, version control, and code-alongs
Riffomonas has great tutorials to get started with R and reproducible research with version control is also a must.
https://riffomonas.org/

##### 🧬 16S rRNA analysis
**Analysis in R**
https://benjjneb.github.io/dada2/tutorial.html

**Analysis in Qiime2**
https://amplicon-docs.qiime2.org/en/stable/tutorials/moving-pictures/
