# WBR-microbiome

"Are there differences between microbial communities on woodchip surfaces and woodchip interiors within woodchip bioreactors?"

Authors: Katie L. Duggan*, James P. Shapleigh, John M. Regan, Matthew C. Reid, M. Todd Walter

*Cornell University Biological & Environmental Engineering; kld93@cornell.edu

Abstract: Excessive amounts of nitrogen (N) and phosphorus (P) can lead to eutrophication in water sources. Woodchip bioreactors have shown success in removing N from agricultural runoff, but less data has been collected regarding P removal. Woodchip bioreactors are subsurface basins filled with woodchips installed downgradient of agricultural land to collect and treat drainage runoff. Microorganisms use the woodchips as a carbon source to transform N and P that are in the runoff. This study aims to explore microbial communities present in the bioreactor, their effect on both N and P cycling, and determine whether milling woodchips to reveal the interior material reveals hidden microbial diversities or activities. Illumina sequencing, metagenomic analyses, and sequence aligning were performed on six woodchip samples (i.e., three unmilled and three milled). All samples had similar DNA purity and yield, read quality, and microbial taxonomy regardless of milling. However, milled woodchips resulted in fewer high-quality genomes, suggesting that milling may cause shorter and more difficult to assemble fragments. In addition, when sequences were aligned against various genes and enzymes, some were more abundant in milled samples versus unmilled.  Our results indicated greater genomic functionality for denitrification and P transformations on the outside of the woodchips (unmilled) while the interior of woodchips (milled) exhibited more functional-gene-abundance for carbohydrate break down.  Thus, it may be important to open woodchip “black box” to ascertain where the functional activity is occurring when studying these types of systems.

Folder structure:
- data = fastq files for each sample and metadata
- analysis = RMarkdown files where analysis takes place
- figures = any figures created during analysis will be saved here

Fastq files (KBase export names):
- sample7G_trim_cut_117450_91_1.FASTQ
- sample8W_trim_cut_117450_100_1.FASTQ
- sample9G_trim_cut_117450_98_1.FASTQ
- sample10W_trim_cut_117450_94_1.FASTQ
- sample11G_trim_cut_117450_102_1.FASTQ
- sample12W_trim_cut_117450_96_1.FASTQ

Sample name key:
- W = "whole" (aka unmilled)
- G = "ground" (aka milled)
- "trim_cut" indicates these are the samples exported from KBase AFTER cleaning/trimming steps



