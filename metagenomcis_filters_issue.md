# Issue with metagenomics filters (context)

This file includes some suggestions for adjustments of features related to the metagenomics filters pipeline from viralunity  v1.3.1 (https://github.com/InstitutoTodosPelaSaude/ViralUnity). Context is: the REVISA table includes results from multiple metagenomic sequencing runs from clinical samples. During the project, we learned, through the usage of multiple negative controls, that the lab was contaminated with libraries from mayaro and HIV, related to other sequencing projects. So this provides an empirical assessment on the usage of such filters. Ocasionally, they won't be able to catch all contaminants, by they should get most. That is what I perceived in my analysis of these results. The adjustments to the filters below are simply a matter of trying to improve the output of the pipeline. 

## Bleed filter 

The bleed threshold used is in RPM, not in RPKM. Whenever RPKM is available, this should be used as the key measure for both the bleed filter and the negative control filters. Apparently RPKM is being used only for the negative control filter, while RPM is used for the bleed. 

## Negative control filters

Beyond what is currently implemented, I also believe it could be potentially useful to have an additional filter with the aggregated signal of a virus in the negative control. This should be based on aggregating read counts across all controls, and treat as a single control, being the basis for an additional (fold-enrichment, log2 filter calculation). Reasons why I think this could be potentially useful: I am seing (empirically) that some viruses cause widespread contamination across the sequencing run, but with great variance (violating normality). For instance, I see that some controls have 2.000> rpkm, while others have 200. I am imagining a scenario where this variance leads to a reduction in thresholds (biasing z-scores), masking contamination signals. The idea is to have a complementary (not substitute) approach that considers the aggregate signal for a virus across the full negative control, without averaging. 


## Fold-enrichment a z-score pass tags

I believe fold enrichment tags (columns) like fold_enrichment_10x_pass and fold_enrichment_100x_pass with boolean values would be useful to flag passing detections. The same goes for assessments of z-score under alternative thresholds neg_pass_5 or neg_pass_10

## Shift from log2 to log10 

At this point, I am not convinced that looking at the log2 ratio is very informative, as it harbors nearly the same information from the fold enrichment. I would pivot it to log10, because at least it brings a bit more information, as log10 is more directly interpretable. In this scenario, the default log ratio threshold should really be a value that marks a 10x difference in fold enrichment in absolute scale.

## Inclusion of final taxonomy

On the right side of nr_correct_species, I believe it would be useful to have a final_species column, which brings information from the column name, and when nr_correct_species contains information, these would be filled in. In this way, we would have a single column with all the correct (confirmed) taxonomy for species across the dataset. 


## Alternative cheap filter

While going through these results I could think of yet another cheap filter to implement, based on the size the largest contig for a taxa. Say I am handling with a virus 10kb long, and the largest assembled contig has 5.000 bp and median sequencing depth of 1000x. This tells me that, for that specific taxa, at least half the genome is covered at a good depth without further heavy computing and processing, as previous outputs on the pipeline already generate said information. So the idea is to include two columns, indicating the size of the largest contig and its median sequencing depth. This should only be possible when --viral-genomes option is provided (as the RPKM calculation). 

## Alternative taxonomic rank assessments

At this point, I am convinced that the in the vast majority of occasions, users will be primarily interested in the species rows. Primarily, I implemented these for users interested in filters operating at other taxonomic ranks. I think it would be interesting to move all summaries to a folder in the outputs within each method (e.g. diamond_contigs, kraken2 contigs) and leave one sub-folder per taxonomic level (family,genus,species). At the same time, information from other taxonomic ranks should be included in lower level tables. For instance, the species table should include new columns with family and genus names, while the genus-level table should include a column for family name. This would make the output more simple to parse and inspect visually. 

## Cleaning up

Remove the source column, it contains only the path for the output, not strictly relevant for users. 

# Making sense out of everything

Inspect the docs and filter implementations and viralunity outputs com the REVISA project. Create a plan with step by step implementation of these adjustments. Work on a separate branch from main and organize the work in logical commits. Use any agents or skills at your disposal that you believe will be useful. Criteria for acceptance: all tests / CI pass. Execution of test runs with the toy dataset (sars-cov-2). Be careful, analysis of REVISA run 5 is still ongoing, so don't launch any analysis until that is securely finished. After it finishes, perform your toy data test. If everything works out correctly, execute the latest summarization steps in the REVISA actual data, for which all the heavy computing has already been performed. Ask as many questions as needed, don't try to figure it out everything on your own. At the end of your work, we are going to have v1.3.2. 

# Useful paths:

* REVISA (data and results): /home/gevop/projects/REVISA
* ViralUnity (repo): /home/gevop/projects/ViralUnity
* Selected REVISA data analysis outputs: /home/gevop/projects/REVISA/selected_sample_outputs