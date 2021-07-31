# Association-Between-Environment-Genome-and-Morphometry-in-marsh-rats
Scripts and input files used in the paper submitted to the Biological Journal of the Linnean Society: 

Prado et al. Biome effect on phenotypic and genomic differences in the South America marsh rats.

(1) Niche Similarity

Environmental variables can be download from WorldClim (Hijmans et al., 2005) and ENVIREM (Title & Bemmels, 2018). 

- 1-Select&Edit_env_variables.R
  Select and edit environmental variables for reduction of environmental self-correlation

- 2-Rarefy_points.R
  Filter occurrence points for reduction of spatial self-correlation

- 3-ModelComplexity.R
  Select the best parameter combination to achieve a balance between goodness‐of‐fit and model complexity

- 4-Cut_model_by_threshold.R
  Construct binary maps based on maxent's threshold
  
- 5-Niche_Similarity.R
  Characterize the ecological niche similarity applying two different approaches: (1) the multivariate environmental niche overlap quantified with the ‘PCA-env’ (Broennimann et al., 2012), and (2) the estimation of n-dimensional environmental hypervolume (Hutchinson, 1957; Blonder et al., 2014). 
 
(2) Genomic data and Bioinformatics

- 6-Create_whitelist_stacks2e.R
  Create a whitelist to remove SNP’s positioned at the end of all loci due to an artificially increased number of SNPs observed at these last positions. Also, it removes loci with high theta values (above the 95 percentiles), given these are suggestive of sequencing and assembly errors. 
  
- 7-sumstat.R
  Check the genetic diversity summary statistics for significant differences, and plots their means values and standard errors.
  
-8-Pop_structure.R
  Perform Principal Component Analysis (PCA), Tracy-Widom test, and Mantel test 
  
-9-moving&best_lhood.R
 Move files among folders to facilitate Fastsimcoal analysis
  
  

