#get the bin burden, and the burden for each accession at different deleterious threshould
Rscript Burden_DifferentCutOff_GERP.R

#plot the distribution of bin burden
F2_BinBurdenDistribution.R

#plot the correlation between burden and phenotype
Plot_Burden_Phenotype_cor.R


#estimated the h by differnet deleterious threshould
Rscript GP.r2_h_GERP.cutoff.R

##genomic prediction
ModelEffect_greml_HomoHeterFixEffect.R

##stat permuattion result
GP_permutation_cat_plot.R