args = commandArgs(trailingOnly=TRUE)

norm_results = args[1]
output = args[2]
th_out = as.numeric(args[3])
name = args[4]

###############################################################################

library(ggplot2)
library(data.table)


norm_info = fread(norm_results)

# Peptide Peptide_type  Tissue  Tissue_type Short_list  median  mean
proS = split(norm_info, list(norm_info$Peptide, norm_info$Peptide_Type), drop = TRUE)

proSplus = lapply(proS, function(x){

  # log on mean
  nbTissue = sum(x$mean > log10(1.0))
  nbTissueShort = sum(x$mean > log10(1.0) & x$Short_list == 'yes')
  superMeanTissue = mean(x$mean)
  superMeanTissueShort = mean(x[x$Short_list == 'yes', mean])

  nbTissue1 = sum(x$mean > log10(1.0 + 1.0) & x$Short_list == 'yes')
  nbTissue5 = sum(x$mean > log10(5.0 + 1.0) & x$Short_list == 'yes')
  nbTissue10 = sum(x$mean > log10(10.0 + 1.0) & x$Short_list == 'yes')
  nbTissue15 = sum(x$mean > log10(15.0 + 1.0) & x$Short_list == 'yes')
  nbTissue20 = sum(x$mean > log10(20.0 + 1.0) & x$Short_list == 'yes')
  nbTissue25 = sum(x$mean > log10(25.0 + 1.0) & x$Short_list == 'yes')

  nbTissueMed = sum(x$median > log10(1.0))
  nbTissueMedShort = sum(x$median > log10(1.0) & x$Short_list == 'yes')
  superMedTissue = median(x$median)
  superMedTissueShort = median(x[x$Short_list == 'yes', median])

  nbTissue1Med = sum(x$median > log10(1.0 + 1.0) & x$Short_list == 'yes')
  nbTissue5Med = sum(x$median > log10(5.0 + 1.0) & x$Short_list == 'yes')
  nbTissue10Med = sum(x$median > log10(10.0 + 1.0) & x$Short_list == 'yes')
  nbTissue15Med = sum(x$median > log10(15.0 + 1.0) & x$Short_list == 'yes')
  nbTissue20Med = sum(x$median > log10(20.0 + 1.0) & x$Short_list == 'yes')
  nbTissue25Med = sum(x$median > log10(25.0 + 1.0) & x$Short_list == 'yes')


  cols = c('nbTissue', 'nbTissueShort',
           'superMeanTissue', 'superMeanTissueShort',
           'nbTissue1short', 'nbTissue5short', 'nbTissue10short',
           'nbTissue15short', 'nbTissue20short', 'nbTissue25short',
           'nbTissueMed', 'nbTissueMedShort',
           'superMedTissue', 'superMedTissueShort',
           'nbTissue1Medshort', 'nbTissue5Medshort', 'nbTissue10Medshort',
           'nbTissue15Medshort', 'nbTissue20Medshort', 'nbTissue25Medshort')

  x[,(cols):= list(nbTissue, nbTissueShort,
                   superMeanTissue, superMeanTissueShort,
                   nbTissue1, nbTissue5, nbTissue10,
                   nbTissue15, nbTissue20, nbTissue25,
                   nbTissueMed, nbTissueMedShort,
                   superMedTissue, superMedTissueShort,
                   nbTissue1Med, nbTissue5Med, nbTissue10Med,
                   nbTissue15Med, nbTissue20Med, nbTissue25Med), by = x]
  return(x)
})

proSplusMerged = Reduce(function(...) rbind(...), proSplus)


# All tissues
filename = sprintf('%s_%s', name, 'all_tissues.pdf') 

label = sprintf('%s%s%s', 'mean > log10(', th_out, '+ 1)' ) 

g = ggplot(proSplusMerged, aes(x = Tissue, y = reorder(Peptide, + nbTissue), fill = mean, color = as.factor(mean > log10(th_out + 1))) )

g = g + labs(col = label) 
g = g + geom_tile(width = 0.75, height = 0.75, size = 0.3)
g = g + scale_color_manual(values=c("grey", "black"))
g = g + scale_fill_gradient(low = 'snow1', high = 'steelblue', na.value = 'white')
#g = g + scale_fill_gradient(limits = c(0,3), low = 'snow1', high = 'steelblue', na.value = 'white')
#g = g + scale_fill_gradient(limits = c(0,2.5), colors=c('snow1', 'steelblue'))
g = g + facet_grid(Peptide_Type ~ Tissue_type, scales = 'free', space = 'free')
g = g + theme_bw() + xlab('') + ylab('')
g = g + guides(fill = guide_colourbar(title = label))
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

total_peptides = length(unique(proSplusMerged$Peptide))


if (total_peptides > 50){
  g = g + theme(axis.text.y = element_text(size = 5))
}
if (total_peptides > 70){
  g = g + theme(axis.text.y = element_text(size = 3))
}
ggsave(sprintf('%s/%s', output, filename), width = 11, height = 11, useDingbats = FALSE)


# Selected tissues

filename = sprintf('%s_%s', name, 'selected_tissues.pdf') 

g = ggplot(proSplusMerged[proSplusMerged$Short_list == 'yes', ],
           aes(x = Tissue, y = reorder(Peptide, + nbTissue15short), fill = mean,
               color = as.factor(mean > log10(th_out + 1))))
g = g + labs(col = label)
g = g + geom_tile(width = 0.75, height = 0.75, size = 0.3)
g = g + scale_color_manual(values=c("grey", "black"))
g = g + scale_fill_gradient(low = 'snow1', high = 'steelblue', na.value = 'white')
#g = g + scale_fill_gradient(limits = c(0,3), low = 'snow1', high = 'steelblue', na.value = 'white')
g = g + facet_grid(Peptide_Type ~ Tissue_type, scales = 'free', space = 'free')
g = g + theme_bw() + xlab('') + ylab('')
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

if (total_peptides > 50){
  g = g + theme(axis.text.y = element_text(size = 5))
}
if (total_peptides > 70){
  g = g + theme(axis.text.y = element_text(size = 3))
}
ggsave(sprintf('%s/%s', output, filename), width = 11, height = 11, useDingbats = FALSE)

