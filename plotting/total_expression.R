args = commandArgs(trailingOnly=TRUE)

norm_results = args[1]
output = args[2]
name = args[3]
label = args[4]
###############################################################################

library(ggplot2)
library(data.table)

tpm = fread(norm_results)

proS = split(tpm, list(tpm$Peptide, tpm$Peptide_Type), drop = TRUE)

proSplus = Reduce(function(...) rbind(...), proS)
g = ggplot(proSplus, aes(x = Sample, y = Peptide, fill = value))
g = g + geom_tile(width = 0.75, height = 0.75, size = 0.3)
g = g + scale_color_manual(values=c("grey", "black"))
g = g + scale_fill_gradient(low = 'snow1', high = 'steelblue', na.value = 'white')
g = g + facet_grid(Peptide_Type ~ ., scales = 'free', space = 'free')
g = g + theme_bw() + xlab('') + ylab('')
g = g + guides(fill = guide_colourbar(title = label))
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

total_peptides = length(unique(proSplus$Peptide))

if (total_peptides > 50){
  g = g + theme(axis.text.y = element_text(size = 5))
}
if (total_peptides > 70){
  g = g + theme(axis.text.y = element_text(size = 3))
}

filename = sprintf('%s.%s', name, 'pdf') 
ggsave(sprintf('%s/%s', output, filename), width = 11, height = 11, useDingbats = FALSE)
file.remove(norm_results)

#https://r-charts.com/correlation/heat-map-ggplot2/
