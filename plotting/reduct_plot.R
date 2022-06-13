args = commandArgs(trailingOnly=TRUE)

norm_results = args[1]
output = args[2]

###############################################################################

library(ggplot2)
library(data.table)


norm_info = fread(norm_results)
filename = 'all_tissues.pdf'

label = 'mean > 0'

g = ggplot(norm_info, aes(x = Sample, y = Peptide, fill = mean, color = as.factor(mean > 0 )) )

g = g + geom_tile(width = 0.75, height = 0.75, size = 0.3)
g = g + scale_color_manual(values=c("FALSE"="grey", "TRUE"="black"))
g = g + guides(color = guide_legend(override.aes = list(fill = "white"))) 

g = g + scale_fill_gradient(low = 'snow1', high = 'steelblue', na.value = 'white')
g = g + theme_bw() + xlab('') + ylab('')
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

total_peptides = length(unique(norm_info$Peptide))

if (total_peptides > 50){
  g = g + theme(axis.text.y = element_text(size = 5))
}
if (total_peptides > 100){
  g = g + theme(axis.text.y = element_text(size = 3))
}
ggsave(sprintf('%s/%s', output, filename), width = 11, height = 11, useDingbats = FALSE)
