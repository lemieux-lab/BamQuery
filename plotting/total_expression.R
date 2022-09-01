args = commandArgs(trailingOnly=TRUE)

norm_results = args[1]
output = args[2]
name = args[3]
label = args[4]
thout = as.numeric(args[5])
###############################################################################

library(ggplot2)
library(data.table)

tpm = fread(norm_results)

proS = split(tpm, list(tpm$Peptide, tpm$Peptide_Type), drop = TRUE)


proSplus = Reduce(function(...) rbind(...), proS)

proSplus <- proSplus[order(Sample),]

if (grepl('Reads_count', label)) {
  g = ggplot(proSplus, aes(x = Sample,  y = Peptide, fill = value, color = as.factor(value > 0)))
  label_2 = sprintf('%s%s', label, ' > 0' ) 
}else if (grepl('TPM', label)) {
  g = ggplot(proSplus, aes(x = Sample,  y = Peptide, fill = value, color = as.factor(value > 0)))
  label_2 = sprintf('%s%s', label, ' > 0' ) 
}else{
  g = ggplot(proSplus, aes(x = Sample,  y = Peptide, fill = value, color = as.factor(value > log10(thout + 1))))
  label_2 = sprintf('%s%s', label, ' > Log10(8.55)' ) 
  #label = sprintf('%s%s%s', 'mean > log10(', th_out, '+ 1)' ) 
}

g = g + labs(col = label_2) 
g = g + geom_tile(width = 0.75, height = 0.75, size = 0.3)
g = g + scale_color_manual(values=c("FALSE"="grey", "TRUE"="black"))
g = g + guides(color = guide_legend(override.aes = list(fill = "white"))) 

g = g + scale_color_manual(values=c("FALSE"="grey", "TRUE"="black"))
g = g + guides(color = guide_legend(override.aes = list(fill = "white"))) 


if (grepl('TPM', label)) {
  g = g + scale_fill_gradient(low = 'snow1', high = 'lightcoral', na.value = 'white')
  }else{
    g = g + scale_fill_gradient(low = 'snow1', high = 'steelblue', na.value = 'white')
  }

g = g + facet_grid(Peptide_Type ~ Tissue, scales = 'free', space = 'free')
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
#file.remove(norm_results)

#https://r-charts.com/correlation/heat-map-ggplot2/
