args = commandArgs(trailingOnly=TRUE)

data = args[1]
output = args[2]

###############################################################################

library(ggplot2)
library(data.table)
library(ggpubr)

my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
"#E16462FF", "#FCA636FF", "#F0F921FF")

data = read.csv(data, header = TRUE, sep = "\t", row.name = 1)
filename = 'positions.pdf'



ggballoonplot(data, fill = "value", show.label = TRUE)+
   scale_fill_gradientn(colors = my_cols)

# g = ggplot(norm_info, aes(x=Positions, y = Peptide, fill=Position))

# g = g + geom_tile(aes(fill = Position))
# #g = g + geom_text(aes(label = round(Position, 1))) 
# g = g + scale_fill_gradient(low = "white", high = "red")

# total_peptides = length(unique(norm_info$Peptide))

# if (total_peptides > 50){
#   g = g + theme(axis.text.y = element_text(size = 5))
# }
# if (total_peptides > 100){
#   g = g + theme(axis.text.y = element_text(size = 3))
# }
ggsave(sprintf('%s/%s', output, filename), width = 11, height = 11, useDingbats = FALSE)
