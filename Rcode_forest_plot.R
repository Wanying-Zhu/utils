# Use this in command line as:
# Rscript Rcode_forest_plot.R mean lower upper label_text fig_title fig_name

library(forestplot)

args = commandArgs(trailingOnly=TRUE)

mean = args[1]
lower = args[2]
upper = args[3]
label_text = args[4]
fig_title = args[5]
fig_name = args[6]

mean = as.numeric(strsplit(mean, ',')[[1]])     # Convert a string to a vector with numeric values
lower = as.numeric(strsplit(lower, ',')[[1]])
upper = as.numeric(strsplit(upper, ',')[[1]])
label_text = strsplit(label_text, ',')[[1]]

number_of_single_analysis = length(mean)-1

jpeg(file=fig_name, width = 800, height = 800, res=150)
par(mai=rep(1,4))

forestplot(labeltext = label_text,
           mean = mean,
           lower = lower,
           upper = upper,
           is.summary=c(rep(FALSE,number_of_single_analysis),TRUE), # Set summary line, so label is be bold and marker is be diamond
           xlog=FALSE, # Set x-axis tick marks to follow a logarithmic scale
           title=fig_title,
           xlab='\nLog odds ratio',
           # txt_gp=fpTxtGp(ticks=gpar(cex=1), xlab=gpar(cex=1), label=gpar(cex=1), title=gpar(cex=2)),
	   shapes_gp=fpShapesGp(zero=(gpar(lty = "longdash"))))

dev.off()
