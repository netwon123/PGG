setwd('G:/博士/泛基因组组/PGG-project/analysis/eval_graph/final/')

list.of.packages <- c("tidyverse", "ggrepel", "svglite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")
require("scales") # For squish


# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name, count
dat <- read.table('results-vg.tsv.gz', header=T)

if (! ("count" %in% names(dat))) {
  # If the count column is not present, add i
  dat$count <- rep(1, nrow(dat))
}



# Determine title
title <- ''

# Determine the order of aligners, based on sorting in a dash-separated tag aware manner
aligner.names <- levels(factor(dat$aligner))
name.lists <- aligner.names %>% (function(name) map(name,  (function(x) as.list(unlist(strsplit(x, "-"))))))
# Transpose name fragments into a list of vectors for each position, with NAs when tag lists end early
max.parts <- max(sapply(name.lists, length))
name.cols <- list()
for (i in 1:max.parts) {
  name.cols[[i]] <- sapply(name.lists, function(x) if (length(x) >= i) { x[[i]] } else { NA })
}
name.order <- do.call(order,name.cols)
aligner.names <- aligner.names[name.order]
dat$aligner <- factor(dat$aligner, levels=aligner.names)
name.lists <- name.lists[name.order]

# Determine colors for aligners
bold.colors <- c("#1f78b4","#e31a1c","#33a02c","#6600cc","#ff8000","#5c415d","#458b74","#698b22","#008b8b")
light.colors <- c("#82B0D2","#FA7F6F","#8ECFC9","#e5ccff","#ffe5cc","#9a7c9b","#76eec6","#b3ee3a","#00eeee")
# We have to go through both lists together when assigning colors, because pe and non-pe versions of a condition need corresponding colors.
cursor <- 1

# This will map from non-pe condition name string to color index.
colors <- c()
for (i in 1:length(name.lists)) {
  # For each name
  name.parts <- unlist(name.lists[[i]])
  if (name.parts[length(name.parts)] == "pe") {
    # Drop the pe tag if present
    name.parts <- name.parts[-c(length(name.parts))]
  }
  if (name.parts[length(name.parts)] == "se") {
    # Drop the se tag if present
    name.parts <- name.parts[-c(length(name.parts))]
  }
  
  # Join up to a string again
  name <- paste(name.parts, collapse='-')
  
  if (! name %in% names(colors)) {
    # No colors assigned for this pair of conditions, so assign them.
    
    if (cursor > length(bold.colors)) {
      write(colors, stderr())
      write(aligner.names, stderr())
      stop('Ran out of colors! Too many conditions!')
    }
    
    # We always assign pe and non-pe colors in lockstep, whichever we see first.
    # We need two entries for -se and no tag which are the same.
    new.colors <- c(bold.colors[cursor], light.colors[cursor], light.colors[cursor])
    names(new.colors) <- c(paste(name, 'pe', sep='-'), paste(name, 'se', sep='-'), name)
    colors <- c(colors, new.colors)
    
    cursor <- cursor + 1
  }
}

# Make colors a vector in the same order as the actually-used aligner names
colors <- colors[aligner.names]

dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))
dat.roc <- dat %>%
  mutate(Positive = (correct == 1) * count, Negative = (correct == 0) * count) %>%
  group_by(aligner, mq) %>%
  summarise(Positive = sum(Positive), Negative = sum(Negative)) %>%
  arrange(-mq) %>%
  mutate(Total=sum(Positive+Negative)) %>%
  mutate(TPR = cumsum(Positive) / Total, FPR = cumsum(Negative) / Total)

# We want smart scales that know how tiny a rate of things we can care about
total.reads <- max(dat.roc$Total)
min.log10 <- floor(log10(1/total.reads))
min.log10=-4
max.log10 <- 0
# Work out a set of bounds to draw the plot on
range.log10 <- min.log10 : max.log10
range.unlogged = 10^range.log10
pdf('roc.pdf',height = 10,width = 9)
ggplot(dat.roc, aes( x= FPR, y = TPR, color = aligner, label=mq)) +
  geom_line(alpha = 0.5,size=1) + geom_text_repel(data = dat.roc, size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5),show.legend = FALSE) +
  geom_point(aes(size=Positive+Negative),alpha = 0.75) +
  scale_color_manual(values=colors, guide=guide_legend(title=NULL)) +
  scale_size_continuous("number", guide=guide_legend(title=NULL),labels = function(x) format(x, scientific = FALSE),range=c(1,10)) +
  scale_x_log10(limits=c(range.unlogged[1],range.unlogged[length(range.unlogged)]), breaks=range.unlogged, oob=squish) +
  theme_bw() + theme(legend.position = c(1, 0),
                     legend.justification = c(1, 0),
                     legend.margin = margin(t = 0, r = 0, b = 10, l = 10),legend.box.background = element_rect(colour="black",size = 0.1),
                     legend.title= element_text( size = 10, color = "black"),
                     legend.background = element_rect(fill = "transparent"))+labs(x = "False-positive rate", y = "True-positive rate")+
  theme(axis.text = element_text( size = 10, color = "black"),axis.title = element_text( size = 12, color = "black"))
dev.off()                                                
A=windowsFont("Arial")





if (title != '') {
  # And a title
  dat.plot + ggtitle(title)
}

