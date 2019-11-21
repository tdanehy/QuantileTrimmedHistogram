## Width calculation function
calcWidth = function(query) { 
  width(query)
}

## plotting function
plotQthist = function(df, EndBarColor = "red", MiddleBarColor = "black" ) {
  n = length(df) # finding number of observations
  if (n > 1000) {n = 1000}
  if (n < 100) {n = 100}
  q = (11 - round(n/100)) # finding quantiles that should be used for end boxes
  b = ((20-(2*q))/q) # finding the number of bins based on the quantiles
  quant = unname(quantile(df, probs = c((q/100), (1-(q/100)))))
  seq_10 = seq(quant[1], quant[2], length = b)
  div = c(-Inf, round(seq_10), Inf)
  
  colors_vect = c( EndBarColor , rep(MiddleBarColor, (length(div)-3)), EndBarColor) # creates a vector for the colors
  
  df = cutDists(df, divisions= div)
  if ("name" %in% names(df)){
    # It has multiple regions
    g = ggplot(df, aes(x=cuts, y=Freq, fill=name)) + 
      facet_grid(. ~name)
  } else {
    g = ggplot(df, aes(x=cuts, y=Freq))
  }
  
  g = g +
    geom_bar(stat="identity", fill = colors_vect) + 
    theme_classic() + 
    theme(aspect.ratio=1) + 
    theme_blank_facet_label() + 
    xlab("Widths of sequences") +
    ylab("Frequency") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
    theme(plot.title = element_text(hjust = 0.5)) + # Center title
    ggtitle("Quantile Trimmed Histogram") +
    theme(legend.position="bottom") +
    geom_text(aes(label= paste(q,"%", sep='')), data=df[c(1,length(df$Freq)),], vjust=-1)
  return(g)
}


plotQthist(y2)

######################################## helper functions from GenomicDistributions


############ LABEL CUTS
### this function doesn't call any of the other internal functions
labelCuts = function(breakPoints, digits=1, collapse="-", infBins=FALSE) {
  labels = 
    apply(round(cbind( breakPoints[-length(breakPoints)],	
                       breakPoints[-1]),digits), 1, paste0, collapse=collapse) 
  
  if (infBins) {
    labels[1] = paste0("<", breakPoints[2])
    labels[length(labels)] = paste0(">", breakPoints[length(breakPoints)-1])
  }
  return(labels)
}

############## CUT DISTANCES
### This function calls label cuts
cutDists = function(dists, divisions = c(-Inf, -1e6, -1e4, -1000, -100, 0, 100, 1000, 10000, 1e6, Inf)) {
  if (is.list(dists)) {
    x = lapply(dists, cutDists)
    
    # To accommodate multiple lists, we'll need to introduce a new 'name'
    # column to distinguish them.
    nameList = names(dists)
    if(is.null(nameList)) {
      nameList = 1:length(query) # Fallback to sequential numbers
    }
    
    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, sapply(x, nrow))
    
    return(xb)
  }
  
  labels = labelCuts(divisions, collapse=" to ", infBins=TRUE)
  cuts = cut(dists, divisions, labels)
  df = as.data.frame(table(cuts))
  return(df)
}

######### CHANGING THEME
theme_blank_facet_label = function() {
  return(theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
  )
}




