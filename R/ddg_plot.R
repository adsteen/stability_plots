# Make stability comparison plots

ddg_plot <- function(df, 
                     with.lines=TRUE,
                     title=NULL, 
                     colors=c("#517C96", "#8D2048"), 
                     ylab=expression(paste("free energy of folding, kJ ", mol^{-1})), 
                     legend.pos="none", 
                     theme=theme_classic(),
                     text.size=NULL) {
  p <- ggplot(df, aes(x=depth, y=total_energy)) +
    geom_boxplot() + 
    geom_point(aes(colour=DDG)) 
  
  if(with.lines) {
    p <- p + geom_point(aes(colour=DDG)) +
      geom_line(aes(x=depth, y=total_energy, colour = DDG, group = paired), alpha = 0.5) 
  }
  
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if(!is.null(colors) & with.lines) {
    p <- p + scale_color_manual(name=NULL, values=colors)
  }
  if(!is.null(ylab)) {
    p <- p + scale_y_continuous(name=ylab)
  }
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  p <- p + theme
  
  if(!is.null(legend)) {
    p <- p + theme(legend.position=legend.pos)
  }
  
  if(!is.null(text.size)) {
    p <- p + theme(text = element_text(size=text.size))
  }
  p
}