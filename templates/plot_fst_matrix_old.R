# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2018/03/04

library(ggplot2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

data <- read.table("AIBST.weighted.fst.estimates.txt")
names(data) <- c("POP1", "POP2", "FST")
upper_tri <- get_upper_tri(data)
# Create a ggheatmap

# Create a ggheatmap
ggheatmap <- ggplot(data, aes(POP2, POP1, fill = FST))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red",
    midpoint = 0, space = "Lab",
    name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
    coord_fixed()
#
# ggheatmap <- ggplot(data, aes(POP2, POP1, fill = FST))+
#     geom_tile(color = "white")+
#     scale_fill_gradient2(low = "red", high = "blue", mid = "white",
#     midpoint = 0, limit = c(-1,1), space = "Lab",
#     name="FST") +
#     theme_minimal()+ # minimal theme
#     theme(axis.text.x = element_text(angle = 45, vjust = 1,
#     size = 12, hjust = 1))+
#     coord_fixed()
# Print the heatmap
# print(ggheatmap)
# ggplot(data = data, aes(x=POP1, y=POP2, fill=FST)) + geom_tile()

ggsave("AIBST.weighted.fst.estimates.png", height=10, width=15.58, units='in', dpi=120)