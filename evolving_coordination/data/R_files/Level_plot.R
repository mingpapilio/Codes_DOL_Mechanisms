### Levelplot
library(lattice)
raw<-read.csv("level_summary.csv",sep=",")
y_grid<- seq(2, 40, by= 2)
x_grid<- seq(0.5, 1.0, length.out= 20)
grid<- expand.grid(X= x_grid, Y= y_grid)
rgb.palette <- colorRampPalette(c("white", "black"), space = "Lab")
rgb.palette_rev <- colorRampPalette(c("black", "white"), space= "Lab")
### Plot
dev.new()
levelplot(raw$avg_s~raw$essentiality*raw$group_size, grid, col.regions=rgb.palette(11), at= seq(0,1,by=0.1), 
          xlab="Trait essentiality, e", ylab="Group size, n")
dev.new()
levelplot(raw$avg_p~raw$essentiality*raw$group_size, grid, col.regions=rgb.palette(11), at= seq(0,1,by=0.1), 
          xlab="Trait essentiality, e", ylab="Group size, n")
dev.new()
levelplot(raw$avg_cord~raw$essentiality*raw$group_size, grid, col.regions=rgb.palette(11), at= seq(0,1,by=0.1), 
          xlab="Trait essentiality, e", ylab="Group size, n")
dev.new()
levelplot(raw$avg_porp~raw$essentiality*raw$group_size, grid, col.regions=rgb.palette(11), at= seq(0,1,by=0.1), 
          xlab="Trait essentiality, e", ylab="Group size, n")
### Proportion data
raw<-read.csv("level_summary.csv",sep=",")
sum(raw$avg_s> 0.1)/ 400
sum(raw$avg_s> 0.9)/sum(raw$avg_s> 0.1)
sum(raw$avg_s> 0.8)/sum(raw$avg_s> 0.1)
