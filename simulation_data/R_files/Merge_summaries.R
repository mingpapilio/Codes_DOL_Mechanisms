library(data.table)
## This function replces double slashes with single slash
file.path2 = function(..., fsep = .Platform$file.sep){
  gsub("//", "/", file.path(..., fsep = fsep))
}

## Please set path before running the rest of the codes
out.file<-""

dirs <- dir(path, pattern="n=[[:alnum:]]+$")
fileName <- file.path(path, paste0(dirs), paste0("summary.txt"))
fileName <- file.path2(fileName)

for (i in seq(along=fileName)){
  file <- fread(fileName[i])
  if(i==1) out.file<- file
  else out.file <- rbind(out.file, file)
}

fwrite(out.file, file = "level_summary.csv",sep=",")
