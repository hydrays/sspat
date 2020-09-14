require(lattice)
##output_file_prefix <- "testold"
output_file_prefix <- "./"

cat("\n ***********run********* ", rep, "\n")
FolderName <- paste("case_", output_file_prefix, "_", as.character(rep), "/out/", sep='')
runcmd <- paste("cd", FolderName, "; screen -d -m Rscript ../../out/plot_debug.r; cd ..")
system(runcmd)
