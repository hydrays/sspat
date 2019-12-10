require(XML)
Cfile <- file("commands", "w")

cmd1 <- paste("cp control.txt control.old")
system(cmd1)

##output_file_prefix <- format(Sys.time(), "%Y%m%d%H%M")
output_file_prefix <- "basecomp1"
output_path = "out"
rep = 0

Lbox = 256
tend = 2000.0
dt = 0.01	
iseed = 1012
tpinc = 5.0

alpha = 0.005
beta = 0.01

moblty = 5.0
k0 = 2.0
phi0 = 0.5

diff = 1.0
lambda = 2.0
gamma = 0.1

R1 = 2
nc = 18

output_dir = "out"
cell_position_file = "cellPos.txt"

##for ( stretching_force in c(seq(15, 15, length.out=10), seq(10, 10, length.out=10) ))
##for ( iseed in seq(1010, 1010+12) )
##for ( nc in c(3, 3, 4, 4, 4, 5, 5, 5, 6, 6) )
##for ( alpha in c(0.003, 0.004, 0.005) )
##for ( moblty in seq(1, 10, by=2) )
##for ( lambda in seq(0, 4, by=0.5) )
for ( gamma in seq(0.1, 5, by=0.5) )
##for ( k0 in seq(6) )
##for ( diff in c(0.1, 0.5, 1, 2, 5, 10) )
{
    ##stretching_force = 4
    iseed = iseed + rep
    rep = rep + 1
    cat("\n ***********run********* ", rep, "\n")
    FolderName <- paste("case_", output_file_prefix, "_", as.character(rep), sep='')
    dir.create(FolderName)

    ## The following part is to produce "control.txt"
    zz <- file("control.txt", "w")
    cat("! Control parameters\n", file=zz, append=T)
    cat("&xdata\n", file=zz, append=T)
    cat(paste("Lbox = ", Lbox, ",\n", sep=''), file=zz, append=T)
    cat(paste("tend = ", tend, ",\n", sep=''), file=zz, append=T)
    cat(paste("dt = ", dt, ",\n", sep=''), file=zz, append=T)
    cat(paste("alpha = ", alpha, ",\n", sep=''), file=zz, append=T)
    cat(paste("diff = ", diff, ",\n", sep=''), file=zz, append=T)
    cat(paste("iseed = ", iseed, ",\n", sep=''), file=zz, append=T)
    cat(paste("tpinc = ", tpinc, ",\n", sep=''), file=zz, append=T)
    cat(paste("R1 = ", R1, ",\n", sep=''), file=zz, append=T)
    cat(paste("beta = ", beta, ",\n", sep=''), file=zz, append=T)
    cat(paste("moblty = ", moblty, ",\n", sep=''), file=zz, append=T)
    cat(paste("k0 = ", k0, ",\n", sep=''), file=zz, append=T)
    cat(paste("phi0 = ", phi0, ",\n", sep=''), file=zz, append=T)    
    cat(paste("lambda = ", lambda, ",\n", sep=''), file=zz, append=T)
    cat(paste("gamma = ", gamma, ",\n", sep=''), file=zz, append=T)
    cat(paste("nc = ", nc, ",\n", sep=''), file=zz, append=T)            
    cat("/\n", file=zz, append=T)

    ## Copy the seqfile.txt into folder
    cmd2 <- paste("cp control.txt", FolderName)
    cmd3 <- paste("mkdir -p ", FolderName, "/out", sep='')
    cmd4 <- paste("cp run", FolderName)
    system(cmd2)
    system(cmd3)
    system(cmd4)

    ##rtime <- runif(1, 20, 200)
    ##runcmd <- paste("sleep ", rtime, "; cd", FolderName, "; screen PhysModel; cd ..")

    ## screen, on PC
    runcmd <- paste("cd", FolderName, "; screen -d -m run; cd ..")
    ## runcmd <- paste("screen -d -m ", FolderName, "/PhysModel", sep='')

    ## ## bsub, on server
    ## runcmd <- paste("bsub ", FolderName, "/PhysModel -o output_", rep, sep='')	
    ## runcmd <- paste("cd ", FolderName, "; bsub PhysModel -o output_", rep, "; cd ..", sep='')	

    writeLines(runcmd, Cfile)

    close(zz)
}

cmd1 <- paste("mv control.old control.txt")
system(cmd1)

close(Cfile)

