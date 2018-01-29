# Author: Alban Ramette
# Description: Given one fasta sequence, creates multifasta files containing mutant fasta sequences
# allows only point mutations. See EMBOSS::msbar (not R) for more complex mutation schemes
# Explanations about filenames generated:
#   e.g. Fasta1_MR0.01_Nmut100.fasta means
# Fasta1 = original name 
# MR0.01 = mutation rate of 0.01 = 1% = (seqlength of 100 bases x 0.01 = 1 point mutation was introduced)
# Nmut100 = 100 mutated sequences. (the last number in each fasta header indicates the mutant id e.g. …._1, …C_2
# Warning message: 
# "Caution: Cannot create non-redundant mutants for Fasta1.txt and MR=0.01sequence= 622!!!"
# indicates that the algorthim could not find a unique mutant (the same mutants occur at least twice in this set)
# Change Ntries to make more search of unique sequences (may take more time then)



#PARAMETERS
FILES=c("Fasta1.txt")
MutationRate=c(0.01,0.02,0.03,0.05,0.10)
Nmutants=100 # number of mutated fasta 
Ntries=10 # how many times to try searching for a new sequence of a variant?
#-----------
for (FASTAFILE in FILES) {
  stub= gsub(pattern=".txt|.fasta|.fna",replacement = "",x=FASTAFILE,perl=TRUE)
  dir.create(stub)
  D=scan(FASTAFILE,what="character",quiet = TRUE)
  Header=D[1]; Header=gsub(pattern = ">","",Header)
  Seq=D[2] ; rm(D)
  Seqlen <- nchar(Seq)
  mutate=function(S=Seq,Nmut=3){
  Pool=c("A","T","G","C")
  #Nmut= number of mutations in the sequence (abs)
  S1 <- seqinr::getSequence(S)
  POSr <- sample(length(S1),Nmut,replace = FALSE) #define random positions
  Smut <- S1
  
  for(i in 1:Nmut){
    Base <- S1[POSr[i]] #actual base
    #BaseLeft <- setdiff(Pool,Base)
    Smut[POSr[i]] <- sample(setdiff(Pool,Base),1)
  }
  return(paste0(Smut,collapse = ""))
} #
  #--MAIN------------------------- 
  for(MRT in MutationRate){ 
        M <- floor(Seqlen*MRT)#Nber mutations
        Mutants <- vector("list",Nmutants)
        names(Mutants) <- paste0(Header,"_",1:Nmutants)
        for(i in 1:Nmutants){
          C=0; flag=0
          Mtemp <- mutate(S=Seq,Nmut=M)
          while(any(Mutants==Mtemp)==TRUE) { # testing to see whether the same seq was already selected
                  Mtemp <- mutate(S=Seq,Nmut=M); #then mutates again
                  C <- C+1  
                  if(C>Ntries){C=0;  #up to Ntries times
                      flag=1
                      break;
                  }
          }
          Mutants[i]  <- Mtemp  # check if already created
          if(flag==1){flag=0;cat(paste0("Caution: Cannot create non-redundant mutants for ",FASTAFILE," and MR=",MRT,", redundant sequence= ",i,"!!!\n"));}
        }
        # write to file
        FILEOUT=paste0(stub,"/",stub,"_MR",MRT,"_Nmut",Nmutants,".fasta")
        seqinr::write.fasta(
          sequences= as.list(Mutants),
          names= names(Mutants),
          file.out=FILEOUT         
        )
  }
}


  

  
  
  
  
