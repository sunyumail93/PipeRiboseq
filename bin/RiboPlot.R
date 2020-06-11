#RiboPlot.R
#This script is used for PipeRiboseq pipeline, and should be kept in $PipelineHomeDir/bin folder
#Y. Sun, 2019-01-04
#Modes: 1) lendis, 2) quadplot, 3) metaAUG, 4) metaSTOP
#Input files:
# 1) lendis mode: InputLendis, NormFactor (to divide), OutputFileName         (Rscript RiboPlot.R lendis sam1.mRNA.mm10.Lendis 26.9496 sam1.mRNA.mm10.Lendis.pdf)
# 2) quadplot mode: Prefix, TranscriptName, AnnotationFile, Gene, NormFactor  (Rscript RiboPlot.R quadplot sam1.mRNA NM_206534 sam1.mm10.mRNAlist.bed12 CHURC1 1.35)
# 3) metaAUG mode: Prefix                                                     (Rscript RiboPlot.R metaAUG sam1.mRNA.AUG)
# 4) metaSTOP mode: Prefix                                                    (Rscript RiboPlot.R metaSTOP sam1.mRNA.STOP)

#This is the default extension length around AUG and STOP. You may change this to a number that can be divided by three. If not divided by three, the frames will be shifted.
DefaultExtension=60

DrawLendis <- function(Input, Norm, Output){
  Data <- read.table(Input, header = F)
  Tablecol <- c("18","19","20",
                "21","22","23","24","25","26","27","28","29","30",
                "31","32","33","34","35","36","37","38","39","40",
                "41","42","43","44","45","46","47","48","49","50")
  opar<- par(lwd = 0.3)
  pdf(Output,width=16,height=9)
  barplot(
    Data$V2/as.numeric(Norm),names.arg=Tablecol,
    xlab="RPF length (nt)",ylab='Normalized counts',
    main=Input,
    col="black",
    beside=TRUE
  )
  mtext(paste0("Normalization factor: ",Norm,sep=""))
  temp <- dev.off()
  print("lendis figure done.")
}

DrawBED12 <- function(BED12Annotation,ypos,barheight){
  Genelength <- BED12Annotation$V3 - BED12Annotation$V2
  VectorPlus <- rep(0,Genelength)
  VectorMinus <- rep(0,Genelength)
  if (BED12Annotation$V10 == 1){
    #This gene contains only one exon
    VectorPlus <- rep(barheight,Genelength)
    VectorMinus <- rep(-barheight,Genelength)
    rect(xleft=1,ybottom=ypos-barheight,xright=Genelength,ytop=ypos+barheight,col = "chocolate3")
  }else{
    #More than 1 exons
    strand <- as.character(Annotation_Curr$V6)
    StartBlocks <- as.character(unlist(BED12Annotation$V12))
    LenBlocks <- as.character(unlist(BED12Annotation$V11))
    StartBlocks_sep <- as.numeric(unlist(strsplit(StartBlocks,",")))
    LenBlocks_sep <- as.numeric(unlist(strsplit(LenBlocks,",")))
    if (strand == "+"){
      for (i in 1:BED12Annotation$V10){
        for (j in (StartBlocks_sep[i]+1):(StartBlocks_sep[i]+LenBlocks_sep[i])){
          VectorPlus[j] = barheight
          VectorMinus[j] = -barheight
        }
        rect(xleft=StartBlocks_sep[i]+1,ybottom=ypos-barheight,xright=StartBlocks_sep[i]+LenBlocks_sep[i],ytop=ypos+barheight,col = "chocolate3")
      }
    }
    else{
      for (i in BED12Annotation$V10:1){
        for (j in (Genelength-StartBlocks_sep[i]-LenBlocks_sep[i]+1):(Genelength-StartBlocks_sep[i])){
          VectorPlus[j] = barheight
          VectorMinus[j] = -barheight
        }
        rect(xleft=Genelength-StartBlocks_sep[i]-LenBlocks_sep[i]+1,ybottom=ypos-barheight,xright=Genelength-StartBlocks_sep[i],ytop=ypos+barheight,col = "chocolate3")
      }
    }
  }
}

ExonLength <- function(BED12Annotation){
  LenBlocks <- as.character(unlist(BED12Annotation$V11))
  LenBlocks_sep <- as.numeric(unlist(strsplit(LenBlocks,",")))
  return(sum(LenBlocks_sep))
}

args <- commandArgs(TRUE)
Mode <- args[1]
if (Mode == "lendis"){
  print("Plotting figures with lendis mode")
  Input=args[2]
  Norm=args[3]
  Output=args[4]
  DrawLendis(Input, Norm, Output)
}else if (Mode == "quadplot"){
  print("Plotting figures with quadplot mode")
  #This mode draws four figures in the same plot!
  Prefix=args[2]
  TranscriptName=args[3]
  AnnotationFile=args[4]
  GeneName=args[5]
  Norm=args[6]
  
  pdf(paste0(Prefix, ".", TranscriptName,".pdf"),width = 16, height = 9)
  par(mfrow=c(2,2))
  
  PlusReadsFile <- paste(Prefix,TranscriptName,"plus.bedGraph.bb",sep = ".")
  MinusReadsFile <- paste(Prefix,TranscriptName,"minus.bedGraph.bb",sep = ".")
  PlusReads <- read.table(PlusReadsFile,col.names = c("chr","start","end","value"))
  MinusReads <- read.table(MinusReadsFile,col.names = c("chr","start","end","value"))
  Annotation <- read.table(AnnotationFile,header = F)
  #Process
  Length <- dim(PlusReads)[1]
  Annotation_Curr <- subset(Annotation,Annotation$V4 == TranscriptName)
  strand <- as.character(Annotation_Curr$V6)
  #Range
  ReadsMaxS <- max(PlusReads$value)
  ReadsMaxA <- min(MinusReads$value)
  Real_ReadsMaxS=ReadsMaxS
  Real_ReadsMaxA=ReadsMaxA
  if (ReadsMaxS == 0){ReadsMaxS=0.0010}
  if (ReadsMaxA == 0){ReadsMaxA=-0.0010}
  
  #Plot 1
  if (strand == "+"){
    Yrange = ReadsMaxS-ReadsMaxA
    Outspace <- Yrange*0.2
    plot(1:Length, PlusReads$value, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="blue", xlab="", ylab="", ann=F, type="h", lwd=2,las=1)
    par(new=T)
    plot(1:Length, MinusReads$value, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="red", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
    par(new=T)
    plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
    title(main=paste0(TranscriptName,", ",GeneName))
    mtext(paste(Prefix,"; ","Max: ",round(Real_ReadsMaxS,2),"; Min: ",round(Real_ReadsMaxA,2),"; Strand: ",strand,"; Exon: ",Annotation_Curr$V10,"; CDS labelled; Norm Factor: ",Norm,sep=""), cex=0.8)
    title(xlab = paste("Gene (",Length," nt)",sep = ""),ylab = "Normalized counts")
    Ts_Start=Annotation_Curr$V7-Annotation_Curr$V2+1
    Ts_End=Annotation_Curr$V8-Annotation_Curr$V2

  }else{
    Yrange = ReadsMaxS-ReadsMaxA
    Outspace <- Yrange*0.2
    plot(1:Length, -rev(MinusReads$value), xlim=c(1,Length), ylim=c(-ReadsMaxS*1.1-Outspace, -ReadsMaxA*1.1), col="blue", xlab="", ylab="", ann=F, type="h", lwd=2,las=1)
    par(new=T)
    plot(1:Length, -rev(PlusReads$value), xlim=c(1,Length), ylim=c(-ReadsMaxS*1.1-Outspace, -ReadsMaxA*1.1), col="red", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
    par(new=T)
    plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(-ReadsMaxS*1.1-Outspace, -ReadsMaxA*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
    title(main=paste0(TranscriptName,", ",GeneName))
    mtext(paste(Prefix,"; ","Max: ",round(-Real_ReadsMaxA,2),"; Min: ",round(-Real_ReadsMaxS,2),"; Strand: ",strand,"; Exon: ",Annotation_Curr$V10,"; CDS labelled; Norm Factor: ",Norm,sep=""), cex=0.8)
    title(xlab = paste("Gene (",Length," nt)",sep = ""),ylab = "Normalized counts")
    Ts_Start=Annotation_Curr$V3-Annotation_Curr$V8+1
    Ts_End=Annotation_Curr$V3-Annotation_Curr$V7
  }
  #Decide bed12 drawing parameters
  yanno <- ReadsMaxA*1.1-Outspace/2
  yheight <- yanno - Outspace/18
  par(new=T)
  plot(1:Length,rep(yanno,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="l", lwd=2,las=1)
  DrawBED12(Annotation_Curr,yanno,Outspace/18)
  rect(xleft = Ts_Start, ybottom = yanno-Outspace/4, xright = Ts_End, ytop = yheight-Outspace/4 - Outspace/12,col = "orchid")
  
  #Plot 2
  DataFile2 <- paste0(Prefix,".Both.",TranscriptName,".sense.bedGraph")
  Data <- read.table(DataFile2,header = F)
  Length <- dim(Data)[1]
  ReadsMaxS <- max(Data$V4)
  ReadsMaxA <- 0
  Real_ReadsMaxS=ReadsMaxS
  Real_ReadsMaxA=ReadsMaxA
  if (ReadsMaxS == 0){ReadsMaxS=0.0010}
  if (ReadsMaxA == 0 & ReadsMaxS>0){ReadsMaxA=-ReadsMaxS*0.1}
  if (ReadsMaxA == 0 & ReadsMaxS==0){ReadsMaxA=-0.001}
  #Decide the plotting color order:
  ColorSeq=c("red","darkgreen","dodgerblue2")
  Yrange = ReadsMaxS-ReadsMaxA
  Outspace <- Yrange*0.2
  plot(1:Length, Data$V4, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col=ColorSeq, xlab="", ylab="", ann=F, type="h", lwd=1,las=1)
  par(new=T)
  plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
  abline(h = 0,col="lightblue")
  title(main=paste0(Prefix," CDS signals"))
  mtext(paste("CDS region, Frame 0: red; Frame 1: green; Frame 2: blue; Frame percentage on CDS",sep=""),cex=0.8)
  title(xlab = paste("Transcript, CDS with 60 nt extension (",Length," nt)",sep = ""),ylab = "Normalized counts")
  yanno <- ReadsMaxA*1.1-Outspace/2
  yheight <- yanno - Outspace/18
  par(new=T)
  rect(xleft = DefaultExtension, ybottom = yheight, xright = Length-DefaultExtension, ytop = yanno, col="grey")
  CDSValues=Data$V4[(DefaultExtension+1):(Length-DefaultExtension)]
  Frame0=0
  Frame1=0
  Frame2=0
  for (i in seq(1,length(CDSValues),3)){Frame0=Frame0+CDSValues[i]}
  for (i in seq(2,length(CDSValues),3)){Frame1=Frame1+CDSValues[i]}
  for (i in seq(3,length(CDSValues),3)){Frame2=Frame2+CDSValues[i]}
  Sum=Frame0+Frame1+Frame2
  if (Sum > 0){
    Frame0p=format(round(Frame0/Sum*100, 2), nsmall = 2)
    Frame1p=format(round(Frame1/Sum*100, 2), nsmall = 2)
    Frame2p=format(round(Frame2/Sum*100, 2), nsmall = 2)
  }else{
    Frame0p=0
    Frame1p=0
    Frame2p=0
  }
  print(paste0("Frame 0 counts: ",Frame0))
  print(paste0("Frame 1 counts: ",Frame1))
  print(paste0("Frame 2 counts: ",Frame2))
  print(paste0("Frame 0 percentage: ",Frame0p," %"))
  print(paste0("Frame 1 percentage: ",Frame1p," %"))
  print(paste0("Frame 2 percentage: ",Frame2p," %"))
  legend(x="topright", fill = c("red","darkgreen","dodgerblue2"),legend = c(paste0(Frame0p,"%"), paste0(Frame1p,"%"), paste0(Frame2p,"%")),bty = "n")
  
  #Plot 3
  DataFile3 <- paste0(Prefix,".AUG.",TranscriptName,".sense.bedGraph")
  GeneFASTAFile <- paste0(Prefix,".AUG.",TranscriptName,".bed12.fa")
  Data <- read.table(DataFile3,header = F,col.names = c("chr","start","end","value"))
  GeneFASTA <- read.table(GeneFASTAFile,header = T,check.names = F)
  Sequence <- as.character(GeneFASTA[1,])
  Middle <- (dim(Data)[1]-1)/2+1
  Length <- dim(Data)[1]
  ReadsMaxS <- max(Data$value)
  ReadsMaxA <- 0
  Real_ReadsMaxS=ReadsMaxS
  Real_ReadsMaxA=ReadsMaxA
  if (ReadsMaxS == 0){ReadsMaxS=0.0010}
  if (ReadsMaxA == 0 & ReadsMaxS>0){ReadsMaxA=-ReadsMaxS*0.1}
  if (ReadsMaxA == 0 & ReadsMaxS==0){ReadsMaxA=-0.001}
  #Decide the plotting color order:
  if (Middle %% 3 ==1){
    ColorSeq=c("red","darkgreen","dodgerblue2")
  }else if (Middle %% 3 ==2){
    ColorSeq=c("dodgerblue2","red","darkgreen")
  }else{
    ColorSeq=c("darkgreen","dodgerblue2","red")
  }
  Yrange = ReadsMaxS-ReadsMaxA
  Outspace <- Yrange*0.2
  plot(1:Length, Data$value, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col=ColorSeq, xlab="", ylab="", ann=F, type="h", lwd=1.5,las=1)
  par(new=T)
  plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
  abline(h = 0,col="lightblue")
  title(main="Signals around AUG")
  mtext(paste("AUG around, Frame 0: red; Frame 1: green; Frame 2: blue; Frame percentage after AUG",sep=""),cex=0.8)
  title(xlab = paste("Transcript, 60 nt extension around AUG (",Length," nt)",sep = ""),ylab = "Normalized counts")
  Letter=1:Length
  yanno <- ReadsMaxA*1.1-Outspace/2
  for (i in 1:Length){
    Letter[i] <- substring(Sequence,i,i)
    if (i == Middle){
      text(x = i,y=yanno,Letter[i],col = "red",cex=0.5)
    }else if(i == Middle+1){
      text(x = i,y=yanno,Letter[i],col = "red",cex=0.5)
    }else if(i == Middle+2){
      text(x = i,y=yanno,Letter[i],col = "red",cex=0.5)
    }else{
      text(x = i,y=yanno,Letter[i],col = "black",cex=0.5)
    }
  }
  CDSValues=Data$value[Middle:Length]
  Frame0=0
  Frame1=0
  Frame2=0
  for (i in seq(1,length(CDSValues),3)){Frame0=Frame0+CDSValues[i]}
  for (i in seq(2,length(CDSValues),3)){Frame1=Frame1+CDSValues[i]}
  for (i in seq(3,length(CDSValues),3)){Frame2=Frame2+CDSValues[i]}
  Sum=Frame0+Frame1+Frame2
  if (Sum > 0){
    Frame0p=format(round(Frame0/Sum*100, 2), nsmall = 2)
    Frame1p=format(round(Frame1/Sum*100, 2), nsmall = 2)
    Frame2p=format(round(Frame2/Sum*100, 2), nsmall = 2)
  }else{
    Frame0p=0
    Frame1p=0
    Frame2p=0
  }
  legend(x="topright", fill = c("red","darkgreen","dodgerblue2"),legend = c(paste0(Frame0p,"%"), paste0(Frame1p,"%"), paste0(Frame2p,"%")),bty = "n")
  
  #Plot 4
  DataFile3 <- paste0(Prefix,".STOP.",TranscriptName,".sense.bedGraph")
  GeneFASTAFile <- paste0(Prefix,".STOP.",TranscriptName,".bed12.fa")
  Data <- read.table(DataFile3,header = F,col.names = c("chr","start","end","value"))
  GeneFASTA <- read.table(GeneFASTAFile,header = T,check.names = F)
  Sequence <- as.character(GeneFASTA[1,])
  Middle <- (dim(Data)[1]-1)/2+1
  Length <- dim(Data)[1]
  ReadsMaxS <- max(Data$value)
  ReadsMaxA <- 0
  Real_ReadsMaxS=ReadsMaxS
  Real_ReadsMaxA=ReadsMaxA
  if (ReadsMaxS == 0){ReadsMaxS=0.0010}
  if (ReadsMaxA == 0 & ReadsMaxS>0){ReadsMaxA=-ReadsMaxS*0.1}
  if (ReadsMaxA == 0 & ReadsMaxS==0){ReadsMaxA=-0.001}
  #Decide the plotting color order:
  if (Middle %% 3 ==1){
    ColorSeq=c("red","darkgreen","dodgerblue2")
  }else if (Middle %% 3 ==2){
    ColorSeq=c("dodgerblue2","red","darkgreen")
  }else{
    ColorSeq=c("darkgreen","dodgerblue2","red")
  }
  Yrange = ReadsMaxS-ReadsMaxA
  Outspace <- Yrange*0.2
  plot(1:Length, Data$value, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col=ColorSeq, xlab="", ylab="", ann=F, type="h", lwd=1.5,las=1)
  par(new=T)
  plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
  abline(h = 0,col="lightblue")
  title(main="Signals around STOP")
  mtext(paste("STOP around, Frame 0: red; Frame 1: green; Frame 2: blue; Frame percentage before STOP",sep=""),cex=0.8)
  title(xlab = paste("Transcript, 60 nt extension around STOP (",Length," nt)",sep = ""),ylab = "Normalized counts")
  Letter=1:Length
  yanno <- ReadsMaxA*1.1-Outspace/2
  for (i in 1:Length){
    Letter[i] <- substring(Sequence,i,i)
    if (i == Middle){
      text(x = i,y=yanno,Letter[i],col = "red",cex=0.5)
    }else if(i == Middle+1){
      text(x = i,y=yanno,Letter[i],col = "red",cex=0.5)
    }else if(i == Middle+2){
      text(x = i,y=yanno,Letter[i],col = "red",cex=0.5)
    }else{
      text(x = i,y=yanno,Letter[i],col = "black",cex=0.5)
    }
  }
  CDSValues=Data$value[1:Middle]
  Frame0=0
  Frame1=0
  Frame2=0
  for (i in seq(1,length(CDSValues),3)){Frame0=Frame0+CDSValues[i]}
  for (i in seq(2,length(CDSValues),3)){Frame1=Frame1+CDSValues[i]}
  for (i in seq(3,length(CDSValues),3)){Frame2=Frame2+CDSValues[i]}
  Sum=Frame0+Frame1+Frame2
  if (Sum > 0){
    Frame0p=format(round(Frame0/Sum*100, 2), nsmall = 2)
    Frame1p=format(round(Frame1/Sum*100, 2), nsmall = 2)
    Frame2p=format(round(Frame2/Sum*100, 2), nsmall = 2)
  }else{
    Frame0p=0
    Frame1p=0
    Frame2p=0
  }
  legend(x="topright", fill = c("red","darkgreen","dodgerblue2"),legend = c(paste0(Frame0p,"%"), paste0(Frame1p,"%"), paste0(Frame2p,"%")),bty = "n")
  #Finishing
  temp <- dev.off()
  print("quadplot figures done.")
  
}else if (Mode == "metaAUG"){
  print("Plotting figures with meta mode, around AUG")
  Prefix <- args[2]
  DataFile <- paste0(Prefix,".bedGraph.summary.txt")
  OutputFile <- paste0(Prefix,".Meta.pdf")
  
  Data <- read.table(DataFile,header = T)
  Data$Mean <- apply(Data[2:ncol(Data)],1,FUN = mean,trim=0.1)
  
  Length <- dim(Data)[1]
  Middle <- (dim(Data)[1]-1)/2+1
  ReadsMaxS <- max(Data$Mean)
  ReadsMaxA <- 0
  Real_ReadsMaxS=ReadsMaxS
  Real_ReadsMaxA=ReadsMaxA
  if (ReadsMaxS == 0){ReadsMaxS=0.0010}
  if (ReadsMaxA == 0 & ReadsMaxS>0){ReadsMaxA=-ReadsMaxS*0.1}
  if (ReadsMaxA == 0 & ReadsMaxS==0){ReadsMaxA=-0.001}
  
  #Decide the plotting color order:
  if (Middle %% 3 ==1){
    ColorSeq=c("red","darkgreen","dodgerblue2")
  }else if (Middle %% 3 ==2){
    ColorSeq=c("dodgerblue2","red","darkgreen")
  }else{
    ColorSeq=c("darkgreen","dodgerblue2","red")
  }
  
  pdf(OutputFile,width = 12,height = 8)
  Yrange = ReadsMaxS-ReadsMaxA
  Outspace <- Yrange*0.2
  plot(1:Length, Data$Mean, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col=ColorSeq, xlab="", ylab="", ann=F, type="h", lwd=2,las=1)
  par(new=T)
  plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
  abline(h = 0,col="lightblue")
  title(main=paste0(Prefix,".subcodon.summary"))
  mtext(paste(Prefix,"; total ",Length," nt; around AUG",sep=""),cex=0.8)
  title(xlab = paste("Transcript (",Length," nt)",sep = ""),ylab = "Normalized counts")
  
  Letter=1:Length
  yanno <- ReadsMaxA*1.1-Outspace/2
  if (Length <= 180){
    CharNum=0.6
  }else if (Length > 180){
    CharNum=0.3
  }
  for (i in 1:Length){
    Letter[i] <- as.character(Data$Nucleotide)[i]
    if (i == Middle){
      text(x = i,y=yanno,Letter[i],col = "red",cex=CharNum)
    }else if(i == Middle+1){
      text(x = i,y=yanno,Letter[i],col = "red",cex=CharNum)
    }else if(i == Middle+2){
      text(x = i,y=yanno,Letter[i],col = "red",cex=CharNum)
    }else{
      text(x = i,y=yanno,Letter[i],col = "black",cex=CharNum)
    }
  }
  CDSValues=Data$Mean[Middle:Length]
  Frame0=0
  Frame1=0
  Frame2=0
  for (i in seq(1,length(CDSValues),3)){Frame0=Frame0+CDSValues[i]}
  for (i in seq(2,length(CDSValues),3)){Frame1=Frame1+CDSValues[i]}
  for (i in seq(3,length(CDSValues),3)){Frame2=Frame2+CDSValues[i]}
  Sum=Frame0+Frame1+Frame2
  if (Sum > 0){
    Frame0p=format(round(Frame0/Sum*100, 2), nsmall = 2)
    Frame1p=format(round(Frame1/Sum*100, 2), nsmall = 2)
    Frame2p=format(round(Frame2/Sum*100, 2), nsmall = 2)
  }else{
    Frame0p=0
    Frame1p=0
    Frame2p=0
  }
  print(paste0("Frame 0 counts: ",Frame0))
  print(paste0("Frame 1 counts: ",Frame1))
  print(paste0("Frame 2 counts: ",Frame2))
  print(paste0("Frame 0 percentage: ",Frame0p," %"))
  print(paste0("Frame 1 percentage: ",Frame1p," %"))
  print(paste0("Frame 2 percentage: ",Frame2p," %"))
  legend(x="topright", fill = c("red","darkgreen","dodgerblue2"),legend = c(paste0(Frame0p,"%"), paste0(Frame1p,"%"), paste0(Frame2p,"%")),bty = "n")
  #Finishing
  temp <- dev.off()
  print("metaAUG figure done.")
  
}else if (Mode == "metaSTOP"){
  print("Plotting figures with meta mode, around STOP")
  Prefix <- args[2]
  DataFile <- paste0(Prefix,".bedGraph.summary.txt")
  OutputFile <- paste0(Prefix,".Meta.pdf")
  
  Data <- read.table(DataFile,header = T)
  Data$Mean <- apply(Data[2:ncol(Data)],1,FUN = mean,trim=0.1)
  
  Length <- dim(Data)[1]
  Middle <- (dim(Data)[1]-1)/2+1
  ReadsMaxS <- max(Data$Mean)
  ReadsMaxA <- 0
  Real_ReadsMaxS=ReadsMaxS
  Real_ReadsMaxA=ReadsMaxA
  if (ReadsMaxS == 0){ReadsMaxS=0.0010}
  if (ReadsMaxA == 0 & ReadsMaxS>0){ReadsMaxA=-ReadsMaxS*0.1}
  if (ReadsMaxA == 0 & ReadsMaxS==0){ReadsMaxA=-0.001}
  
  #Decide the plotting color order:
  if (Middle %% 3 ==1){
    ColorSeq=c("red","darkgreen","dodgerblue2")
  }else if (Middle %% 3 ==2){
    ColorSeq=c("dodgerblue2","red","darkgreen")
  }else{
    ColorSeq=c("darkgreen","dodgerblue2","red")
  }
  
  pdf(OutputFile,width = 12,height = 8)
  Yrange = ReadsMaxS-ReadsMaxA
  Outspace <- Yrange*0.2
  plot(1:Length, Data$Mean, xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col=ColorSeq, xlab="", ylab="", ann=F, type="h", lwd=2,las=1)
  par(new=T)
  plot(1:Length,rep(0,Length), xlim=c(1,Length), ylim=c(ReadsMaxA*1.1-Outspace, ReadsMaxS*1.1), col="black", xlab="", ylab="", axes=F, type="h", lwd=2,las=1)
  abline(h = 0,col="lightblue")
  title(main=paste0(Prefix,".subcodon.summary"))
  mtext(paste(Prefix,"; total ",Length," nt; around STOP",sep=""),cex=0.8)
  title(xlab = paste("Region (",Length," nt)",sep = ""),ylab = "Normalized counts")
  
  Letter=1:Length
  yanno <- ReadsMaxA*1.1-Outspace/2
  if (Length <= 180){
    CharNum=0.6
  }else if (Length > 180){
    CharNum=0.3
  }
  for (i in 1:Length){
    Letter[i] <- as.character(Data$Nucleotide)[i]
    if (i == Middle){
      text(x = i,y=yanno,Letter[i],col = "red",cex=CharNum)
    }else if(i == Middle+1){
      text(x = i,y=yanno,Letter[i],col = "red",cex=CharNum)
    }else if(i == Middle+2){
      text(x = i,y=yanno,Letter[i],col = "red",cex=CharNum)
    }else{
      text(x = i,y=yanno,Letter[i],col = "black",cex=CharNum)
    }
  }
  CDSValues=Data$Mean[1:Middle]
  Frame0=0
  Frame1=0
  Frame2=0
  for (i in seq(1,length(CDSValues),3)){Frame0=Frame0+CDSValues[i]}
  for (i in seq(2,length(CDSValues),3)){Frame1=Frame1+CDSValues[i]}
  for (i in seq(3,length(CDSValues),3)){Frame2=Frame2+CDSValues[i]}
  Sum=Frame0+Frame1+Frame2
  if (Sum > 0){
    Frame0p=format(round(Frame0/Sum*100, 2), nsmall = 2)
    Frame1p=format(round(Frame1/Sum*100, 2), nsmall = 2)
    Frame2p=format(round(Frame2/Sum*100, 2), nsmall = 2)
  }else{
    Frame0p=0
    Frame1p=0
    Frame2p=0
  }
  print(paste0("Frame 0 counts: ",Frame0))
  print(paste0("Frame 1 counts: ",Frame1))
  print(paste0("Frame 2 counts: ",Frame2))
  print(paste0("Frame 0 percentage: ",Frame0p," %"))
  print(paste0("Frame 1 percentage: ",Frame1p," %"))
  print(paste0("Frame 2 percentage: ",Frame2p," %"))
  legend(x="topright", fill = c("red","darkgreen","dodgerblue2"),legend = c(paste0(Frame0p,"%"), paste0(Frame1p,"%"), paste0(Frame2p,"%")),bty = "n")
  #Finishing
  temp <- dev.off()
  print("metaSTOP figure done.")
}else{
  print("Input mode incorrect...")
}

###End of the script