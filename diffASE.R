require("data.table")
require("metap")
require("limma")
require("biomaRt")

#mapping between gene names
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
geneMap = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), mart = ensembl)

#purity file with sample names in first column and purity in second
allPurity = fread("purity.csv")
#for differential expression, cpm matrices for tumor and normal where columns are samples and rows are genes, we assume they are ensembl transcipt ids (as produced by TCGA)
cpmTumor = read.table("cpmTumor.csv", row.names=1)
cpmNormal = read.table("cpmNormal.csv", row.names=1)

#ASE info has the first 7 columns as follows:
#Gene name, sample name,  variant location, normal major read count, normal total read count, tumor major read count, tumor total read count
computeASEbaseline = function(ASEinfo){
	#if multiple samples are being tested together
	samples = unique(ASEinfo$V2)

	#trim sites with no tumor reads
	ASEinfo = ASEinfo[ASEinfo$V7!=0,]

	#compute pvalues per heterozygous site
	ASEraw_pval = c()
	for(x in 1:dim(ASEinfo)[1]){
		downsampledCount = min(ASEinfo$V7[x],100)
		cancer_reads_raw = downsampledCount
		pNormal = ASEinfo$V4[x]/ASEinfo$V5[x]
	    pCancer_raw = ASEinfo$V6[x]/ASEinfo$V7[x]
		ASEraw_pval = c(ASEraw_pval, 2*min(pnorm(pCancer_raw*cancer_reads_raw, pNormal*cancer_reads_raw, sqrt(cancer_reads_raw*pNormal*(1-pNormal))),pnorm(pCancer_raw*cancer_reads_raw, pNormal*cancer_reads_raw, sqrt(cancer_reads_raw*pNormal*(1-pNormal)), lower.tail = FALSE)))
	}

	#Combine p-values, downsampling as needed
	SampSigdown10 = c()
	allTestGenes = c()
	allSamples = c()
	for(sample in samples){
		testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
		allTestGenes = c(allTestGenes, testGenes)
		allSamples = c(allSamples, rep(sample, length(testGenes)))
		for(gene in testGenes){
	        if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))==0){SampSigdown10=c(SampSigdown10,NA)}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))==1){SampSigdown10=c(SampSigdown10,min(ASEraw_pval[ASEinfo$V1==gene & ASEinfo$V2==sample],1))}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))>10){
	            SampSigdown10 = c(SampSigdown10,median(sapply(1:10, function(q) sumlog(sample(sapply(ASEraw_pval[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,1)), 10))[['p']])))
	        }
	        else{SampSigdown10 = c(SampSigdown10,sumlog(sapply(ASEraw_pval[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,1)))[['p']])}
    	}
	}	
    #Compute weighted ASE
    ASEsW100 = c()
    for(sample in samples){
    	testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
    	ASEsW100 = c(ASEsW100, sapply(testGenes, function(gene) weighted.median(abs((ASEinfo$V6/ASEinfo$V7-ASEinfo$V4/ASEinfo$V5)[ASEinfo$V1==gene & ASEinfo$V2==sample]), sapply(ASEinfo$V7[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,100)))))
	}
    data.frame(sample=allSamples, gene=allTestGenes, pval=SampSigdown10, ASE=ASEsW100)
}
computeASEpurity = function(ASEinfo){
	#if multiple samples are being tested together
	samples = unique(ASEinfo$V2)

	#trim sites with no tumor reads
	ASEinfo = ASEinfo[ASEinfo$V7!=0,]

	#get purity scores for each ASE
	asePurity = allPurity$V2[match(ASEinfo$V2,allPurity$V1)]

	#compute purity corrected Rc
	Rc = (ASEinfo$V6/ASEinfo$V7-(1-asePurity)*ASEinfo$V4/ASEinfo$V5)/asePurity
	Rc = sapply(Rc, function(x) max(0,min(x,1)))
	#compute pvalues per heterozygous site
	ASErawPurity_pval = c()
	for(x in 1:dim(ASEinfo)[1]){
		downsampledCount = min(ASEinfo$V7[x],100)
		cancer_reads_rawPurity = asePurity[x]*downsampledCount
		pNormal = ASEinfo$V4[x]/ASEinfo$V5[x]
	    pCancer_rawPurity = Rc[x]
	    ASErawPurity_pval = c(ASErawPurity_pval, 2*min(pnorm(pCancer_rawPurity*cancer_reads_rawPurity, pNormal*cancer_reads_rawPurity, sqrt(cancer_reads_rawPurity*pNormal*(1-pNormal))),pnorm(pCancer_rawPurity*cancer_reads_rawPurity, pNormal*cancer_reads_rawPurity, sqrt(cancer_reads_rawPurity*pNormal*(1-pNormal)), lower.tail = FALSE)))
	}

	#Combine p-values, downsampling as needed
	SampSigdown10 = c()
	allTestGenes = c()
	allSamples = c()
	for(sample in samples){
		testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
		allTestGenes = c(allTestGenes, testGenes)
		allSamples = c(allSamples, rep(sample, length(testGenes)))
		for(gene in testGenes){
	        if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))==0){SampSigdown10=c(SampSigdown10,NA)}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))==1){SampSigdown10=c(SampSigdown10,min(ASErawPurity_pval[ASEinfo$V1==gene & ASEinfo$V2==sample],1))}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))>10){
	            SampSigdown10 = c(SampSigdown10,median(sapply(1:10, function(q) sumlog(sample(sapply(ASErawPurity_pval[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,1)), 10))[['p']])))
	        }
	        else{SampSigdown10 = c(SampSigdown10,sumlog(sapply(ASErawPurity_pval[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,1)))[['p']])}
    	}
	}	
    #Compute weighted ASE
    ASEsW100 = c()
    for(sample in samples){
    	testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
    	ASEsW100 = c(ASEsW100, sapply(testGenes, function(gene) weighted.median(abs((Rc-ASEinfo$V4/ASEinfo$V5)[ASEinfo$V1==gene & ASEinfo$V2==sample]), sapply(ASEinfo$V7[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,100)))))
	}

    data.frame(sample=allSamples, gene=allTestGenes, pval=SampSigdown10, ASE=ASEsW100)
}
computeASEexp = function(ASEinfo){
	#if multiple samples are being tested together
	samples = unique(ASEinfo$V2)

	#trim sites with no tumor reads
	ASEinfo = ASEinfo[ASEinfo$V7!=0,]

	#get purity scores for each site
	asePurity = allPurity$V2[match(ASEinfo$V2,allPurity$V1)]

	#get expression levels for each site
	aseExp = c()
	for(i in 1:length(ASEinfo$V1)){
		aseExp = rbind(aseExp, c(sum(cpmNormal[rownames(cpmNormal) %in% geneMap$ensembl_gene_id[match(ASEinfo$V1[i],geneMap$entrezgene)],colnames(cpmNormal)==ASEinfo$V2[i]]),sum(cpmTumor[rownames(cpmTumor) %in% geneMap$ensembl_gene_id[match(ASEinfo$V1[i],geneMap$entrezgene)],colnames(cpmTumor)==ASEinfo$V2[i]])))
	}

	#compute fraction of transcripts 
	ft = (aseExp[,2] - (1 - asePurity)*aseExp[,1])/(aseExp[,2])
	ft = sapply(ft, function(x) max(0,min(x,1)))

	#compute purity corrected Rc
	Rc = (ASEinfo$V6/ASEinfo$V7-(1-ft)*ASEinfo$V4/ASEinfo$V5)/ft
	Rc = sapply(Rc, function(x) max(0,min(x,1)))
	#compute pvalues per heterozygous site
	ASErawExp_pval = c()
	for(x in 1:dim(ASEinfo)[1]){
		downsampledCount = min(ASEinfo$V7[x],100)
		cancer_reads_rawExp = ft[x]*downsampledCount
		pNormal = ASEinfo$V4[x]/ASEinfo$V5[x]
	    pCancer_rawExp = Rc[x]
	    ASErawExp_pval = c(ASErawExp_pval, 2*min(pnorm(pCancer_rawExp*cancer_reads_rawExp, pNormal*cancer_reads_rawExp, sqrt(cancer_reads_rawExp*pNormal*(1-pNormal))),pnorm(pCancer_rawExp*cancer_reads_rawExp, pNormal*cancer_reads_rawExp, sqrt(cancer_reads_rawExp*pNormal*(1-pNormal)), lower.tail = FALSE)))
	}

	#Combine p-values, downsampling as needed
	SampSigdown10 = c()
	allTestGenes = c()
	allSamples = c()
	for(sample in samples){
		testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
		allTestGenes = c(allTestGenes, testGenes)
		allSamples = c(allSamples, rep(sample, length(testGenes)))
		for(gene in testGenes){
	        if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample & !is.na(ASErawExp_pval)))==0){SampSigdown10=c(SampSigdown10,NA)}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample & !is.na(ASErawExp_pval)))==1){SampSigdown10=c(SampSigdown10,min(ASErawExp_pval[ASEinfo$V1==gene & ASEinfo$V2==sample & !is.na(ASErawExp_pval)],1))}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample & !is.na(ASErawExp_pval)))>10){
	            SampSigdown10 = c(SampSigdown10,median(sapply(1:10, function(q) sumlog(sample(sapply(ASErawExp_pval[ASEinfo$V1==gene & ASEinfo$V2==sample & !is.na(ASErawExp_pval)], function(x) min(x,1)), 10))[['p']])))
	        }
	        else{SampSigdown10 = c(SampSigdown10,sumlog(sapply(ASErawExp_pval[ASEinfo$V1==gene & ASEinfo$V2==sample & !is.na(ASErawExp_pval)], function(x) min(x,1)))[['p']])}
    	}
	}	
    #Compute weighted ASE
    ASEsW100 = c()
    for(sample in samples){
    	testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
    	ASEsW100 = c(ASEsW100, sapply(testGenes, function(gene) weighted.median(abs((Rc-ASEinfo$V4/ASEinfo$V5)[ASEinfo$V1==gene & ASEinfo$V2==sample]), sapply(ASEinfo$V7[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,100)))))
	}
    data.frame(sample=allSamples, gene=allTestGenes, pval=SampSigdown10, ASE=ASEsW100)
}
computeTumorSampleASE = function(ASEinfo){
	#if multiple samples are being tested together
	samples = unique(ASEinfo$V2)

	#trim sites with no tumor reads
	ASEinfo = ASEinfo[ASEinfo$V7!=0,]

	#compute pvalues per heterozygous site
	ASEraw_pval = c()
	for(x in 1:dim(ASEinfo)[1]){
		downsampledCount = min(ASEinfo$V7[x],100)
		cancer_reads_raw = downsampledCount
		pNormal = .5
	    pCancer_raw = ASEinfo$V6[x]/ASEinfo$V7[x]
		ASEraw_pval = c(ASEraw_pval, 2*min(pnorm(pCancer_raw*cancer_reads_raw, pNormal*cancer_reads_raw, sqrt(cancer_reads_raw*pNormal*(1-pNormal))),pnorm(pCancer_raw*cancer_reads_raw, pNormal*cancer_reads_raw, sqrt(cancer_reads_raw*pNormal*(1-pNormal)), lower.tail = FALSE)))
	}

	#Combine p-values, downsampling as needed
	SampSigdown10 = c()
	allTestGenes = c()
	allSamples = c()
	for(sample in samples){
		testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
		allTestGenes = c(allTestGenes, testGenes)
		allSamples = c(allSamples, rep(sample, length(testGenes)))
		for(gene in testGenes){
	        if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))==0){SampSigdown10=c(SampSigdown10,NA)}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))==1){SampSigdown10=c(SampSigdown10,min(ASEraw_pval[ASEinfo$V1==gene & ASEinfo$V2==sample],1))}
	        else if(length(which(ASEinfo$V1==gene & ASEinfo$V2==sample))>10){
	            SampSigdown10 = c(SampSigdown10,median(sapply(1:10, function(q) sumlog(sample(sapply(ASEraw_pval[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,1)), 10))[['p']])))
	        }
	        else{SampSigdown10 = c(SampSigdown10,sumlog(sapply(ASEraw_pval[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,1)))[['p']])}
    	}
	}	
    #Compute weighted ASE
    ASEsW100 = c()
    for(sample in samples){
    	testGenes = unique(ASEinfo$V1[ASEinfo$V2==sample])
    	ASEsW100 = c(ASEsW100, sapply(testGenes, function(gene) weighted.median(abs((ASEinfo$V6/ASEinfo$V7-.5)[ASEinfo$V1==gene & ASEinfo$V2==sample]), sapply(ASEinfo$V7[ASEinfo$V1==gene & ASEinfo$V2==sample], function(x) min(x,100)))))
	}
    data.frame(sample=allSamples, gene=allTestGenes, pval=SampSigdown10, ASE=ASEsW100)
}

#Simulate ASE reads at a heterozyogus site
generateReads = function(Rc, Rn=.5, p=.8, en=1, ec=1, AF=1, numReads=20, numCells=10000, errorRate=1/10000){
	trueASE = abs(Rc-Rn)
	cancerCells = rbinom(1, numCells, p)
	mutatedCells = rbinom(1, cancerCells, min(AF/p,1))
	unmutatedCells = cancerCells-mutatedCells
	normalCells = numCells-cancerCells
	mutatedRna = mutatedCells*ec
	unmutatedRna = unmutatedCells*en
	normalRna = normalCells*en
	mutatedRnaMajor = rbinom(1,round(mutatedRna),Rc)
	unmutatedRnaMajor = rbinom(1,round(unmutatedRna),Rn)
	normalRnaMajor = rbinom(1,round(normalRna),Rn)
	mutatedRnaMajor = rbinom(1,mutatedRnaMajor,1-errorRate)+rbinom(1,mutatedRna-mutatedRnaMajor,errorRate)
	unmutatedRnaMajor = rbinom(1,unmutatedRnaMajor,1-errorRate)+rbinom(1,unmutatedRna-unmutatedRnaMajor,errorRate)
	normalRnaMajor = rbinom(1,normalRnaMajor,1-errorRate)+rbinom(1,normalRna-normalRnaMajor,errorRate)
	Rs = (mutatedRnaMajor+unmutatedRnaMajor+normalRnaMajor)/(mutatedRna+unmutatedRna+normalRna)
	majorReads = rbinom(1,numReads,Rs)
	normalMajorReads = rbinom(1,round(Rn*numReads),1-errorRate)+rbinom(1,numReads-round(Rn*numReads),errorRate)
	c(normalMajorReads, numReads, majorReads, numReads, trueASE)
}

#simple test
testASE = data.frame(1,"sample1", "mutLoc", t(generateReads(.7, Rn=.6, p=.8, en=2)[1:4]))
colnames(testASE) = paste0("V",1:7)
aseBaseline = computeASEbaseline(testASE)
allPurity = data.frame(V1=c("sample1"), V2=c(.8))
asePurity = computeASEpurity(testASE)
cpmTumor = data.frame("sample1"=1)
rownames(cpmTumor) = "ENSG00000121410"
cpmNormal = data.frame("sample1"=2)
rownames(cpmNormal) = "ENSG00000121410"
aseExp = computeASEexp(testASE)
aseTumorSample = computeTumorSampleASE(testASE)