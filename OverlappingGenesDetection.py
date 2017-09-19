# -*- coding: utf-8 -*-

import time

def LambdaModule(pos, variable):
	return lambda x: x[pos] == variable

def Accession2GeneNamesParser(infileName = "Human_GeneNames.tsv", infilesPath = "InputFiles"):
	infile = open("{}/{}".format(infilesPath, infileName),'r')
	result = {}
	for row in infile:
		tmp = row.strip().split('\t')
		result[tmp[0]] = tmp[1]
	infile.close()
	return result

def RefSeqTranscriptsParser(Accession2GeneNames, infileName = "Human_RefSeqTranscripts.tsv", infilesPath = "InputFiles"):
	infile = open("{}/{}".format(infilesPath, infileName),'r')
	RefSeqTranscripts = {}
	for row in infile:
		tmp = row.strip().split('\t')
		RefSeqTranscripts[(Accession2GeneNames[tmp[3]], tmp[3], tmp[0], tmp[5], int(tmp[1]), int(tmp[2]))] = 1 #GeneName, AccesionNumber, Chromosome, strand, TranscriptStart, TranscriptEnd
	infile.close()
	return RefSeqTranscripts

def PrepareRepresentativeGeneCoords(GeneNames, RefSeqTranscripts):
	TmpAllGenes = {}
	for transcript in RefSeqTranscripts: #example "transcript": ('C16orf57', 'NM_001195302', 'chr16', '+', 58035276, 58055527)
		try:
			if (TmpAllGenes[transcript[0]][3] > int(transcript[4])): #check if new 'gene start' is smaller then the old one
				TmpAllGenes[transcript[0]] = TmpAllGenes[transcript[0]][0], TmpAllGenes[transcript[0]][1], TmpAllGenes[transcript[0]][2], int(transcript[4]), TmpAllGenes[transcript[0]][4]
			if (TmpAllGenes[transcript[0]][4] < int(transcript[5])): #check if new 'gene end' is smaller then the old one
				TmpAllGenes[transcript[0]] = TmpAllGenes[transcript[0]][0], TmpAllGenes[transcript[0]][1], TmpAllGenes[transcript[0]][2], TmpAllGenes[transcript[0]][3], int(transcript[5])
		except KeyError:
			TmpAllGenes[transcript[0]] = transcript[0], transcript[2], transcript[3], int(transcript[4]), int(transcript[5])
	AllGenes = {}
	for key in TmpAllGenes.keys():
		AllGenes[TmpAllGenes[key]] = 1
	return AllGenes

def DBTSS_DataParser(Species, line, Accession2GeneNames, ReprGeneCoords, infilesPath = "InputFiles", reportsPath = "reports"):
	infile = open("{}/DBTSS_{}_{}.tab".format(infilesPath, Species, line), 'r')
	reportOutfile = open("{}/DBTSS_DataParser_REPORT_{}_{}.tab".format(reportsPath, Species, line),'w')
	TSSs = {}
	for row in infile:
		tmp = row.strip().split('\t')
		AccNums = tmp[4].split(',')
		if (tmp[1] != ''):
			if ((int(tmp[1]) == 0) and (tmp[2] == "NM_ correspond")): #here we choose only TSS with confident status (DBTSS status == 0), that were assigned to mRNAs (Type == "NM_ correspond")
				for elem in AccNums:
					if (elem != ''):
						switch = 1
						DBTSS_ID = int(tmp[0])
						TagPattern = int(tmp[1])
						TagType = tmp[2]
						GroupID = int(tmp[3])
						
						try:
							GeneName = Accession2GeneNames[elem]
						except KeyError: #Key error in this place means that there is no accesion number in the data taken from UCSC - these are rare cases
							switch = 0
							reportOutfile.write(elem + '\t"Accession2GeneNames[elem]" KeyError -> called Accession number from DBTSS data is not present in our list of genes/transcripts\n')
						
						try:
							check = filter(LambdaModule(0,Accession2GeneNames[elem]), ReprGeneCoords.keys())[0][4]
						except IndexError:
							switch = 0
							reportOutfile.write(elem + '\t"GeneStart = filter(LambdaModule(0,Accession2GeneNames[elem]), AllGenes.keys())[0][3]" IndexError -> RefSeqTranscripts are missing AccNumber for a gene\n')
						except KeyError:
							switch = 0
							reportOutfile.write(elem + '\t"GeneStart = filter(LambdaModule(0,Accession2GeneNames[elem]), AllGenes.keys())[0][3]" KeyError -> RefSeqTranscripts are missing AccNumber for a gene\n')
						
						try: ### I had to add this exception, because there are empty fields of unknown origin in the DBTSS files. This is only in one cell line - DLD1NormoxiaHIF2minus in the following row: 376174	0	NM_ correspond	11096	NM_001191005,NM_001191006,NM_001191007,NM_001191009,NM_006625,NM_054016		2		chr1	R	24021556	7	8349807	24021556	24021566					
							AlterPromoterID = int(tmp[7])
						except ValueError:
							switch = 0
							reportOutfile.write(elem + '\t\t\t VelueError -> There is gap ('') instead of alternative promoter id column (column 7 (by python counting))\n')
						
						if switch:
							RefseqName = elem
							AntiSenseRefseqName = tmp[5]
							PromoterID = int(tmp[6])
							AlterPromoterID = int(tmp[7])
							Chromosome = tmp[8]
							if (Species == "Human"):
								if (tmp[9] == 'R'): # in previous versions reverse strand was designated as 0 (as it is described now in mouse datasets)
									Strand = '-'
								elif (tmp[9] == 'F'): # in previous versions positive was described as 1 (as it is described now in mouse datasets)
									Strand = '+'
							elif (Species =="Mouse"):
								if (int(tmp[9]) == 0):
									Strand = '-'
								elif (int(tmp[9]) == 1):
									Strand = '+'
							RepresentTSS = int(tmp[10])
							TagNum = int(tmp[11])
							TotalTagNum = int(tmp[12])
							TagsPPM = ((TagNum * 1000000)/(TotalTagNum * 1.0))
							try:
								TSSStart = int(tmp[13])
							except IndexError:
								TSSStart = 0
							try:
								TSSEnd = int(tmp[14])
							except IndexError:
								TSSEnd = 0
							GeneStart = filter(LambdaModule(0,Accession2GeneNames[elem]), ReprGeneCoords.keys())[0][3]
							GeneEnd = filter(LambdaModule(0,Accession2GeneNames[elem]), ReprGeneCoords.keys())[0][4]
						
							InternalSwitch = 0
							if (Strand == '+'):
								if ((RepresentTSS  >= (GeneStart - 5000)) and (RepresentTSS <= GeneEnd)):
									InternalSwitch = 1
								else:
									reportOutfile.write(row.strip() + '\t"InternalSwitch" Error -> Representant TSS position was rejected as being further then 5 kb from the assigned gene\n')
							elif (Strand == '-'):
								if ((RepresentTSS <= (GeneEnd + 5000)) and (RepresentTSS >= GeneStart)):
									InternalSwitch = 1
								else:
									reportOutfile.write(row.strip() + '\t"InternalSwitch" Error -> Representant TSS position was rejected as being further then 5 kb from the assigned gene\n')
							if InternalSwitch:
								TSSs[(DBTSS_ID, TagPattern, TagType, GroupID, GeneName, RefseqName, AntiSenseRefseqName, PromoterID, AlterPromoterID, Chromosome, Strand, RepresentTSS, TagNum, TotalTagNum, TagsPPM, TSSStart, TSSEnd, GeneStart, GeneEnd)] = 1
	reportOutfile.close()
	infile.close()
	return TSSs

def UnitParser(Species, ReprGeneCoords, TSSs, CUTOFF = 5):
#	"Unit" is a gene with its TSSs. It is created by taking the annotated END of the gene and the "first" TSS.
#	These units will be later on searched for the overlapping pars in an easy way.
#	each Unit will have the following fields:
#		Name		:	GeneName
#		Chromosome	:
#		Strand	:
#		UnitStart	:	Dependently on the strand:
#					-	Positive strand	- 	The first TSS (with the smallest coordinate number)
#					-	Negative strand	- 	Gene start (as in fact this is the gene end -> on the genome, 
#										genes are annotated FROM (start) TO (END), what is later on interpreted 
#										by the gene start or the gene end (dependetly on the strand))
#		UnitEnd	:	Dependently on the strand:
#					-	Positive strand	- 	Gene end
#					-	Negative strand	- 	TSS with the highest coordinate number -> the first TSS
#		TSSsPos	:	All alternative start sites annotated positions (ordered from the smallest to the highest "position" number)
#					These parameters are stored in tuple
#		TSSExpLevel	:	Expression levels assigned to each TSS. These parameters are stored in tuple 
#					Order is the same as of the TSSs in the above field.
#		UnitExpLvl	:	Total Expression level of the Unit
#		CUTOFF	:	Cutoff for the minimal expression level of the TSSs
#		GeneStart	:
#		GeneEnd	:
#	
#		Example output of "Units" variable:
#			
#
#	Warning:
#		If some gene has many alternative transcripts, TSSs with the same representative position, 
#		and the same expression level are assigned multiple times, so if we want to compute total 
#		expression level we cannot simply summ up all the TSSs expression levels, we must check which 
#		one are repeating and sum up only unique TSS positions
	
	Units = {}
	for gene in ReprGeneCoords:
		GeneData = filter(LambdaModule(4,gene[0]), TSSs.keys())
		TSSPositions = {} #this dictionary will store the TSS position (key), and corresponding expression level (value)
		for elem in GeneData:
			if (elem[14] >= CUTOFF):	# Here we filter out all TSSs with the exp level lower then CUTOFF
				TSSPositions[elem[11]] = elem[14]
		
		if (len(TSSPositions) > 0):
			TSSPositionsList = TSSPositions.keys()
			TSSPositionsList.sort()
			
			Name = gene[0]
			Chromosome = gene[1]
			Strand = gene[2]
			
			if (Strand == '+'):
				UnitStart = TSSPositionsList[0]
			elif (Strand == '-'):
				UnitStart = gene[3]
			
			if (Strand == '+'):
				UnitEnd = gene[4]
			elif (Strand == '-'):
				UnitEnd = TSSPositionsList[-1]
			GeneStart = gene[3]
			GeneEnd = gene[4]
			TSSsPos = tuple(TSSPositionsList)
			
			TSSExpLevelList = []
			for TSS in TSSPositionsList:
				TSSExpLevelList.append(TSSPositions[TSS])
			TSSExpLevel	= tuple(TSSExpLevelList)
			
			UnitExpLvl = 0
			for TSS in TSSPositionsList:
				UnitExpLvl += TSSPositions[TSS]
			
			Units[(Name, Chromosome, Strand, UnitStart, UnitEnd, TSSsPos, TSSExpLevel, UnitExpLvl, CUTOFF, GeneStart, GeneEnd)] = 1

	return Units

def Pos2NegOverlap_LambdaModule(StrandField, ChromosomeField, GOI_ChromosomeVar, UnitEndField, GOI_UnitStartVar, GeneStartField, GOI_GeneStart, GeneEndField, GOI_GeneEnd): # GOI stands for "Gene Of Interest"
	return lambda x: x[StrandField] == '-' and x[ChromosomeField] == GOI_ChromosomeVar and x[UnitEndField] > GOI_UnitStartVar and x[GeneStartField] < GOI_GeneStart and x[GeneEndField] < GOI_GeneEnd

def Neg2PosOverlap_LambdaModule(StrandField, ChromosomeField, GOI_ChromosomeVar, UnitStartField, GOI_UnitEndVar, GeneStartField, GOI_GeneStart, GeneEndField, GOI_GeneEnd): # GOI stands for "Gene Of Interest"
	return lambda x: x[StrandField] == '+' and x[ChromosomeField] == GOI_ChromosomeVar and x[UnitStartField] < GOI_UnitEndVar and x[GeneStartField] > GOI_GeneStart and x[GeneEndField] > GOI_GeneEnd

def OverlapParser(Species, RefSeqTranscripts, ReprGeneCoords, Units, CUTOFF = 5):
#	This function will find overlapping regions in "Units". 
#	Searching is made gene by gene (unit by unit), finding all genes with the particular gene overlap --> simple FOR loop is used
#	To consider two genes overlapping, the following conditions must be fulfilled:
#		I.	Genes must be located on the opposite STRANDS
#		II.	Genes must be located on the same CHROMOSOME
#		III.	START coordinate of the gene (unit) on positive strand must be bigger than START of the gene on negative DNA strand
#		IV.	END coordinate of the gene (unit) on positive strand must be bigger than END of the gene on negative DNA strand
#		V.	START coordinate of the gene (unit) on positive strand must be smaller then the END of the gene on negative DNA strand
#	
#	While going gene by gene, we will alternately have 'positive' or 'negative' genes, 
#	so the above conditions will be adjusted accordingly to the situation.
#	
#			          >>>>>>		Gene (unit) on positive DNA strand	
#			         /\v.  /iv.
#			        /  \  /
#			   iii./    \/						
#			       <<<<<<		Gene (unit) on negative DNA strand
#
#	Parameterss description:
#		PairNumber		:	integer	
#		GeneName		:	string			--> 	eg. 'NUDC'
#		AccNums		:	tuple	with strings	--> 	eg. ('NM_000213','NM_732134')
#		Chromosome		:	string			--> 	eg. 'chrX'
#		Strand		:	string			--> 	'+' or '-'
#		TSSPositions	:	tuple with integers	--> 	eg. (12534424,12537853,12540981)
#		TSSExprLvls		:	tuple with floats 	--> 	eg. (432.6,10.3)
#		GeneStart		:	integer
#		GeneEnd		:	integer
#		UnitStart		:	integer
#		UnitEnd		:	integer
#		ExpressionLvl	:	float
#		OvrGenes		:	tuple with strings	--> 	eg. ('SPOCK1','SPOCK2')
#		TOR			:	float				--> 	stands for 'Theoretical Overlap Ratio', to count it, one check 
#											how many TSSs overlap with the other gene's TSSs, and divide this 
#											number by total number of TSSs of the particular gene.
#		OR			:	float				--> 	stands for 'Overlap Ratio', it is computed similarily to TOR, but 
#											instead number of TSSs, one sum up expression levels assigned to 
#											the TSSs overlapping with the TSSs from the opposite gene, and 
#											divide this number by the total expression level of the gene
#		PTO			:	float				--> 	stands for 'Preference To Overlap', to compute it, one must divide 
#											OR by TOR
#		CUTOFF		:	integer			-->	Minimal expression level assigned for each TSS

	
	OverlappingGenes = {}
	OverlappingGeneNumber = 0
	for unit in Units:
		Switch = 0		#"switch" is activated when there is an overlapping gene matching our criteria
		if (unit[2] == "+"):
			OvrGene = filter(Pos2NegOverlap_LambdaModule(2, 1, unit[1], 4, unit[3], 9, unit[9], 10, unit[10]), Units.keys())
			if (OvrGene != []):
				Switch = 1
		elif (unit[2] == "-"):
			OvrGene = filter(Neg2PosOverlap_LambdaModule(2, 1, unit[1], 3, unit[4], 9, unit[9], 10, unit[10]), Units.keys())			
			if (OvrGene != []):
				Switch = 1
		if Switch:
			OverlappingGeneNumber += 1
			GeneName = unit[0]
			
			AccNumsData = filter(LambdaModule(0,unit[0]), RefSeqTranscripts.keys())
			AccNumsDict = {}
			for elem in AccNumsData:
				AccNumsDict[elem[1]] = 1
			AccNums = tuple(AccNumsDict.keys())
			
			Chromosome = unit[1]
			Strand = unit[2]
			TSSPositions = unit[5]
			TSSExprLvls = unit[6]
			GeneStart = filter(LambdaModule(0,unit[0]), ReprGeneCoords.keys())[0][3]
			GeneEnd = filter(LambdaModule(0,unit[0]), ReprGeneCoords.keys())[0][4]
			UnitStart = unit[3]
			UnitEnd = unit[4]
			ExpressionLvl = unit[7]
			
			OvrGenesDict = {}
			for elem in OvrGene:
				OvrGenesDict[elem[0]] = 1
			OvrGenes = tuple(OvrGenesDict.keys())
			
			StartSites = {}
			OvrTSS = 0
			OvrTSSExpLvl = 0
			for position, ExpLvl in zip(unit[5], unit[6]):
				StartSites[position] = ExpLvl
			if (unit[2] == "+"):
				MaxBorder = OvrGene[0][4]
				for elem in OvrGene:
					if (MaxBorder < elem[4]):
						MaxBorder = elem[4]
				for tss in unit[5]:
					if (tss < MaxBorder):
						OvrTSS += 1
						OvrTSSExpLvl += StartSites[tss]
			elif (unit[2] == "-"):
				MinBorder = OvrGene[0][3]
				for elem in OvrGene:
					if (MinBorder > elem[3]):
						MinBorder = elem[3]
				for tss in unit[5]:
					if (tss > MinBorder):
						OvrTSS += 1
						OvrTSSExpLvl += StartSites[tss]
			TOR = OvrTSS / (len(unit[5]) * 1.0)
			OR = OvrTSSExpLvl / (ExpressionLvl * 1.0)
			PTO = OR / TOR
			CUTOFF = unit[8]
			OverlappingGenes[(OverlappingGeneNumber, GeneName, AccNums, Chromosome, Strand, TSSPositions, TSSExprLvls, GeneStart, GeneEnd, UnitStart, UnitEnd, ExpressionLvl, OvrGenes, TOR, OR, PTO, CUTOFF)] = 1
	return OverlappingGenes

def OverlappingPairsDump(Species, line, OverlappingPairs, CUTOFF):
	dump = open("OverlappingGenePairs_{}_{}_{}ppmCUTOFF.tsv".format(Species, line, CUTOFF), 'w')
	dump.write("Gene pair number	Gene name	Chromosome	Strand	Representative gene start	Representative gene end	Positions of TSSs (comma-separated)	Expression level of individual TSSs (comma-separated)	Gene expression level	Overlap ratio (OR)	\n")
	KeysList = OverlappingPairs.keys()
	KeysList.sort()
	for key in KeysList:
		for elem in OverlappingPairs[key]:
			dump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, elem[1], elem[3], elem[4], elem[9], elem[10], elem[5], elem[6], elem[11], elem[14]))
	dump.close()

def OverlappingPairsParser(Species, line, OverlappingGenes, CUTOFF = 5):
#	Overlapping genes parser is here to 'join' overlapping genes generated by the 'Overlap' function into pairs
	DictOfOverlaps = {} # {GeneName : (Original entry (row) from OverlappingGenes dict)}
	for gene in OverlappingGenes:
		DictOfOverlaps[gene[1]] = gene
	
	TMPOverlappingPairs = {}
	for gene in DictOfOverlaps:
		GOI = DictOfOverlaps[gene] #GOI stands for "Gene Of Interest"
		for elem in GOI[12]:
			if (GOI[4] == '+'):
				TMPOverlappingPairs[(GOI, DictOfOverlaps[elem])] = 1
			elif (GOI[4] == '-'):
				TMPOverlappingPairs[(DictOfOverlaps[elem], GOI)] = 1
			else:
				print "#!#\tOverlappingPairsParser Error (for %s %s with CUTOFF = %i) --> This message should never be here! something's wrong..:/" % (Species, tissue, CUTOFF)

	counter = 0
	OverlappingPairs = {}
	for pair in TMPOverlappingPairs:
		counter += 1
		OverlappingPairs[counter] = pair
	
	OverlappingPairsDump(Species, line, OverlappingPairs, CUTOFF)

def main(Species, line, CUTOFF):
	start = time.time()
	print "Overlapping genes will be detected for {} {} cell line with the minimal expression level set to {} ppm.".format(Species, line, CUTOFF)
	Accession2GeneNames = Accession2GeneNamesParser()
	RefSeqTranscripts = RefSeqTranscriptsParser(Accession2GeneNames)
	ReprGeneCoords = PrepareRepresentativeGeneCoords(Accession2GeneNames, RefSeqTranscripts)
	print "The next step is parsing of the TSS details file from DBTSS, which is the longest step that may take up to few hours."
	TSSs = DBTSS_DataParser(Species, line, Accession2GeneNames, ReprGeneCoords)
	Units = UnitParser(Species, ReprGeneCoords, TSSs)
	OverlappingGenes = OverlapParser(Species, RefSeqTranscripts, ReprGeneCoords, Units)
	OverlappingPairs = OverlappingPairsParser(Species, line, OverlappingGenes, CUTOFF)
	print "Overlapping gene pairs detection was sucessfully completed for {} {} cell line with the minimal expression level set to {} ppm.\nIt took {} minutes.".format(Species, line, CUTOFF, (time.time()-start)/60)
	
Species = "Human" # Either "Human" or "Mouse"
line = "AdultTestis" # Name of the studied TSS-Seq library / tissue / cell line. 
CUTOFF = 5 # Minimal TSS-Seq expression level in parts per million (PPM)

main(Species, line, CUTOFF)