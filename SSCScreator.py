#!/usr/bin/python3 -tt
'''

SSCScreator
By Oddmund Nordgård, Stavanger University Hospital. 
Including patches from the DuplexSequencing project by Kennedy and colleagues. See file LICENCE for more information about licensing.

See SSCScreator.CHANGELOG for more information about changes

Written for Python 3
Required modules: Pysam, Samtools

Inputs: 
    A position-sorted single-end BAM file containing reads with a molecular tag in the header. All reads must be oriented in the same direction as the reference genome. The reverse reads must be aligned to a reverse copy of the reference genome. Can also read from stdin, if - is given as filename.

Outputs:
    1: A single-end FASTQ file containing SSCSs
    2: A single-end BAM file containing reads with bad tag, unmapped (or softclipped if s filter)
    3: A tagcounts file
  
The program starts at the position of the first good read, determined by the type of read specified on startup.  It then goes through the file until it finds a new position, saving all reads as it goes, in families according to molecular tag.   When it finds a new position, it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are written to the output FASTQ file.  After making consensuses with the reads from the first position, it continues on through the origional file until it finds another new position, sends those reads to the consensus maker, and so on until the end of the file. The SSCScreator function is only called if there is more than minmem reads in a family. It also handles the read lenght and stops calling consensus bases when there is less than minmem sequences long enough.This version is made only or single end reads. In contrast to the ConsensusMaker of Kennedy, this program collects information for all reads in a family, not only those with the most common cigar string. Both a specific cigar code and a base has to be present in more than cutoff fraction of the reads to be included in the consensus.

Outputs both a FASTQ file and a _filtereout.bam BAM file with all filtered reads. Mainly because they have not aligned. But also softclipped, if we apply the softclipping filter. In addition, if we use the --tagprefix variable, we output some tag statics files. The sequence name of the reads in the FASTQ file consist of the position, the tag and the family size

usage: SSCScreator.py [-h] --infile INFILE [--tagprefix TAGFILE]
                             --outfile OUTFILE [--rep_filt REP_FILT]
                             [--minmem MINMEM] [--maxmem MAXMEM]
                             [--cutoff CUTOFF] [--Ncutoff NCUTOFF]
                             [--read_out ROUT] [--filt FILT]

optional arguments:
  -h, --help           show this help message and exit
  --infile INFILE      input BAM file. Stdin if -.
  --tagprefix TAGFILE  Prefix for output tagcounts (.tags) and tag family size
                       (.familysizes) files
  --outfile OUTFILE    output FASTQ file. Stdout if -.
  --rep_filt REP_FILT  Remove tags with homomeric runs of nucleotides of
                       length x. [9]
  --minmem MINMEM      Minimum number of reads allowed to comprise a
                       consensus. [3]
  --maxmem MAXMEM      Maximum number of reads allowed to comprise a
                       consensus. [1000]
  --cutoff CUTOFF      Percentage of nucleotides at a given position in a read
                       that must be identical in order for a consensus to be
                       called at that position. [0.7]
  --Ncutoff NCUTOFF    With --filt 'n', maximum fraction of Ns allowed in a
                       consensus [1.0]
  --read_out ROUT      How often you want to be told what the program is
                       doing. [1000000]
  --filt FILT          A string indicating which filters should be
                       implemented. Filters: s: Softclipping filter. n: N
                       filter. None: No filter applied ['sn']
  --remove_false       A filter that compares molecular barcodes and removes SSCS 
                       families that is caused by seq errors in mol barcode



Details of different arguments:
    --minmem and --maxmem set the range of family sizes (constrained by cigar score) that can be used to make a consensus sequence.  Examples use --minmem of 3 and --maxmem of 1000
        Example 3: 
            A family with over 1000 members exists.  A random sample of 1000 reads from that family is used to make a SSCS.
    --cutoff sets the strictness of the consensus making.    
        Example (--cutoff = 0.7):
            Four reads (readlength = 10) are as follows:
                Read 1: ACTGATACTT
                Read 2: ACTGAAACCT
                Read 3: ACTGATACCT
                Read 4: ACTGATACTT
            The resulting SSCS is:
                ACTGATACNT
    --Ncutoff, with --filt n enabled, sets the maximum percentage of Ns allowed in a SSCS.  
        Example (--Ncutoff = .1, --readlength = 20):
            Two SSCSs are generated as follows:
                SSCS 1: ACGTGANCTAGTNCTNTACC
                SSCS 2: GATCTAGTNCATGACCGATA
SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, while SSCS 1 does not with 3/20 = 15% Ns.
    --filt sets which filters are used.  Options are: 		
		s: Softclipping filter.  Filters out any reads which have been soft-clipped in alignment.  This avoids later problems with hard-clipping.  
		n: N filter. Filters out consensus sequences with a higher percentage of Ns than the threshold imposed by --Ncutoff.  Without this option, --Ncutoff doesn't do anything.  

'''

import sys
import pysam
import re
import random
from collections import defaultdict
from argparse import ArgumentParser
import pdb
import array
import numpy
import editdistance
import time

###################################################
# CLASS DEFINITIONS 
###################################################

class Cigar:
    cigartuples = []  #  A list of (c,n) tuples storing the code c and the number of bases n
    baseindex = -1 # Current base position in the sequence the cigar belongs to. First base is 0
    cigarpos = 0 # Current cigar position (code). First code is 0
    cigarcounter = 0 #To count down the number of bases within a cigarpos
    previous_cigarpos = 0  # To keep the previous cigarpos
    
    # Init function that creates the new cigar based on the input cigartuples
    def __init__(self,cigtup):
        self.cigartuples = cigtup
        self.cigarcounter = self.cigartuples[0][1]  # initalize the cigarcounter
        self.cigarpos = 0 # The position in the cigar string (tuple number)
        self.baseindex=-1
        self.previous_cigarpos = 0 

    # Function that returns the next cigarcode tuple, i.e. the cigar code for the next base
    # in the corresonding sequence. Return False when no more bases are available. E.g. (M,1)
    def nextCigarPos(self):
        # We start by storing the present position in previous_position, in case we must step
        # back later
        self.previous_cigarpos = self.cigarpos
        self.previous_cigarcounter = self.cigarcounter
        if self.cigarpos < (len(self.cigartuples)):  # Check that we still have more cigar codes
            if self.cigarcounter>0:   #As long as we still have more bases in a cig. code
                if self.cigartuples[self.cigarpos][0] != 2: #if not a D
                    self.baseindex += 1   # Just increase the base position (in sequence) if we don't have a deletion
                self.cigarcounter -= 1
#                print("Cigarcounter in if: ",self.cigarcounter)
                return((self.cigartuples[self.cigarpos][0],self.baseindex))
            else:  #Jump to the next cigar code
                self.cigarpos += 1
                if self.cigarpos < len(self.cigartuples):
                    self.cigarcounter = self.cigartuples[self.cigarpos][1]
                    self.cigarcounter -= 1   #Count down the counter when used
                    if self.cigartuples[self.cigarpos][0]!=2: # Only if not deletion
                        self.baseindex += 1
                    return((self.cigartuples[self.cigarpos][0],self.baseindex))
                else: # If we have reached the end of the cigar
                    return((-1,-1))
        else:
            return((-1,-1))  # Return (-1,-1) if we reach the end of the cigar

    #Function that computes the length of the aligned part of the cigar (Ms and Is)    
    def getAlignedLength(self):
        length=0
        for ct in self.cigartuples:
            if (ct[0]==0 or ct[0]==1):
                length += ct[1]
        return(length)

    # Function that steps the counters for the current cigar one step back
    def oneStepBack(self):
        self.cigarpos = self.previous_cigarpos
        self.cigarcounter = self.previous_cigarcounter
        self.baseindex -= 1
        
                    





###################################################
# FUNCTION DEFINITIONS
###################################################


def consensusMaker (cigar_list,sequence_list,quality_list,fam_cutoff,cutoff) :
    '''The consensus maker uses a simple "majority rules" algorithm to qmake a consensus at each base position.  If no nucleotide majority reaches above the minimum theshold (--cutoff), the position is considered undefined and an 'N' is placed at that position in the read.'''
    nucIdentityList=[0, 0, 0, 0, 0, 0] # In the order of T, C, G, A, N, Total
    nucKeyDict = {0:'T', 1:'C', 2:'G', 3:'A', 4:'N'}
    consensusRead = ''
    consensusQual = ''  #Store the consensus base qualities, which is the mean of the matched bases' qualities in a posision
    consensusNoBases = []  # Store the number of bases supporting a posistion in a SSCS
    consensusNoA = [] # Store the number of As supporting a position in a SSCS
    consensusNoC = [] # Store the number of Cs supporting a position in a SSCS
    consensusNoG = [] # Store the number of Gs supporting a position in a SSCS
    consensusNoT = [] # Store the number of Ts supporting a position in a SSCS
    consensusNoN = [] # Store the number of Ns supporting a position in a SSCS

    consensuslength=0
    seqlengths = []
    seqnum = len(cigar_list)

    

    #Determine the concensusread length = the family_cutoff highest length
    for cig in cigar_list:
        seqlengths.append(cig.getAlignedLength())

    # Set consensuslength to the cutoff longest read
    for i in range(fam_cutoff):
       maxlength = max(seqlengths)
       seqlengths.remove(maxlength)
    consensuslength = maxlength

    
    
    currentcigars = []   # List of tuples (code,position). Initialize for every posision
    codeDict = {0:0,1:0,2:0,4:0}

    # Loop through the cigars for the current family and fetch the first cigar not an S
    for cig in cigar_list:
        n = cig.nextCigarPos()
        while n[0] != -1 and n[0] == 4:    # Loop until we get the first cigar that is not an S
            n = cig.nextCigarPos()
            
            # All cigars, also -1, are collected for the position
        currentcigars.append(n)
            
        # Keep control of the number of remaining cigars
    remainingcigs = len(currentcigars)-currentcigars.count((-1,-1))


        #Loop through the posisions of the sequences until we don't have enough remaning bases,
    #according to the fam_cutoff
    
    while remainingcigs >= fam_cutoff:  # remainingcigs is adjusted before making consensus

        # If we still have more than cutoff sequences with a base in this position,
        # check which cigar code is predominant and use this as a consensus
        currentcodes=[x[0] for x in currentcigars]  # The codes (M/I/D) of the cigar tups
        currentpositions = [x[1] for x in currentcigars]  #the number of bases in the tups
        Ms = currentcodes.count(0)  # Number of Matches in current cigar
        Is = currentcodes.count(1)
        Ds = currentcodes.count(2)
        Ss = currentcodes.count(4)
        determinators = []  # List of bases that is used to decide consensus base
        qualdeterminators = [] # List of base qualitites that will be used to decide consensus base qualitites

        #If M is the consensus cigar,
        #Pick the bases that will be used to decide this position in the consensus
        if Ms/remainingcigs > cutoff:
            # Scan the cigar codes and collect the bases and base qualities that correponds to all Ms
            for j  in range(len(currentcigars)):

                if currentcodes[j] == 0:  # an M
                    try:
                        determinators.append(sequence_list[j][currentpositions[j]])
                        qualdeterminators.append(quality_list[j][currentpositions[j]])
                    except:
                        pdb.set_trace()                            
                #If another sequence has an inserted base here, just skip that one
                # this function does not handle longer inserts than 2 in a nonconcensus readp
                elif currentcodes[j] == 1: #If we have an I
                    n  = cigar_list[j].nextCigarPos()
                    while n[0]==1:   #Scan the sequnce to the next base not an I, expect it to be an M. 
                        n = cigar_list[j].nextCigarPos()
                    if n[0] == 0:  #Assure that we retrieve an M before we add base/qual
                        determinators.append(sequence_list[j][n[1]])  
                        qualdeterminators.append(quality_list[j][n[1]])

                # Determine which base is the most frequent -> consensus
            for base in determinators:
                    
                try:
                    if base == 'T' :
                        nucIdentityList[0] += 1
                    elif base == 'C':
                        nucIdentityList[1] += 1
                    elif base == 'G':
                        nucIdentityList[2] += 1
                    elif base == 'A':
                        nucIdentityList[3] += 1
                    elif base == 'N':
                        nucIdentityList[4] += 1
                    else:
                        nucIdentityList[4] += 1
                    nucIdentityList[5] += 1     # Keep track of the number of bases in total

                except:
                    break
            consensusNoT.append(nucIdentityList[0])   #Store the number of As in this position of the SSCS family
            consensusNoC.append(nucIdentityList[1])
            consensusNoG.append(nucIdentityList[2])
            consensusNoA.append(nucIdentityList[3])
            consensusNoN.append(nucIdentityList[4])
            consensusNoBases.append(nucIdentityList[5])  #Store the total number of bases in this pos. of the fam
            try:
                for j in [0, 1, 2, 3, 4] :
                    
                    if float(nucIdentityList[j])/float(nucIdentityList[5]) > cutoff :
                        consensusRead += nucKeyDict[j]
                        consensusQual += chr(33+int(round(numpy.mean([int(qualdeterminators[i]) for i in range(len(qualdeterminators)) if determinators[i] == nucKeyDict[j]]))))
                        
                        break
                    elif j==4:
                        consensusRead += 'N'
                        consensusQual += chr(33)
                        #pdb.set_trace()
            except:
                consensusRead += 'N'
                consensusQual += chr(33)
            nucIdentityList=[0, 0, 0, 0, 0, 0] # Reset for the next nucleotide position

            # 1f I (insert) is the predominant cigar code, insert the most frequent I base
        elif Is/remainingcigs > cutoff:
            # Scan the cigar codes and collect the bases that correponds to all Is
            for j  in range(len(currentcigars)):
                if currentcodes[j] == 1:
                    determinators.append(sequence_list[j][currentpositions[j]])
                    qualdeterminators.append(quality_list[j][currentpositions[j]])
                elif currentcodes[j] == 0 or currentcodes[j] == 2:
                    # The cigars/sequences without an I must be reversed one step
                    # before next cycle
                    cigar_list[j].oneStepBack()

            # Determine which base is the most frequent -> consensus
            for base in determinators:
                
                try:
                    if base == 'T' :
                        nucIdentityList[0] += 1
                    elif base == 'C':
                        nucIdentityList[1] += 1
                    elif base == 'G':
                        nucIdentityList[2] += 1
                    elif base == 'A':
                        nucIdentityList[3] += 1
                    elif base == 'N':
                        nucIdentityList[4] += 1
                    else:
                        nucIdentityList[4] += 1
                    nucIdentityList[5] += 1
                  
                except:
                    break
                
            consensusNoT.append(nucIdentityList[0])   #Store the number of As in this position of the SSCS family
            consensusNoC.append(nucIdentityList[1])
            consensusNoG.append(nucIdentityList[2])
            consensusNoA.append(nucIdentityList[3])
            consensusNoN.append(nucIdentityList[4])
            consensusNoBases.append(nucIdentityList[5])  #Store the total number of bases in this pos. of the fam.
            try:
                for j in [0, 1, 2, 3, 4] :
                    if float(nucIdentityList[j])/float(nucIdentityList[5]) > cutoff :
                        consensusRead += nucKeyDict[j]
                        consensusQual += chr(33+int(round(numpy.mean([int(qualdeterminators[i]) for i in range(len(qualdeterminators)) if determinators[i] == nucKeyDict[j]]))))
                        break
                    elif j==4:
                        consensusRead += 'N'
                        consensusQual += chr(33)
            except:
                consensusRead += 'N'
                consensusQual += chr(33)
            nucIdentityList=[0, 0, 0, 0, 0, 0] # Reset for the next nucleotide position

            # 1f D (deletion) is the predominant cigar code, move all other cigars one
            # position ahead (skip that base)
        elif Ds/remainingcigs > cutoff:
            # Scan the cigar codes and collect the bases that correponds to all Ds
            for j  in range(len(currentcigars)):
                # If another sequence has an M, just skip that base
                #if currentcodes[j] == 0:
                #    cigar_list[j].nextCigarPos() Not neccessary
                # If another sequence has an I, skip both that base and the next
                if currentcodes[j] == 1:
                    cigar_list[j].nextCigarPos()

        elif Ss/remainingcigs > cutoff:  # If we just have softclipped sequences left, break the for loop
            break
                        
                
            # In case no cigar is more frequent than the cutoff
        else:
            consensusRead += 'N'
            consensusQual += chr(33)
            consensusNoA.append(0)   #If we don't have a trustable cigar, just append 0s, to make all lists synchronized
            consensusNoC.append(0)
            consensusNoG.append(0)
            consensusNoT.append(0)
            consensusNoN.append(0)
            consensusNoBases.append(0)

            for j  in range(len(currentcigars)):
                if currentcodes[j] == 1: #If we have an I, move counting one base out
                    cigar_list[j].nextCigarPos()
                                

        currentcigars = []   # List of tuples (code,position). Initialize for every posision
        codeDict = {0:0,1:0,2:0,4:0}
        for cig in cigar_list:
            n = cig.nextCigarPos()
            currentcigars.append(n)
        # Keep control of the number of remaining cigars
        remainingcigs = len(currentcigars)-currentcigars.count((-1,-1))
        
            
    return (consensusRead,consensusQual,consensusNoA,consensusNoC,consensusNoG,consensusNoT,consensusNoN,consensusNoBases)


#Function write_debuginfo writes information about every read in a SSCS family

def write_debuginfo(consensusdup):
    counter = 1
    for base in range(len(consensusdup[0])):
        debugfile.write(str(readDict[dictTag][1]+1)+":"+str(readDict[dictTag][2]+counter)+":" + str(readDict[dictTag][2])+":"+str(dictTag) +":")
        debugfile.write(str(consensusdup[0][base])+":")  #Output consensusbase
        # debugfile.write(str(consensusdup[1][base]) +":")  #Quality
        debugfile.write(str(consensusdup[2][base]) +":")  #Number of bases
        debugfile.write(str(consensusdup[3][base]) +":")  #
        debugfile.write(str(consensusdup[4][base]) +":")  #
        debugfile.write(str(consensusdup[5][base]) +":")  #
        debugfile.write(str(consensusdup[6][base]) +":")  #
        debugfile.write(str(consensusdup[7][base]) +"\n")  #
        counter = counter + 1
        



def tagStats(tagCountsFile):
    familySizeCounts=defaultdict( lambda: 0 )

    fIn = open(tagCountsFile, 'r')
    fOut = open("tagstats.txt", 'w')
    for line in fIn:
        familySizeCounts[int(line.strip().split()[1].split(":")[0])] += 1

    # Process readDict1 first, compare within    
    totals = 0
    for size in familySizeCounts.keys():
        familySizeCounts[size] *= int(size)
        totals += int(familySizeCounts[size])
    
    for size in sorted(familySizeCounts.keys()):
        fOut.write("%s\t%s\n" % (size, float(familySizeCounts[size])/float(totals)))
    
    fOut.close()
    return(True)




# CLEAN_READ_DICT FUNCTION
# This function compares the tags in the two input dictionaries, if they are adjacent og removes families that seems to be caused by sequencing errors in the tags

def cleanReadDict(readDict1,readDict2,falsefamfile):

# The function removes potential false SSCS only if the number of SSCS with the same start point is low enough to keep the expected number of random barcode collisisons to be less than 2%. That means, we can loose up to 2% of all SSCS because random barcode collisons.
    # When the number of SSCS starting in a position is low enough for cleaning, we removes SSCS with edit distance (number of differences) up to 2 if the SSCS have the same starting position. If they have different starting positions(neightbours), we accept up to 1 difference.
    #For example, the limit number of SSCS for clieaning if we consider only 1 difference between tags, is RandomBarcodeMatchLimitED1 = 22. That means if we have 22 SSCS starting in this position, the expected number of barcode collisions by random is 1,08, which is 5% of 22. That's why we use this limit for activation of cleaning. If the difference between SSCS families in famsize is > 10, we double the limits. That means accepts more reads before not using cleaning. 
# The reason why we accept higher expectation values if the size ratio is so high, is that we expect that random barcode matches will have family size of similar order of magnitude. If we observe extreme ratios, it's more likely due to sequencing/PCR errors in the barcode.

    # COMPARE THE TWO READDICTS AND RD2 WITH ITSELF

    falsefamcounter = 0 # Count the number of tags that are removed

    for run in [1,2]:  #The search for matches will be repeated with taglist1 and 2 interchanged
       
        matchlistlist = [] # To store the matching barcodes, list of lists of tuples(tag,1 or 2).
        famsizelistlist = [] # To store the famsizes of the matching SSCS, list of lists
        edlistlist = []   #To store the edit distances, list of lists
        maxratiolist = [] # To store the list of max ratios internally in each barcode group
        tags_to_remove1 = []  # The tags to remove for readD1
        tags_to_remove2 = [] # The tags to remove for readD2
        tags_to_remove1_sizes = [] # To store the sizes of the SSCS families corresponding to the tags that will be removed
        tags_to_remove2_sizes = []
        taglist1 = [tag for tag in readDict1.keys()]  #Get the actual tags (barcodes)
        taglist2 = [tag for tag in readDict2.keys()]
        max_template_no1 = len(taglist1) #Get the number of SSCSs with the current start position
        max_template_no2 = len(taglist2) # Get the  number off SSCS in the previous position
        
            
        startpos1=0
        startpos2=0
        if max_template_no1 > 0:   #In case the cleaning has wiped all remaining tags from dict
            startpos1 = readDict1[taglist1[0]][2]
        if max_template_no2 > 0:
            startpos2 = readDict2[taglist2[0]][2]

        


        #DO THE PAIRWISE COMPARISONS AND STORE BARCODE MATCHES

        for i in range(len(taglist2)):   #Go through the tags in readDict2
            matchlist = [(taglist2[i],2)] # To store the temporary matches in readD1 and readD2, begin with the comparator from ReadDict2. The number 2 is to indicate which ReadDict
            famsizelist = [(len(readDict2[taglist2[i]][7]))]  #Store the first famsize
            edlist = [] # To store the inner loops edit dists.

            
            # SCAN ALL TAGS IN READDICT1 FOR MATCHES, ONLY IF RD1 and RD2 are ADJACENT
            # IF two similar tags belongs to readDicts with different starting positions (neighbour), we only store tags with 0 or 1 difference (ED)
            if abs(startpos2-startpos1)==1 and max_template_no1 > 0:  # Check that the readDics are really neighbours and that RD is not empty
                for j in range(len(taglist1)):
                    ed = editdistance.eval(taglist2[i],taglist1[j])
                    if (ed <2): # If the number of differences is 0 or 1, store tags
                        #Store distance, tag and famsize for each match
             
                        matchlist.append((taglist1[j],1)) # store the tag from readDict 1 in case of matches, 1 to tell that it's from a neighbouring position.
                        famsizelist.append(len(readDict1[taglist1[j]][7])) #Store famsize of the match from readDict1
                        edlist.append(ed) # Store the number of edit differences between the tags that match
                        # End of inner for loop

            #SCAN ALL REMAINING TAGS IN READDICT2 FOR MATCHES
            for j in range(i+1,len(taglist2)):
                ed = editdistance.eval(taglist2[i],taglist2[j])
                if (ed<3):
                    matchlist.append((taglist2[j],2))
                    famsizelist.append(len(readDict2[taglist2[j]][7])) #Store famsize of the first
                    edlist.append(ed)
            # End of for loop scanning tags to compare with taglist2[i]


                                       
            # STORE MATCHES FOR THE CURRENT TAG IN readDict2
            if len(matchlist)>1:  #Only store matches, not single barcodes
                matchlistlist.append(matchlist)
                famsizelistlist.append(famsizelist)
                edlistlist.append(edlist)
                mx = max(famsizelist)
                mn = min(famsizelist)
                famsizeratio = mx/mn
                maxratiolist.append(famsizeratio)  # Needed to estimate how much difference in family sizes

        #End of outer for loop, that means, all tags in readDict2 and 1 have been compared and matches stored

        #SEARCH MATCHLISTLIST FOR TAGS TO DELETE AND STORE THEM
        
        for i in range(len(matchlistlist)):  #Loop through the matches/groups of SSCS in this position to decide barcode families to delete

            # Set default limits for false family removal, threshold give 2% erroneous removal of SSCS by random
            RandomBarcodeMatchLimitED0 = 160  # Estimated by simulation with simulate_barcodes.py
            RandomBarcodeMatchLimitED1 = 9
            RandomBarcodeMatchLimitED2 = 2 


            
            # Adjust limits for cleaning if the size differnce between families is high
            if (maxratiolist[i] > 10 ):  #If we have a very big SSCS family in a SSCS group, we accept higher SSCS numbers in a position before not removing potential false family
            
                RandomBarcodeMatchLimitED0 *= 2  
                RandomBarcodeMatchLimitED1 *= 2
                RandomBarcodeMatchLimitED2 *= 2 

                        
            #Adjust limits for cleaning even more, if we have more than 2 matching tags
            if len(matchlistlist[i])>2:  #If we have more than 2 matching barcodes, it's unlikely to be random. Increase limits. Multiplication is the right model for two independent events, if we have three matching barcodes. But as we get really many SSCS in a position, also a triple barcode collision is getting likely. 
            
                RandomBarcodeMatchLimitED0 *=  2
                RandomBarcodeMatchLimitED1 *=  2
                RandomBarcodeMatchLimitED2 *=  2

            #Determine whether the matchlist is from both rD1 and rD2 or only from rD2    
            meanReadDictNo = numpy.mean([matchlistlist[i][j][1] for j in range(len(matchlistlist[i]))])
            mx = max(famsizelistlist[i])  # Maximum SSCS family size for the match in question

            #DO THE ACTUAL SELECTION OF TAGS TO REMOVE
            
            if meanReadDictNo != 2:  # Reads from both read dictionarires
                if max(edlistlist[i]) == 0 and (max_template_no1 + max_template_no2) <= RandomBarcodeMatchLimitED0:
                     tags_to_remove1 = tags_to_remove1 + [matchlistlist[i][j][0] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==1]  #Removes all but the biggest SSCS from each match group
                     tags_to_remove1_sizes = tags_to_remove1_sizes + [famsizelistlist[i][j] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==1]  #Removes all but the biggest SSCS from each match group
         
                     tags_to_remove2 = tags_to_remove2 + [matchlistlist[i][j][0] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==2]
                     tags_to_remove2_sizes = tags_to_remove2_sizes + [famsizelistlist[i][j] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==2]  #Removes all but the biggest SSCS from each match group
                     
                elif max(edlistlist[i]) == 1 and (max_template_no1 + max_template_no2) <= RandomBarcodeMatchLimitED1:
                     tags_to_remove1 = tags_to_remove1 + [matchlistlist[i][j][0] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==1]
                     tags_to_remove1_sizes = tags_to_remove1_sizes + [famsizelistlist[i][j] for j in range(len(famsizelistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==1]  #Removes all but the biggest SSCS from each match group
                     tags_to_remove2 = tags_to_remove2 + [matchlistlist[i][j][0] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==2]
                     tags_to_remove2_sizes = tags_to_remove2_sizes + [famsizelistlist[i][j] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx and matchlistlist[i][j][1]==1]  #Removes all but the biggest SSCS from each match group

            
            elif meanReadDictNo ==2:   #Both from readD2, ED can then only be 1 or 2
                if (max(edlistlist[i]) == 1 and max_template_no2 <= RandomBarcodeMatchLimitED1): #Least stringent test if the distance in this within_rD match is 1                
                    tags_to_remove2 = tags_to_remove2 + [matchlistlist[i][j][0] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx]  # Find the tags for the smaller families, probably mistakes.
                    tags_to_remove2_sizes = tags_to_remove2_sizes + [famsizelistlist[i][j] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx]  # Find the tags for the smaller families, probably mistakes.
                elif max(edlistlist[i]) == 2 and max_template_no2 <= RandomBarcodeMatchLimitED2:
                    tags_to_remove2 = tags_to_remove2 + [matchlistlist[i][j][0] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx]
                    tags_to_remove2_sizes = tags_to_remove2_sizes + [famsizelistlist[i][j] for j in range(len(matchlistlist[i])) if famsizelistlist[i][j] < mx]  # Find the tags for the smaller families, probably mistakes.
                    

        #Remove unique tags from the readDict2. 2 first because this is the first position
        if len(tags_to_remove2)>0:
            unique_tags_to_remove2 = list(set(tags_to_remove2))
            falsefamcounter = falsefamcounter + len(unique_tags_to_remove2)
            falsefamfile.write(str(readDict2[tags_to_remove2[0]][1])  + "\t" + str(readDict2[tags_to_remove2[0]][2]) +"\t" + str(tags_to_remove2) + "\t" + str(tags_to_remove2_sizes) +  str("\n"))
            for tag2 in unique_tags_to_remove2:
                del(readDict2[tag2])   #Deletes the SSCSs that are false. By converting to set, we get only unique tags

        #REMOVE UNIQE TAGS FROM DICT 1
        if len(tags_to_remove1)>0:
            unique_tags_to_remove1 = list(set(tags_to_remove1))
            falsefamcounter = len(unique_tags_to_remove1)
            falsefamfile.write(str(readDict1[tags_to_remove1[0]][1])  + "\t" + str(readDict1[tags_to_remove1[0]][2]) + "\t" + str(tags_to_remove1) + "\t" + str(tags_to_remove1_sizes) + '\n')
            for tag1 in unique_tags_to_remove1:

                del(readDict1[tag1])   #Deletes the SSCSs that are false. By converting to set, we get only unique tags
                

    #We are repeating the run with the roles of the two read Dicts changed. The reason is that we would like to compare all tags with all, both ways, to find all kinds of matches
        tempDict = readDict2  # In the frist Run, 1 becomes 2 and vica versa. 
        readDict2 = readDict1 # In the second rund 1 becomes 1 again
        readDict1 = tempDict # Then its OK to return them
               
                
    return (readDict1,readDict2,falsefamcounter)
                    
                                 
                                

########################################################
# MAIN FUNCTION
########################################################


def main():
    ########################
    #Parameters to be input and variables to initatilize
    ########################
    overallStartTime = time.process_time()
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file. Stdin if -.", required=True)
    parser.add_argument("--tagprefix",  action="store",  dest="tagfile", help="Prefix for output tagcounts (.tags) and tag family size (.familysizes) files",  default='', required=False)
    parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file. Stdout: -.", required=True)
    parser.add_argument("--rep_filt", action="store",  type=int, dest='rep_filt', help="Remove tags with homomeric runs of nucleotides of length x. [9]", default=9 )
    parser.add_argument('--minmem', type=int, default=3, dest='minmem', help="Minimum number of reads allowed to comprise a consensus. [3] ")
    parser.add_argument('--maxmem', type=int, default=1000, dest='maxmem', help="Maximum number of reads allowed to comprise a consensus. [1000]")
    parser.add_argument('--cutoff', type=float, default=.7, dest='cutoff', help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [0.7]")
    parser.add_argument('--Ncutoff', type=float, default=1, dest='Ncutoff', help="With --filt 'n', maximum fraction of Ns allowed in a consensus [1.0]")
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut', help = 'How often you want to be told what the program is doing. [1000000]')
    parser.add_argument('--filt', type=str, default='n', dest='filt', help="A string indicating which filters should be implemented.  Filters: s: Softclipping 5' filter. n: N filter. 0: No filter applied ['n']")
    parser.add_argument('--reportfile', type=str, default='', dest='reportfile', help="A short report of results and call will be written to this file")
    parser.add_argument('--debugfile', type=str, default='', action="store",dest='debugfilename', help="Ouput a diagnostics file with SSCS sequences, qualities and frequencies of bases in each position")
    parser.add_argument('--remove_false', action="store", type=str, dest='remove_false_file', help="Remove SSCS families caused by sequencing errors in molecular barcode and output barcode to logfile")
    parser.add_argument('--start_pos_file', action="store",type=str,dest='start_pos_file', help="Count number of SSCS starting in every position of the genome covered by the sequencing and store to this file")
    parser.add_argument('--expand_SSCSs',action='store_true',default=False,help="Expand all SSCSs to the number of reads they originated from, all identical to the consensus sequence determined.")
    
    o = parser.parse_args()

    # Initialization of all global variables, main input/output files, and main iterator and dictionaries.
    goodFlag=[0,16]  # Which flags are accepted as good flags (only single end reads)

    # If the user want to read from stdin (Unix piping) instead, use -
    # Pysam handles that automatically
    inBam = pysam.AlignmentFile( o.infile, "rb" ) # Open the input BAM file

        #If the user want to write to stdout instead of a file

    if o.outfile == '-':
        outfile = sys.stdout
    else:
        outfile = open(o.outfile,'w')  # Output FASTQ file

    nonMap = pysam.AlignmentFile( o.outfile + "_filtered.bam", "wb", template = inBam )
    # File for reads with strange flags or bad tags
    outfile.Nfiltered = open(o.outfile + ".Nfiltered","w")
    # File for SSCS reads filtered because of too many Ns

    readNum = 0
    nM = 0   # Nonmapped files, also with bad tag
    bF = 0   # Number of reads with bad flag
    sC = 0   # No of soft clipped reads
    rT = 0   # Number of reads with repetetive tags
    nC = 0   # Number of consensuses with too many Ns
    
    LCC = 0
    ConMade = 0

    fileDone=False # Initialize end of file bool
    finished=False
    fam_size_list = [] #To store family sizes, family = reads with same genome position and tag
    readOne=True  # True when the current read is the first read in a new position, otherwise false


    bamEntry = inBam.fetch( until_eof = True ) # Initialize the iterator

    #Returns all reads in the whole file

    # The readWin is a list of two reads, odd numbered reads are being stored in the 1 position
    # Even numbered reads are being stored in the 0 position
    readWin = [next(bamEntry), '']  # Fetching the first readto 0 position
    winPos = 0    # Simple sequence
    falsefamcounter = 0 # To count SSCSs that are deleted because they seems to be genereated be tag sequencing error

    consensustimer = 0 # Count the time spent in the consensusmaker
    cleantimer = 0 # measure the time spent in the clean_read_dict
    readtimer = 0 # To measure the time spent reading BAM files into readDict
    writetimer = 0  # To measure the time spent on writing to FQ file
    
    readDict = {} # Initialize the main read dictionary
    newReadDict = {}  # To store the newly read read dictionary
    tagDict = defaultdict( lambda: 0 ) # Initialize the tag dictionary

    consensusDict={}

    positionlist = []  # Count the number of SSCS reads that starts in a certain position

    if o.debugfilename != '': #True if the debugfile option has been set
        debugfile = open(o.debugfilename,"w")   # Open the debug file for writing validation data
        debugfile.write("Chr:Pos:Start:Tag:Consensus:A:C:G:T:N:Total\n")

    if o.remove_false_file:  #Open file to log false barcodes (seq error in barcode)
        remove_false_file = open(o.remove_false_file,"w")
        remove_false_file.write("Chr\tPos\tRemovedTags\tFamsizes\n")
            

################################    
#Start going through the input BAM file, one position at a time, storing reads in families
################################
    for line in bamEntry:   # Fetching the next read
        winPos += 1
        readWin[winPos%2] = line #Place the next read in the right position in readWin

        if readOne==True:
            winPos -= 1

            
        # Loop in for loop until we find the first two sequences with the same start point in the alignment (reference_start)
        # Do while loop as long as we are reading sequences starting at the same alignment point. If it is the first read in a position, readOne 
          
        while (readWin[winPos%2].reference_start == readWin[(winPos-1)%2].reference_start and fileDone==False and readOne==False) or readOne == True:
            starttime = time.process_time()
            if readNum % o.rOut == 0 and readNum != 0:
                sys.stderr.write("[SSCScreator] Reads processed:" + str(readNum) + "\n")
            
            try:
                #read the tag 
                tag = readWin[winPos%2].query_name.split(':')[3]
                tagDict[tag] += 1


            except:
                print(readNum)
                raise
            

            readNum +=1

           
            # Softclip filter: filters out 5' softclipped reads (with --filt s)
            softClip=False
            if 's' in o.filt:
                if readWin[winPos%2].cigartuples != None:
                    if readWin[winPos%2].cigartuples[0][0]==4:  #Only consider 5' end softclipping now
#                    for tupple in readWin[winPos%2].cigar:
#                        if tupple[0]==4:
                        softClip=True
                        
        
            # Check if the given read is good data

            if int( readWin[winPos%2].flag ) in goodFlag  and softClip==False:
                
                 #Check if barcode has too long homopolymer stretches, according to rep_filt
                if ('A'*o.rep_filt in tag) or ('C'*o.rep_filt in tag) or ('G'*o.rep_filt in tag) or ('T'*o.rep_filt in tag): 
                    # Store reads with bad tag in nonMap file
                    nM += 1
                    nonMap.write(readWin[winPos%2])
                    rT += 1
                else :
                    
                    # Add the sequence to the read dictionary
                    # If the tag is new, store a new tag in newReadDict
                    # The read dictionary is a dictionary of lists, but the cigar string
                    # is stored as a list of Cigar objects
                    if tag not in newReadDict:
                        newReadDict[tag] = [readWin[winPos%2].flag, readWin[winPos%2].reference_name, readWin[winPos%2].reference_start, readWin[winPos%2].next_reference_id, readWin[winPos%2].next_reference_start, readWin[winPos%2].template_length,[Cigar(readWin[winPos%2].cigartuples)],[readWin[winPos%2].query_sequence],[readWin[winPos%2].query_qualities]] #Store data for a certain SSCS family, including qualitites (query_qualities)

                        #newReadDict is a dictionary of lists, one list for each SSCS.
                        #List structure:
                        # [0]: flag
                        # [1]: reference_id (kromosomnavn)
                        # [2]: reference_start (posisjon i kromosom)
                        # [3]: next_reference_id (neste sekvens sitt kromosomnavn)
                        # [4]: next_reference_start
                        # [5]: template_length
                        # [6]: cigars (as list of tuples)
                        # [7]: The sequences, as list of strings
                        # [8]: The base qualities, as list of strings
                        
                    else:  #If tag exist, just append new cigar and sequence
                        newReadDict[tag][6].append(Cigar(readWin[winPos%2].cigartuples)) # The Cigar strings as Cigar tuples
                        newReadDict[tag][7].append(readWin[winPos%2].query_sequence) # the sequence
                        newReadDict[tag][8].append(readWin[winPos%2].query_qualities) #Base qualitites

                        
                                         

                     
            else:
                nM += 1
                # Write unmapped reads (and softclipped) to nonMap file
                nonMap.write(readWin[winPos%2])
                if int(readWin[winPos%2].flag) not in goodFlag:
                    bF += 1
                elif softClip == True:
                    sC += 1
            
            winPos += 1

            if readOne == False:
                try: # Keep StopIteration error from happening at the end of a file
                    # The call to next belongs to the while loop (while at same genome pos)
                    readWin[winPos%2] = next(bamEntry) # Iterate the line
                except:
                    fileDone = True # Tell the program that it has reached the end of the file
            else:
                readOne = False

            stoptime = time.process_time()
            readtimer = readtimer + (stoptime - starttime)
        
        else:  # This else belongs to the while: Is being done only when the while clause fails
            # That means before we move to a new genome position 

            # Send reads to consensusMaker
            readOne=True     #To allow the while loop to start again in a new position

           
            #Remove too small readDictionaries first, to save computation time
            
            if newReadDict:
                taglist = [tag for tag in newReadDict.keys()]
                famsizelist = [len(newReadDict[t][7]) for t in taglist]
                for i in range(len(taglist)):
                    if famsizelist[i]<o.minmem:
                           del(newReadDict[taglist[i]])

                # If remove_false is defined and not the first position with SSCSs, remove false families from both read dictionaries
            if (o.remove_false_file)  :
                starttime = time.process_time()
                #Clean the two dictionaries against each other                   
                (newReadDict,readDict,falsefams) = cleanReadDict(newReadDict,readDict,remove_false_file)
                stoptime = time.process_time()
                cleantimer = cleantimer + (stoptime-starttime)
                falsefamcounter = falsefamcounter + falsefams


                
            ##############################
            # Send sequences in the previous readDict to the consensus maker
            ##############################

            if readDict:    #Check that we have a dictionary to analyze

                taglist = [tag for tag in readDict]
                chrom = readDict[taglist[0]][1]
                position = readDict[taglist[0]][2]
                number = len(readDict)  # Store the size of each Dict
                positionlist.append((chrom,position,number)) #Store number of sequences starting in each position
                for dictTag in readDict.keys():
                    starttime = time.process_time()  # measure computation time
                    fam_size = len(readDict[dictTag][6])
                     
                    if fam_size <= o.maxmem:   # If the number of reads in a family is below mx
                        consensusdup = consensusMaker(readDict[dictTag][6],readDict[dictTag][7],readDict[dictTag][8],o.minmem,o.cutoff)
                        if o.debugfilename != '':
                            write_debuginfo(consensusdup)
                        
                    else:
                        rsample = random.sample(range(fam_size),o.maxmem)
                        consensusdup = consensusMaker([readDict[dictTag][6][i] for i in rsample],[readDict[dictTag][7][i] for i in rsample],[readDict[dictTag][8][i] for i in rsample],o.minmem,o.cutoff)
                        
                    stoptime = time.process_time()
                    consensustimer = consensustimer + (stoptime-starttime)
                        
                                
                    # PRINT SSCS
                    # The first and optional SSCS writing routine expands the SSCS to a number of reads equal to the number of family members in base for the SSCS (Yves idea)
                    if o.expand_SSCSs == True:
                        if (consensusdup[0].count("N" )/ len(consensusdup[0]) <= o.Ncutoff and 'n' in o.filt) or ('n' not in o.filt):
                            starttime = time.process_time()
                            for nr in range(fam_size):
                                seqname = str(readDict[dictTag][2]) + ":" + dictTag + ":" + str(fam_size) + '-' + str(nr+1) # Sekvensnavnet består av posisjon, tag, running number
                                fam_size_list.append(fam_size)
                                outfile.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))  # Writing one read for every founding member of the SSCS
                            ConMade += 1    #Count written SSCS
                            stoptime = time.process_time()
                            writetimer = writetimer + (stoptime - starttime)

                        else:
                            seqname = str(readDict[dictTag][2]) + ":" + dictTag + ":" + str(fam_size) # Sekvensnavnet består av posisjon, tag, familiestørrelse
                            nC += 1
                            outfile.Nfiltered.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))
                        
                    else:
                    #Write single consensus sequences to FASTQ file
                        seqname = str(readDict[dictTag][2]) + ":" + dictTag + ":" + str(fam_size) # Sekvensnavnet består av posisjon, tag, familiestørrelse
                        # Filter out consensuses with too many Ns in them (the if keeps only the others for the loop)
                        if (consensusdup[0].count("N" )/ len(consensusdup[0]) <= o.Ncutoff and 'n' in o.filt) or ('n' not in o.filt):
                            fam_size_list.append(fam_size)
                            starttime = time.process_time()
                            outfile.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))
                            ConMade += 1    #Count written SSCS
                            stoptime = time.process_time()
                            writetimer = writetimer + (stoptime - starttime)
                        else:
                            nC += 1
                            outfile.Nfiltered.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))
        
                        
             
        readDict=newReadDict
        newReadDict = {} # Reset the read dictionary
    # For loop ends, recursing through genome positions

    ############################################################
    # Send sequences in the last readDict to the consensus maker
    ############################################################

    
    #Make consensuses and write for the last readDict
    if readDict:    #Check that we have a dictionary to analyze
        for dictTag in readDict.keys(): 
            fam_size = len(readDict[dictTag][6])
                
            if fam_size >= o.minmem:  #Only process a family if it is big enough

                if fam_size <= o.maxmem:   # If the number of reads in a family is below mx
                    consensusdup = consensusMaker(readDict[dictTag][6],readDict[dictTag][7],readDict[dictTag][8],o.minmem,o.cutoff) #Making the SSCS
                    if o.debugfilename != '':
                        write_debuginfo(consensusdup)
                        
                else:
                    rsample = random.sample(range(fam_size),o.maxmem)
                    consensusdup = consensusMaker([readDict[dictTag][6][i] for i in rsample],[readDict[dictTag][7][i] for i in rsample],[readDict[dictTag][8][i] for i in rsample],o.minmem,o.cutoff)


                # Write expanded SSCS from the last readDict to file, optional
                if o.expand_SSCSs == True:
                    if (consensusdup[0].count("N" )/ len(consensusdup[0]) <= o.Ncutoff and 'n' in o.filt) or ('n' not in o.filt):
                        starttime = time.process_time()
                        for nr in range(fam_size):
                            seqname = str(readDict[dictTag][2]) + ":" + dictTag + str(fam_size) + '-' + str(nr+1) # Sekvensnavnet består av posisjon, tag, running number
                            fam_size_list.append(fam_size)
                            outfile.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))  # Writing one read for every founding member of the SSCS
                        ConMade += 1    #Count written SSCS
                        stoptime = time.process_time()
                        writetimer = writetimer + (stoptime - starttime)

                    else:  # If we have too many Ns
                        seqname = str(readDict[dictTag][2]) + ":" + dictTag + ":" + str(fam_size) # Sekvensnavnet består av posisjon, tag, familiestørrelse
                        nC += 1
                        outfile.Nfiltered.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))
                else:
                    #Write single consensus sequences to FASTQ file
                
                    seqname = str(readDict[dictTag][2]) + ":" + dictTag + ":" + str(fam_size) # Sekvensnavnet består av posisjon, tag, familiestørrelse
                    if (consensusdup[0].count("N" )/ len(consensusdup[0]) <= o.Ncutoff and 'n' in o.filt) or ('n' not in o.filt):
                        fam_size_list.append(fam_size)
                        starttime = time.process_time()
                        outfile.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))
                        ConMade += 1    #Count written SSCS
                        stoptime = time.process_time()
                        writetimer = writetimer + (stoptime - starttime)
                    else:
                        nC += 1
                        outfile.Nfiltered.write(">%s\n%s\n+\n%s\n"  % (seqname, consensusdup[0],consensusdup[1]))
            
                                
        
        
# Close BAM files

    if o.debugfilename != '':
        debugfile.close()

    if o.remove_false_file:
        remove_false_file.close()
        

    inBam.close()
    outfile.close()
    outfile.Nfiltered.close()
    nonMap.close()

    overallStopTime = time.process_time()  # Final timepoint

    reportfile = sys.stderr
    
    if o.reportfile != '':
        reportfile = open(o.reportfile,"w")
    
    # Write summary statistics
    reportfile.write("REPORT \n")
    reportfile.write("********\n")
    reportfile.write("Call: " + ' '.join(sys.argv) + "\n")
    reportfile.write("Reads processed:" + str(readNum) + "\n")
    reportfile.write("Bad reads: %s\n" % nM)
    reportfile.write("\tReads with Bad Flags: %s\n" % bF)
    reportfile.write("\t5' softclipped Reads: %s\n" %sC)
    reportfile.write("\tRepetitive Tag: %s\n" % rT)
    reportfile.write("Consensuses Made: %s\n" % ConMade)
    reportfile.write("Consensuses with Too Many Ns: %s\n" % nC)
    reportfile.write("Mean number of reads per SSCS family: %f\n" % (sum(fam_size_list)/len(fam_size_list)))
    reportfile.write("Max number of reads per SSCS family: %f\n" % max(fam_size_list))
    reportfile.write("Number of false potential SSCS families deleted: " + str(falsefamcounter) + "\n")
    reportfile.write("\nOverall processing time (s): " + str(overallStopTime-overallStartTime) + "\n")
    reportfile.write("Time spent on reading BAM file (s): "+ str(readtimer) + "\n")
    reportfile.write("Time spent on SSCScreator function (s): "+ str(consensustimer) + "\n")
    reportfile.write("Time spent on cleaning read dictionaries (s): "+ str(cleantimer) + "\n")
    reportfile.write("Time spent on writing to Fastq file (s): "+ str(writetimer) + "\n")
    if o.reportfile != '':
        reportfile.close()
    

    # Write the tag counts and the tag family size files.
    if o.tagfile != '':
        tagFile = open( o.tagfile+".tags", "w" )
        tagFile.write ( "\n".join( [ "%s\t%d" % ( SMI, tagDict[SMI] ) for SMI in sorted( tagDict.keys(), key=lambda x: tagDict[x], reverse=True ) ] ))
        tagFile.close()
        tagStats(o.tagfile+".tags")

        famsizefile = open(o.tagfile+"-famsizes.csv","w")
        famsizefile.write("SSCS read family sizes\n")
        for i in fam_size_list:
            famsizefile.write("%s\n" % i)
        famsizefile.close()

        #write position count file
    if o.start_pos_file:
        posfile = open(o.start_pos_file,"w")
        posfile.write("Chr\tPos\tCount\n")
        for tup in positionlist:
            posfile.write(str(tup[0]) +"\t" + str(tup[1]) + "\t" + str(tup[2]) + "\n")
        posfile.close()

if __name__ == "__main__":
    main()
