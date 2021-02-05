# SSCScreator
Software to collapse single strand consensus reads from molecular tags in sequencing data

Inputs: 
    A position-sorted single-end BAM file containing reads with a molecular tag 
in the header. All reads must be oriented in the same direction as the reference
 genome. The reverse reads must be aligned to a reverse copy of the reference ge
nome. Can also read from stdin, if - is given as filename.

Outputs:
    1: A single-end FASTQ file containing SSCSs
    2: A single-end BAM file containing reads with bad tag, unmapped (or softcli
pped if s filter)
    3: A tagcounts file
  
The program starts at the position of the first good read, determined by the type of read specified on startup.  It then goes through the file until it finds a new position, saving all reads as it goes, in families according to molecular tag.   When it finds a new position, it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are written to the output FASTQ file.  After making consensuses with the reads from the first position, it continues on through the origional file until it finds another new position, sends those reads to the consensus maker, and so on until the end of the file. The SSCScreator function is only called if there is more than minmem reads in a family. It also handles the read lenght and stops calling consensus bases when there is less than minmem sequences long enough.This version is made only or single end reads. In contrast to the ConsensusMaker of Kennedy, this program collects information for all reads in a family, not only those with the most common cigar string. Both a specific cigar code and a base has to be present in more than cutoff fraction of the reads to be included in the consensus.

Outputs both a FASTQ file and a _filtereout.bam BAM file with all filtered reads
. Mainly because they have not aligned. But also softclipped, if we apply the so
ftclipping filter. In addition, if we use the --tagprefix variable, we output so
me tag statics files. The sequence name of the reads in the FASTQ file consist o
f the position, the tag and the family size

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
    --Ncutoff, with --filt n enabled, sets the maximum percentage of Ns allowed 
in a SSCS.  
        Example (--Ncutoff = .1, --readlength = 20):
            Two SSCSs are generated as follows:
                SSCS 1: ACGTGANCTAGTNCTNTACC
                SSCS 2: GATCTAGTNCATGACCGATA
SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, while SSCS 1 does not with 3/20 = 15% Ns.
    --filt sets which filters are used.  Options are:           
                s: Softclipping filter.  Filters out any reads which have been soft-clipped in alignment.  This avoids later problems with hard-clipping.  
                n: N filter. Filters out consensus sequences with a higher percentage of Ns than the threshold imposed by --Ncutoff.  Without this option, --Ncutoff doesn't do anything.  


