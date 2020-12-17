#NOTES ON RUNNING MiXCR
###############################################################

#"MiXCR is a universal framework that processes big immunome data from raw sequences to quantitated clonotypes. MiXCR efficiently handles paired- and single-end
#reads, considers sequence quality, corrects PCR errors and identifies germline hypermutations. The software supports both partial- and full-length profiling 
#and employs all available RNA or DNA information, including sequences upstream of V and downstream of J gene segments."

###############################################################
#GENERALIZED STEPS

#1. Align - align reads to reference genes
#2. Assemble - assemble clonotypes using alignments generated in #1 step
#3. Export -  either alignments or clonotype library

#RECOMMENDED FOR RNA-SEQ:
  #include "assemblePartial" and "extend" (optional otherwise)
  #recommend to use "analyze" command: packs a complicated execution pipeline (alignment, assembly, exporting) into a single command. See more below.

#Main challenges for RNA-Seq: 
  #1)extraction and alignment of target molecule fragments. Alingnment parameters (-p rna-seq) was specifically optimized for false positive rate.
  #2)assembly of overlapping fragmented sequencing reads into long-enough CDR3 containing contigs. Special action to perform such an assembly of reads, 
      #partially covering CDR3 is "assemblePartial." This Performs an overlap of already aligned reads from *.vdjca file, realigns resulting contig, 
      #and checks if initial overlap has covered enough part of a non-template N region.

#All default parameters should be suitable for RNA-Seq

###############################################################
#BREAKDOWN OF STEPS / TYPICAL ANALYSIS WORKFLOW

#
- align step: 
  - Species: Using -s .. paramater. hsa = HomoSapiens, can use either
  - Data source origin: Genomic or transcriptomic. Affects which part of the 
  references V gene sequences will be used for alignment. transcriptomic source is 
  assumed by default. 



# 1. align sequencing reads against reference V, D, J and C genes.
     
     mixcr align -p rna-seq -s hsa -OallowPartialAlignments=true data_R1.fastq.gz data_R2.fastq.gz alignments.vdjca
                                                                  #[input file 1]  [input file 2]  [output]
     
     #NEED to specify -s parameter : indicates Species. "-s hsa" (or "-s HomoSapiens")
     #NEED to specify -p parameter : indicated data source origin, genomic or transcriptomic. Affects which part of the references V gene sequences will be used
      #for alignment. Transcriptomic source is assumed by default. 
     #"OallowPartialAlignments=true" : is default, successfully aligned reads - number of reads with at least one of V or J alignments, that passed all alignment
        #thresholds and cover at least one nucleotide of CDR3. Needed to prevent MiXCR from filtering out partial alignments, that don’t fully cover CDR3 


# 2. perform two rounds of contig assembly

#a)
      mixcr assemblePartial alignments.vdjca alignments_rescued_1.vdjca
                          #[alignment output from #1] [output]
#b)
      mixcr assemblePartial alignments_rescued_1.vdjca alignments_rescued_2.vdjca
                          #[alignment output from #2a] [output]

# 3. (optional) Perform extension of incomplete TCR CDR3s with uniquely determined V and J genes using germline sequences
    #RECOMMENDED FOR RNA-SEQ

      mixcr extend alignments_rescued_2.vdjca alignments_rescued_2_extended.vdjca
                           #[output from #2b] [output]

# 4. assebmle clonotypes

      mixcr assemble alignments_rescued_2_extended.vdjca clones.clns
                           #[output from #3 or #2b] [output]      

# 5. export all clonotypes

      mixcr exportClones clones.clns clones.txt
                           #[output from #4] [output]

#or clonotypes for a specific immunological chain:

#mixcr exportClones -c TRB clones.clns clones.TRB.txt
#mixcr exportClones -c IGH clones.clns clones.IGH.txt
#...

#The resulting *.txt files will contain clonotypes along with comprehansive biological information like V, D, J and C genes, clone abundances, etc…

###############################################################
#USE OF ANALYZE COMMAND (RECOMMENDED FOR RNA-SEQ), OPTIONAL

#ANALYZE command
  # - two types of library prep : 1) analyze amplicon vs. 2) analyze shotgun : random fragments (RNA-Seq)
  # - shotgun: implements pipeline for the analysis of non-enriched RNA-seq. Includes alignment of raw sequencin using "align", assembly of overlapping 
      #fragmented reads using "assemblePartial", imputing good TCR alignments using "extend", assembly of aligned sequences into clonotypes using "assemble",
      #exporting the resulting clonotypes into tab-delimited file w/ "export"
  # - optional: can also assembles full receptor sequences with "assembleContigs"

#USAGE: 

mixcr analyze shotgun

#  -s <species> \  : must define.
#  --starting-material <startingMaterial> \   : must define. type of starting material, ("rna").
#  [OPTIONS] input_file1 [input_file2] analysis_name

#The pipeline above is equivalent to execution of the following MiXCR actions:

# align raw reads
mixcr align -s <species> -p <aligner> \
    -OvParameters.geneFeatureToAlign=<vFeatureToAlign> \
    -OvParameters.parameters.floatingLeftBound=false \
    -OvParameters.parameters.floatingRightBound=false \
    -OvParameters.parameters.floatingRightBound=false \
    [align options] input_R1.fastq [input_R2.fastq] my_analysis.vdjca

# assemble overlapping fragmented sequencing reads
mixcr assemblePartial [assemblePartial options] my_analysis.vdjca my_analysis.rescued_1.clna
mixcr assemblePartial [assemblePartial options] my_analysis.rescued_1.vdjca my_analysis.rescued_2.clna

# impute germline sequences for good TCR alignments
mixcr extend [extend options] my_analysis.rescued_2.vdjca my_analysis.rescued_2.extended.vdjca

# assemble CDR3 clonotypes
mixcr assemble --write-alignments [assemble options] my_analysis.rescued_2.extended.vdjca my_analysis.clna

# assemble contigs: execute only if --assembleContigs is specified
mixcr assembleContigs [assembleContigs options] my_analysis.clna my_analysis.clns

# export to tsv
mixcr exportClones [export options] my_analysis.clns my_analysis.txt

#As in the case of analyze amplicon, required option --starting-material affects the choice of V gene region which will be used as 
#target in align step (vParameters.geneFeatureToAlign, see align documentation): rna corresponds to the VTranscriptWithout5UTRWithP and dna to 
#VGeneWithP (see Gene features and anchor points for details).
