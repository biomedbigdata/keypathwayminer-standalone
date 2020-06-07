
KeyPathwayMiner 5.0
============
Given a biological network and a set of case-control studies, KeyPathwayMiner efficiently extracts all maximal connected sub-networks. These sub-networks contain the genes that are mainly dysregulated, e.g., differentially expressed, in most cases studied.

For more informaiton please visit our website [(Key Pathway Miner website)](https://keypathwayminer.compbio.sdu.dk/keypathwayminer/).

## ABOUT
    KeyPathwayMiner version 5.0
    Copyright 2016 by
    Nicolas Alcaraz: nalcaraz@mpi-inf.mpg.de and
    Jan Baumbach: jan.baumbach@imada.sdu.dk

## USAGE
   General structure:
      
      java -jar [jvm options] KPM-5.0.jar [-KEY1=VAL1] .... [-KEYN=VALN] 

   Execution examples:

      java -jar -Xmx2G KPM-5.0.jar -strategy=INES -algo=GREEDY -K=2 -L1=5
      java -jar KPM-5.0.jar -strategy=GLONE -algo=GREEDY -L1=5

      java -jar KPM-5.0.jar -strategy=GLONE -algo=GREEDY -L1=5 -matrix1=resources/datasets/colon-gene-expression-DOWN-p0.05.txt -L2=6 -matrix2=resources/datasets/colon-gene-expression-UP-p0.05.txt

      java -jar KPM-5.0.jar -strategy=GLONE -algo=ACO -datasetsFile=resources/datasets_file.txt

   Note: If the input is large and/or complex then the java virtual machine options must be set.

### PARAMETERS

    All parameters can be defined in the "kpm.properties" file
    and most through command line arguments. In case the same parameter
    is defined both in the properties file and through a command
    line argument, then the command line argument will have preference
    and override the value in the properties file. 


#### INPUT FILES

    KPM takes as input several files

    -- GRAPH file:  The file containing all the interactions
    of the protein interaction network.  This must be in sif format:

    NODE1   INTERACTION_TYPE   NODE2
    112 pp 342
    12  pp 342
    ...
    ...

    -- MATRIX file(s):  The file containing and indicator matrix
    for expression studies. KPM can take as input several of these
    files, which can be defined either through the command line 
    or in the DATASETS file (see below). 

    GENE_ID   CASE1	  CASE2   CASE3
    10203         1	      0       1
    3232          0           0       1
    ...
    ...		

    -- DATASETS file:  This file contains the paths to each individual
    indicator matrix file and it's corresponding L parameter. The format
    should be the following:

    ID      L       PATH
    1       10      path/to/matrix1.txt
    2       15      path/to/matrix2.txt
    ...
    ...

    -- POSITIVE/NEGATIVE list file: these are optional files that contain
    a list of genes that will be given high (POSITIVE list) or low (NEGATIVE list)
    priority when searching for pathways.  The format is:

    NODE_ID
    11423
    1213
    ...
    ...

#### BASIC OPTIONS

      -K {'integer'}                
      Gene exceptions (only used for INES)
  
      -graphFile {graph.sif}   
      The path to the graph file

      -matrixN
      The path to matrix N 
      
      -LN {'integer'} The case (L) exceptions for matrix N

      -datasetsFile  
      The matrix files and their L parameters can also 
      be specified in the datasets file. 

      -positiveFile  
      The path to the file with the positive gene list
  
      -negativeFile  
      The path to the file with the negative gene list

      -strategy {'INES','GLONE'}  
      The strategy that will be used to extract pathways:

      INES: Will extract maximal pathways where all but at most K nodes
            are NOT active/diff. expressed in at most L cases.

      GLONE: Will extract maximal pathways where the total sum of 
             NOT-active/diff. exp. cases is at most L.
              
      -algo {'ACO','GREEDY','OPTIMAL'} 
      The algorithm that will be used to extract the pathways:

      GREEDY: A simple greedy algorithm, performs fast

      ACO: Ant-colony-optimization algorithm. Convergence times
               can vary depending on input size and parameters 
               (see advanced options section).

      OPTIMAL: An exact fixed-parameter-tracktability algorithm. 
               Running times increase exponentially with input size(GLONE)
               and parameter K (INES). Note: due to large running times,
               this will only extract the best pathway. 

	  
#### ADVANCED OPTIONS
      -numProc  {'integer'}
      Number of threads to use (for parallel computing)

      -combineOp {OR,AND,CUSTOM}
      How to combine multiple matrices using boolean operators.
      If CUSTOM is chosen then the logical predicate
      defined in the kpm.properties file will be used. 

      -maxSolutions {'integer'}
       Maximum number of reported pathways 
      (default is 20). Ignored if
      OPTIMAL algorithm is selected.

##### ADVANCED OPTIONS (ACO)
      -alpha {'double'}
      Parameter to control the importance given 
      to the pheromone.

      -beta {'double'}
      Parameter to control the importance given to the
      heuristic value of the node

      -rho {'double'}
      Parameter that controls the pheromone decay rate

      -iterations {'integer'} 
      Maximum number of iterations

      -maxrunswithoutchange {'integer'}
      Maximum number of iterations allowed without
      improvement in the best solution.

      -iterationbased {'true','false'}
      If iteration or global best ACO strategy should be used 
      (only for GLONE). 
             

 