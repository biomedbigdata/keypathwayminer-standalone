#!/bin/bash
      #this is just a small script wit same calls for testing. It does not include alle edge cases!
      java -jar -Xmx2G KPM-5.0.jar -strategy=INES -algo=GREEDY -K=2 -L1=5 -matrix1=resources/indicator_matrix_hd.txt
      java -jar KPM-5.0.jar -strategy=GLONE -algo=GREEDY -L1=5 -matrix1=resources/indicator_matrix_hd.txt
      java -jar KPM-5.0.jar -strategy=GLONE -algo=GREEDY -L1=5 -matrix1=resources/datasets/colon-gene-expression-DOWN-p0.05.txt -L2=6 -matrix2=resources/datasets/colon-gene-expression-UP-p0.05.txt
      java -jar KPM-5.0.jar -strategy=GLONE -algo=ACO -datasetsFile=resources/datasets_file.txt
	java -jar KPM-5.0.jar -strategy=INES -algo=OPTIMAL -datasetsFile=resources/datasets_file.txt -K=3

	  