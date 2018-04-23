package de.mpg.mpiinf.ag1.kpm;

import java.util.HashMap;

import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.graph.KPMGraph;

public class KPMParameterConfigurator {
	private volatile KPMGraph graph;
	
	private volatile KPMAlgorithm algorithm;
	
	private volatile KPMStrategy strategy;
	
	private String prefix = "Configuration is incorrect:\t";

    private volatile KPMSettings kpmSettings;
    
	public KPMParameterConfigurator(Parameters params, KPMSettings settings){
		this.kpmSettings = settings;
		this.graph = kpmSettings.MAIN_GRAPH;
		this.algorithm = params.ALGORITHM;
		this.strategy = params.STRATEGY;
	}
	
	public boolean areSettingsValid(){
		kpmSettings.STARTING_TIME = System.nanoTime();
		
		if(graph == null){
			System.out.println(prefix + "No graph loaded.");
			return false;
		}
		
		int maxNumberOfGeneExceptions = graph.getNodeIdSet().size();
		
		if (kpmSettings.IS_BATCH_RUN) {
			if(strategy == KPMStrategy.INES){
				if(kpmSettings.MAX_K <= 0){
					System.out.println(prefix + "Invalid MAX_K.");
					return false;
				}
				
				if(kpmSettings.MIN_K > kpmSettings.MAX_K){
					System.out.println(prefix + "Invalid MIN_K. It is larger than MAX_K.");
					return false;
				}
				
				if(kpmSettings.INC_K == 0){
					System.out.println(prefix + "Invalid INC_K. Must be larger than 0.");
					return false;
				}

			}else{
				kpmSettings.MIN_K = kpmSettings.MAX_K = kpmSettings.INC_K = 0;
				kpmSettings.GENE_EXCEPTIONS = 0;
			}
			
			// Checking L's

			//sanityCheckAndSetMinCaseExceptionsMap();
			//sanityCheckAndSetStepSizeCaseExceptionsMap();
			
		}else{
			if(strategy == KPMStrategy.INES){
				if(maxNumberOfGeneExceptions < kpmSettings.GENE_EXCEPTIONS){
					System.out.println(prefix + "Invalid number of GENE_EXCEPTIONS. Must be less than the number of nodes.");
					return false;
				}
			}
		}

		return true;
	}
	
	private void sanityCheckAndSetMaxCaseExceptionsMap(){
		if(strategy != KPMStrategy.INES){
			return;
		}
	}
	
	private void sanityCheckAndSetMinCaseExceptionsMap() {
		HashMap<String, Integer> checkedMinCaseExceptionsMap = new HashMap<String, Integer>();
		for (String dataSetID : kpmSettings.MIN_L.keySet()) {
			System.out.println("Sanity = " + dataSetID);
				int lMinTemp = kpmSettings.MIN_L.get(dataSetID);
				int maxLMin = kpmSettings.MAX_L.get(dataSetID) - 1;
				if (lMinTemp > maxLMin) {
					checkedMinCaseExceptionsMap.put(dataSetID, maxLMin);
					System.out.println("Min L for data set "
							+ dataSetID
							+ " was set to " + maxLMin);
				} else {
					checkedMinCaseExceptionsMap.put(dataSetID, lMinTemp);
				}
		}
		kpmSettings.MIN_L = checkedMinCaseExceptionsMap;
	}
	
	private void sanityCheckAndSetStepSizeCaseExceptionsMap() {
        HashMap<String, Integer> checkedStepSizeCaseExceptionsMap = new HashMap<String, Integer>();
		for (String dataSetID : kpmSettings.INC_L .keySet()) {
				int lStepTemp = kpmSettings.INC_L .get(dataSetID);
				int maxStepSize = kpmSettings.MAX_L.get(dataSetID)
						- kpmSettings.MIN_L.get(dataSetID);
				if (lStepTemp > maxStepSize) {
					checkedStepSizeCaseExceptionsMap
					.put(dataSetID, maxStepSize);
					System.out.println("Step L for data set "
							+ dataSetID
							+ " was set to " + maxStepSize);
				} else {
					checkedStepSizeCaseExceptionsMap.put(dataSetID, lStepTemp);
				}
		}
		kpmSettings.INC_L = checkedStepSizeCaseExceptionsMap;
	}
}
