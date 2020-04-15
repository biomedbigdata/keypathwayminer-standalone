/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.main;

import java.io.IOException;

import de.mpg.mpiinf.ag1.kpm.KPMRunHandler;
import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.Program;
import de.mpg.mpiinf.ag1.kpm.parsers.ArgsParametersParser;
import de.mpg.mpiinf.ag1.kpm.parsers.InputFileParser;
import de.mpg.mpiinf.ag1.kpm.parsers.PropertiesFileParser;
import de.mpg.mpiinf.ag1.kpm.shortestpath.ShortestPathAlgorithms;
import dk.sdu.kpm.KPMSettings;

/**
 *
 * @author nalcaraz
 *
 */

public class Main {

    private volatile KPMSettings kpmSettings;
    
	public static void main(String[] args) throws Exception {
		Main main = new Main();
		main.run(args);
	}

	public void run(String[] args){
        int mb = 1024*1024;
        Runtime runtime = Runtime.getRuntime();
        System.out.println("Used Memory:" + (runtime.freeMemory()) / mb);

		kpmSettings = new KPMSettings();
		
		
		// parse parameters from the kpm.properties file
		//PropertiesFileParser pfp = new PropertiesFileParser(kpmSettings);
		//Parameters params = pfp.parse("src/main/resources/kpm.properties");

		// Setting folder prefix for testing:
		String prefix = "resources/";
		Parameters params = new Parameters();
		params.DATASETS_FILE = prefix + params.DATASETS_FILE;
		params.POSITIVE_FILE = prefix + params.POSITIVE_FILE;
		params.NEGATIVE_FILE = prefix + params.NEGATIVE_FILE;
		params.RESULTS_FOLDER = prefix + params.RESULTS_FOLDER;
		params.VALIDATION_FILE = prefix + params.VALIDATION_FILE;
		params.GRAPH_FILE = prefix+ params.GRAPH_FILE;

		try {
			start(args, prefix + "kpm.properties", params);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	

	public void start(String[] args, String propertiesFile, Parameters params) throws Exception {

		ArgsParametersParser argspp = new ArgsParametersParser(kpmSettings);
		params = argspp.parse(args, params);

		// check if all input files are correct
		InputFileParser ifp = new InputFileParser(kpmSettings);
		params = ifp.parse(params);

		if (params.PROGRAM == Program.SP) {
			System.out.println("PROGRAM SELECTED: Shortest Paths");
			ShortestPathAlgorithms.shortestPathways(params.GRAPH_FILE, params);
		} else {
			KPMRunHandler kpmHandler = new KPMRunHandler(kpmSettings);
			kpmHandler.runBatch(params);
		}
	}
}
