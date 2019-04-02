package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import au.com.bytecode.opencsv.CSVReader;
import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.main.Main;
import dk.sdu.kpm.KPMSettings;

public class DatasetsFileParser {

	private volatile KPMSettings kpmSettings;

	public DatasetsFileParser(KPMSettings settings){
		this.kpmSettings = settings;
	}

	public Parameters parse(char separator, boolean hasHeader, Parameters params) {
		try {

            HashMap<String, String> id2path = new HashMap<String, String>();
            HashMap<String, Integer> id2param = new HashMap<String, Integer>();
			CSVReader cr =
					new CSVReader(new FileReader(params.DATASETS_FILE), separator);
			String[] nextLine = null;


			if (hasHeader) {
				nextLine=cr.readNext();
			}

			int row =1;
			nextLine = cr.readNext();
			while (nextLine != null) {
				String id;
				String filePath;
				int lParam;
				if (nextLine.length < 2) {
					break;
				}
				if (nextLine.length == 2) {
					id = String.valueOf(row);
					lParam = Integer.parseInt(nextLine[0]);
					filePath = nextLine[1];
				} else {
					id = nextLine[0];
					lParam = Integer.parseInt(nextLine[1]);
					filePath = nextLine[2];
				}
				if (!new File(filePath).isFile()) {
					System.out.println("Matrix file " + filePath + " does not exist !");
					System.exit(-1);
				}
				// TODO: Delete this code. It is evil. Min max must be set properly. Question is whether this datasetfile parser makes sense??
				kpmSettings.MIN_L.put(id, lParam);
				kpmSettings.MAX_L.put(id, lParam);
				kpmSettings.INC_L.put(id,lParam);
				// end
				id2param.put(id, lParam);
				id2path.put(id, filePath);
				nextLine =cr.readNext();
			}

			kpmSettings.CASE_EXCEPTIONS_MAP = id2param;
			kpmSettings.MATRIX_FILES_MAP = id2path;

			cr.close();
		} catch (IOException ex) {
			Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
		}

		return params;
	}
}
