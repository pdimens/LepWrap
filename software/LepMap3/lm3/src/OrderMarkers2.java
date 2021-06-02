/**
    This file is part of Lep-MAP3.

    Lep-MAP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Lep-MAP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Lep-MAP.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2013-2017 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/

//Errorcodes 3xx
//Command line interface for marker ordering
import java.util.ArrayList;
public class OrderMarkers2 {

	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(301);
		String extraParameters = "";
		for (int i = 0; i < args.length; ++i) {
			extraParameters +=  " " + args[i];
		}
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(301);
		
		pp.warning(new String[]{"refineParentOrder", "maskIgnoreParentOrder", "calculateIntervals", "data", "informativeMask", "useKosambi", "useMorgan", "improveOrder", "evaluateOrder", "map", "numThreads", "numMergeIterations", "chromosome", "scale", "scaleMode", "minError", "outputPhasedData", "removeMarkers", "sexAveraged", "phasedData", "recombination1", "recombination2", "interference1", "interference2", "identicalLimit", "computeLODScores", "randomPhase", "randSeed", "grandparentPhase", "hyperPhaser", "phasingIterations", "selfingPhase", "families", "usePhysical"});

		if (pp.getValueAsString("randSeed", null) != null)
			Misc.setSeed(Long.parseLong(pp.getValueAsString("randSeed", null)));
		else
			extraParameters += " randSeed=" + Misc.getSeed();
			
		System.out.println("#java OrderMarkers2" + extraParameters); 
		
		Order o3 = new Order(); 
		Data2 data = new Data2();
		LGMap2 lgm = null;
		
		if (pp.getValueAsString("useKosambi", "0").equals("1"))
			o3.setKosambi();
		else if (pp.getValueAsString("useMorgan", "0").equals("1"))
			o3.setMorgan();
		o3.setNumThreads(Integer.parseInt(pp.getValueAsString("numThreads", "1")));
		o3.setNumMergeIterations(Integer.parseInt(pp.getValueAsString("numMergeIterations", "6")));
		o3.setSexAveraged(pp.getValueAsString("sexAveraged", "0").equals("1"));

		double r1 = Double.parseDouble(pp.getValueAsString("recombination1", "0.001"));
		double r2 = Double.parseDouble(pp.getValueAsString("recombination2", "0.001"));
		double i1 = Double.parseDouble(pp.getValueAsString("interference1", "0.001"));
		double i2 = Double.parseDouble(pp.getValueAsString("interference2", "0.001"));
		o3.setRecombination(r1, r2, i1, i2);

		o3.setIdenticalLimit(Double.parseDouble(pp.getValueAsString("identicalLimit", "0.01")));

		o3.setUsePhysical(Integer.parseInt(pp.getValueAsString("usePhysical", 0, "0")), Double.parseDouble(pp.getValueAsString("usePhysical", 1, "0.01")));

		System.err.println(pp.getValueAsString("usePhysical", 0, "0"));
		System.err.println(pp.getValueAsString("usePhysical", 1, "0.01"));
		
		String chr = pp.getValueAsString("chromosome", null);
		String mapFile = pp.getValueAsString("map", null);
		boolean mask[] = null;
		if (mapFile != null) {  
			lgm = new LGMap2();
			if (!lgm.loadLGMap(mapFile))
				Error.error(302);

			int numMarkers = lgm.getNumMarkers();  
			mask = new boolean[numMarkers];// make a mask to skip unnecessary markers
			int chromosome = (chr == null) ? -1 : Integer.parseInt(chr);
			for (int i = 0; i < numMarkers; ++i) {
				int lg = lgm.getLGName(i); 
				if (lg != 0 && (chromosome < 0 || chromosome == lg))
					mask[i] = true;
			}
			
		}
		
		String filename = pp.getValueAsString("data", null);
		if (filename == null)
			Error.error(301);
		
		boolean gpPhase = false;
		if (pp.getValueAsString("grandparentPhase", "0").equals("1")) {
			gpPhase = true;
			//System.err.println("Error:grandparentPhase not yet implemented");
			//System.exit(-1);
		}
		
		data.addFamilyFromFile(filename, pp.getValueAsString("informativeMask", "0123"), mask, gpPhase); //somehow input gpPhase
		
		if (lgm == null)
			lgm = new LGMap2(data.getNumMarkers());

		o3.init(data, lgm);

		
		String dataScale = pp.getValueAsString("scale", 0, "M/N");

		double maxDataScale = Double.parseDouble(pp.getValueAsString("scale", 1, "2"));
		double capDataScale = Double.parseDouble(pp.getValueAsString("scale", 2, "" + maxDataScale));
		//data.scale(dataScale);
		double minError = Double.parseDouble(pp.getValueAsString("minError", "0.001"));
		data.setMinError(minError);

		data.setPhasedData(pp.getValueAsString("phasedData", "0").equals("1") || gpPhase);
		
		if (pp.getValueAsString("randomPhase", "0").equals("1"))
			data.setRandomPhase();

		if (pp.getValueAsString("selfingPhase", "0").equals("1"))
			o3.setSelfingPhase(true);		
		
		o3.setDataScale(dataScale, maxDataScale, capDataScale, Integer.parseInt(pp.getValueAsString("scaleMode", "2")));
		
		o3.setHyperPhaser(pp.getValueAsString("hyperPhaser", "0").equals("1"));
		o3.setPhasingIterations(Integer.parseInt(pp.getValueAsString("phasingIterations", "1")));		

		for (int m = 0; m < pp.getNumberOfValues("removeMarkers"); ++m) 
			data.removeMarker(Integer.parseInt(pp.getValueAsString("removeMarkers", m, "0")) - 1);

		if (pp.getNumberOfValues("families") > 0) {
			ArrayList<String> f = pp.getValue("families");
			data.keepFamilies(f);
		}
		
		boolean maskParentOrder = pp.getValueAsString("maskIgnoreParentOrder", "0").equals("1"); 
		if (maskParentOrder)
			data.maskIgnoreParentOrder();
	
		
		String orderFile = pp.getValueAsString("evaluateOrder", null);
		boolean refineParentOrder = pp.getValueAsString("refineParentOrder", "0").equals("1"); 
		if (refineParentOrder && (maskParentOrder || gpPhase))
			Error.error(305);
		
		o3.setRefineParentOrder(refineParentOrder);
		
		o3.setPrintPhased(Integer.parseInt(pp.getValueAsString("outputPhasedData", "1")));
		if (orderFile != null) {
			o3.evaluateOrder(orderFile, pp.getValueAsString("improveOrder", "1").equals("1"));
		}
		else {
			if (mapFile == null)
				Error.error(303);	
			if (chr == null)
				o3.orderChromosomes(true, pp.getValueAsString("improveOrder", "1").equals("1"), pp.getValueAsString("evaluateIndividualErrors", "0").equals("1"));
			else
				o3.orderChromosomes(Integer.parseInt(chr), true, pp.getValueAsString("improveOrder", "1").equals("1"), pp.getValueAsString("evaluateIndividualErrors", "0").equals("1"));
		}
		
		if (pp.getValueAsString("computeLODScores", null) != null) {
			o3.computeLODScores(pp.getValueAsString("computeLODScores", null));		
		}

		if (pp.getValueAsString("calculateIntervals", 0, null) != null) {
			double limit = Double.parseDouble(pp.getValueAsString("calculateIntervals", 1, "1.0"));
			o3.calculateIntervals(pp.getValueAsString("calculateIntervals", 0, null), limit);		
		}
		
	}

}
