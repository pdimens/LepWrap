import java.util.ArrayList;

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

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/

public class SeparateChromosomes2 {
	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(1401);
		String extraParameters = "";
		System.out.print("#java SeparateChromosomes2 ");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(1401);
		pp.warning(new String[]{"data", "lodLimit", "informativeMask", "femaleTheta", "maleTheta", "sizeLimit", "theta", "numThreads", "lod3Mode", "subsample", "distortionLod", "map", "lg", "renameLGs", "families", "samplePairs", "phasedData", "grandparentPhase", "minLod"});
	
		double theta = Double.parseDouble(pp.getValueAsString("theta", "0.03"));
		double theta1 = Double.parseDouble(pp.getValueAsString("maleTheta", "" + theta));
		double theta2 = Double.parseDouble(pp.getValueAsString("femaleTheta", "" + theta));
		
		boolean gpPhase = false;
		if (pp.getValueAsString("grandparentPhase", "0").equals("1")) {
			gpPhase = true;
		}
		
		Separate2 sc = new Separate2(); 
		sc.setTheta1(theta1);
		sc.setTheta2(theta2);
		sc.setMinLod(Double.parseDouble(pp.getValueAsString("minLod", "" + Double.NEGATIVE_INFINITY)));

		String mapFile = pp.getValueAsString("map", null); 
		LGMap2 map2 = new LGMap2();
		int lg = Integer.parseInt(pp.getValueAsString("lg", "1"));
		boolean mask[] = null;
		if (mapFile != null) { //Refine previous lg assignment
			if (!(map2.loadLGMap(mapFile)))
					Error.error(1011);
			int nv = pp.getNumberOfValues("lg");
			if (nv > 1)
				Error.error(1012);

			int numMarkers = map2.getNumMarkers();  
			mask = new boolean[numMarkers];// make a mask to skip unnecessary markers
			for (int i = 0; i < numMarkers; ++i) {
				if (map2.getLGName(i) == lg) 
					mask[i] = true;
			}
		}
		

		String dataFile = pp.getValueAsString("data", null);
		if (dataFile == null)
			Error.error(1402);
		
		sc.addFamilyFromFile(dataFile, pp.getValueAsString("informativeMask", "0123"), mask, gpPhase);

		if (mapFile != null) {
			if (sc.getNumMarkers() != map2.getNumMarkers())
				Error.error(1013);
			sc.maskLG(map2, lg);
		}
		sc.setPhasedData(pp.getValueAsString("phasedData", "0").equals("1") || gpPhase);

		int nv = pp.getNumberOfValues("lodLimit");
		if (nv > 1)
			Error.error(1010);
		double lodLimit = Double.parseDouble(pp.getValueAsString("lodLimit", "10.0"));
		
		sc.maskInformative(pp.getValueAsString("informativeMask", "0123"));
		sc.setLod3Mode(Integer.parseInt(pp.getValueAsString("lod3Mode", "1")));

		double subFraction = Double.parseDouble(pp.getValueAsString("subsample", "1.0"));
		sc.subsample(subFraction);
		
		if (pp.getValueAsString("distortionLod", "0").equals("1"))
			sc.setDistortionLodMode();
		
		
		if (pp.getNumberOfValues("families") > 0) {
			ArrayList<String> f = pp.getValue("families");
			sc.keepFamilies(f);
		}
		
		
		LGMap2 map = sc.separateChromosomes(lodLimit, Integer.parseInt(pp.getValueAsString("numThreads", "1")), Double.parseDouble(pp.getValueAsString("samplePairs", "1.0")));
		
		if (mapFile != null) { //Refine previous lg assignment
			int numLGs = map2.getNumLGs();
			for (int i = 0; i < map.getNumMarkers(); ++i) {
				int lgMap = map.getLGName(i);
				int lgMap2 = map2.getLGName(i);
				if (lgMap == 0) {
					if (lgMap2 != lg)
						map.setLGName(i, lgMap2);
				} else { // lgMap > 0
					if (lgMap == 1) //largest group keeps name lg
						map.setLGName(i, lg);
					else
						map.setLGName(i, lgMap + numLGs - 1); //others go to the end
				}
			}
			if (pp.getValueAsString("renameLGs", "1").equals("1"))
				map.renameLGs();
		} else
			map.renameLGs();
		
		map.removeSmallLGs(Integer.parseInt(pp.getValueAsString("sizeLimit", "1")));
		map.printLGAssignment(sc.getNumMarkers());
	
		//data.maskInformative(pp.getValueAsString("informativeMask", "0123"));
		
		//sep.separateChromosomesLOD2(lodLimit);
		//LGMap cm = sep.getLGMap();
		//cm.removeSmallLGs(Integer.parseInt(pp.getValueAsString("sizeLimit", "0")));
		//cm.printLGAssignment(data);
	}
}
