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
//Should this work only for a single family?
public class SeparateIdenticals {
	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(401);
		String extraParameters = "";
		System.out.print("#java SeparateIdenticals ");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(401);
		pp.warning(new String[]{"data", "lodLimit", "informativeMask", "femaleTheta", "maleTheta", "sizeLimit", "theta", "numThreads", "numParts", "removeSingles", "lod3Mode", "keepRate", "betweenSameType", "families"});
	
		double theta = Double.parseDouble(pp.getValueAsString("theta", "0.0"));
		double theta1 = Double.parseDouble(pp.getValueAsString("maleTheta", "" + theta));
		double theta2 = Double.parseDouble(pp.getValueAsString("femaleTheta", "" + theta));
		
		Separate2 sc = new Separate2(); 
		sc.setTheta1(theta1);
		sc.setTheta2(theta2);
	
		int numDataFiles = pp.getNumberOfValues("data");
		if (numDataFiles == 0 || numDataFiles > 1)
			Error.error(402);
		for (int i = 0; i < numDataFiles; ++i) {
			String filename = pp.getValueAsString("data", i, null);
			if (filename == null)
				Error.error(403);
			sc.addFamilyFromFile(filename, pp.getValueAsString("informativeMask", "0123"));
		}
		
		double lodLimit[];
		int nv = pp.getNumberOfValues("lodLimit");
		if (nv >= 2) {
			lodLimit = new double[nv];
			for (int i = 0; i < nv; ++i)
				lodLimit[i] = Double.parseDouble(pp.getValueAsString("lodLimit", i, "10.0"));
		} else
			lodLimit = new double[]{Double.parseDouble(pp.getValueAsString("lodLimit", "10.0"))};
		
		sc.maskInformative(pp.getValueAsString("informativeMask", "0123"));
		int numThreads = Integer.parseInt(pp.getValueAsString("numThreads", "1"));
		
		sc.setLod3Mode(Integer.parseInt(pp.getValueAsString("lod3Mode", "1")));

		boolean sameType = pp.getValueAsString("betweenSameType", "0").equals("1");
		if (sameType) {
			if (lodLimit.length != 3) // either 3 values...
				Error.error(4499);
		}

		if (pp.getNumberOfValues("families") > 0) {
			ArrayList<String> f = pp.getValue("families");
			sc.keepFamilies(f);
		}
		
		LGMap2 map = sc.separateIdenticals(lodLimit, Integer.parseInt(pp.getValueAsString("numParts", "" + numThreads)), numThreads, pp.getValueAsString("removeSingles", "1").equals("1"), Double.parseDouble(pp.getValueAsString("keepRate", "1.0")), sameType);
		///if (pp.getValueAsString("outputData", "0").equals("1")) {
		//	sc.outputData(Integer.parseInt(pp.getValueAsString("sizeLimit", "0")));
		//} else {
			map.renameLGs();
			map.removeSmallLGs(Integer.parseInt(pp.getValueAsString("sizeLimit", "1")));
			map.printLGAssignment(sc.getNumMarkers());
		//}
	
		//data.maskInformative(pp.getValueAsString("informativeMask", "0123"));
		
		//sep.separateChromosomesLOD2(lodLimit);
		//LGMap cm = sep.getLGMap();
		//cm.removeSmallLGs(Integer.parseInt(pp.getValueAsString("sizeLimit", "0")));
		//cm.printLGAssignment(data);
	}
}
