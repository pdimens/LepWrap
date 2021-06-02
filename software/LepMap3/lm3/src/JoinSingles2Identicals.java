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
// Joins singular markers to LGs
import java.util.ArrayList;
//Errorcodes 7xx
public class JoinSingles2Identicals {
	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(701);
		String extraParameters = "";
		System.out.print("#java JoinSingles2 ");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(701);

		pp.warning(new String[]{"map", "data", "lodLimit", "informativeMask", "femaleTheta", "maleTheta", "sizeLimit", "theta", "alleleLimit", "lodDifference", "lod3Mode", "betweenSameType", "numThreads"});
		
		double theta = Double.parseDouble(pp.getValueAsString("theta", "0.0"));
		double theta1 = Double.parseDouble(pp.getValueAsString("maleTheta", "" + theta));
		double theta2 = Double.parseDouble(pp.getValueAsString("femaleTheta", "" + theta));
		
		Separate2 sc = new Separate2(); 
		sc.setTheta1(theta1);
		sc.setTheta2(theta2);

		int numDataFiles = pp.getNumberOfValues("data");
		if (numDataFiles == 0)
			Error.error(4402);
		for (int i = 0; i < numDataFiles; ++i) {
			String filename = pp.getValueAsString("data", i, null);
			if (filename == null)
				Error.error(4403);
			sc.addFamilyFromFile(filename, pp.getValueAsString("informativeMask", "0123"));
		}
		
		double lodLimit[];
		int nv = pp.getNumberOfValues("lodLimit");
		if (nv >= 1) {
			lodLimit = new double[nv];
			for (int i = 0; i < nv; ++i)
				lodLimit[i] = Double.parseDouble(pp.getValueAsString("lodLimit", i, "10.0"));
		} else
			lodLimit = new double[]{Double.parseDouble(pp.getValueAsString("lodLimit", "10.0"))};
		
		LGMap2 map = new LGMap2();
		if (!map.loadLGMap(pp.getValueAsString("map", null)))
			Error.error(4406);
		
		sc.maskInformative(pp.getValueAsString("informativeMask", "0123"));
		
		sc.setLod3Mode(Integer.parseInt(pp.getValueAsString("lod3Mode", "1")));
		
		sc.joinSingles2Identicals(lodLimit, Double.parseDouble(pp.getValueAsString("lodDifference", "0.0")), map, pp.getValueAsString("betweenSameType", "0").equals("1"), Integer.parseInt(pp.getValueAsString("numThreads", "1")));

		//sc.SeparateChromosomes(lodLimit);
		//map.renameLGs();
		//map.removeSmallLGs(Integer.parseInt(pp.getValueAsString("sizeLimit", "1")));
		//map.printLGAssignment(sc.getNumMarkers());
		
	}
	
}
