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

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@helsinki.fi, University of Helsinki, University of Cambridge
*/

//TODO: Test that String builder works
public class OutputData {
	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(901);
		String extraParameters = "";
		System.out.print("#java OutputData ");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(901);
		pp.warning(new String[]{"data", "map", "sizeLimit", "lod3Mode", "informativeMask"});
	
		double theta = Double.parseDouble(pp.getValueAsString("theta", "0.0")); // not possible to set...
		double theta1 = Double.parseDouble(pp.getValueAsString("maleTheta", "" + theta)); // not possible to set...
		double theta2 = Double.parseDouble(pp.getValueAsString("femaleTheta", "" + theta)); // not possible to set...
		
		Separate2 sc = new Separate2(); 
		sc.setTheta1(theta1);
		sc.setTheta2(theta2);
	
		int numDataFiles = pp.getNumberOfValues("data");
		if (numDataFiles == 0 || numDataFiles > 1)
			Error.error(4402);
		for (int i = 0; i < numDataFiles; ++i) {
			String filename = pp.getValueAsString("data", i, null);
			if (filename == null)
				Error.error(4403);
			sc.addFamilyFromFile(filename, pp.getValueAsString("informativeMask", "0123"));
		}

		//sc.maskInformative(pp.getValueAsString("informativeMask", "0123"));
		
		sc.setLod3Mode(Integer.parseInt(pp.getValueAsString("lod3Mode", "1")));

		
		LGMap2 map = new LGMap2();
		map.loadLGMap(pp.getValueAsString("map", null));

		sc.maskInformative(pp.getValueAsString("informativeMask", "0123"));

		
		sc.joinIdenticalsFast(map);// check how works...
		
		sc.outputData(Integer.parseInt(pp.getValueAsString("sizeLimit", "0")));
		
		
		//sep.separateChromosomesLOD2(lodLimit);
		//LGMap cm = sep.getLGMap();
		//cm.removeSmallLGs(Integer.parseInt(pp.getValueAsString("sizeLimit", "0")));
		//cm.printLGAssignment(data);
	}
}
