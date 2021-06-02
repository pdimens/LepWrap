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
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
//Parses the pedigree and posterior data from a file
//takes into account grandparents (and removes them from families)
//getParentIndex gives the indexes for parents and grandparents...
public class DataParser {

	private String names[] = null;
	private String fathers[] = null;
	private String mothers[] = null;
	
	private int numIndividuals;
	
	private String markerName = "";
	
	private FileParser fp = null;
	private BufferedReader br = null;
	private String fn2 = null;
	private int lineNumber = 0;
	
	private ArrayList<ArrayList<Integer>> families = new ArrayList<ArrayList<Integer>>();
	private ArrayList<int[]> parentIndex = new ArrayList<int[]>();
	private ArrayList<String> familyNames = new ArrayList<String>(); 
	private ArrayList<Integer> sex = new ArrayList<Integer>(); 
	private HashMap<String, Integer> familyHash = new HashMap<String, Integer>();

	public DataParser() {
 	}
	
	public String getIndividualName(int index) {
		return names[index + 2];
	}

	public ArrayList<String> getFamilyNames(){
		 ArrayList<String> ret = new ArrayList<String>();
		 ret.addAll(familyNames);
		 return ret;
	}
	public String getFamilyName(int index){
		 return familyNames.get(index);
	}

	public String getFName() {
		return markerName;
	}
	
	public String getMarkerName() {
		return markerName;
	}
	public ArrayList<ArrayList<Integer>> getFamilies() {
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();
		for (ArrayList<Integer> al : families) {
			ArrayList<Integer> l = new ArrayList<Integer>();
			l.addAll(al);
			ret.add(l);
		}
		return ret;		
	}

	public ArrayList<int[]> getParentIndex() {
		ArrayList<int []> ret = new ArrayList<int[]>();
		for (int[] al : parentIndex) {
			int r[] = Arrays.copyOf(al, al.length);
			ret.add(r);
		}
		return ret;		
	}
	public ArrayList<Integer> getSex() {
		ArrayList<Integer> ret = new ArrayList<Integer>();
		ret.addAll(sex);
		return ret;		
	}
	public int getNumIndividuals() {
		return numIndividuals; 
	}
	
	private int[] getParentsAndGrandParents(String names[], String fathers[], String mothers[])
	{		
		names = Arrays.copyOfRange(names, 2, names.length);
		fathers = Arrays.copyOfRange(fathers, 2, fathers.length);
		mothers = Arrays.copyOfRange(mothers, 2, mothers.length);
		int fi = 0;
		for (ArrayList<Integer> familyIndex : families)
			parentIndex.add(getParentsAndGrandParents(names, fathers, mothers, familyIndex, familyNames.get(fi++)));
		return null;
	}

	private int[] getParentsAndGrandParents(String names[], String fathers[], String mothers[], ArrayList<Integer> familyIndex, String familyName)
	{		
		
		ArrayList<Integer> founders = new ArrayList<Integer>();
		ArrayList<Integer> children = new ArrayList<Integer>();

		HashMap<String, Integer> nameHash = new HashMap<String, Integer>();
		
		HashMap<String, Boolean> parents = new HashMap<String, Boolean>();

		for (int i: familyIndex) {
			if (fathers[i].equals("0") || mothers[i].equals("0"))
				founders.add(i);
			nameHash.put(names[i], i);
			if (!fathers[i].equals("0"))
				parents.put(fathers[i], true);
			if (!mothers[i].equals("0"))
				parents.put(mothers[i], true);
			//System.err.println(names[i]);
			
			//if (motherHash.containsKey(names[i]))
			//	Error.error(-99999);
			//if (fatherHash.containsKey(names[i]))
			//	Error.error(-99999);
			//fatherHash.put(names[i], fathers[i]);
			//motherHash.put(names[i], mothers[i]);
		}
		//System.err.println(parents);
		for (int i: familyIndex) {
			if (!parents.containsKey(names[i])) {
				children.add(i);
				//System.err.println(names[i]);
			}
		}
		
		int father = -1;
		int mother = -1; 
		for (int i: children) {
			if (!nameHash.containsKey(fathers[i]))
				Error.error(520, familyName + ":" + fathers[i]);

			if (father >= 0 && (father != nameHash.get(fathers[i]))) //two fathers
				Error.error(512, familyName + ":" + fathers[i]);
			//System.err.println(names[i] + " " + nameHash.containsKey(fathers[i]) + " " + fathers[i]);
			
			father = nameHash.get(fathers[i]);

			if (!nameHash.containsKey(mothers[i]))
				Error.error(520, familyName + ":" + mothers[i]);

			if (mother >= 0 && mother != nameHash.get(mothers[i])) //two mothers
				Error.error(512, familyName + ":" + mothers[i]);

			mother = nameHash.get(mothers[i]);
			
			if (sex.get(nameHash.get(mothers[i])) != 2) //mother not female 			
				Error.error(511, familyName + ":" + mothers[i]);
			if (sex.get(nameHash.get(fathers[i])) != 1) //father not male
				Error.error(511, familyName + ":" + fathers[i]);
		}
		int ret[] = new int[]{father, mother, -1, -1, -1, -1};
		
		if (!familyIndex.remove(new Integer(father)))
			Error.error(502);
		familyIndex.add(father);

		if (!familyIndex.remove(new Integer(mother)))
			Error.error(502);
		familyIndex.add(mother);
		
		//System.err.println(names[father]);
		//System.err.println(names[mother]);
		
		if (familyIndex.size() == 2 + children.size()) { // all done
			System.err.println("No grandparents present in family " + familyName);
		} else {
			ArrayList<Integer> gps = new ArrayList<Integer>(); 
			for (int i : founders) {
				if (i != father && i != mother) {
					gps.add(i);
				} else {
					if (mothers[father].equals(names[mother])) {
							System.err.println("Backcross to mother" + familyName);
							gps.add(i);
					}
					if (fathers[mother].equals(names[father])) {
							System.err.println("Backcross to father in family " + familyName);
							gps.add(i);
					}
				}
			}
			int numGps = 0;
			if (!fathers[father].equals("0")) {
				if (!nameHash.containsKey(fathers[father]))
					Error.error(513, familyName + ":" + fathers[father]);
				ret[2] = nameHash.get(fathers[father]);
				if (sex.get(ret[2]) != 1)
					Error.error(511,  familyName + ":" + fathers[father]);
				++numGps;
				if (ret[2] != father && ret[2] != mother) // backcross
					familyIndex.remove(new Integer(ret[2]));
				
			}
			if (!mothers[father].equals("0")) {
				if (!nameHash.containsKey(mothers[father]))
					Error.error(513, familyName + ":" + mothers[father]);
				ret[3] = nameHash.get(mothers[father]);
				if (sex.get(ret[3]) != 2)
					Error.error(511,  familyName + ":" + mothers[father]);
				++numGps;
				if (ret[3] != father && ret[3] != mother) // backcross
					familyIndex.remove(new Integer(ret[3]));			
				
			}
			if (!fathers[mother].equals("0")) {
				if (!nameHash.containsKey(fathers[mother]))
					Error.error(513, familyName + ":" + fathers[mother]);
				
				ret[4] = nameHash.get(fathers[mother]);
				if (sex.get(ret[4]) != 1)
					Error.error(511, familyName + ":" + fathers[mother]);
				++numGps;
				if (ret[4] != father && ret[4] != mother) // backcross
					familyIndex.remove(new Integer(ret[4]));	
			}
			if (!mothers[mother].equals("0")) {
				if (!nameHash.containsKey(mothers[mother]))
					Error.error(513, familyName + ":" + mothers[mother]);
				ret[5] = nameHash.get(mothers[mother]);
				if (sex.get(ret[5]) != 2)
					Error.error(511, familyName + "," + mothers[mother]);
				++numGps;
				if (ret[5] != father && ret[5] != mother) // backcross
					familyIndex.remove(new Integer(ret[5]));			
			}
			
			if (numGps != gps.size()) {
				//clone grandparents...
				//System.err.println(gps);
				//System.err.println(numGps);
				//Error.error(513);
				System.err.println("Warning: Different number of grandparents (" + numGps + " and " + gps.size() +  ") in family " + familyName);
			}

			System.err.println("Found " + gps.size() + " grandparents in family " + familyName);
		}
		
		return ret;
	}
	
	private interface FileParser{
		public String[] parseNext(String line);
	}
	
	
	private class VCFParser implements FileParser{
		private boolean firstLine = true;
		private int mapping[] = null;

		private double pow[] = new double[4001];
		
		private HashMap<String, String> gtHash = new HashMap<String, String>();   
		
		public VCFParser() {		
			for (int i = 0; i < pow.length; ++i)
				pow[i] = Misc.exp10(-i * 0.1);
			gtHash.put("./.", "1 1 1 1 1 1 1 1 1 1");
			gtHash.put("0/0", "1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("0/1", "0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("1/0", "0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("0/2", "0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("2/0", "0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("0/3", "0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("3/0", "0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("1/1", "0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("1/2", "0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001");
			gtHash.put("2/1", "0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001");
			gtHash.put("1/3", "0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001");
			gtHash.put("3/1", "0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001");
			gtHash.put("2/2", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001");
			gtHash.put("2/3", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001");
			gtHash.put("3/2", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001");
			gtHash.put("3/3", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000");
			
			gtHash.put(".|.", "1 1 1 1 1 1 1 1 1 1");
			gtHash.put("0|0", "1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("0|1", "0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("1|0", "0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("0|2", "0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("2|0", "0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("0|3", "0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("3|0", "0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("1|1", "0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001 0.001");
			gtHash.put("1|2", "0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001");
			gtHash.put("2|1", "0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001 0.001");
			gtHash.put("1|3", "0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001");
			gtHash.put("3|1", "0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001 0.001");
			gtHash.put("2|2", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001 0.001");
			gtHash.put("2|3", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001");
			gtHash.put("3|2", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000 0.001");
			gtHash.put("3|3", "0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 1.000");			
		}
		
		private double pow(int value) {
			return pow[Math.min(value, pow.length - 1)]; 
		}

		private String getField(String fields, int index)
		{
			int prevPos = 0;
		    int pos = fields.indexOf(':');

		    while (index > 0 && pos >= 0) {
		    	prevPos = pos + 1;
		        pos = fields.indexOf(':', pos + 1);
		    	--index;
		    }
		    if (index > 0) {// field not found
		    	//System.err.println(fields + "\t" + index);
		    	//Error.error(518);
		    	return "0"; // support for dropped trailing fields
		    }
		    	
		    int end = (pos < 0) ? fields.length() : pos;
		    return fields.substring(prevPos, end);
		}

		private String fieldPL2Posterior(String field)
		{
			int separators[] = new int[10];
		    int pos = field.indexOf(',');

		    int numSeparators = 0;
		    while (numSeparators < 10 && pos >= 0) {
		    	separators[numSeparators++] = pos; 
		    	pos = field.indexOf(',', pos + 1);
		    }
		    if (numSeparators == 2 || numSeparators == 5 || numSeparators == 9) {
			    if (numSeparators == 2) {
			    	return "" + pow(Integer.parseInt(field.substring(0              , separators[0]))) + " "            //AA
			                  + pow(Integer.parseInt(field.substring(separators[0]+1, separators[1]))) + " 0 0 "        //AC
			    			  + pow(Integer.parseInt(field.substring(separators[1]+1, field.length()))) + " 0 0 0 0 0"; //CC
			    }
			    if (numSeparators == 5) {
			    	//Order is AA,AC,CC,AG,CG,GG
			    	//         0  1  2  3  4  5
			    	return "" + pow(Integer.parseInt(field.substring(0              , separators[0]))) + " "     //AA
			                  + pow(Integer.parseInt(field.substring(separators[0]+1, separators[1]))) + " "     //AC
			                  + pow(Integer.parseInt(field.substring(separators[2]+1, separators[3]))) + " 0 "   //AG
			                  + pow(Integer.parseInt(field.substring(separators[1]+1, separators[2]))) + " "     //CC
			                  + pow(Integer.parseInt(field.substring(separators[3]+1, separators[4]))) + " 0 "   //CG
      		  			  	  + pow(Integer.parseInt(field.substring(separators[4]+1, field.length())))+ " 0 0"; //GG
			    }
			    if (numSeparators == 9) {
			    	//Order is AA,AC,CC,AG,CG,GG,AT,CT,GT,TT
			    	//         0  1  2  3  4  5  6  7  8  9
			    	return "" + pow(Integer.parseInt(field.substring(0              , separators[0]))) + " " //AA
			                  + pow(Integer.parseInt(field.substring(separators[0]+1, separators[1]))) + " " //AC
			                  + pow(Integer.parseInt(field.substring(separators[2]+1, separators[3]))) + " " //AG
			                  + pow(Integer.parseInt(field.substring(separators[5]+1, separators[6]))) + " " //AT
                		  	  + pow(Integer.parseInt(field.substring(separators[1]+1, separators[2]))) + " " //CC
        		  			  + pow(Integer.parseInt(field.substring(separators[3]+1, separators[4]))) + " " //CG
				              + pow(Integer.parseInt(field.substring(separators[6]+1, separators[7]))) + " " //CT
				              + pow(Integer.parseInt(field.substring(separators[4]+1, separators[5]))) + " " //GG
				              + pow(Integer.parseInt(field.substring(separators[7]+1, separators[8]))) + " " //GT
    		                  + pow(Integer.parseInt(field.substring(separators[8]+1, field.length())));     //TT
			    }
		    }
	    	return "1 1 1 1 1 1 1 1 1 1"; // unknown variant
	    	//Error.error(?); 
		}
		
		private String fieldGL2Posterior(String field)
		{
			int separators[] = new int[10];
		    int pos = field.indexOf(',');

		    int numSeparators = 0;
		    while (numSeparators < 10 && pos >= 0) {
		    	separators[numSeparators++] = pos; 
		    	pos = field.indexOf(',', pos + 1);
		    }
		    if (numSeparators == 2 || numSeparators == 5 || numSeparators == 9) {
			    if (numSeparators == 2) {
			    	return "" + Misc.exp10(Double.parseDouble(field.substring(0              , separators[0]))) + " "            //AA
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[0]+1, separators[1]))) + " 0 0 "        //AC
			    			  + Misc.exp10(Double.parseDouble(field.substring(separators[1]+1, field.length()))) + " 0 0 0 0 0"; //CC
			    }
			    if (numSeparators == 5) {
			    	return "" + Misc.exp10(Double.parseDouble(field.substring(0              , separators[0]))) + " "     //AA
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[0]+1, separators[1]))) + " "     //AC
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[2]+1, separators[3]))) + " 0 "   //AG
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[1]+1, separators[2]))) + " "     //CC
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[3]+1, separators[4]))) + " 0 "   //CG
      		  			  	  + Misc.exp10(Double.parseDouble(field.substring(separators[4]+1, field.length())))+ " 0 0"; //GG
			    }
			    if (numSeparators == 9) {
			    	return "" + Misc.exp10(Double.parseDouble(field.substring(0              , separators[0]))) + " " //AA
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[0]+1, separators[1]))) + " " //AC
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[2]+1, separators[3]))) + " " //AG
			                  + Misc.exp10(Double.parseDouble(field.substring(separators[5]+1, separators[6]))) + " " //AT
                		  	  + Misc.exp10(Double.parseDouble(field.substring(separators[1]+1, separators[2]))) + " " //CC
        		  			  + Misc.exp10(Double.parseDouble(field.substring(separators[3]+1, separators[4]))) + " " //CG
				              + Misc.exp10(Double.parseDouble(field.substring(separators[6]+1, separators[7]))) + " " //CT
				              + Misc.exp10(Double.parseDouble(field.substring(separators[4]+1, separators[5]))) + " " //GG
				              + Misc.exp10(Double.parseDouble(field.substring(separators[7]+1, separators[8]))) + " " //GT
    		                  + Misc.exp10(Double.parseDouble(field.substring(separators[8]+1, field.length())));     //TT
			    }
		    }
	    	return "1 1 1 1 1 1 1 1 1 1"; // unknown variant
	    	//Error.error(?); 
		}
		
		private String fieldGP2Posterior(String field)
		{
			int separators[] = new int[10];
		    int pos = field.indexOf(',');

		    int numSeparators = 0;
		    while (numSeparators < 10 && pos >= 0) {
		    	separators[numSeparators++] = pos; 
		    	pos = field.indexOf(',', pos + 1);
		    }
		    if (numSeparators == 2 || numSeparators == 5 || numSeparators == 9) {
			    if (numSeparators == 2) {
			    	return "" + (Double.parseDouble(field.substring(0              , separators[0]))) + " "            //AA
			                  + (Double.parseDouble(field.substring(separators[0]+1, separators[1]))) + " 0 0 "        //AC
			    			  + (Double.parseDouble(field.substring(separators[1]+1, field.length()))) + " 0 0 0 0 0"; //CC
			    }
			    if (numSeparators == 5) {
			    	return "" + (Double.parseDouble(field.substring(0              , separators[0]))) + " "     //AA
			                  + (Double.parseDouble(field.substring(separators[0]+1, separators[1]))) + " "     //AC
			                  + (Double.parseDouble(field.substring(separators[2]+1, separators[3]))) + " 0 "   //AG
			                  + (Double.parseDouble(field.substring(separators[1]+1, separators[2]))) + " "     //CC
			                  + (Double.parseDouble(field.substring(separators[3]+1, separators[4]))) + " 0 "   //CG
      		  			  	  + (Double.parseDouble(field.substring(separators[4]+1, field.length())))+ " 0 0"; //GG
			    }
			    if (numSeparators == 9) {
			    	return "" + (Double.parseDouble(field.substring(0              , separators[0]))) + " " //AA
			                  + (Double.parseDouble(field.substring(separators[0]+1, separators[1]))) + " " //AC
			                  + (Double.parseDouble(field.substring(separators[2]+1, separators[3]))) + " " //AG
			                  + (Double.parseDouble(field.substring(separators[5]+1, separators[6]))) + " " //AT
                		  	  + (Double.parseDouble(field.substring(separators[1]+1, separators[2]))) + " " //CC
        		  			  + (Double.parseDouble(field.substring(separators[3]+1, separators[4]))) + " " //CG
				              + (Double.parseDouble(field.substring(separators[6]+1, separators[7]))) + " " //CT
				              + (Double.parseDouble(field.substring(separators[4]+1, separators[5]))) + " " //GG
				              + (Double.parseDouble(field.substring(separators[7]+1, separators[8]))) + " " //GT
    		                  + (Double.parseDouble(field.substring(separators[8]+1, field.length())));     //TT
			    }
		    }
	    	return "1 1 1 1 1 1 1 1 1 1"; // unknown variant
	    	//Error.error(?); 
		}	
		
		
		private String fieldGT2Posterior(String field) {
			if (gtHash.containsKey(field))
				return gtHash.get(field);
			return "1 1 1 1 1 1 1 1 1 1"; // unknown variant
	    	//Error.error(?); 
		}				
		
		public String[] parseNext(String line) {
			if (firstLine) {
				if (line.length() < 2 || line.substring(0, 2).equals("##")) // comment in vcf
					return null;

				if (line.charAt(0) != '#') //individual names, first not starting with ##
					Error.error(517);
				
				String postNames[] = line.split("\t");
				HashMap<String, Integer> nameHash = new HashMap<String, Integer>();
				for (int i = 9; i < postNames.length; ++i) {
					if (nameHash.containsKey(postNames[i])) // same individual twice in posteriors...
						Error.error(515);
					else
						nameHash.put(postNames[i], i);
				}
				if (names == null) {
					names = new String[postNames.length - 7];
					names[0] = postNames[0];
					names[1] = postNames[1];
					for (int i = 2; i < names.length; ++i)
						names[i] = postNames[i + 7];
					numIndividuals = names.length - 2;
				}
				mapping = new int[names.length];
				for (int i = 2; i < names.length; ++i) {
					if (!nameHash.containsKey(names[i])) {
						System.err.println("Warning: Individual " + names[i] + " not contained in the data, set to all missing"); //posteriors not found
						mapping[i] = -1;
					} else
						mapping[i] = nameHash.get(names[i]);
				}
				firstLine = false;
			} else {
				String postData[] = line.split("\t");
				
				
				int fieldType = 0;
				String fields[] = postData[8].split(":"); //find PL field
				int index = 0;
				for (String f : fields) {
					if (f.equals("PL"))
						break;
					++index;
				}
				if (index == fields.length) {// PL field not found, try finding GL
					++fieldType; // fieldType = 1
					index = 0;
					for (String f : fields) {
						if (f.equals("GL"))
							break;
						++index;
					}
					if (index == fields.length) {// GL field not found, try finding GP
						++fieldType;  // fieldType = 2
						index = 0;
						for (String f : fields) {
							if (f.equals("GP"))
								break;
							++index;
						}
						if (index == fields.length) { //finally try to find a GT field
							++fieldType;  // fieldType = 3
							index = 0;
							for (String f : fields) {
								if (f.equals("GT"))
									break;
								++index;
							}
							if (index == fields.length) // not even GT not found
								Error.error(519);
						}
					}
				}
				
				boolean warning = false;
				String ret[] = new String[names.length];
				ret[0] = postData[0];
				ret[1] = postData[1];
				for (int i = 2; i < names.length; ++i)
					if (mapping[i] >= 0) {
						String field = getField(postData[mapping[i]], index);
						try {
							switch (fieldType) {
							case 0: ret[i] = fieldPL2Posterior(field); break;
							case 1: ret[i] = fieldGL2Posterior(field); break;
							case 2: ret[i] = fieldGP2Posterior(field); break;
							case 3: ret[i] = fieldGT2Posterior(field); break;
							}
						} catch (NumberFormatException e) { // types of .,.,. or 20,.,.
							warning = true;
							ret[i] = "1 1 1 1 1 1 1 1 1 1";
						}
					} else
						ret[i] = "1 1 1 1 1 1 1 1 1 1"; // set missing
				if (warning) {
					System.err.println("Warning: variant " + ret[0] + " " + ret[1] + " has strange value(s), set to missing");
				}
				return ret;
			}
			return null;
		}
	}
	
	private class PosteriorParser implements FileParser{
		private boolean firstLine = true;
		private int mapping[] = null;
		
		public String[] parseNext(String line) {
			if (firstLine) {
				String postNames[] = line.split("\t");
				HashMap<String, Integer> nameHash = new HashMap<String, Integer>();
				for (int i = 2; i < postNames.length; ++i) {
					if (nameHash.containsKey(postNames[i])) // same individual twice in posteriors...
						Error.error(515);
					else
						nameHash.put(postNames[i], i);
				}
				if (names == null) {
					names = postNames;
				  	numIndividuals = names.length - 2;		
				}
						
				mapping = new int[names.length];
				mapping[0] = 0;
				mapping[1] = 1;
				for (int i = 2; i < names.length; ++i) {
					if (!nameHash.containsKey(names[i])) { //posteriors not found
						System.err.println("Warning: Individual " + names[i] + " not contained in the data, set to all missing"); 
						mapping[i] = -1;
					} else
						mapping[i] = nameHash.get(names[i]);
				}
				firstLine = false;
			} else {
				String postData[] = line.split("\t");
				String ret[] = new String[names.length];
				for (int i = 0; i < names.length; ++i)
					if (mapping[i] >= 0)
						ret[i] = postData[mapping[i]];
					else
						ret[i] = "1 1 1 1 1 1 1 1 1 1"; // set missing
				return ret;
			}
			return null;
		}
	
	}	

	public boolean loadFile(String filename, String vcfFile, String posteriorFile) {
		fp = null;
		fn2 = posteriorFile; 
		if (posteriorFile != null || vcfFile != null) {
			if (posteriorFile != null && vcfFile != null)
				Error.error(514);
			if (posteriorFile != null)
				fp = new PosteriorParser();
			if (vcfFile != null) {
				fp = new VCFParser();
				fn2 = vcfFile; 
			}
		}
		try {
			if (filename != null) {
				if (filename.equals("-"))
					br = new BufferedReader(new InputStreamReader(System.in));
				else
					br = new BufferedReader(new FileReader(filename));
			}
		}  catch (Exception e) {
			e.printStackTrace();
			Error.error(504);
		}
		return true;
	}

	public ArrayList<Double> getNextLine(boolean printPedigree){
		return getNextLine(printPedigree, false, false); 
	}

	public ArrayList<Double> getNextLine(boolean printPedigree, boolean skip){
		return getNextLine(printPedigree, skip, false); 
	}

	public ArrayList<Double> getNextLine(boolean printPedigree, boolean skip, boolean genotypesOnly) {
		try {
			String line = "";

			String split[] = null;

			boolean firstLine = true;

			//allow data without pedigree... 
			if (genotypesOnly && lineNumber < 6) { 
				lineNumber = 6;
				if (fn2.equals("-"))
					br = new BufferedReader(new InputStreamReader(System.in));
				else
					br = new BufferedReader(new FileReader(fn2));
			}
			
out:		while (true) {
				do {
					line = br.readLine();
					
					if (line != null) {
						if (lineNumber < 6 || fp == null) {
							int index = line.indexOf('#');
							if (index >= 0)
								line = line.substring(0, index);
							if (lineNumber < 6 || !skip)
								split = line.split("\t");
						} else {
							split = fp.parseNext(line);
							if (split == null)
								line = "";
						}
						
					}
				} while (line != null && line.length() == 0);
				
				if (line == null) {
					if (fp != null && lineNumber == 6 && firstLine) {
						firstLine = false;
						if (fn2.equals("-"))
							br = new BufferedReader(new InputStreamReader(System.in));
						else
							br = new BufferedReader(new FileReader(fn2));
						continue out;
					}
					break;
				}
				++lineNumber;
				
				if (lineNumber < 7 && printPedigree)
					System.out.println(line);
					
				int end = line.indexOf('\t', line.indexOf('\t') + 1);
				markerName = line.substring(0, end);
				
				switch (lineNumber) {
				case 1:
					for (int i = 2; i < split.length; ++i) { 
						String fid = split[i];
						if (!familyHash.containsKey(fid)) {
							familyHash.put(fid, families.size());
							families.add(new ArrayList<Integer>()); 
							familyNames.add(fid); 
						}
						int index = familyHash.get(fid);
					  	families.get(index).add(i - 2);
					}
				  	numIndividuals = split.length - 2;
					break;
				case 2:
					names = split;
				  	break;
				case 3:
					fathers = split;
				  	break;
				case 4:
					mothers = split;
				  	break;
				case 5:
					for (int i = 2; i < split.length; ++i) 
						sex.add(Integer.parseInt(split[i]));
					break;
				case 6:
					//parsePedigree(names, fathers, mothers);
					getParentsAndGrandParents(names, fathers, mothers);
					System.err.println("Number of individuals = " + numIndividuals);
					System.err.println("Number of families = " + families.size());
					break;
				}
				if ((lineNumber < 7 || !skip) && split.length - 2 != numIndividuals && (lineNumber < 7 || split.length - 2 != 10 * numIndividuals)) {
					//System.err.println(columnIndex + "vs" + numMarkers);
					Error.error(503);
				}
				if (lineNumber >= 7) {
					if (skip)
						return null;
					ArrayList<Double> doubleLine = new ArrayList<Double>();
					for (int i = 2; i < split.length; ++i) {
						int start = 0;
						for (int j = 0; j < 10; ++j) {
							int stop = split[i].indexOf(" ", start);
							if (stop < 0)
								doubleLine.add(Double.parseDouble(split[i].substring(start)));
							else {
								doubleLine.add(Double.parseDouble(split[i].substring(start, stop)));
								start = stop + 1;
							}
							
						}
					}
					return doubleLine;
				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
			Error.error(504);
		}
		return null;
	}
	
}
