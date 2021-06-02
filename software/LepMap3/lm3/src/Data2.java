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

	Copyright (C) 2013-2014 Pasi Rastas, pasi.rastas@helsinki.fi, University of Helsinki
	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/

import java.io.BufferedReader;
import java.io.InputStreamReader;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class Data2 {
	// indexes in LINKAGE file
	static final int FAMILY_ID = 0;
	static final int INDIVIDUAL_ID = 1;
	static final int FATHER_ID = 2;
	static final int MOTHER_ID = 3;
	static final int SEX = 4;
	static final int TRAIT = 5;
	static final int FIRST_GENOTYPE = 6;
	
	private boolean phasedData = false; 
	
	private ArrayList<Family2> families = new ArrayList<Family2>();
	private ArrayList<String> markerNames = new ArrayList<String>(); 
	private ArrayList<Long> markerNameNumbers = new ArrayList<Long>(); 

	public String getMarkerName(int marker)
	{
		return markerNames.get(marker);
	}
	
	private boolean numberInit = false;
	HashMap<String, Integer> scaffoldHash = new HashMap<String, Integer>();
	
	
	public long getMarkerNameNumber(int marker)
	{
		if (!numberInit) {
			int n = 0;
			long mi = 0;
			for (String m : markerNames) {
				String s[] = m.split("\t");
				if (!scaffoldHash.containsKey(s[0]))
					scaffoldHash.put(s[0], n++);
				long pos = mi;
				try {
					pos = Long.parseLong(s[1]);
				} catch (Exception e) {
				}
				markerNameNumbers.add(((long) (scaffoldHash.get(s[0])) << 32) + pos); // max chromosome pos = 2**32 - 1
				++mi;
			}
			numberInit = true;
		}
		return markerNameNumbers.get(marker);
	}
	
	
	public String getIndividualName(int index)
	{
		return individualHash.get(index);
	}
	
	public boolean isPhased(){
		return phasedData;
	}

	public void setPhasedData(boolean phased){
		phasedData = phased;
	}
	public void setRandomPhase(){
		if (!phasedData)
			for (Family2 f : families)
				f.randomPhase();

	}
	
	public ArrayList<Family2> getFamilies()
	{
		return families;
	}
	public int getNumMarkers()
	{
		if (families.size() == 0)
			return 0;
		return getFamilies().get(0).getNumMarkers();
	}
	
	public void subSampleData(double rate)
	{
		if (rate < 1.0) {
			int numM = getNumMarkers();
			for (int m = 0; m < numM; ++m)
				if (Math.random() > rate)
					for (Family2 f: families)
						f.removeMarker(m);
		}
		
	}
	
	public int getNumFamilies()
	{
		return families.size();
	}
	
	private float[] addToArray(float values[], float array[], int start)
	{
		if (start + values.length > array.length)
			array = Misc.resizeArray(array, (int) (1 + 1.2 * (start + values.length)));
		for (float f : values)
			array[start++] = f;
		return array;
	}

	public void scale(double value) {
		for (Family2 f: families) {
			f.scale(value);
		}
	}

	public void setMinError(double value) {
		for (Family2 f: families) {
			f.minError((float) value);
		}
	}
	
	public void removeMarker(int m)
	{
		for (Family2 f: families)
			f.removeMarker(m);			
		
	}
	public void keepFamilies(ArrayList<String> famNames)
	{
		HashMap<String, Integer> nameHash = new HashMap<String, Integer>();
		for (String name : famNames)
			nameHash.put(name, 1);
		for (Family2 f: families)
			if (!nameHash.containsKey(f.getName()))
				for (int m = 0; m < f.getNumMarkers(); ++m)
					f.removeMarker(m);			
	}
	public void addFamilyFromFile(String filename, String informativeMask)
	{
		addFamilyFromFile(filename, informativeMask, null, false);
	}
	
	HashMap<Integer, String> individualHash = new HashMap<Integer, String>();
	
	public void addFamilyFromFile(String filename, String informativeMask, boolean mask[], boolean grandparentPhase)
	{
		// TODO: check consistency of the input pedigree file...

		int numIndividuals = 0;

		System.err.println("Loading file");
		
		DataParser dp = new DataParser();
		dp.loadFile(filename, null, null);
	
		ArrayList<Double> dl = dp.getNextLine(false, false);

		ArrayList<ArrayList<Integer>> familyIndex = dp.getFamilies(); 
		ArrayList<int[]> parentIndex = dp.getParentIndex();
		int individualIndex = 0;
		for (int fi = 0; fi < familyIndex.size(); ++fi ) {
			String familyName = dp.getFamilyName(fi);
			families.add(new Family2(familyName));
			int fis = familyIndex.get(fi).size();
			for (int i = 0; i < fis - 2; ++i)
				individualHash.put(individualIndex++, familyName + "\t" + dp.getIndividualName(familyIndex.get(fi).get(i)));
			numIndividuals += fis;

		}
		ArrayList<Integer> sex = dp.getSex();
		
		int lineNumber = 0;
		
		while (dl != null) {
			markerNames.add(dp.getMarkerName());
			
			if (mask != null && mask.length > lineNumber && !mask[lineNumber]) {
				for (int fi = 0; fi < familyIndex.size(); ++fi )
					families.get(fi).addMarker(); // add dummy marker
			}
			else {
				for (int fi = 0; fi < familyIndex.size(); ++fi ) {
					families.get(fi).addMarker(dl, familyIndex.get(fi), sex, informativeMask, grandparentPhase, parentIndex.get(fi));
				}
			}

			++lineNumber;
			if (lineNumber % 100000 == 0)
				System.err.println("processed " + lineNumber + " lines" );
			
			if (mask != null && mask.length > lineNumber && !mask[lineNumber])
				dp.getNextLine(false, true); //do not process line
			else
				dl = dp.getNextLine(false, false); // no skip
		}
		System.err.println("File loaded with " +  lineNumber + " SNPs" );

		System.err.println("Number of individuals = " + numIndividuals + " excluding grandparents");
		System.err.println("Number of families = " + familyIndex.size());
	}

	public char[] getIgnoreParentOrder(int marker) {
		int numF = getNumFamilies();
		char ret[] = new char[numF];
		Arrays.fill(ret, '_');
		String m = markerNames.get(marker);
		if (m.length() > numF) {
			int start = m.length() - numF;
			//only '+' and '_' allowed for the last characters...
			for (int f = 0; f < numF; ++f)
				switch (m.charAt(f + start)) {
					case '+' : break;
					case '_' : break;
					default: return ret;
				}
			for (int f = 0; f < numF; ++f)
				ret[f] = m.charAt(f + start);
		}
		return ret;
		
	}

	public void maskIgnoreParentOrder() {
		int numF = getNumFamilies();
		int numM = getNumMarkers();
out:	for (int mi = 0; mi < numM; ++mi) {
			String m = markerNames.get(mi);
			if (m.length() > numF) {
				int start = m.length() - numF;
				//only '+' and '_' allowed for the last characters...
				for (int i = 0; i < numF; ++i)
					switch (m.charAt(i + start)) {
						case '+' : break;
						case '_' : break;
						default: continue out;
					}
				for (int i = 0; i < numF; ++i)
					if (m.charAt(i + start) == '+')
						families.get(i).removeMarker(mi);
			}
		}
	}		

}
