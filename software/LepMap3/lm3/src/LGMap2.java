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
// Stores a linkage map
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class LGMap2 {
	private ArrayList<Integer> lGNames = new ArrayList<Integer>();
	
	LGMap2()
	{
		lGNames.clear();

	}
	LGMap2(int numMarkers)
	{
		lGNames.clear();
		for (int i = 0; i < numMarkers; ++i)
			lGNames.add(0);
	}
	
	public int getNumLGs()
	{
		int numLGs = 0;
		for (int cn : lGNames)
			numLGs = Math.max(numLGs, cn);
		return numLGs;
	}
	
	
	private class intArrayComparator implements Comparator<int[]>{
		public int compare(int a1[], int a2[])
		{
			int r = 0;
			for (int i = 0; r == 0 && i < Math.min(a1.length, a2.length); ++i)
				r = new Integer(a1[i]).compareTo(a2[i]);
			return r;
		}
				
	}
	
	public void renameLGs()
	{
		int sizes[] = new int[getNumLGs() + 1];
		for (int m = 0; m < getNumMarkers(); ++m) {
			int name = lGNames.get(m); 
			if (name > 0)
				++sizes[name];
		}
		ArrayList<int[]> al = new ArrayList<int[]>();
		for (int m = 0; m < getNumMarkers(); ++m) {
			int name = lGNames.get(m); 
			if (name > 0 && sizes[name] > 0) {
				al.add(new int[]{-sizes[name], name});
				sizes[name] = 0;
			}
		}
		Collections.sort(al, new intArrayComparator());
		int names[] = new int[getNumLGs() + 1];
		int j = 0;
		for (int sn[] : al)
			names[sn[1]] = ++j;

		for (int m = 0; m < getNumMarkers(); ++m) {
			int name = lGNames.get(m); 
			if (name > 0)
				lGNames.set(m, names[name]);
		}
		
	}
	// returns markers of each LG 
	public ArrayList<ArrayList<Integer>> getMarkersOfLGs(int numMarkers)
	{
		if (numMarkers != getNumMarkers()) {
			System.err.println("Warning: The map has different number of markers as the data, ");
			System.err.println(numMarkers + " and " + getNumMarkers());			
		}
		//renameLGs();
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();
		int numLGs = getNumLGs();
		for (int i = 0; i <= numLGs; ++i) 
			ret.add(new ArrayList<Integer>());
		for (int m = 0; m < numMarkers; ++m) {
			int name = lGNames.get(m);
			ret.get(name).add(m);
		}
		return ret;
	}
	
	public int getLGName(int marker)
	{
		return lGNames.get(marker);
	}
	public void setLGName(int marker, int name)
	{
		lGNames.set(marker, name);
	}
	public ArrayList<Integer> getLGNames()
	{
		return lGNames;
	}
	public int getNumMarkers()
	{
		return lGNames.size();
	}

/*	public void printLGAssignment(Data d)
	{
		int singles = 0;
		int numChr = 0;
		int numMarkers = d.getNumMarkers();
		int names[] = new int[numMarkers];

		for (int m = 0; m < numMarkers; ++m) {
			String chrType = "";
			int name = lGNames.get(m);
			if (d.getMarkerType(m) == Family.Z_CHR && name != 0)
				chrType = "\tZ";
			System.out.println(name + chrType);
			if (name == 0)
				++singles;
			else
				if (names[name] == 0)
					names[name] = ++numChr;
		}
		System.err.println("Number of LGs = " + numChr + ", markers in LGs = " + (numMarkers - singles) + ", singles = " + singles);
	}*/
	public void printLGAssignment(int numMarkers)
	{
		int singles = 0;
		int numChr = 0;
		int names[] = new int[numMarkers];

		StringBuilder sb = new StringBuilder(); 
		for (int m = 0; m < numMarkers; ++m) {
			int name = lGNames.get(m);
			//System.out.println(name);
			sb.append(name);
			sb.append('\n');
			if (sb.length() > 10000 || (m == numMarkers - 1)) {
				System.out.print(sb);
				sb.setLength(0);
			}
				
			if (name == 0)
				++singles;
			else
				if (names[name] == 0)
					names[name] = ++numChr;
		}
		System.err.println("Number of LGs = " + numChr + ", markers in LGs = " + (numMarkers - singles) + ", singles = " + singles);
	}

	
	
	public void removeSmallLGs(int sizeLimit)
	{
		int numNames = getNumLGs();
		int sizes[] = new int[numNames + 1];
		int numMarkers = getNumMarkers();
		for (int i = 0; i < numMarkers; ++i)
			++sizes[lGNames.get(i)];
		for (int i = 0; i < numMarkers; ++i)
			if (sizes[lGNames.get(i)] < sizeLimit)
				lGNames.set(i, 0);
		for (int i = 1; i <= numNames; ++i)
			if (sizes[i] > 0 && sizes[i] < sizeLimit)
				System.err.println("Removing LG " + i + " (as too small) with size " + sizes[i] );
	}
	public boolean loadLGMap(String filename)
	{
		ArrayList<ArrayList<String>> table = Input.loadTable(filename, " \t");
		if (table == null)
			return false;

		lGNames.clear();

		for (ArrayList<String> line : table) {
			//if (line.size() != 1 && line.size() != 2)
			//	Error2.error(104);
			lGNames.add(Integer.parseInt(line.get(0)));
		}
		return true;
	}

}
