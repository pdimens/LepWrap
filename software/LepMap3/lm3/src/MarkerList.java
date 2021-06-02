/**
    This file is part of Lep-MAP.

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

	Copyright (C) 2013 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki
*/
// Stores a list of markers
import java.util.ArrayList;
public class MarkerList {

	private UnionFind uf = null;
	private ArrayList<ArrayList<Integer>> markers=new ArrayList<ArrayList<Integer>>();
	
	public MarkerList(int numMarkers)
	{
		uf = new UnionFind(numMarkers);
		for (int i = 0; i < numMarkers; ++i) {
			ArrayList<Integer> al = new ArrayList<Integer>();
			al.add(i);
			markers.add(al);
		}
			
	}
	// This could be made more efficient by using (modified) double linked lists
	// However, this is probably as efficient as is needed
	public int union(int m1, int m2)
	{
		int c1 = uf.find(m1);
		int c2 = uf.find(m2);
		if (c1 == c2)
			return c1;
		int c = uf.union(c1, c2);
		if (c == c1) {
			markers.get(c1).addAll(markers.get(c2));
			markers.set(c2, null);
		} 
		else {
			markers.get(c2).addAll(markers.get(c1));
			markers.set(c1, null);
		}
		return c;
	}

	public int unionBetween(int m1, int m2, int markerToAdd)
	{
		int c1 = uf.find(m1);
		int c2 = uf.find(m2);
		if (c1 == c2) {
			int c = uf.union(c1, markerToAdd);
			for (int i = 0; i < markers.get(c2).size(); ++i) {
				int m = markers.get(c2).get(i);
				if (m == m1 || m == m2) {
					markers.get(c2).add(i, markerToAdd);
					return c;
				}
			}
			Error.error(504);
		}	else
			Error.error(503);
		return 0;
	}

	public int getFirstMarker(int c)
	{
		return markers.get(c).get(0);
	}
	public int getLastMarker(int c)
	{
		return markers.get(c).get(markers.get(c).size() - 1);
	}
	public int numMarkers(int c)
	{
		return markers.get(c).size();
	}
	
	
	// This could be made more efficient by using (modified) double linked lists
	// However, this is probably as efficient as is needed
	public ArrayList<Integer> unionMarkersWithOrder(int m1, int m2)
	{
		int c1 = uf.find(m1);
		int c2 = uf.find(m2);
		ArrayList<Integer> list1 = markers.get(c1); 
		ArrayList<Integer> list2 = markers.get(c2); 
		ArrayList<Integer> al = new ArrayList<Integer>();
		int s1 = list1.size();
		int s2 = list2.size();
		
		if (list1.get(s1 - 1) == m1)
			al.addAll(list1);  
		else 
			if (list1.get(0) == m1)
				for (int i = s1 - 1; i >= 0; --i)
					al.add(list1.get(i));
			else {
				System.err.println(list1);
				Error.error(501);
			}

		if (list2.get(0) == m2)
			al.addAll(list2);
		else 
			if (list2.get(s2 - 1) == m2)
				for (int i = s2 - 1; i >= 0; --i)
					al.add(list2.get(i));
			else {
				System.err.println(list2); 
				System.err.println(m2); 
				Error.error(502);
			}
		return al;
	}
	
	public void unionWithOrder(int m1, int m2)
	{
		ArrayList<Integer> al = unionMarkersWithOrder(m1, m2);
			
		int c1 = uf.find(m1);
		int c2 = uf.find(m2);
		int c = uf.union(c1, c2);
		markers.set(c1, null);
		markers.set(c2, null);
		markers.set(c, al);
	}
	public int find(int marker)
	{
		return uf.find(marker);
	}
	public boolean isEmpty(int m)
	{
		return (markers.get(m) == null);
	}
	public void removeMarkers(int c)
	{
		markers.get(c).clear();
	}
	
	public ArrayList<Integer> getMarkers(int c)
	{
		return markers.get(c);
	}
	public ArrayList<Integer> getMarkersCopy(int c)
	{
		return new ArrayList<Integer>(getMarkers(c)); // if we should copy
	}	


}
