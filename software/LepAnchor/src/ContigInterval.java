/**
    This file is part of Lep-Anchor.

    Lep-Anchor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Lep-Anchor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Lep-Anchor.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2019 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki
	
*/
import java.util.*;

public class ContigInterval implements Comparable<ContigInterval>{
	private String name;
	private long start;
	private long end;
	
	private long startI[];
	private long endI[];
	
	private int chr;
	
	private int currentOrientation = 0;
	private int rank = 0;
	
	public ArrayList<ArrayList<Marker>> markers = new ArrayList<ArrayList<Marker>>();
	
	public void setRank(int rank) {
		this.rank = rank;
	}
	
	public String getContig()
	{
		return name;
	}

	public int getChromosome()
	{
		return chr;
	}
	public void setChromosome(int value)
	{
		chr = value;
	}

	public long getMinStart() {
		if (startI != null)
			return startI[0];
		return start;
	}
	public long getMaxStart() {
		if (startI != null)
			return startI[1];
		return start;
	}

	public long getMaxEnd() {
		if (endI != null)
			return endI[1];
		return end;
	}
	public long getMinEnd() {
		if (endI != null)
			return endI[0];
		return end;
	}
	
	public long getStart() {
		return start;
	}
	
	public long[] getStartI()
	{
		long[] ret = new long[2];
		ret[0] = getMinStart();
		ret[1] = getMaxStart();
		return ret;
	}

	public long[] getEndI()
	{
		if (endI == null)
			return new long[]{end, end};

		long[] ret = new long[endI.length];
		ret[0] = getMinEnd();
		ret[1] = getMaxEnd();
		return ret;
	}


	public void setStart(long value) {
		start = value;
	}

	public long getEnd() {
		return end;
	}

	public void setEnd(long value) {
		end = value;
	}

	public int getRank(){
		return rank;
	}
	public void setOrientation(int o) {
		this.currentOrientation = o;
	}
	public int getOrientation(){
		return currentOrientation;
	}
	public void flipOrientation()
	{
		if (currentOrientation == 0)
			currentOrientation = -1;
		else 
			currentOrientation = -currentOrientation;
	}
	//copy ContigInterval
	public ContigInterval(ContigInterval c) {
		this.name = c.name;
		this.start = c.start;
		this.end = c.end;

		//TODO: consider Arrays.copyOf
		this.startI = c.startI;
		this.endI = c.endI;
		
		this.chr = c.chr;
		this.currentOrientation = c.currentOrientation;
		this.rank = c.rank;

		this.markers = new ArrayList<ArrayList<Marker>>();
		for (ArrayList<Marker> mm: c.markers) {
			ArrayList<Marker> mm_new = new ArrayList<Marker>(); 
			this.markers.add(mm_new);
			for (Marker m : mm)
				mm_new.add(new Marker(m, this));
		}
	}
	
	public ContigInterval(String name, long start, long end)
	{
		this.name = name;
		this.start = start;
		this.end = end;
	}

	public ContigInterval(String name, long start[], long end[])
	{
		this.name = name;
		if (start[0] == 1)
			this.start = 1;
		else
			this.start = (start[0] + start[1]) / 2;
		
		if (end.length > 2)
			this.end = end[1];
		else
			this.end = (end[0] + end[1]) / 2;
		startI = start;
		endI = end;
	}
	
	public boolean inside(long p)
	{
        if (p >= start && p <= end)
            return true;
        return false;
	}
    @Override
    public int hashCode() {
        return name.hashCode() + Long.hashCode(start);
    }
    @Override
    public boolean equals(Object o) {
        if (this == o) 
        	return true;
        if (o == null) 
        	return false;
        if (this.getClass() != o.getClass()) 
        	return false;
        
        ContigInterval cio = (ContigInterval) o;
        return (name.equals(cio.name) && start == cio.start);
    }
	@Override
	public int compareTo(ContigInterval o) {
		int strcmp = name.compareTo(o.name);
		if (strcmp != 0)
			return strcmp; 
		
		if (start < o.start)
			return -1;
		if (start > o.start)
			return 1;
		return 0;
	}
	@Override
	public String toString()
	{
		return name + "\t" + start + "\t" + end;
	}
	public void addMap()
	{
		markers.add(new ArrayList<Marker>());
	}
	public void addMarker(Marker m, int map)
	{
		while (markers.size() <= map)
			addMap();
		markers.get(map).add(m);
		
		if (chr > 0 && chr != m.getChromosome()) {
			System.err.println(this);
			System.err.println("Error: Multiple chromosomes in one interval in bed");
			System.exit(-1);
		}
		chr = m.getChromosome();
		m.ci = this;
	}
	public int getNumMarkers()
	{
		int ret = 0;
		for (ArrayList<Marker> m : markers)
			ret += m.size();
		return ret;
	}
	
	public ArrayList<Marker> getMarkers(int map) {
		if (markers.size() > map) {
			if (currentOrientation >= 0) {			
				ArrayList<Marker> ret = new ArrayList<Marker>();
				ret.addAll(markers.get(map));
				return ret;
			}
			else {
				ArrayList<Marker> ret = new ArrayList<Marker>();
				ret.addAll(markers.get(map));
				Collections.reverse(ret);
				return ret;
			}
		} else
			return new ArrayList<Marker>();
	}
	public ArrayList<Marker> getMarkersReverse(int map) {
		return getMarkers(map, (currentOrientation >= 0) ? -1 : 0);
	}

	
	public ArrayList<Marker> getMarkers(int map, int orientation) {
		if (markers.size() > map) {
			if (orientation == 0) {		
				ArrayList<Marker> ret = new ArrayList<Marker>();
				ret.addAll(markers.get(map));
				return ret;
			}
			else {
				ArrayList<Marker> ret = new ArrayList<Marker>();
				ret.addAll(markers.get(map));
				Collections.reverse(ret);
				return ret;
			}
		} else
			return new ArrayList<Marker>();
	}

	public Marker getMarker(int map, int index) {
		return markers.get(map).get(index);
	}

	//calculate gap length for
	//  ci in orientation -> c2 with alignment [astart to aend]
	public long calculateGapLength(int orientation, long astart, long aend){
		long l = 0;
		if (orientation == 0) { // + orientation
			if (astart > getEnd())
				l = astart - getEnd() - 1;
			else if (aend < getEnd())
				l = getEnd() - aend;
		} else { // - orientation
			if (aend < getStart())
				l = getStart() - aend - 1;
			else if (astart > getStart())
				l = astart - getStart();
		}
		return l;
	}

	//ci is the haplotype, calculate how much of it does not align...
	public long calculateGapLengthHaplotype(long astart, long aend){

		long g1 = Math.max(0,  astart - getStart()); // getStart()
		long g2 = Math.max(0,  getEnd() - aend); // getEnd()
		
		return g1 + g2; 
	}
	
	//split ContigInterval (this) to start..position and position+1,...end 
	public ContigInterval[] splitContigInterval(long position[])
	{
		if (position[0] < start || position[1] >= end) // no need to split
			return new ContigInterval[]{this};

		long pos = (position[0] + position[1]) / 2;
		
		ContigInterval c1 = new ContigInterval(name, new long[]{getMinStart(), getMaxStart()}, position);
		ContigInterval c2 = new ContigInterval(name, new long[]{position[0] + 1, position[1] + 1}, new long[]{getMinEnd(), getMaxEnd()});
		int numMaps = markers.size();
		
		for (int map = 0; map < numMaps; ++map) {
			ArrayList<Marker> list = markers.get(map);
			for (Marker m : list) {
				if (m.getPosition() > pos) {
					Marker mNew = new Marker(m);
					mNew.ci = c2;
					c2.addMarker(mNew, map);
				} else {
					Marker mNew = new Marker(m);
					mNew.ci = c1;
					c1.addMarker(mNew, map);
				}
			}
		}
		if (c1.getNumMarkers() == 0 || c2.getNumMarkers() == 0)
			return new ContigInterval[]{this};
		
		return new ContigInterval[]{c1, c2};
	}
	public ContigInterval[] splitContigInterval(long position1[], long position2[], long position3[])
	{
		if (position1[0] > position2[0]) {
			long tmp[] = position1;
			position1 = position2;
			position2 = tmp;
		}
		if (position1[0] > position3[0]) {
			long tmp[] = position1;
			position1 = position3;
			position3 = tmp;
		}
		ContigInterval sp1[] = splitContigInterval(position2, position3);
		ContigInterval sp2[] = sp1[0].splitContigInterval(position1);
		ArrayList<ContigInterval> list = new ArrayList<ContigInterval>();
		for (int i = 0; i < sp2.length; ++i)
			list.add(sp2[i]);
		for (int i = 1; i < sp1.length; ++i)
			list.add(sp1[i]);
		ContigInterval s[] = new ContigInterval[list.size()]; 
		return list.toArray(s);
			
	}

	public ContigInterval[] splitContigInterval(long position1[], long position2[])
	{
		if (position1[0] > position2[0]) {
			long tmp[] = position1;
			position1 = position2;
			position2 = tmp;
		}
		ContigInterval sp1[] = splitContigInterval(position1);
		if (sp1.length == 2) {
			ContigInterval sp2[] = sp1[1].splitContigInterval(position2);
			if (sp2.length == 2)
				return new ContigInterval[]{sp1[0], sp2[0], sp2[1]};
			else
				return sp1;
		} else
			return splitContigInterval(position2);
	}
	
}
