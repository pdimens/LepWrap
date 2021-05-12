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
public class Marker implements Comparable<Marker>{
	String contig;
	int chromosome;
	long position;
	int intervals[];
	ContigInterval ci;
	
	int pPlus;
	int pMinus;


	public Marker(String contig, long position, int chromosome, int intervals[])
	{
		this.contig = contig;
		this.position = position;
		this.chromosome = chromosome;
		//TODO: consider cloning intervals...
		this.intervals = intervals;
	}

	public Marker(Marker m, long position)
	{
		this(m.contig, position, m.chromosome, m.intervals);
	}

	public Marker(Marker m)
	{
		this(m.contig, m.position, m.chromosome, m.intervals);
	}

	public Marker(Marker m, ContigInterval ci)
	{
		this(m.contig, m.position, m.chromosome, m.intervals);
		this.ci = ci;
	}

	//return next position inside from p or p + 1 if p is inside
	public int nextInside(int p) {
		int min = Integer.MAX_VALUE;
        for (int j = 0; j < intervals.length; j+=2)
            if (p <= intervals[j + 1] && min > intervals[j])
            	min = intervals[j];
        if (min <= p)
        	return p + 1;
        return min;
	}
	
	public int inside(int p)
	{
        for (int j = 0; j < intervals.length; j+=2)
            if (p >= intervals[j] && p <= intervals[j + 1])
                return 1;
        return 0;
	}
	@Override
	public String toString()
	{
		return contig + "\t" + position + "\t" + chromosome;
	}
	//orders first by contig and then by position...
	@Override
	public int compareTo(Marker o) {
		int ret = getContig().compareTo(o.getContig());
		if (ret != 0)
			return ret;
		if (position < o.position)
			return -1;
		if (position > o.position)
			return 1;
		return 0;
	}
	public String getName(){
		return contig + '\t' +  position;
	}
	public String getContig(){
		return contig;
	}
	public long getPosition(){
		return position;
	}
	
	public int getChromosome(){
		return chromosome;
	}
	public int[] getIntervals() {
		return intervals;
	}
	
	public int maxBin() // get maximum bin value...
	{
		int ret = 0;
        for (int j = 1; j < intervals.length; j+=2)
            ret = Math.max(ret, intervals[j]);
        return ret;
	}
	public int minBin() // get minimum bin value...
	{
		int ret = Integer.MAX_VALUE;
        for (int j = 0; j < intervals.length; j+=2)
            ret = Math.min(ret, intervals[j]);
        return ret;
	}
	public void flip(int maxBin)
	{
		int i[] = new int[intervals.length];
        for (int j = 0; j < intervals.length; ++j)
        	i[j] = maxBin - intervals[j ^ 1];
        intervals = i;
	}

	public ContigInterval getContigInterval() {
		return ci;
	}
	public static void main(String args[]){
		Marker m = new Marker("Contig", 1, 1, new int[]{10,20,40,60});
		for (int i = 0; i <= 61; ++i)
			System.out.println(i + ":" + m.nextInside(i));
	}

}
