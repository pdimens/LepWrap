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

	Copyright (C) 2013 Pasi Rastas, pasi.rastas@helsinki.fi, University of Helsinki
*/

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class Misc {
	
	//compare based on first element of long[]
	private static class LongComparator0 implements Comparator<long[]>
	{
		@Override
		public int compare(long[] l1, long[] l2){
			return Long.compare(l1[0], l2[0]);
		}
	}

	//compare based on second element of long[]
	private static class LongComparator1 implements Comparator<long[]>
	{
		@Override
		public int compare(long[] l1, long[] l2){
			return Long.compare(l1[1], l2[1]);
		}
	}
	//list elements contains 3 values, start, end and weight as long[3]
	//dynamic programming for finding coverage of multiple weighted intervals
	public static ArrayList<Long> cov(ArrayList<long[]> list_)
	{
		ArrayList<long[]> listStop = new ArrayList<long[]>();
		listStop.addAll(list_);
		ArrayList<long[]> list = new ArrayList<long[]>();
		list.addAll(list_);
		Collections.sort(list, new LongComparator0());
		Collections.sort(listStop, new LongComparator1());
		
		ArrayList<Long> ret = new ArrayList<Long>(); 
		int n = list.size();
		int end = 0;
		int start = 0;
		long pos = 0;

		long cov = 0;
		long ps = list.get(start)[0]; // start
		long pe = listStop.get(end)[1]; // end

		long prevP = 0;
		long prevC = 0;
		while (end < n) {
			if (start < n && ps <= pe) {
				pos = ps;
				while (start < n && pos == list.get(start)[0]) { // take all intervals with same start
					cov += list.get(start)[2];
					++start;
					if (start < n)
						ps = list.get(start)[0];
				}
			} else {
				pos = pe;
				while (end < n && pos == listStop.get(end)[1]) { // teka all intervals with same end
					cov -= listStop.get(end)[2];
					++end;
					if (end < n)
						pe = listStop.get(end)[1];
				}
				++pos;
			}
			if (pos != prevP) {
				//System.err.println(prevP + "\t" + prevC);
				ret.add(prevP);
				ret.add(prevC);
			}
			prevP = pos;
			prevC = cov;
		}
		ret.add(pos);
		ret.add(0L);
		//System.err.println(pos + "\t0");
		return ret;
	}

	
	public static boolean intersectIntervals(long b1, long e1, long b2, long e2){
		if (e1 < b2 || e2 < b1)
			return false;
		return true;
	}

	public static long intersectIntervalsLength(long b1, long e1, long b2, long e2){
		if (e1 <= b2 || e2 <= b1)
			return 0;
		return Math.min(e1 - b2, e2 - b1);
	}
	
	private static final double LN10 = Math.log(10.0);

	public static double exp10(double a)
	{
		return Math.exp(a * LN10);
	}
	
    public static class ArrayIndexComparator<T extends Comparable<T>> implements Comparator<Integer>
    {
            private final T[] array;

            public ArrayIndexComparator(T[] array)
        {
            this.array = array;
        }
        @SuppressWarnings("unchecked")
            public ArrayIndexComparator(ArrayList<T> alArray)
        {
                this.array = (T[]) alArray.toArray(new Comparable[0]);
        }       

        public Integer[] createIndexArray()
        {
            Integer[] indexes = new Integer[array.length];
            for (int i = 0; i < array.length; i++)
            {
                indexes[i] = i;
            }
            return indexes;
        }

            public int compare(Integer index1, Integer index2)
        {
            return array[index1].compareTo(array[index2]);
        }
    }       
    //Misc.ArrayIndexComparator<Double> comparator = new Misc.ArrayIndexComparator<Double>(table);
    //Integer[] indexes = comparator.createIndexArray();
    //Arrays.sort(indexes, comparator);     

    public static class ArrayIndexComparator2<T> implements Comparator<Integer>
    {
        private final T[] array;
        private Comparator<T> C;

        public ArrayIndexComparator2(T[] array, Comparator<T> C)
        {
            this.array = array;
            this.C = C;
        }
        @SuppressWarnings("unchecked")
        public ArrayIndexComparator2(ArrayList<T> alArray, Comparator<T> C)
        {
        	this.array = (T[]) alArray.toArray(new Comparable[0]);
            this.C = C;
        }       

        public Integer[] createIndexArray()
        {
            Integer[] indexes = new Integer[array.length];
            for (int i = 0; i < array.length; i++)
            {
                indexes[i] = i;
            }
            return indexes;
        }

        public int compare(Integer index1, Integer index2)
        {
            return C.compare(array[index1], array[index2]);
        }
    }       
    //Misc.ArrayIndexComparator2<Double> comparator = new Misc.ArrayIndexComparator2<Double>(table, customDoubleComparator);
    //Integer[] indexes = comparator.createIndexArray();
    //Arrays.sort(indexes, comparator);     
    
    //returns random merge of two lists...
    public static <T> ArrayList<T> randomMerge(ArrayList<T> l1, ArrayList<T> l2)
    {
    	ArrayList<T> ret = new ArrayList<T>();
		int n = l1.size();
		int m = l2.size();
		int j = 0;
		int k = 0;
		for (int i = 0; i < n + m; ++i) {
			if (j < n && (k == m || Math.random() * (n - j) > (m - k) * 0.5))
				ret.add(l1.get(j++));
			else
				ret.add(l2.get(k++));
		}
    	return ret;
    }
    
}
