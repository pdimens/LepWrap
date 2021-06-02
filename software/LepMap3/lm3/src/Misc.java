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
//Miscellaneous functions 
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;
public class Misc {
	
	
	private static Random seedGenerator = new Random();
	private static long seed = seedGenerator.nextLong();
	private static Random rand = new Random(seed);
	
	public static void setSeed(long seed) {
		Misc.seed = seed;
		rand.setSeed(seed);
	}

	public static long getSeed() {
		return seed;
	}
	
	public static double random() {
		return rand.nextDouble();
	}
	
	//exponential distribution capped to [0,1[.
	public static double exponential01(double lambda) {
		double ret = 2.0;
		while (ret > 1.0)
			ret = -Math.log(Misc.random()) / lambda;
		return ret;
	}
	public static double entropy(double p[])
	{
		double sum = 0.0;
		for (double x : p)
			sum += x;
		
		assert(sum > 0);
		
		double iSum = 1.0 / sum; 
		double ret = 0;
		for (double x : p)
			if (x > 0) {
				double y = x * iSum;
				ret += y * Math.log(y);
			}
		return - ret * Constants.invLN2;
	}

	public static double entropy(float p[])
	{
		double sum = 0.0;
		for (double x : p)
			sum += x;
		
		assert(sum > 0);
		
		double iSum = 1.0 / sum; 
		double ret = 0;
		for (double x : p)
			if (x > 0) {
				double y = x * iSum;
				ret += y * Math.log(y);
			}
		return - ret * Constants.invLN2;
	}
	
	public static int[][] permutation(int n){
		if (n > 12)
			return null;
		int m = 1;
		ArrayList<Integer> elements = new ArrayList<Integer>(); 
		for (int i = 1; i <= n; ++i) {
			m *= i;
			elements.add(i - 1);
		}
			
		ArrayList<ArrayList<Integer>> perm = permutation(elements);
		int i = 0;
		int ret[][] = new int[m][n];
		for (ArrayList<Integer> p : perm) {
			for (int j = 0; j < n; ++j)
				ret[i][j] = p.get(j);
			++i;
		}
		return ret;
	}
	
	private static ArrayList<ArrayList<Integer>> permutation(ArrayList<Integer> elements){
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();
		if (elements.size() <= 1) {
			ArrayList<Integer> tmp = new ArrayList<Integer>();
			tmp.addAll(elements);
			ret.add(tmp);
			return ret;
		}
		for (int i = 0; i < elements.size(); ++i) {
			Integer element = elements.remove(i);
			ArrayList<ArrayList<Integer>> perm0 = permutation(elements);
			for (ArrayList<Integer> p : perm0) {
				p.add(0, element);
				ret.add(p);
			}
			elements.add(i, element);
		}
		return ret;
	}
	
	//no permutations that have been tried earlier
	public static int[][] prunePermutation1(int perm[][])
	{
		int k = perm[0].length;
		int ret[][] = new int[perm.length - perm.length / k + 1][];
		int i = 0;
		for (int p[] : perm) {
			if (p[0] != 0 || i == 0)
				ret[i++] = p;
		}
		return ret;
	}

	//no symmetric permutations
	public static int[][] prunePermutation2(int perm[][])
	{
		int k = perm[0].length;
		int ret[][] = new int[perm.length / 2][];
		int i = 0;
		for (int p[] : perm) {
			if (p[0] < p[k - 1])
				ret[i++] = p;
		}
		return ret;
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
	
	private static final double LN10 = Math.log(10);
	private static final double invLN10 = 1.0 / Math.log(10);

	public static double logSumFast(double logX, double logY)
	{
		if (logX >= logY) {
			if (logX == Double.NEGATIVE_INFINITY)
				return Double.NEGATIVE_INFINITY;
			return logX + log1pExpFast(logY - logX);
		}
    	else
			return logY + log1pExpFast(logX - logY);
	}
	
	public static double getLookup(double x) {
		if (x >= lookup.length - 1)
			return lookup[lookup.length - 1][0];
		int ix = (int) x;
		return lookup[ix][(int)(lookup[ix].length * (x - ix))];
	}
	
	private static double lookup[][] = null;
	static {
		int bins = 10;
		int size = 1 << bins;
		lookup = new double[bins + 1][];
		for (int i = 0; i <= bins; ++i) {
			lookup[i] = new double[size];
			for (int j = 0; j < size; ++j)
				lookup[i][j] = Math.log1p(Math.exp(-i - (j + 0.5) / size)); //add 0.5 to get mid point of each interval
			size >>= 1;
		}
		lookup[bins][0] = 0.0;
	}
	
	//x <= 0
	private static double log1pExpFast(double x) {
		return getLookup(-x);
	}
	
	
	public static double logSum(double logX, double logY)
	{
		if (logX >= logY) {
			if (logX == Double.NEGATIVE_INFINITY)
				return Double.NEGATIVE_INFINITY;
			return logX + Math.log1p(Math.exp(logY - logX));
		}
    	else
			return logY + Math.log1p(Math.exp(logX - logY));
		
	}
	public static double logSub(double logX, double logY)
	{
		return logX + Math.log1p(-Math.exp(logY - logX));
	}

	public static double exp10(double a)
	{
		return Math.exp(a * LN10);
	}
	private static class KthSmallest{
		double t[];
		public KthSmallest(double t[])
		{
			this.t = t;
		}
		public double getKthSmallest(int k){ // k = 0, 1, ... t.length - 1
			getKthSmallest(k, 0, t.length - 1);
			return t[k];
		}
		private void getKthSmallest(int k, int low, int high) {
			int i = low, j = high;
			double pivot = t[low + (high-low)/2];
			while (i <= j) {
				while (t[i] < pivot)
					i++;
				while (t[j] > pivot)
					j--;
				if (i <= j) {
					double tmp = t[i];
					t[i] = t[j];
					t[j] = tmp;
					i++;
					j--;
				}
			}
		    if (low < j && k <= j)
		    	getKthSmallest(k, low, j);
		    else if (i < high && k >= i)
		    	getKthSmallest(k, i, high);
		  }
	};
	public static double getKthSmallest(int k, double t[])
	{
		KthSmallest ks = new KthSmallest(t);  
		return ks.getKthSmallest(k);
	}
	public static int[] randomPermutation(int n){
		int perm[] = new int[n];
		randomPermutation(perm, n);
		return perm;
	}
	public static void randomPermutation(int perm[], int n){
		randomPermutation(perm, n, n - 1);
	}
	public static void randomPermutation(int perm[], int n, int k){
		for (int i = 0; i < n; ++i)
			perm[i] = i;
		for (int i = 0; i < k; ++i) {
			int j = i + (int) (Misc.random() * (n - i));
			int tmp = perm[i];
			perm[i] = perm[j];
			perm[j] = tmp;
		}
	}
	static float[] resizeArray(float array[], int newSize) {
		float[] na = new float[newSize];
		System.arraycopy(array, 0, na, 0, Math.min(newSize, array.length));
		return na;
	}
	static int[] resizeArray(int array[], int newSize) {
		int[] na = new int[newSize];
		System.arraycopy(array, 0, na, 0, Math.min(newSize, array.length));
		return na;
	}
	
	public static byte[][] transpose(byte array[][])
	{
		if (array.length == 0)
			return new byte[0][0]; //???		
		byte[][] ret = new byte[array[0].length][array.length];
		for (int i = 0; i < ret.length; ++i)
			for (int j = 0; j < ret[i].length; ++j)
				ret[i][j] = array[j][i];
		return ret;
	}
	
	public static void printArray(int array[][])
	{
		for (int a[] : array) {
			for (int ai : a) 
				System.out.print(ai + "\t");
			System.out.println();
		}
	}
	public static void printArray(double array[][])
	{
		for (double a[] : array) {
			for (double ai : a) 
				System.out.print(ai + "\t");
			System.out.println();
		}
	}

	public static void printList(double array[])
	{
		for (double ai : array) 
			System.out.print(ai + "\t");
		System.out.println();
	}
	
	public static void main(String args[])
	{
		for (BigDecimal x = new BigDecimal("0.0"); x.doubleValue() > -11; x=x.subtract(new BigDecimal("0.01"))) {
			double v1 = log1pExpFast(x.doubleValue());
			double v2 = Math.log1p(Math.exp(x.doubleValue()));
			System.err.println(x + " " + v1 + " " + v2 + " " + (v1 - v2));
		}
		
		
		//printArray(permutation(4));
		printArray(prunePermutation1(permutation(5)));
		
		
		
		
		for (int i = 0; i < 10; ++i)
			System.out.println(exponential01(10));
		
		int n = 10;
		double t[] = new double[n];
		for (int k = 0; k < n; ++k)
			t[k] = Math.random();
		double old = 0;
		for (int k = 0; k < n; ++k) {
			double sk = getKthSmallest(k, t);
			if (sk < old)
				System.err.println("error");
			System.err.println(sk + "\t" + 1.0 * k / n );
			old = sk;
		}
	}
	

	
}
