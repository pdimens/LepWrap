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
    along with Lep-MAP3.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
*/

import java.util.ArrayList;
import java.util.Arrays;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReentrantLock;
import java.util.HashMap;

// Posterior aware (subset of) Filtering 
public class IBD {
//	private final ReentrantLock lock = new ReentrantLock();
	
	private final double SMALL = 1e-20; // to avoid problems with zero probabilities
	private int numIndividuals = 0;
	
	private double mafLimit=0.1;

	private double gLimit=0.1;
	
	ArrayList<float[][]> data = new ArrayList<float[][]>();
	ArrayList<double[]> frequencies = new ArrayList<double[]>();
	
	private double missingLimit = Double.MAX_VALUE;
	private double missingTolerance = 0.1;
	
	private final double LIKELIHOOD_TOLERANCE = 1e-12; // used in EM
	private final double LIKELIHOOD_TOLERANCE2 = 1e-6; // used in IBD
	
	private int numThreads = 4;

	String names[] = null;
	
	double biAllelicLimit = 1;
	
	public void setBiAllelic(double value){
		biAllelicLimit = value;
	}

	public void setMAFLimit(double value){
		mafLimit = value;
	}
	public void setGLimit(double value){
		gLimit = value;
	}	
	public void setMissingLimit(double value){
		missingLimit = value;
	}
	public void setMissingTolerance(double value){
		missingTolerance = value;
	}

	public void setNumThreads(int value){
		numThreads = value;
	}

	//print also...
	public void filter(DataParser dp) {
		int lineNumber = 0;
		while (true) {
			ArrayList<Double> doubleLine = dp.getNextLine(false, false, true);
			if (numIndividuals == 0)
				numIndividuals = dp.getNumIndividuals();
			if (doubleLine == null)
				break;

			if (lineNumber++ % 100000 == 99999)
				System.err.println("processed " + lineNumber + " lines");
			
			filter(dp.getMarkerName(), doubleLine);
		}
		System.err.println("Number of markers = " + lineNumber + " of which " + data.size() + " pass filtering");
		System.err.println("Number of individuals = " + numIndividuals);
		
		names = new String[numIndividuals + 2];
		for (int i = 0; i < numIndividuals; ++i)
			names[i + 2] = dp.getIndividualName(i);
	}
	
	public boolean filter(String prefix, ArrayList<Double> posteriors)
	{
		if (numMissingGenotypes(posteriors) > missingLimit)
			return false;
		
		int  numInd = posteriors.size() / 10;
		float p[][] = new float[numInd][10];
		for (int i = 0; i < numInd; ++i)
			for (int j = 0; j < 10; ++j)
				p[i][j]	= posteriors.get(10 * i + j).floatValue();
		double afreq[] = new double[4];
		double freq[] = new double[10];
		em(p, freq, afreq);
		
		double gfreq[] = new double[10];
		em_10(p, gfreq);
		
		int max = 0;
		for (int j = 1; j < 4; ++j)
			if (afreq[j] > afreq[max])
				max = j;
		int max2 = max ^ 1;
		for (int j = 0; j < 4; ++j)
			if (j != max)
				if (afreq[j] > afreq[max2])
					max2 = j;
	
		int maxG = 0;
		for (int j = 1; j < 10; ++j)
			if (gfreq[j] > gfreq[maxG])
				maxG = j;
		
		if (afreq[max] <= 1.0 - mafLimit && gfreq[maxG] <= 1.0 - gLimit && afreq[max] + afreq[max2] >= 1 - biAllelicLimit) {
			//lock.lock();
			//try {
				data.add(p);
				frequencies.add(freq);
			//} finally {
			//	lock.unlock();
			//}
			//frequencies.add(freq);
			//System.err.println(prefix + "\t" + freq[max] );
				return true;
		}
		return false;
	}
	
	private class FilterThread implements Runnable{
		String markers[] = null;
		ArrayList<Double> doubleLines[] = null;
		int index1 = 0;
		int index2 = 0;

		public FilterThread(int index1, int index2, ArrayList<Double> doubleLines[], String markers[]) {
			this.index1 = index1; 
			this.index2 = index2; 
			this.markers = markers; 
			this.doubleLines = doubleLines; 
		}
		public void run() {
			for (int i = index1; i < index2; ++i)
				filter(markers[i], doubleLines[i]);
		}
	}	
	
	private String s(double v)
	{
		int vi =  (int)(1000 * v + 0.5);
		if (vi == 0)
			return "0";
		if (vi == 1000)
			return "1";
			
		return "" + (vi / 1000.0);
	}
	
	public double em_10(float posteriors[][], double finalP[])
	{
		double p[] = new double[]{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
		double finalL = Double.NEGATIVE_INFINITY;
			
		for (int iteration = 0; iteration < 100; ++iteration) {
			double q[] = new double[10];
			double logL = 0.0;
			double ll = 1.0;
			for (float post[] : posteriors) {
				double l = 0.0;
				for (int j = 0; j < 10; ++j)
					l += (p[j] * post[j]); 
				double il = 1.0 / l;
				for (int j = 0; j < 10; ++j) 
					q[j] += (p[j] * post[j] * il);
				ll *= l;
				if (ll < 1e-200) { // log when needed
					logL += Math.log(ll);
					ll = 1.0;
				}
				
			}
			logL += Math.log(ll);

			double sum = 0.0; 
			for (int j = 0; j < 10; ++j)
				sum += q[j]; 
			double iSum = 1.0 / sum;
			for (int j = 0; j < 10; ++j)
				p[j] = q[j] * iSum;
			//System.err.println("ll = " + ll);
			double oldL = finalL;  
			finalL = logL;
			if (finalL - oldL < LIKELIHOOD_TOLERANCE && iteration >= 10)
				break;
			
		}
		//System.err.println("p = " + p[0] + "," + s(p[1]) + "," + s(p[2]) + "," + s(p[3]) + "," + s(p[4]) + "," + s(p[5]) + "," + s(p[6])  + "," + s(p[7]) + "," + s(p[8]) + "," + s(p[9]));
		for (int i = 0; i < 10; ++i)
			finalP[i] = p[i];
		
		return finalL;
	}

	public double em(float posteriors[][], double finalP[], double aFreq[])
	{
		double p[] = new double[]{0.25, 0.25, 0.25, 0.25};
		double p2[] = new double[]{0.25, 0.25, 0.25, 0.25};
		double finalL = Double.NEGATIVE_INFINITY;
			
		for (int iteration = 0; iteration < 100; ++iteration) {
			double q[] = new double[4];
			double logL = 0.0;
			double ll = 1.0;
			for (float post[] : posteriors) {
				double l = 0.0;
				for (int a = 0; a < 4; ++a)
					for (int b = 0; b < 4; ++b)
						l += (p[a] * p[b] * post[mapAlleles(a, b)]);
				
				double iL = 1.0 / l; 
				
				for (int a = 0; a < 4; ++a)
					p2[a] = p[a] * iL;
					
				for (int a = 0; a < 4; ++a)
					for (int b = 0; b < 4; ++b)
						q[a] += (p2[a] * p[b] * post[mapAlleles(a, b)]);

				ll *= l;
				if (ll < 1e-200) { // log when needed
					logL += Math.log(ll);
					ll = 1.0;
				}
			}
			logL += Math.log(ll);

			double sum = 0.0; 
			for (int j = 0; j < 4; ++j)
				sum += q[j]; 
			for (int j = 0; j < 4; ++j)
				p[j] = q[j] / sum;
			//System.err.println("ll = " + ll);
			double oldL = finalL;  
			finalL = logL;
			if (finalL - oldL < LIKELIHOOD_TOLERANCE && iteration >= 10)
				break;
			
		}
		
		//System.err.println("p = " + p[0] + "," + s(p[1]) + "," + s(p[2]) + "," + s(p[3]));
		
		Arrays.fill(finalP, 0.0);
		for (int a = 0; a < 4; ++a)
			for (int b = 0; b < 4; ++b)
				finalP[mapAlleles(a, b)] += p[a] * p[b];
		
		for (int a = 0; a < 4; ++a)
			aFreq[a] = p[a];
		
		return finalL;
	}
	
	public int numMissingGenotypes(ArrayList<Double> posteriors) {
		ArrayList<Integer> familyIndex = new ArrayList<Integer>();
		for (int i = 0; i <= numIndividuals; ++i)
			familyIndex.add(i);
		return numMissingGenotypes(posteriors, familyIndex);
	}
	
	public int numMissingGenotypes(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex)
	{
		int ret = 0;
		for (int i = 0; i < familyIndex.size() - 2; ++i) {
			int index = 10 * familyIndex.get(i);		

			double sum = 0.0;
			for (int j = 0; j < 10; ++j)
				sum += posteriors.get(index + j);

			if (sum >= 1.0 + missingTolerance)
				++ret;
		}
		return ret;
	}
	

	private int mapAlleles(int a1, int a2)
	{
		if (a1 > a2)
			return mapAlleles(a2, a1);	// a1 <= a2
		if (a2 > 3)
			Error.error(1005);
		if (a1 == 0)
			return a2;
		if (a1 == 1)
			return 3 + a2;
		if (a1 == 2)
			return 5 + a2;
		if (a1 == 3)
			return 9;
		Error.error(1005);
		return -1;
	}

	private int allele1(int g)
	{
		if (g < 4)
			return 0;
		if (g < 7)
			return 1;
		if (g < 9)
			return 2;
		if (g == 9)		
			return 3;
		Error.error(1005);
		return -1;
	}	
	
	private int allele2(int g)
	{
		if (g < 4)
			return g;
		if (g < 7)
			return g - 3;
		if (g < 9)
			return g - 5;
		if (g == 9)		
			return 3;
		Error.error(1005);
		return -1;
	}
	
	private int ibs(int g1, int g2){
		if (g1 == g2)
			return 2;
		
		int g1a1 = allele1(g1); 
		int g1a2 = allele2(g1); 
		int g2a1 = allele1(g2); 
		int g2a2 = allele2(g2);

		if (g1a1 == g2a1 || g1a1 == g2a2 || g1a2 == g2a1 || g1a2 == g2a2)
			return 1;
		
		return 0;
	}

	private final int PRINT_SIZE = 10000;
	
	
	private class PairwiseThread implements Runnable{
		private AtomicInteger index = null;
		StringBuilder sb = null;
		
		public PairwiseThread(AtomicInteger index) {
			this.index = index;
		}
		private void printResult() {
			synchronized(System.out) {
				System.out.print(sb.toString());
			}
			sb.setLength(0);
		}
		
		public void run() {
			sb = new StringBuilder(); 
			for (int ind1 = index.getAndIncrement(); ind1 < numIndividuals; ind1 = index.getAndIncrement())
				for (int ind2 = ind1 + 1; ind2 < numIndividuals; ++ind2) {
					sb.append(pairwise(ind1, ind2));
					sb.append('\n');
					if (sb.length() >= PRINT_SIZE)
						printResult();
						
				}
			printResult();
		}
	}

	
	
	public void pairwise()
	{
/*		String s[]={"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
		for (int g1 = 0; g1 < 10; ++g1) {
			for (int g2 = 0; g2 < 10; ++g2) {
				System.err.println(s[g1] + "\t" + s[g2] + "\t" + ibs(g1, g2));
			}
			System.err.println();
		}*/
		PairwiseThread pairwiseThreads[] = new PairwiseThread[numThreads];
		Thread threads[] = new Thread[numThreads];
		AtomicInteger index = new AtomicInteger();  
		for (int i = 0; i < numThreads; ++i) {
			threads[i] = new Thread(new PairwiseThread(index));
			threads[i].start();
		}
		for (int i = 0; i < numThreads; ++i) {
			while (threads[i].isAlive())
				try {
					threads[i].join();
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(-1);
				}
		}
	}

	private boolean validTrio(int p1, int p2, int c)
	{
		int p11 = allele1(p1);
		int p12 = allele2(p1);
		int p21 = allele1(p2);
		int p22 = allele2(p2);
		return (c == mapAlleles(p11, p21) || c == mapAlleles(p11, p22) || c == mapAlleles(p12, p21) || c == mapAlleles(p12, p22)); 
	}
	private boolean validDuo(int p1, int c)
	{
		int p11 = allele1(p1);
		int p12 = allele2(p1);
		int c1 = allele1(c);
		int c2 = allele2(c);
		return (c1 == p11 || c2 == p11 || c1 == p12 || c2 == p12); 
	}
	
	private class ParentCheckThread implements Runnable{
		int ind1, ind2;
		StringBuilder ret = new StringBuilder();
		
		public String getOutput(){
			return ret.toString();
		}
		public ParentCheckThread(int ind1, int ind2) {
			this.ind1 = ind1; 
			this.ind2 = ind2; 
		}
		public void run() {
			int numMarkers = data.size();
			if (ind2 < 0) { // check single parent case
				String prefix = names[ind1 + 2] + "\tNA\t";
				for (int ci = 0; ci < numChildren; ++ci) {
					int child = childIndex[ci];
					double err = 0.0;
					double n[] = new double[2];
					double p[] = new double[2];
					for (int m = 0; m < numMarkers; ++m) {
						Arrays.fill(p, 0.0);
						float posterior[][] = data.get(m);
			
						for (int g1 = 0; g1 < 10; ++g1)
							if (posterior[ind1][g1] > 1e-10)
								for (int g2 = 0; g2 < 10; ++g2)
									if (posterior[child][g2] > 1e-10) {
										int vIndex = validDuo(g1, g2) ? 0 : 1;
										p[vIndex] = Math.max(p[vIndex], posterior[ind1][g1] * posterior[child][g2]);
									}
						if (p[0] > p[1]) {
							n[0] += (p[0] - p[1]);
							err += (p[0] - p[1]) * p[1] / (p[0] + p[1]);
						}
						else {
							n[1] += (p[1] - p[0]);
							err += (p[1] - p[0]) * p[0] / (p[0] + p[1]);
						}
					}
					ret.append(prefix + names[child + 2] + '\t' + (Math.max(0, n[1] - err) / (n[0] + n[1] - err)) + '\t' + n[1] + '\t' + n[0] + '\t' + err + '\n');
				}
				
			} else { // two parents...
				String prefix = names[ind2 + 2] + "\t" + names[ind1 + 2] + "\t";
				
				for (int ci = 0; ci < numChildren; ++ci) {
					int child = childIndex[ci];
					double n[] = new double[2];
					double p[] = new double[2];
					double err = 0.0;
					for (int m = 0; m < numMarkers; ++m) {
						Arrays.fill(p, 0.0);
						float posterior[][] = data.get(m);
		
						for (int g1 = 0; g1 < 10; ++g1)
							if (posterior[ind1][g1] > 1e-10)
								for (int g2 = 0; g2 < 10; ++g2)
										if (posterior[ind2][g2] > 1e-10) {
											double p12 = posterior[ind1][g1] * posterior[ind2][g2];
											for (int g3 = 0; g3 < 10; ++g3)
												if (posterior[child][g3] > 1e-10) {
													int vindex = validTrio(g1, g2, g3) ? 0 : 1; 
													p[vindex] = Math.max(p[vindex], p12 * posterior[child][g3]);
												}
										}
							if (p[0] > p[1]) {
								n[0] += (p[0] - p[1]);
								err += (p[0] - p[1]) * p[1] / (p[0] + p[1]); 
							}
							else {
								n[1] += (p[1] - p[0]);
								err += (p[1] - p[0]) * p[0] / (p[0] + p[1]); 
							}
					}
					ret.append(prefix + names[child + 2] + '\t' + (Math.max(0, n[1] - err) / (n[0] + n[1] - err)) + '\t' + n[1] + '\t' + n[0] + '\t' + err + '\n');
				}
			}
		}
	}	
	int numChildren = 0;
	int numParents = 0;
	int parentIndex[] = null;
	int childIndex[] = null;
	
	public void checkParents(ArrayList<String> parents, boolean allpairs)
	{
		
		boolean isParent[] = new boolean[names.length];
		HashMap<String, Integer> hm = new HashMap<String, Integer>(); 
		
		int index = -2;
		for (String name : names) {
			if (index >= 0)
				hm.put(name, index);
			++index;
		}
		for (String parent : parents) {
			if (!hm.containsKey(parent))
				System.err.println("Warning: parent " + parent + " is missing");
			else
				isParent[hm.get(parent)] = true;
		}
		parentIndex = new int[numIndividuals];
		childIndex = new int[numIndividuals];
		for (int i = 0; i < numIndividuals; ++i) {
			if (isParent[i])
				parentIndex[numParents++] = i;
			if (!isParent[i] || allpairs)
				childIndex[numChildren++] = i;		
		}
		
		int activeThreads = 0;
		Thread threads[] = new Thread[numThreads];
		ParentCheckThread pthreads[] = new ParentCheckThread[numThreads];
		activeThreads = 0;
		for (int pi1i = -1; pi1i < numParents; ++pi1i) {
			int pi1 = (pi1i < 0) ? pi1i: parentIndex[pi1i];
			for (int pi2i = pi1i + 1; pi2i < numParents; ++pi2i) {
				int pi2 = parentIndex[pi2i];
				if (activeThreads < numThreads) {
					pthreads[activeThreads] = new ParentCheckThread(pi2, pi1);
					(threads[activeThreads] = new Thread(pthreads[activeThreads])).start();
					++activeThreads;
				}
				if (activeThreads == numThreads || pi1i == numParents - 2) {
					for (int i = 0; i < activeThreads; ++i) {
						if (threads[i].isAlive())
							try {
								threads[i].join();
							} catch (Exception e) {
								e.printStackTrace();
								System.exit(-1);
							}
						System.out.print(pthreads[i].getOutput());
					}
					activeThreads = 0;
				}
			}
		}
	}
	
	public String pairwise(int ind1, int ind2)
	{
		int numMarkers = data.size();
		
		double p[] = new double[3];
		for (int j = 0; j <= 2; ++j)
            p[j] = 1.0 / 3.0;
		double oldL = -1e100;
		
		for (int iteration = 0; iteration < 100; ++iteration) {
			double logL = 0.0;
			double q[] = new double[3];
			double dibs[] = new double[3];
			double pibs[] = new double[3];

			double ll = 1.0;
			for (int m = 0; m < numMarkers; ++m) {
				Arrays.fill(dibs, 0.0);
				Arrays.fill(pibs, 0.0);
				
				float posterior[][] = data.get(m);
				double freq[] = frequencies.get(m);
				
				for (int g1 = 0; g1 < 10; ++g1)
					for (int g2 = 0; g2 < 10; ++g2) {
                        int ibs = ibs(g1, g2);
                        //System.err.println(ibs);
                        dibs[ibs] = Math.max(dibs[ibs], posterior[ind1][g1] * posterior[ind2][g2]);
                        pibs[ibs] += freq[g1] * freq[g2];
	                }
				//for (int a = 0; a <= 2; ++a)
				//	System.err.println(pibs[a] + " vs " + dibs[a]);
				//System.err.println("*");

				double c2 = dibs[2];
				double c1 = (pibs[2] * dibs[2] + pibs[1] * dibs[1]) / (1 - pibs[0]);
				double c0 = pibs[2] * dibs[2] + pibs[1] * dibs[1] + pibs[0] * dibs[0];
				
				double c[] = {c0, c1, c2};
				double l = 0.0;
				for (int j = 0; j <= 2; ++j)
					l += p[j] * c[j];
				double iL = 1.0 / l;
				
				for (int j = 0; j <= 2; ++j)
					q[j] += p[j] * c[j] * iL;
				ll *= l;
				if (ll < 1e-200) { // log when needed
					logL += Math.log(ll);
					ll = 1.0;
				}
			}
			logL += Math.log(ll);
			double sum = 0.0;
			for (int j = 0; j <= 2; ++j)
				sum += q[j];
			
			double iSum = 1.0 / sum;
			for (int j = 0; j <= 2; ++j)
				p[j] = q[j] * iSum;
			
			//System.err.println("logL =" + ll);
			if (logL - oldL < LIKELIHOOD_TOLERANCE2)
				break;
			oldL = logL;
		}
		//System.err.println("**");
		return names[ind1 + 2] + "\t" + names[ind2 + 2]+ "\t" + s(p[2] + 0.5 * p[1]) + "\t" +  s(p[0]) + "\t" + s(p[1]) + "\t" + s(p[2]);
		
	}

	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(2001);
		String extraParameters = "";
		System.out.print("#java IBD");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(2001);
		pp.warning(new String[]{"data", "posteriorFile", "vcfFile", "MAFLimit", "missingLimit", "genotypeLimit", "numThreads", "parents", "allParentChildPairs"});
		
		String filename = pp.getValueAsString("data", null);
		DataParser dp = new DataParser();
		if (filename == null) {
			//only vcf or posterior input?
			String pfn = pp.getValueAsString("posteriorFile", null);
			String vfn = pp.getValueAsString("vcfFile", null);
			if (pfn != null || vfn != null) {
				dp.loadFile(null, vfn, pfn);
			}
			else
				Error.error(2001);
		} else
			dp.loadFile(filename, pp.getValueAsString("posteriorFile", null), pp.getValueAsString("vcfFile", null));

		IBD ibd = new IBD();

		ibd.setMAFLimit(Double.parseDouble(pp.getValueAsString("MAFLimit", "0.05")));
		ibd.setGLimit(Double.parseDouble(pp.getValueAsString("genotypeLimit", "0.1")));

		ibd.setMissingLimit(Double.parseDouble(pp.getValueAsString("missingLimit", 0, "" + Double.MAX_VALUE)));
		ibd.setMissingTolerance(Double.parseDouble(pp.getValueAsString("missingLimit", 1, "0.1")));

		ibd.setNumThreads(Integer.parseInt(pp.getValueAsString("numThreads", "4")));
		
		String ba = pp.getValueAsString("biAllelicLimit", null);
		if (ba != null)
			ibd.setBiAllelic(Double.parseDouble(ba));
		
		System.err.println("filtering markers...");
		ibd.filter(dp);
		
		if (!pp.getValueAsString("onlyFilter", "0").equals("1"))
			if (pp.getNumberOfValues("parents") > 0) {
				ArrayList<String> parents = new ArrayList<String>();  
				parents = pp.getValuesAsList("parents");
				//int n = parents.size();
				System.err.println("Starting computations to check parenthood...");
				//for (int p1 = 0; p1 < n ; ++p1)
				//	ibd.checkParent(parents.get(p1));
	
				//for (int p1 = 0; p1 < n ; ++p1)
				//	for (int p2 = p1 + 1; p2 < n; ++p2)
				//		ibd.checkParents(parents.get(p1), parents.get(p2));
				ibd.checkParents(parents, pp.getValueAsString("allParentChildPairs", "0").equals("1"));
			} else {
				System.err.println("Starting pairwise computations...");
				ibd.pairwise();
			}
	}
}
