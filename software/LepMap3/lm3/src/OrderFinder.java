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
    along with Lep-MAP2.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/
//TODO: superPhaserer
import java.util.ArrayList;

import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.*;
//import java.util.concurrent.atomic.DoubleAdder; in java 1.8
import java.io.PrintStream;

//TODO: Scaling to recombination parameters... To reduce end effects 
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;

public class OrderFinder {
    ArrayList<float[][][]> originalHaplotypes;
    ArrayList<int[]> originalInformative;
    ArrayList<float[]> originalInformation;
    ArrayList<char[]> originalParentOrder;

	private int numMergeIterations = 6;
	private int numOriginalMarkers;
	
	private double RECOMBINATION1 = 0.01;
	private double RECOMBINATION2 = 0.01;
	private double INTERFERENCE1  = 0.01;
	private double INTERFERENCE2  = 0.01;
	
	private int phasingIterations = 1;
	
	private String dataScaleString = "1.0";
	private double fullDataScale = 1.0;
	private double dataScale = 1.0;
	private double MAX_SCALE = 1.0;
	private double CAP_SCALE = 1.0;
	private int dataScalingMode = 1;
	
	private boolean selfingPhase = false;
	
	private boolean hyperPhasing = false;
	
	private double DISTANCE_LIMIT = 0.01;
	
	private boolean phasedData = false;
	
	private int usePhysicalPositions = 0;
	private double logPhysicalPenalty = 0.0;

	private int NUM_THREADS = 1;

	float qTables[][][] = null;
	double __tmpRec[][][] = new double[2][2][];

	//to avoid new calls in merge class 
	double forwardM1[] = null;
	double forwardM2[] = null;
	double backwardM1[] = null;
	double backwardM2[] = null;
	double scoresM[] = null;
	
	ArrayList<Integer> markerNames = new ArrayList<Integer>();
	
	long originalPhysicalPositions[] = null; 
	
	String individualNames[] = null;
	
	private SingleFamily[] sfs;

	private PhysicalFamily[] pfs;
	
	OrderFinder()
	{
	}

	// for future parallel implementation
	OrderFinder(OrderFinder other)
	{
		numMergeIterations = other.numMergeIterations;
		RECOMBINATION1 = other.RECOMBINATION1; 
		RECOMBINATION2 = other.RECOMBINATION2; 
		INTERFERENCE1 = other.INTERFERENCE1; 
		INTERFERENCE2 = other.INTERFERENCE2;
	
		phasingIterations = other.phasingIterations;
		
		dataScaleString = other.dataScaleString;
		fullDataScale = other.fullDataScale;
		dataScale = other.dataScale;
		MAX_SCALE = other.MAX_SCALE;
		CAP_SCALE = other.CAP_SCALE;
		dataScalingMode = other.dataScalingMode;
		hyperPhasing = other.hyperPhasing;
		DISTANCE_LIMIT = other.DISTANCE_LIMIT;
		phasedData = other.phasedData;
		selfingPhase = other.selfingPhase;
		NUM_THREADS = other.NUM_THREADS;
	}
	
	public void setDataScale(String s, double maxS, double capS, int mode)
	{
		MAX_SCALE = maxS;
		CAP_SCALE = 1.0 / capS;
		dataScaleString = s;
		dataScalingMode = mode;
	}
	public void setPhasingIterations(int value)
	{
		phasingIterations = value;
	}

	public void setSelfingPhase(boolean value)
	{
		selfingPhase = value;
	}

	public void setUsePhysical(int value, double prob)
	{
		usePhysicalPositions = value;
		logPhysicalPenalty = Math.log(prob);
	}
	
	private void initDataScale()
	{
		int numIndividuals = 0;
		for (SingleFamily sf : sfs)
			numIndividuals += sf.numIndividuals;
		if (dataScaleString.equals("M/N")) {
			dataScale = (double)numIndividuals / numOriginalMarkers;
		} else {
			int mpn = dataScaleString.indexOf("M/N");
			if (mpn >= 0) {
				if (mpn == 0 || mpn != dataScaleString.length() - 3) {
					Error.error(304);
					System.exit(-1);
				}
				dataScale = Double.parseDouble(dataScaleString.substring(0, mpn)) * numIndividuals / numOriginalMarkers;
			}  else {
				int pn = dataScaleString.indexOf("/N");
				if (pn >= 0) {
					if (pn == 0 || pn != dataScaleString.length() - 2) {
						//System.err.println("|" + dataScaleString + "|");
						Error.error(303);
						System.exit(-1);
					}
					dataScale = Double.parseDouble(dataScaleString.substring(0, pn)) / numOriginalMarkers;
				} else {
					dataScale = Double.parseDouble(dataScaleString);
				}
			}
		}
		if (dataScale > 1.0)
			dataScale = 1.0;
		fullDataScale = dataScale;
		System.err.println("Setting dataScale to " + dataScale + " all likelihoods will be multiplied by " + (1.0 / dataScale));
	}
	
	public void setRecombination(double rec1, double rec2, double inter1, double inter2) {
		RECOMBINATION1 = rec1;
		RECOMBINATION2 = rec2;
		INTERFERENCE1 = inter1;
		INTERFERENCE2 = inter2;
		initTmpRec();
	}
	
	private final double WORK_PER_THREAD = 300.0 * 300.0 * 100.0; // at least this amount of individual score evaluations should be done in a single thread...  
	private final int SPACE_PER_THREAD = 2000000; // space for each thread to store scores...  
	
	private int calcNumThreadsPolish(double numMarkers, double numPositions, double numIndividuals) {
		double n = numMarkers;
		double m = numPositions;
		double k = numIndividuals;
		double threads = n * m * k / WORK_PER_THREAD;
		if (threads >= 0.25 * numMarkers) 
			threads = 0.25 * numMarkers; // at least 4 jobs per thread, overhead <= 25%
		return (int) Math.min(NUM_THREADS, threads); 
	}

	private int calcNumThreadsScore(int numMarkers, int numPositions, int numIndividuals, int maxJobs) {
		double n = numMarkers;
		double m = numPositions;
		double k = numIndividuals;
		double threads = Math.min(maxJobs, n * m) * k / WORK_PER_THREAD;
		if (threads >= 0.25 * numMarkers * numPositions) 
			threads = 0.25 * numMarkers * numPositions; // at least 4 jobs per thread, overhead <= 25%
		return (int) Math.min(NUM_THREADS, threads); 
	}
	
	public void initTmpRec(){ //double scale) {

		double r1 = RECOMBINATION1;
		double r2 = RECOMBINATION2;
		double i1 = INTERFERENCE1;
		double i2 = INTERFERENCE2;
		
		//tmpRec[0][0] = new double[]{(1 - r1) * (1 - r2), i1 * (1 - r2), (1 - r1) * i2,  i1 * i2};
		//tmpRec[0][1] = new double[]{(1 - r1) * (1 - r2), i1 * (1 - r2), (1 - r1) * r2,  i1 * r2};
		//tmpRec[1][0] = new double[]{(1 - r1) * (1 - r2), r1 * (1 - r2), (1 - r1) * i2,  r1 * i2};
		//tmpRec[1][1] = new double[]{(1 - r1) * (1 - r2), r1 * (1 - r2), (1 - r1) * r2,  r1 * r2}; //phew...
		
		//assumed 1-r1 and 1-r2 are close to 1
		__tmpRec[0][0] = new double[]{1.0, i1, i2,  i1 * i2};
		__tmpRec[0][1] = new double[]{1.0, i1, r2,  i1 * r2};
		__tmpRec[1][0] = new double[]{1.0, r1, i2,  r1 * i2};
		__tmpRec[1][1] = new double[]{1.0, r1, r2,  r1 * r2}; //phew...
		
		for (double t2[][] : __tmpRec)
			for (double t[] : t2)
				for (int i = 0; i < t.length; ++i)
					t[i] = Math.log(t[i]);

		//for first version of recombination scaling
/*		double tmpRec1[][][] = new double[2][2][];
		double tmpRec2[][][] = new double[2][2][];
		
		double logI1 = Math.log(i1);
		tmpRec1[0][0] = new double[]{0.0, logI1, 0.0, logI1};
		tmpRec1[0][1] = new double[]{0.0, logI1, 0.0, logI1};
		double logR1 = Math.log(r1);
		tmpRec1[1][0] = new double[]{0.0, logR1, 0.0, logR1};
		tmpRec1[1][1] = new double[]{0.0, logR1, 0.0, logR1}; //phew...
		
		double logI2 = Math.log(i2);
		tmpRec2[0][0] = new double[]{0.0, 0.0, logI2, logI2};
		double logR2 = Math.log(r2);
		tmpRec2[0][1] = new double[]{0.0, 0.0, logR2, logR2};
		tmpRec2[1][0] = new double[]{0.0, 0.0, logI2, logI2};
		tmpRec2[1][1] = new double[]{0.0, 0.0, logR2, logR2}; //phew...*/		
	}

	public void setIdenticalLimit(double limit) {
		DISTANCE_LIMIT = limit;
	}
	
	
	public void setNumThreads(int numThreads)
	{
		NUM_THREADS = numThreads;
	}
	
	public void setNumMergeIterations(int numIterations)
	{
		numMergeIterations = numIterations;
	}
	
	public void setHyperPhasing(boolean value)
	{
		hyperPhasing = value;
	}

	public boolean isHyperPhasing()
	{
		return hyperPhasing;
	}
	

	public ArrayList<int[][]> getPhasedHaplotypes(ArrayList<Integer> markers, boolean silent, boolean maskHaplotypes)
	{
		numOriginalMarkers = originalHaplotypes.get(0).length;
		ArrayList<int[][]> ret = new ArrayList<int[][]>();
		ArrayList<int[]> retI = new ArrayList<int[]>();
		
		for (float h[][][] : originalHaplotypes) {
			int numIndividuals =  h[0].length;
			int ph[][] =  new int[numOriginalMarkers][numIndividuals * 2 + 2];
			ret.add(ph);
			retI.add(new int[numIndividuals]);
		}
		double c[] = new double[2];
		for (SingleFamily sf : sfs) {
			double tmp[] = sf.Viterbi2(markers, ret, retI, silent, maskHaplotypes);
			c[0] += tmp[0];
			c[1] += tmp[1];
		}
		if (!silent) {
			System.err.println("number of recombinations = " + c[0] + " logL = " + c[1]);
			int i = 0;
			for (int r[] : retI) {
				for (int j = 0; j < r.length; ++j)
					if (r[j] >= 1) {
						System.err.println("Individual\t" + individualNames[i + j] + "\trecombines\t" + r[j] + "\ttimes");
					}
				i += r.length;
			}
		}
		return ret;
	}

	public ArrayList<int[][]> getPhasedHaplotypes(ArrayList<Integer> markers)
	{
		return getPhasedHaplotypes(markers, false, false);
	}
	
	public ArrayList<float[][]> getPhasedQualityData(ArrayList<Integer> markers)
	{
		numOriginalMarkers = originalHaplotypes.get(0).length;
		ArrayList<float[][]> quality = new ArrayList<float[][]>();
		
		for (float h[][][] : originalHaplotypes) {
			int numIndividuals =  h[0].length;
			float q[][] =  new float[numOriginalMarkers][numIndividuals * 4];
			quality.add(q);
		}
		for (SingleFamily sf : sfs) {
			sf.Viterbi3(markers, quality);
		}
		return quality;
	}	

	//checks whether posterior is informative for paternal (parent=0) or maternal (parent=1) side
	private boolean nonMissing(int parent, float prob[]){
		final double MISSING_LIMIT = Math.log(0.8); // second posterior of 0.8 or higher is considered missing...
		if (parent == 0)
			return (prob[0] + prob[2] < MISSING_LIMIT && prob[1] + prob[3] < MISSING_LIMIT); // + is product but max(a,b) = 0 so a + b is min...
		else {
			assert(parent == 1);
			return (prob[0] + prob[1] < MISSING_LIMIT && prob[2] + prob[3] < MISSING_LIMIT); // + is product but max(a,b) = 0 so a + b is min...
		}
	}
	
	
	
	public void setHaplotypes(ArrayList<float[][][]> data, ArrayList<int[]> dataInf, ArrayList<char[]> dataParentOrder, Data2 data2, ArrayList<Integer> order)
	{
		if (order.size() == 0) // nothing to do...
			return;
		
		phasedData = data2.isPhased();

		originalHaplotypes = new ArrayList<float[][][]>();

		numOriginalMarkers = data.get(0).length;
		
//		haplotypes = new ArrayList<float[][][]>();
		originalInformative = new ArrayList<int []>();
		//informative = new ArrayList<int []>();

		originalInformation = new ArrayList<float []>();
		//information = new ArrayList<float []>();

		originalParentOrder = new ArrayList<char []>();
		
		int numIndividuals = 0;
		int fam = 0;
		for (float h[][][] : data) {
			char parentOrder[] = null;
			if (refineParentOrder)
				parentOrder = dataParentOrder.get(fam);
			
			int inf[] = new int[numOriginalMarkers];
			float inf2[] = new float[numOriginalMarkers];
			int numInd = (h.length > 0) ? h[0].length : 0;
			//float ohn[][][] = new float[numOriginalMarkers][numInd][5];
			float ohn[][][] = new float[numOriginalMarkers][numInd][4];
			float hn[][][] = new float[numOriginalMarkers][][];
			for (int m = 0; m < numOriginalMarkers; ++m) {
				double sume = 0.0; 
				for (int ind = 0; ind < numInd; ++ind) {
					for (int haplotype = 0; haplotype < 4; ++haplotype) 
						ohn[m][ind][haplotype] = (float) Math.log(h[m][ind][haplotype]);
					
					double e = 2.0 - Misc.entropy(h[m][ind]);
					sume += e;
					//ohn[m][ind][4] = (float) (e); //store information content for element 5... for future use ...
				}
				hn[m] = ohn[m];
				inf2[m] = (float) (sume / numInd);
				//System.err.println(m + "\t" + (sume / numInd));

				if (dataInf == null) {
					boolean allMissing1 = true;
					for (float hap[] : ohn[m])
						if (nonMissing(0, hap))
							allMissing1 = false;
		
					boolean allMissing2 = true;
					for (float hap[] : ohn[m])
						if (nonMissing(1, hap))
							allMissing2 = false;
		
					if (!allMissing1)
						inf[m] += 1;
					if (!allMissing2)
						inf[m] += 2;
				} else
					inf[m] = dataInf.get(fam)[m];
			}
			++fam;
			//haplotypes.add(hn);
			originalHaplotypes.add(ohn);

			originalInformative.add(inf);
			//informative.add(new int[numOriginalMarkers]);

			originalInformation.add(inf2);
			//information.add(new float[numOriginalMarkers]);
			
			originalParentOrder.add(parentOrder);
			
		//	
		//	
			numIndividuals += hn[0].length;
		}
		
		markerNames = order;
		
		if (usePhysicalPositions > 0) {
			//physicalPositions = new long[numOriginalMarkers];
			originalPhysicalPositions = new long[numOriginalMarkers];
			//order.size == numOriginalMarkers?
			for (int i = 0; i < order.size(); ++i) {
				originalPhysicalPositions[i] = data2.getMarkerNameNumber(order.get(i));
				//System.err.println(order.get(i) + "\t" + (data2.getMarkerNameNumber(order.get(i)) >>> 32));
			}

			pfs = new PhysicalFamily[1];
			pfs[0] = new PhysicalFamily();
		}
		
		individualNames = new String[numIndividuals];
		for (int i = 0; i < numIndividuals; ++i)
			individualNames[i] = data2.getIndividualName(i);
		
		qTables = new float[numIndividuals][numOriginalMarkers + 1][4];
		//mergePath = new boolean[numOriginalMarkers / 2 + 1][numOriginalMarkers / 2 + 2];	
		
		int numFamilies = originalHaplotypes.size();
		
		sfs = new SingleFamily[numFamilies];

		for (int i = 0; i < numFamilies; ++i) {
			sfs[i] = new SingleFamily(i);
		}
		
		forwardM1 = new double[numOriginalMarkers];
		forwardM2 = new double[numOriginalMarkers];
		backwardM1 = new double[numOriginalMarkers];
		backwardM2 = new double[numOriginalMarkers];
		scoresM = new double[SPACE_PER_THREAD * NUM_THREADS]; //space to store enough values to justify creating threads...
		
		initDataScale();
	}

	public void initMarkers(ArrayList<Integer> markers)
	{
	}

		
	public double scoreBetter(ArrayList<Integer> markers, double limit)
	{
		return score(markers);			
	}
	
	public double score(ArrayList<Integer> markers)
	{
		return likelihood(markers); 
	}

	public double likelihood(ArrayList<Integer> markers)
	{
		return likelihood(markers, false);
	}

	public double likelihood(ArrayList<Integer> markers, boolean keepTables)
	{
		initMarkers(markers);
		double ll = 0.0;

		/*if (NUM_THREADS > 1 && originalHaplotypes.size() > 1 && markers.size() >= 32) {
	        ExecutorService ftb = Executors.newFixedThreadPool(NUM_THREADS);
			CompletionService<Double> ecs = new ExecutorCompletionService<Double>(ftb);
			
			for (SingleFamily sf : sfs)
				ecs.submit(sf);
			
			try {
				for (int t = 0; t < sfs.length; ++t)
					ll += ecs.take().get();
			} catch (Exception e){
				Error.error(-999999);
			}
			if (keepTables)
				for (SingleFamily sf : sfs)
					ftb.execute(sf);
			
			ftb.shutdown();
		}  else*/ {
			for (SingleFamily sf : sfs)
				ll += sf.phase(markers);
			if (keepTables) {
				for (SingleFamily sf : sfs)
					sf.computeTables(markers);
			}
		}
		//physical info...
		if (pfs != null) {
			for (PhysicalFamily pf : pfs)
				ll += pf.likelihood(markers, keepTables);
		}
		
		return ll;
	}

	//TODO: DISTANCE_LIMIT
	private class PhysicalFamily{
		//private double logScore = Math.log(0.01);
		private ArrayList<Integer> initMarkers = null;
		
		private double scoreLong(long c1, long c2) {
			if ((c1 | 0xFFFFFFFFl) == (c2 | 0xFFFFFFFFl)) // same contig
				return 0;
			return logPhysicalPenalty; //TODO: position within contigs...
			// positions (c1 & 0xFFFFFFFFl) and (c2 & 0xFFFFFFFFl)
		}
		
		
		private double likelihood(ArrayList<Integer> markers_, boolean keepTables)
		{
			int numMarkers = markers_.size();
			double ll = 0.0;
			for (int m = 1; m < numMarkers; ++m)
				ll += scoreLong(originalPhysicalPositions[markers_.get(m - 1)], originalPhysicalPositions[markers_.get(m)]);
			//use physicalPosition
			if (keepTables)
				initMarkers = markers_;
			return ll;
		}
		
		private double polishScore(int perm[], int pos, ArrayList<Integer> markers)
		{
			double ll = 0.0;
			int k = perm.length;
			for (int i = 1; i < k; ++i)
				ll += scoreLong(originalPhysicalPositions[perm[i - 1]], originalPhysicalPositions[perm[i]]);
			
			if (pos > 0)
				ll += scoreLong(originalPhysicalPositions[markers.get(pos - 1)], originalPhysicalPositions[perm[0]]);

			if (pos + k < markers.size())
				ll += scoreLong(originalPhysicalPositions[perm[k - 1]], originalPhysicalPositions[markers.get(pos + k)]);
			
			return ll;
		}
		private double scoreDistance(int pos1, int pos2)
		{
			return scoreLong(originalPhysicalPositions[initMarkers.get(pos1)], originalPhysicalPositions[pos2]);
		}

		// assumes likelihood(markers, true) has been called
		private double score(int marker, int pos)
		{
			//if (true)
			//	return 0.0;
			int prev = pos - 1;
			int next = pos;
			if (prev >= 0 && initMarkers.get(prev).equals(marker))
				--prev;
			if (next < initMarkers.size() && initMarkers.get(next).equals(marker))
				++next;
			double ll = 0.0;
			if (prev >= 0)
				ll += scoreLong(originalPhysicalPositions[initMarkers.get(prev)], originalPhysicalPositions[marker]);
			if (next < initMarkers.size())
				ll += scoreLong(originalPhysicalPositions[marker], originalPhysicalPositions[initMarkers.get(next)]);
			
			return ll;
		}
	}

	private class SingleFamily implements Callable<Double>, Runnable{
		public int fam;
		
		public int numIndividuals;
		
		//private int hyperMarkers[] = new int[numOriginalMarkers];
		double forwardSP[][][];
		double backwardSP[][][];
		
		double logLs[];
		
		double tmpB[] = new double[16];
		
		double forwardV[][];
		double forwardV2[][];
		int path[][];
		
		RecombinationScale recScale = (dataScalingMode == 2) ? new RecombinationScale2() : new RecombinationScale1();

		public int getFamilyIndex()
		{
			return fam;
		}
		
		//index 4 binary indexes into one array index
		private int index4(int i1, int i2, int i3, int i4)
		{
			return 8 * i1 + 4 * i2 + 2 * i3 + i4; 
			//return (i1 << 3) + (i2 << 2) + (i3 << 1) + i4;
//			return ((2 * i1 + i2) * 2 + i3) * 2 + i4; 
//			int ret = i1 + i1 + i2;
//			ret += ret + i3;
//			ret += ret + i4;
//			return ret;
		}
		
		public SingleFamily(int fam)
		{
			this.fam = fam;
			
			numIndividuals = originalHaplotypes.get(fam)[0].length;
			
			forwardSP = new double[numIndividuals][numOriginalMarkers + 1][16];
			backwardSP = new double[numIndividuals][2][16];
			logLs = new double[numIndividuals];
			
			for (int hi = 0; hi < numIndividuals; ++hi) {
				for (int r1 = 0; r1 < 2; ++r1) 
					for (int r2 = 0; r2 < 2; ++r2)
						if (r1 > 0 || r2 > 0) {
							forwardSP[hi][0][index4(0,0,r1,r2)] = Double.NEGATIVE_INFINITY;
							forwardSP[hi][0][index4(0,1,r1,r2)] = Double.NEGATIVE_INFINITY;
							forwardSP[hi][0][index4(1,0,r1,r2)] = Double.NEGATIVE_INFINITY;
							forwardSP[hi][0][index4(1,1,r1,r2)] = Double.NEGATIVE_INFINITY;
						}
			}		
			
			forwardV = new double[numOriginalMarkers + 1][16];
			for (int r1 = 0; r1 < 2; ++r1) 
				for (int r2 = 0; r2 < 2; ++r2)
					if (r1 > 0 || r2 > 0) {
						forwardV[0][index4(0,0,r1,r2)] = Double.NEGATIVE_INFINITY;
						forwardV[0][index4(0,1,r1,r2)] = Double.NEGATIVE_INFINITY;
						forwardV[0][index4(1,0,r1,r2)] = Double.NEGATIVE_INFINITY;
						forwardV[0][index4(1,1,r1,r2)] = Double.NEGATIVE_INFINITY;
					}

			forwardV2 = new double[numOriginalMarkers + 1][16];

			
			path =  new int[numOriginalMarkers + 1][16];

			
		}

		//swaps t[i] and t[j]
		private void swap(int i, int j, float t[])
        {
                float tmp = t[i];
                t[i] = t[j];
                t[j] = tmp;                     
        }
		

		public double phase(ArrayList<Integer> markers)
		{
			double ret = 0.0;
			//System.err.println("phasing...");
			
//			int numInfMarkers = 0;
			//numInfMarkers1 = 0;
			//numInfMarkers2 = 0;

//			int inf[] = originalInformative.get(fam);
//			int numMarkers = markers_.size();
			//int markers[] = new int[numMarkers];
			
			
			
			//for (int mi = 0; mi < numMarkers; ++mi) {
			//	int m = markers_.get(mi);
			//	if (inf[m] > 0) {
			//		markers[numInfMarkers++] = m;
			//		//if ((inf[m] & 1) == 1)
					//	++numInfMarkers1;
					//if ((inf[m] & 2) == 2)
					//	++numInfMarkers2;			

			//	}
			//}

			//TODO: Check if this is really needed (or is markers_ <=> markers enough)
			int inf[] = originalInformative.get(fam);
			ArrayList<Integer> markers_ = new ArrayList<Integer>(); 
			for (int m : markers)
				if (inf[m] > 0)
					markers_.add(m);
			recScale.init(markers_);
			
			//System.err.println(numMarkers1 + "\t" + numMarkers2 + "\t" + numMarkersAll);
			
			if (phasedData) {
				ret += Viterbi2Fast(markers_);
			} else {
				if (hyperPhasing) {
					double ll = hyperPhaser(markers_); 
					//iterate
					if (phasingIterations > 1) {
						//System.err.println(ll);
						float phap[][][] = getPhase(markers_);
						for (int i = 1; i < phasingIterations; ++i) {
							randomPhase(markers_);
							double ll2 = hyperPhaser(markers_);
							//System.err.println(ll2);
							if (ll2 > ll) {
								ll = ll2;
								phap = getPhase(markers_);
							}
						}
						setPhase(markers_, phap);						
					}
					ret += ll;
				} else {
					//double l[];
					//while ((l = superPhaser(numInfMarkers, markers))[0] > 0){
					//	n += l[0];
					//}
					//ret += l[1];
					double ll = superPhaser(markers_);
					//iterate
					if (phasingIterations > 1) {
						float phap[][][] = getPhase(markers_);
						//System.err.println(ll);
						for (int i = 1; i < phasingIterations; ++i) {
							randomPhase(markers_);
							double ll2 = superPhaser(markers_);
							//System.err.println(ll2);
							if (ll2 > ll) {
								ll = ll2;
								phap = getPhase(markers_);
							}
						}
						setPhase(markers_, phap);						
					}
					ret += ll;
				}
			}
			//System.err.println("running superPhaser!\n" + n + " changes");
			return ret;
		}
		
		private float[][][] getPhase(ArrayList<Integer> markers)
		{
			int numMarkers = markers.size();
            
			float[][][] hap = originalHaplotypes.get(fam);
            float[][][] ret = new float[numMarkers][numIndividuals][4];
            for (int mi = 0; mi < numMarkers; ++mi) {
				int m = markers.get(mi);
				for (int i = 0; i < numIndividuals; ++i) {
					for (int j = 0; j < 4; ++j)
						ret[mi][i][j] = hap[m][i][j];
				}
            }
            
            return ret;
		}
		
		private void setPhase(ArrayList<Integer> markers, float newhap[][][])
		{
			int numMarkers = markers.size();

			float[][][] hap = originalHaplotypes.get(fam);
            for (int mi = 0; mi < numMarkers; ++mi) {
				int m = markers.get(mi);
				for (int i = 0; i < numIndividuals; ++i) {
					for (int j = 0; j < 4; ++j)
						hap[m][i][j] = newhap[mi][i][j];
				}
            }
		}
		
		private void randomPhase(ArrayList<Integer> markers)
		{
            float[][][] hap = originalHaplotypes.get(fam);

            for (int m : markers) {
            	if (!selfingPhase) {
		            if (Misc.random() < 0.5)
						for (int i = 0; i < numIndividuals; ++i) {
							swap(0, 2, hap[m][i]);
							swap(1, 3, hap[m][i]);
						}
					if (Misc.random() < 0.5)
						for (int i = 0; i < numIndividuals; ++i) {
							swap(0, 1, hap[m][i]);
							swap(2, 3, hap[m][i]);									
						}
            	} else { // do not mix phase with selfingPhase
		            if (Misc.random() < 0.5)
						for (int i = 0; i < numIndividuals; ++i) {
							swap(0, 2, hap[m][i]);
							swap(1, 3, hap[m][i]);
							swap(0, 1, hap[m][i]);
							swap(2, 3, hap[m][i]);									
						}
            	}
            }
		}
		
		
		
		public double hyperPhaser(ArrayList<Integer> markers)
		{
			//totalMarkers = markers.size();
			return hyperPhaser_(1, markers); 
		}

		public double hyperPhaser_(int adder, ArrayList<Integer> markers) 
		{
			int totalMarkers = markers.size();
			if (totalMarkers >= 6 * adder)
				hyperPhaser_(2 * adder, markers);
			ArrayList<Integer> hyperMarkers = new ArrayList<Integer>(); 
			for (int i = 0; i < totalMarkers; i+=adder)
				hyperMarkers.add(markers.get(i));
			dataScale = Math.min(1.0, fullDataScale * totalMarkers / hyperMarkers.size());
			return superPhaser(hyperMarkers);
		}
		public double superPhaser(ArrayList<Integer> markers)
		{
			double ret = 0.0;
			int n = 0;
			double l[];
			while ((l = superPhaser_(markers))[0] > 0){
				n += l[0];
			}
			ret += l[1];
			//System.err.println("running superPhaser!\n" + n + " changes");
			return ret;
		}
		
		//phase data by maximizing its likelihood
		public double[] superPhaser_(ArrayList<Integer> markers)
		{
			double logL = 0.0;

            float[][][] hap = originalHaplotypes.get(fam);
            char[] parentOrder = originalParentOrder.get(fam);
            
            int numPhases = 7;
            if (refineParentOrder)
            	numPhases = 14;
            double pSP[] = new double[numPhases];
            double QPhaseSP[] = new double[numPhases];
            
            int numParentChanges = 0;
			
			// forwardSP and iscaleSP have been initialized already...

			for (int hi = 0; hi < numIndividuals; ++hi) {
				logL += ViterbiForwardIndividual(hi, forwardSP[hi], markers, null);
			}
			//System.err.println("logL = " + logL);
			
			int numMarkers = markers.size();

			initBackward(numMarkers & 1);
			//for (int hi = 0; hi < numIndividuals; ++hi)
			//	for (int b1 = 0; b1 < 2; ++b1)
			//		for (int b2 = 0; b2 < 2; ++b2)
			//			for (int r1 = 0; r1 < 2; ++r1)
			//				for (int r2 = 0; r2 < 2; ++r2)
			//					backwardSP[hi][numMarkers2 & 1][b1][b2][r1][r2] = 0.0;			

			int s1 = 0;
			int s2 = 0;
			int phasesChanged = 0;
			for (int mi = numMarkers - 1; mi >= 0; --mi) {
				int m = markers.get(mi);
				
				Arrays.fill(QPhaseSP, 0.0);
				for (int hi = 0; hi < numIndividuals; ++hi) {
					double fTable[] = forwardSP[hi][mi + 1];
					double bTable[] = backwardSP[hi][(mi + 1) & 1];
					Arrays.fill(pSP, Double.NEGATIVE_INFINITY);
					for (int r1 = 0; r1 < 2; ++r1)
						for (int r2 = 0; r2 < 2; ++r2)
							for (int a1 = 0; a1 < 2; ++a1)
								for (int a2 = 0; a2 < 2; ++a2) {
									int b1 = a1 ^ s1;
									int b2 = a2 ^ s2;
		
									double commonTerm =       fTable[index4(b1,b2,r1,r2)] + bTable[index4(a1,a2,r1,r2)];
									pSP[0] = Math.max(pSP[0], commonTerm); // no change
									if (!selfingPhase) {
										pSP[1] = Math.max(pSP[1], fTable[index4(b1 ^ 1,b2    ,r1,r2)] + bTable[index4(a1,a2,r1,r2)]); // change phase of forward for male
										pSP[2] = Math.max(pSP[2], fTable[index4(b1    ,b2 ^ 1,r1,r2)] + bTable[index4(a1,a2,r1,r2)]); // change phase of forward for female
									}
									pSP[3] = Math.max(pSP[3], fTable[index4(b1 ^ 1,b2 ^ 1,r1,r2)] + bTable[index4(a1,a2,r1,r2)]); // change phase of forward for both
									commonTerm -= emission(hap[m][hi], b1, b2);
									if (!selfingPhase) {
										pSP[4] = Math.max(pSP[4], commonTerm + emission(hap[m][hi], b1 ^ 1, b2    )); // change phase of marker m for male
										pSP[5] = Math.max(pSP[5], commonTerm + emission(hap[m][hi], b1    , b2 ^ 1)); // change phase of marker m for female
									}
									pSP[6] = Math.max(pSP[6], commonTerm + emission(hap[m][hi], b1 ^ 1, b2 ^ 1)); // change phase of marker m for both
									
									if (refineParentOrder && parentOrder[m] == '+') { // change parental genotypes for marker m...
										commonTerm = fTable[index4(b1,b2,r1,r2)] + bTable[index4(a1,a2,r1,r2)] - emission(hap[m][hi], b1, b2);
										pSP[7] = Math.max(pSP[7], commonTerm + emission(hap[m][hi], b2, b1)); // no change, but parents
										if (!selfingPhase) {
											pSP[8] = Math.max(pSP[8], fTable[index4(b1 ^ 1,b2    ,r1,r2)] + bTable[index4(a1,a2,r1,r2)] - emission(hap[m][hi], b1 ^ 1, b2) + emission(hap[m][hi], b2, b1 ^ 1)); // change phase of forward for male
											pSP[9] = Math.max(pSP[9], fTable[index4(b1    ,b2 ^ 1,r1,r2)] + bTable[index4(a1,a2,r1,r2)] - emission(hap[m][hi], b1, b2 ^ 1) + emission(hap[m][hi], b2 ^ 1, b1)) ; // change phase of forward for female
										}
										pSP[10] = Math.max(pSP[10], fTable[index4(b1 ^ 1,b2 ^ 1,r1,r2)] + bTable[index4(a1,a2,r1,r2)]  - emission(hap[m][hi], b1 ^ 1, b2 ^ 1) + emission(hap[m][hi], b2 ^ 1, b1 ^ 1)); // change phase of forward for both
										if (!selfingPhase) {
											pSP[11] = Math.max(pSP[11], commonTerm + emission(hap[m][hi], b2, b1 ^ 1)); // change phase of marker m for male
											pSP[12] = Math.max(pSP[12], commonTerm + emission(hap[m][hi], b2 ^ 1, b1)); // change phase of marker m for female
										}
										pSP[13] = Math.max(pSP[13], commonTerm + emission(hap[m][hi], b2 ^ 1, b1 ^ 1)); // change phase of marker m for both
									}

									//double em = emission(hap[m][hi], b1 ^ 1, b2 ^ 1);
									//pSP[7] = Math.max(pSP[7], fTable[b1 ^ 1][b2    ][r1][r2] + bTable[a1][a2][r1][r2] - emission(hap[m][hi], b1 ^ 1, b2) + em); // change forward for male and female for marker m 
									//pSP[8] = Math.max(pSP[8], fTable[b1    ][b2 ^ 1][r1][r2] + bTable[a1][a2][r1][r2] - emission(hap[m][hi], b1, b2 ^ 1) + em); // change forward for female and male for marker m
									//pSP[9] = Math.max(pSP[9], fTable[b1 ^ 1][b2    ][r1][r2] + bTable[a1][a2][r1][r2] - emission(hap[m][hi], b1 ^ 1, b2) + emission(hap[m][hi], b1, b2 ^ 1)); // change forward for male and female for marker m 
									//pSP[10] = Math.max(pSP[10], fTable[b1    ][b2 ^ 1][r1][r2] + bTable[a1][a2][r1][r2] - emission(hap[m][hi], b1, b2 ^ 1)  + emission(hap[m][hi], b1 ^ 1, b2)); // change forward for female and male for marker m
								}
					
						for (int d=0; d < numPhases; ++d)
							QPhaseSP[d] += pSP[d];
					
				}
				
				//System.err.println(QPhaseSP[0] + "\t" + QPhaseSP[1] + "\t"  + QPhaseSP[2] + "\t"  + QPhaseSP[3] + "\t"  + QPhaseSP[4] + "\t"  + QPhaseSP[5] + "\t"  + QPhaseSP[6]);
				int e1 = 0;
				int e2 = 0;
				
				int max = 0;
				for (int j = 1; j < numPhases; ++j)
					if (QPhaseSP[j] > QPhaseSP[max] + 1e-9)
						max = j;
				
				if (max >= 7) { //refineParentOrder
					++numParentChanges;
					for (int i = 0; i < numIndividuals; ++i)
						swap(1, 2, hap[m][i]); // swap mother and father...
					max -= 7; // continue with phase
				}
				
				if (max < 4) {
					s1 = s1 ^ (max & 1);
					s2 = s2 ^ (max >> 1);
				}
				else if (max < 7){ 
					e1 = (max - 3) & 1;
					e2 = (max - 3) >> 1;
				} else { // unreachable code...
					e1 = (max & 1) ^ 1;
					e2 = (max & 1);
					s1 = s1 ^ e2;
					s2 = s2 ^ e1;
				}
				//System.err.println("" + s1 + s2);
				
				
				if ((s1 ^ e1) == 1 || (s2 ^ e2) == 1) {
					++phasesChanged;
					if ((s1 ^ e1) == 1)
						for (int i = 0; i < numIndividuals; ++i) {
							swap(0, 2, hap[m][i]);
							swap(1, 3, hap[m][i]);
						}
					if ((s2 ^ e2) == 1) {
						for (int i = 0; i < numIndividuals; ++i) {
							swap(0, 1, hap[m][i]);
							swap(2, 3, hap[m][i]);									
						}
					}
				}
				if (max > 0) { 
					//System.err.print(max + "\t");
					//System.err.print((QPhaseSP[max] -  QPhaseSP[0]) + "\t");
					//System.err.println(QPhaseSP[0] + "\t" + QPhaseSP[1] + "\t"  + QPhaseSP[2] + "\t"  + QPhaseSP[3] + "\t"  + QPhaseSP[4] + "\t"  + QPhaseSP[5] + "\t"  + QPhaseSP[6]);
					//System.err.println(QPhaseSP[0] + "\t" + QPhaseSP[1] + "\t"  + QPhaseSP[2] + "\t"  + QPhaseSP[3] + "\t"  + QPhaseSP[4] + "\t"  + QPhaseSP[5] + "\t"  + QPhaseSP[6] + "\t"  + QPhaseSP[7] + "\t"  + QPhaseSP[8]);
					
					//scale = Math.exp((QPhaseSP[0] - QPhaseSP[max]) * onePerNumIndividuals); //Math.pow(QPhase[0] / QPhase[max], 1.0 / numIndividuals);
					//System.err.println(Math.log10(QPhaseSP[max] / QPhaseSP[0]) + "\t" + max + "\t" + mi);
				}

				updateBackward(mi, 1, markers);
			}
			//if (numParentChanges > 0) {
			//	System.err.println("Changing " + numParentChanges + " parents");
			//}
			//if (phasesChanged > 0)
			//	System.err.println("Family " + fam + " changed " + phasesChanged + " marker phases");
			//System.err.println("logL = " + logL);
			return new double[]{phasesChanged + numParentChanges, logL};
		}
		
		public void computeTables(ArrayList<Integer> markers_)
		{
			
            int fi = 0; //index for each individual
            for (int i = 0; i < fam; ++i)
                fi += originalHaplotypes.get(i)[0].length;

            int numMarkers = markers_.size();
			
			recScale.init(markers_);
			
			// forwardSP and iscaleSP have been initialized already...

			for (int hi = 0; hi < numIndividuals; ++hi) {
				logLs[hi] = ViterbiForwardIndividual(hi, forwardSP[hi], markers_, null);
			}
			
			
			
			//System.err.println(iscaleSP[0][numMarkers]);

			//for (int hi = 0; hi < numIndividuals; ++hi)
			//	for (int b1 = 0; b1 < 2; ++b1)
			//		for (int b2 = 0; b2 < 2; ++b2)
			//			for (int r1 = 0; r1 < 2; ++r1)
			//				for (int r2 = 0; r2 < 2; ++r2)
			//					backwardSP[hi][numMarkers & 1][b1][b2][r1][r2] = 0.0;
			initBackward(numMarkers & 1);
			

			for (int m = numMarkers - 1; m >= -1; --m) {
				double _tmpRec[][][] = recScale.getRecMatrix(m + 1); // TODO: check
				
				for (int hi = 0; hi < numIndividuals; ++hi) {
					double fTable[] = forwardSP[hi][m + 1];
					double bTable[] = backwardSP[hi][(m + 1) & 1];

					float qt[][] = qTables[fi + hi];
					Arrays.fill(qt[m + 1], Float.NEGATIVE_INFINITY);
					
					for (int b1 = 0; b1 < 2; ++b1)
						for (int b2 = 0; b2 < 2; ++b2)
							for (int r1 = 0; r1 < 2; ++r1)
								for (int r2 = 0; r2 < 2; ++r2) {
									double f[] = {_tmpRec[r1][r2][0] + fTable[index4(b1     ,b2     ,r1    ,r2    )] + bTable[index4(b1,b2,r1,r2)],
							                      _tmpRec[r1][r2][1] + fTable[index4(b1 ^ r1,b2     ,r1 ^ 1,r2    )] + bTable[index4(b1,b2,r1,r2)],
									              _tmpRec[r1][r2][2] + fTable[index4(b1     ,b2 ^ r2,r1    ,r2 ^ 1)] + bTable[index4(b1,b2,r1,r2)],   
									              _tmpRec[r1][r2][3] + fTable[index4(b1 ^ r1,b2 ^ r2,r1 ^ 1,r2 ^ 1)] + bTable[index4(b1,b2,r1,r2)]};
									int max = 0;
									for (int p = 1; p < 4; ++p)
										if (f[p] > f[max])
											max = p;
									qt[m + 1][2 * b1 + b2] = Math.max(qt[m + 1][2 * b1 + b2], (float) (f[max] - logLs[hi]));
								}
				}
				if (m >= 0)
					updateBackward(m, 1, markers_);
			}
		}

		//init backward and forward tables...
		//to evaluate local changes fast
		private void calcForward(ArrayList<Integer> markers){

			recScale.init(markers);			
			// forwardSP and iscaleSP have been initialized already...

			for (int hi = 0; hi < numIndividuals; ++hi) {
				ViterbiForwardIndividual(hi, forwardSP[hi], markers, null);
			}
		}

        private final void initBackward(int position){
            for (int hi = 0; hi < numIndividuals; ++hi)
                for (int b1 = 0; b1 < 2; ++b1)
                    for (int b2 = 0; b2 < 2; ++b2)
                        for (int r1 = 0; r1 < 2; ++r1)
                            for (int r2 = 0; r2 < 2; ++r2)
                                backwardSP[hi][position][index4(b1,b2,r1,r2)] = 0.0;                 
        }

		//	s1 and s2 are the two states
		private double emission(float posterior[], int s1, int s2)
		{
			//return dataScale * posterior[2 * s1 + s2];
			return posterior[2 * s1 + s2];
		}

		//TODO: finish RecombinationScale
		private abstract class RecombinationScale{
			public abstract void update(int perm[], int start); 
			public abstract void init(ArrayList<Integer> markers);
			public abstract double[][][] getRecMatrix(int position);
		}
		private class RecombinationScale1 extends RecombinationScale{
			private int inf1[] = null;
			private int inf2[] = null;
			private int infLimitLow = 0;
			
			private double scalePerMax = 1.0;
			private double tmpRec2[][][] = new double[2][2][4];
			
			private double prevScale1 = 1.0;
			private double prevScale2 = 1.0;
			
			public void init(ArrayList<Integer> markers) {
				scalePerMax = 1.0 / MAX_SCALE;
				infLimitLow = (scalePerMax >= 1.0 / dataScale) ? 0 : (int) (1.0 / (dataScale * scalePerMax));
	
				if (infLimitLow == 0)
					return;
				int numMarkers = markers.size();
				
				inf1 = new int[numMarkers + 1];
				inf2 = new int[numMarkers + 1];
				int inf[] = originalInformative.get(fam);//originalInformative.get(fam);
				int numInf1 = 0;
				int numInf2 = 0;
				for (int mi = 0; mi < numMarkers; ++mi) {
					int m = markers.get(mi);
					if ((inf[m] & 1) == 1)
						++numInf1;
					if ((inf[m] & 2) == 2)
						++numInf2;
					inf1[mi + 1] = numInf1;
					inf2[mi + 1] = numInf2;
				}
				
			}
			public void update(int perm[], int start) { // update for polishScore
				if (infLimitLow == 0)
					return;
				int k = perm.length;
				int inf[] = originalInformative.get(fam);//originalInformative.get(fam);

				int numInf1 = inf1[start]; 
				int numInf2 = inf2[start]; 

				for (int mi = 0; mi < k; ++mi) {
					int m = perm[mi];
					if ((inf[m] & 1) == 1)
						++numInf1;
					if ((inf[m] & 2) == 2)
						++numInf2;
					inf1[start + mi + 1] = numInf1;
					inf2[start + mi + 1] = numInf2;
				}
				
			}
			private void initMatrix(double scale1, double scale2) {
				if (infLimitLow == 0 || (prevScale1 == scale1 && prevScale2 == scale2))
					return;
				double logI1 = scale1 * __tmpRec[0][0][1];
				double logI2 = scale2 * __tmpRec[0][0][2];
				double logR1 = scale1 * __tmpRec[1][1][1];
				double logR2 = scale2 * __tmpRec[1][1][2];

//				tmpRec[0][0] = new double[]{1.0, i1, i2,  i1 * i2};
//				tmpRec[0][1] = new double[]{1.0, i1, r2,  i1 * r2};
//				tmpRec[1][0] = new double[]{1.0, r1, i2,  r1 * i2};
//				tmpRec[1][1] = new double[]{1.0, r1, r2,  r1 * r2}; //phew...

				tmpRec2[0][0][1] = logI1;
				tmpRec2[0][0][2] = logI2;
				tmpRec2[0][0][3] = logI1 + logI2;
				
				tmpRec2[0][1][1] = logI1;
				tmpRec2[0][1][2] = logR2;
				tmpRec2[0][1][3] = logI1 + logR2;

				tmpRec2[1][0][1] = logR1;
				tmpRec2[1][0][2] = logI2;
				tmpRec2[1][0][3] = logR1 + logI2;
	
				tmpRec2[1][1][1] = logR1;
				tmpRec2[1][1][2] = logR2;
				tmpRec2[1][1][3] = logR1 + logR2;
				prevScale1 = scale1;
				prevScale2 = scale2;
			}
			
			public double[][][] getRecMatrix(int position)
			{
				if (infLimitLow == 0)
					return __tmpRec;
			
				double scale1 = 1.0 / dataScale;
				double scale2 = scale1;
				double cap = Math.min(CAP_SCALE, scale1);
				
				int infTotal1 = inf1[inf1.length - 1];
				int infTotal2 = inf2[inf2.length - 1];
				int numInf1 = inf1[position];
				int numInf2 = inf2[position];
				
				if (numInf1 < infLimitLow || infTotal1 - numInf1 < infLimitLow) {
					int k = Math.min(numInf1, infTotal1 - numInf1);
					scale1 = (k <= 1) ? scalePerMax : k * scalePerMax;
				}
				if (numInf2 < infLimitLow || infTotal2 - numInf2 < infLimitLow) {
					int k = Math.min(numInf2, infTotal2 - numInf2);
					scale2 = (k <= 1) ? scalePerMax : k * scalePerMax;
				}
				//System.err.println("getRecMatrix " + position + " scale " + scale1 + " " + scale2 + " informative " + (inf1[position + 1] - inf1[position] + 2 *(inf2[position + 1] - inf2[position])) );

				scale1 = Math.max(scale1, cap);
				scale2 = Math.max(scale2, cap);
				
				if (scale1 == 1.0 && scale2 == 1.0) {
					//System.err.println(scale1 + "\t" + scale2);// + "\t" + informative.get(fam)[position]);
					return __tmpRec;					
				}
				else {
					//System.err.println(scale1 + "\t" + scale2 + "\t" + informative.get(fam)[position]);
					initMatrix(scale1, scale2);
					return tmpRec2;
				}
			}
		}
		private class RecombinationScale2 extends RecombinationScale{
			private double inf1[] = null;
			private double inf2[] = null;
			private int infLimitLow = 0;
			
			private double scalePerMax = 1.0;
			private double tmpRec2[][][] = new double[2][2][4];
			
			private double prevScale1 = 1.0;
			private double prevScale2 = 1.0;
			
			private double scale1 = 1.0;
			private double scale2 = 1.0;
			
			public void init(ArrayList<Integer> markers) {
				scalePerMax = 1.0 / MAX_SCALE;
				infLimitLow = (scalePerMax >= 1.0 / dataScale) ? 0 : (int) (1.0 / (dataScale * scalePerMax));
				//infLimitLow = (1.0 / (dataScale * scalePerMax));
	
				if (infLimitLow == 0)
					return;
				int numMarkers = markers.size();
				
				inf1 = new double[numMarkers + 1];
				inf2 = new double[numMarkers + 1];
				int inf[] = originalInformative.get(fam);//originalInformative.get(fam);
				float info[] = originalInformation.get(fam);//originalInformative.get(fam);
				double numInf1 = 0.0;
				double numInf2 = 0.0;
				for (int mi = 0; mi < numMarkers; ++mi) {
					int m = markers.get(mi);
					if (inf[m] > 0) {
						if (inf[m] == 3) {
							numInf1 += 0.5 * info[m];
							numInf2 += 0.5 * info[m];
						}
						else if (inf[m] == 1) {
							numInf1 += info[m];
						}
						else if (inf[m] == 2) {
							numInf2 += info[m];
						}
					}
					inf1[mi + 1] = numInf1;
					inf2[mi + 1] = numInf2;
				}
				
				//for (int mi = 0; mi < numMarkers; ++mi) {
				//	getRecMatrix(mi);
				//	System.err.println(mi + "\t" + scale1 + "\t" + scale2);
				//}
				
				
			}
			public void update(int perm[], int start) { // update for polishScore
				if (infLimitLow == 0)
					return;
				int k = perm.length;
				int inf[] = originalInformative.get(fam);//originalInformative.get(fam);
				float info[] = originalInformation.get(fam);

				double numInf1 = inf1[start]; 
				double numInf2 = inf2[start]; 

				for (int mi = 0; mi < k; ++mi) {
					int m = perm[mi];
					if (inf[m] > 0) {
						if (inf[m] == 3) {
							numInf1 += 0.5 * info[m];
							numInf2 += 0.5 * info[m];
						}
						else if (inf[m] == 1) {
							numInf1 += info[m];
						}
						else if (inf[m] == 2) {
							numInf2 += info[m];
						}
					}
					inf1[start + mi + 1] = numInf1;
					inf2[start + mi + 1] = numInf2;
				}
				
			}
			private void initMatrix() {
				if (infLimitLow == 0 || (prevScale1 == scale1 && prevScale2 == scale2))
					return;
				double logI1 = scale1 * __tmpRec[0][0][1];
				double logI2 = scale2 * __tmpRec[0][0][2];
				double logR1 = scale1 * __tmpRec[1][1][1];
				double logR2 = scale2 * __tmpRec[1][1][2];

//				tmpRec[0][0] = new double[]{1.0, i1, i2,  i1 * i2};
//				tmpRec[0][1] = new double[]{1.0, i1, r2,  i1 * r2};
//				tmpRec[1][0] = new double[]{1.0, r1, i2,  r1 * i2};
//				tmpRec[1][1] = new double[]{1.0, r1, r2,  r1 * r2}; //phew...

				tmpRec2[0][0][1] = logI1;
				tmpRec2[0][0][2] = logI2;
				tmpRec2[0][0][3] = logI1 + logI2;
				
				tmpRec2[0][1][1] = logI1;
				tmpRec2[0][1][2] = logR2;
				tmpRec2[0][1][3] = logI1 + logR2;

				tmpRec2[1][0][1] = logR1;
				tmpRec2[1][0][2] = logI2;
				tmpRec2[1][0][3] = logR1 + logI2;
	
				tmpRec2[1][1][1] = logR1;
				tmpRec2[1][1][2] = logR2;
				tmpRec2[1][1][3] = logR1 + logR2;
				prevScale1 = scale1;
				prevScale2 = scale2;
			}
			
			public double[][][] getRecMatrix(int position)
			{
				if (infLimitLow == 0)
					return __tmpRec;
			
				scale1 = 1.0 / dataScale;
				scale2 = scale1;
				double cap = Math.min(CAP_SCALE, scale1);
				
				double infTotal1 = inf1[inf1.length - 1];
				double infTotal2 = inf2[inf2.length - 1];
				double numInf1 = inf1[position];
				double numInf2 = inf2[position];
				
				if (numInf1 < infLimitLow || infTotal1 - numInf1 < infLimitLow) {
					double k = Math.min(numInf1, infTotal1 - numInf1);
					scale1 = (k <= 1) ? scalePerMax : k * scalePerMax;
				}
				if (numInf2 < infLimitLow || infTotal2 - numInf2 < infLimitLow) {
					double k = Math.min(numInf2, infTotal2 - numInf2);
					scale2 = (k <= 1) ? scalePerMax : k * scalePerMax;
				}
				//System.err.println("getRecMatrix " + position + " scale " + scale1 + " " + scale2 + " informative " + (inf1[position + 1] - inf1[position] + 2 *(inf2[position + 1] - inf2[position])) );

				scale1 = Math.max(scale1, cap);
				scale2 = Math.max(scale2, cap);
				
				if (scale1 == 1.0 && scale2 == 1.0) {
					//System.err.println(scale1 + "\t" + scale2);// + "\t" + informative.get(fam)[position]);
					return __tmpRec;					
				}
				else {
					//System.err.println(scale1 + "\t" + scale2 + "\t" + informative.get(fam)[position]);
					initMatrix();
					return tmpRec2;
				}
			}
		}
		

		//forward Viterbi calculations for single individual
		private double ViterbiForwardIndividual(int hi, double forward[][], ArrayList<Integer> markers, int path[][]){
			return ViterbiForwardIndividual(hi, forward, markers, path, false);
		}
		
		final static double VITERBI_SMALL = 1e-9;
		
		private double ViterbiForwardIndividual(int hi, double forward[][], ArrayList<Integer> markers, int path[][], boolean allSolutions){
			float[][][] hap = originalHaplotypes.get(fam);
			int numMarkers = markers.size();
			
			for (int mi = 0; mi < numMarkers; ++mi) {
				int m = markers.get(mi);
				double _tmpRec[][][] = recScale.getRecMatrix(mi);
				
				double fTable[] = forward[mi]; 
				double fTable2[] = forward[mi + 1]; 
				for (int b1 = 0; b1 < 2; ++b1)
					for (int b2 = 0; b2 < 2; ++b2) {
						double emis = emission(hap[m][hi], b1, b2); 
						//double fTable22[][] = fTable2[b1][b2];
						for (int r1 = 0; r1 < 2; ++r1) {
							for (int r2 = 0; r2 < 2; ++r2) {
								double f[] = {_tmpRec[r1][r2][0] + fTable[index4(b1     ,b2     ,r1    ,r2    )],
					                          _tmpRec[r1][r2][1] + fTable[index4(b1 ^ r1,b2     ,r1 ^ 1,r2    )],
								              _tmpRec[r1][r2][2] + fTable[index4(b1     ,b2 ^ r2,r1    ,r2 ^ 1)],  
								              _tmpRec[r1][r2][3] + fTable[index4(b1 ^ r1,b2 ^ r2,r1 ^ 1,r2 ^ 1)]};
								int max = 0;
								for (int p = 1; p < 4; ++p)
									if (f[p] > f[max])
										max = p;
								
								double fMax = f[max] + emis;
								//fTable22[r1][r2] = fMax;
								fTable2[index4(b1, b2, r1, r2)] = fMax;
								if (path != null)
									switch (max) {
										case 0: path[m + 1][index4(b1,b2,r1,r2)] = b1        + 2 * b2           + 4 * r1       + 8 * r2;break;
										case 1: path[m + 1][index4(b1,b2,r1,r2)] = (b1 ^ r1) + 2 * b2           + 4 * (r1 ^ 1) + 8 * r2;break;
										case 2: path[m + 1][index4(b1,b2,r1,r2)] = b1        + 2 * (b2 ^ r2)    + 4 * r1       + 8 * (r2 ^ 1);break;
										case 3: path[m + 1][index4(b1,b2,r1,r2)] = (b1 ^ r1) + 2 * (b2 ^ r2)    + 4 * (r1 ^ 1) + 8 * (r2 ^ 1);break; 
									}
								if (allSolutions) {
									int states = 0;
									for (int p = 0; p < 4; ++p)
										if (f[p] >= f[max] - VITERBI_SMALL)
											switch (p) {
												case 0: 
													states += 16 << (b1        + 2 * b2           + 4 * r1       + 8 * r2);
													break;
												case 1: 
													states += 16 << ((b1 ^ r1) + 2 * b2           + 4 * (r1 ^ 1) + 8 * r2);
													break;
												case 2: 
													states += 16 << (b1        + 2 * (b2 ^ r2)    + 4 * r1       + 8 * (r2 ^ 1));
													break;
												case 3: 
													states += 16 << ((b1 ^ r1) + 2 * (b2 ^ r2)    + 4 * (r1 ^ 1) + 8 * (r2 ^ 1));
													break;
											}
									path[m + 1][index4(b1,b2,r1,r2)] += states;
								} 								
							}
						}
					}

			}
			double logL = Double.NEGATIVE_INFINITY;
			double fTable[] = forward[numMarkers]; 
			for (int b1 = 0; b1 < 2; ++b1)
				for (int b2 = 0; b2 < 2; ++b2)
					for (int r1 = 0; r1 < 2; ++r1)
						for (int r2 = 0; r2 < 2; ++r2)
							logL = Math.max(logL, fTable[index4(b1,b2,r1,r2)]);
			return logL;
		}
		
		private void ViterbiForwardIndividualFaster(int hi, double forward[][], int numMarkers2, int markers[], int start){ // for polishScore
			float[][][] hap = originalHaplotypes.get(fam);
			
			double f[] = new double[5]; // save new...?
			for (int mi = 0; mi < numMarkers2; ++mi) {
				double _tmpRec[][][] = recScale.getRecMatrix(mi + start);
				
				int m = markers[mi];

				double fTable[] = forward[mi];
				for (int b1 = 0; b1 < 2; ++b1)
					for (int b2 = 0; b2 < 2; ++b2) {
						//double fTable2[][] = forward[mi + 1][b1][b2];
						double fTable2[] = forward[mi + 1];
						double emis = emission(hap[m][hi], b1, b2); 
						for (int r1 = 0; r1 < 2; ++r1) {
							for (int r2 = 0; r2 < 2; ++r2) {
								f[0] = _tmpRec[r1][r2][0] + fTable[index4(b1     ,b2     ,r1    ,r2    )];
								f[1] = _tmpRec[r1][r2][1] + fTable[index4(b1 ^ r1,b2     ,r1 ^ 1,r2    )];
								f[2] = _tmpRec[r1][r2][2] + fTable[index4(b1     ,b2 ^ r2,r1    ,r2 ^ 1)];  
								f[3] = _tmpRec[r1][r2][3] + fTable[index4(b1 ^ r1,b2 ^ r2,r1 ^ 1,r2 ^ 1)];
								
								int max = 0;
								for (int p = 1; p < 4; ++p)
									if (f[p] > f[max])
										max = p;
								fTable2[index4(b1, b2, r1, r2)] = f[max] + emis;
							}
						}
					}
			}
		}
		
		
        private final void updateBackward(int mi, int mask, ArrayList<Integer> markers){
            float hap[][][] = originalHaplotypes.get(fam);
            int m = markers.get(mi);
			double _tmpRec[][][] = recScale.getRecMatrix(mi); // TODO: check if + 1 is needed
            
            for (int hi = 0; hi < numIndividuals; ++hi) {

                    double bTable[] = backwardSP[hi][(mi + 1) & mask];
                    double bTable2[] = backwardSP[hi][mi &  mask];

                    for (int b1 = 0; b1 < 2; ++b1)
                        for (int b2 = 0; b2 < 2; ++b2) {
                            double emis = emission(hap[m][hi], b1, b2);
                                for (int r1 = 0; r1 < 2; ++r1)
                                    for (int r2 = 0; r2 < 2; ++r2)
                                    tmpB[index4(b1,b2,r1,r2)] =  bTable[index4(b1,b2,r1,r2)] + emis;
                        }
                    
                    
                    for (int r1 = 0; r1 < 2; ++r1)
                        for (int r2 = 0; r2 < 2; ++r2) {
                            for (int a1 = 0; a1 < 2; ++a1)
                                for (int a2 = 0; a2 < 2; ++a2) {
                                    double f[] =   {_tmpRec[r1    ][r2    ][0]  + tmpB[index4(a1         ,a2         ,r1    ,r2    )],
                            						_tmpRec[r1 ^ 1][r2    ][1]  + tmpB[index4(a1 ^ r1 ^ 1,a2         ,r1 ^ 1,r2    )],
                            						_tmpRec[r1    ][r2 ^ 1][2]  + tmpB[index4(a1         ,a2 ^ r2 ^ 1,r1    ,r2 ^ 1)],  
                            						_tmpRec[r1 ^ 1][r2 ^ 1][3]  + tmpB[index4(a1 ^ r1 ^ 1,a2 ^ r2 ^ 1,r1 ^ 1,r2 ^ 1)]};
                                                    
                                    int max = 0;
                                    for (int p = 1; p < 4; ++p)
                                    	if (f[p] > f[max])
                                    		max = p;
                            		bTable2[index4(a1,a2,r1,r2)] = f[max];
                                }
                        }
            }
        }

	    private double polishScore(int perm[], int pos) {
	
	        int k = perm.length;
	        
	        double logL = 0.0;
	        for (int hi = 0; hi < numIndividuals; ++hi) {
	            for (int b1 = 0; b1 < 2; ++b1)
	                for (int b2 = 0; b2 < 2; ++b2)
	                    for (int r1 = 0; r1 < 2; ++r1)
	                        for (int r2 = 0; r2 < 2; ++r2)
	                            forwardV2[0][index4(b1,b2,r1,r2)] = forwardSP[hi][pos][index4(b1,b2,r1,r2)];

	        recScale.update(perm, pos);
	        ViterbiForwardIndividualFaster(hi, forwardV2, k, perm, pos);
	
	        double max = Double.NEGATIVE_INFINITY;
	        for (int b1 = 0; b1 < 2; ++b1)
	            for (int b2 = 0; b2 < 2; ++b2)
	                for (int r1 = 0; r1 < 2; ++r1)
	                    for (int r2 = 0; r2 < 2; ++r2)
	                        max = Math.max(max, forwardV2[k][index4(b1,b2,r1,r2)] + backwardSP[hi][(pos + k) & 1][index4(b1,b2,r1,r2)]);
	            //System.err.
	            
	        	logL += max;
	        }
	        return logL;
	    }
		
		public double Viterbi2Fast(ArrayList<Integer> markers) {
			double logL = 0.0;
			
		 	//forwardV and tmpRec are initialized 
			
			for (int hi = 0; hi < numIndividuals; ++hi) {
				logL += ViterbiForwardIndividual(hi, forwardV, markers, null);
			}
			return logL;
		}
		
		
		private int[] getBits16(int number){
			int n = 0;
			for (int bit = 0; bit < 16; ++bit)
				if ((number & (1 << bit)) != 0)
					++n;

			int ret[] = new int[n];
			n = 0;
			for (int bit = 0; bit < 16; ++bit)
				if ((number & (1 << bit)) != 0)
					ret[n++] = bit;

			return ret;
		}		
		
		public void Viterbi3(ArrayList<Integer> markers_, ArrayList<float[][]> quality_ret)
		{
            float quality[][] = quality_ret.get(fam);
    		
            int numMarkers = markers_.size();
			
            double logL = 0.0;
			for (int hi = 0; hi < numIndividuals; ++hi) {
				logL += ViterbiForwardIndividual(hi, forwardSP[hi], markers_, null);
			}
						
    		initBackward(numMarkers & 1);   
			double q[][] = new double[2][2];
			
			for (int m = numMarkers - 1; m >= 0; --m) {
				for (int hi = 0; hi < numIndividuals; ++hi) {
					double fTable[] = forwardSP[hi][m + 1];
					double bTable[] = backwardSP[hi][(m + 1) & 1];
					q[0][0] = Double.NEGATIVE_INFINITY;
					q[0][1] = Double.NEGATIVE_INFINITY;
					q[1][0] = Double.NEGATIVE_INFINITY;
					q[1][1] = Double.NEGATIVE_INFINITY;
					
					for (int r1 = 0; r1 < 2; ++r1)
						for (int r2 = 0; r2 < 2; ++r2)
							for (int a1 = 0; a1 < 2; ++a1)
								for (int a2 = 0; a2 < 2; ++a2) {
									q[a1][a2] = Misc.logSum(q[a1][a2], fTable[index4(a1,a2,r1,r2)] + bTable[index4(a1,a2,r1,r2)]);
								}
					double max = Double.NEGATIVE_INFINITY;
					for (int a1 = 0; a1 < 2; ++a1)
						for (int a2 = 0; a2 < 2; ++a2)
							max = Math.max(max, q[a1][a2]);

					quality[m][4 * hi + 0] = (float)Math.exp(q[0][0] - max);
					quality[m][4 * hi + 1] = (float)Math.exp(q[0][1] - max);
					quality[m][4 * hi + 2] = (float)Math.exp(q[1][0] - max);
					quality[m][4 * hi + 3] = (float)Math.exp(q[1][1] - max);
				}
				updateBackward(m, 1, markers_);
			}
		}
		
		public double[] Viterbi2(ArrayList<Integer> markers_, ArrayList<int[][]> ret, ArrayList<int[]> ret2, boolean silent, boolean maskHaplotypes)
		{
			
			double logL = 0.0;
			int count = 0;
			
			int numMarkers = markers_.size();
			
			//forwardV already initialized
			//tmpRec already initialized

			recScale.init(markers_);
			
            float[][][] hap = originalHaplotypes.get(fam);
            int [][] retF = ret.get(fam);
            int [] retI = ret2.get(fam); //individual recombination rate
			
			
			for (int hi = 0; hi < numIndividuals; ++hi) {
				logL += ViterbiForwardIndividual(hi, forwardV, markers_, path, maskHaplotypes);

				double max = Double.NEGATIVE_INFINITY;
				int maxB1 = 0;
				int maxB2 = 0;
				int maxR1 = 0;
				int maxR2 = 0;

				for (int r1 = 0; r1 < 2; ++r1)
					for (int r2 = 0; r2 < 2; ++r2)
						for (int b1 = 0; b1 < 2; ++b1)
							for (int b2 = 0; b2 < 2; ++b2) {
								double f = forwardV[numMarkers][index4(b1,b2,r1,r2)]; 
								if (f > max) {
									max = f;
									maxB1 = b1;
									maxB2 = b2;
									maxR1 = r1;
									maxR2 = r2;
								}
							}
				
				
				for (int m = numMarkers - 1; m >= 0; --m) { //backtrack most likely path
					retF[m][hi] = maxB1; 
					retF[m][hi + hap[0].length] = maxB2;
					
					int p = path[m + 1][index4(maxB1,maxB2,maxR1,maxR2)];
					int tmpB1 = (p & 1);
					int tmpB2 = (p & 2) >> 1;
					int tmpR1 = (p & 4) >> 2;
					int tmpR2 = (p & 8) >> 3;
						
					if (tmpB1 != maxB1) {
						++count;
						++retF[m][2 * hap[0].length]; 
						++retI[hi]; 
					}
					if (tmpB2 != maxB2) {
						++count;
						++retF[m][2 * hap[0].length + 1]; 
						++retI[hi]; 
					}
					maxB1 = tmpB1;
					maxB2 = tmpB2;
					maxR1 = tmpR1;
					maxR2 = tmpR2;
				}
				//mask uncertain haplotypes
				if (maskHaplotypes) {
					int states = 0; //set of possible states
					for (int r1 = 0; r1 < 2; ++r1)
						for (int r2 = 0; r2 < 2; ++r2)
							for (int b1 = 0; b1 < 2; ++b1)
								for (int b2 = 0; b2 < 2; ++b2)
									if (forwardV[numMarkers][index4(b1,b2,r1,r2)] >= max - VITERBI_SMALL)
										states += 1 << (b1 + 2 * b2 + 4 * r1 + 8 * r2);
					
					//ArrayList<Integer> tmp = new ArrayList<Integer>();
					//for (int s : getBits16(states))
					//	tmp.add(s);
					//System.err.println(states + " " + tmp);
					
					
					for (int m = numMarkers - 1; m >= 0; --m) { //backtrack most likely paths
						int nextStates = 0;
						int bits[] = getBits16(states);
						int b1Mask = 0;
						int b2Mask = 0;
						for (int p: bits) {
							int b1 = (p & 1);
							int b2 = (p & 2) >> 1;
							int r1 = (p & 4) >> 2;
							int r2 = (p & 8) >> 3;
							b1Mask |= (1 << b1);
							b2Mask |= (1 << b2);
							nextStates |= (path[m + 1][index4(b1,b2,r1,r2)] >> 4);
						}
						states = nextStates;
						if (b1Mask == 3)
							retF[m][hi] = 2; 
						if (b2Mask == 3)
							retF[m][hi + hap[0].length] = 2;
					}
					
				}
				
			}
			
			if (!silent)
				System.err.println("logL = " + logL);
			return new double[]{count, logL};
		}
		
		@Override
		public Double call()
		{
			return 0.0;//phase();
		}
		@Override
		public void run()
		{
			//computeTables();
		}
		
	}
	
	// assumes likelihood(markers, true) has been called
	// and then initData(markers2)
	//TODO: handle phased data (maxPhase should do it)
	//TODO: handle selfing data (phaseAdder should do it)
	public double score(int marker, int position)
	{
		double ret = 0.0;

		int phaseAdder = 1;
		int maxPhase = 4;
		if (phasedData)
			maxPhase = 1;
		else if (selfingPhase) // only 0 and 3 as phase for selfing data
			phaseAdder = 3;
			
		
		int numF = originalHaplotypes.size();
		int fi = 0;
		for (int fam = 0; fam < numF; ++fam) {
			float[][][] hap = originalHaplotypes.get(fam); // TODO: Check ... (haplotypes.get does not work...)
			double retPhase[] = new double[4];
			int numIndividuals = hap[0].length;
			for (int hi = 0; hi < numIndividuals; ++hi) {
				float Q[] = qTables[fi][position];
				++fi;
				for (int phase = 0; phase < maxPhase; phase+=phaseAdder) {
					//double f = Q[phase] + dataScale * hap[marker][hi][0];
					double f = Q[phase] + hap[marker][hi][0];
					for (int haplotype = 1; haplotype < 4; ++haplotype)
						f = Math.max(f, Q[haplotype ^ phase] + hap[marker][hi][haplotype]); 
						//f = Math.max(f, Q[haplotype ^ phase] + dataScale * hap[marker][hi][haplotype]); 

					/* //possible improvement, take the sum of two best scores 
					double f = Q[phase] + hap[marker][hi][0];
					double f2 = Double.NEGATIVE_INFINITY;
					for (int haplotype = 1; haplotype < 4; ++haplotype) {
						double nf = Q[haplotype ^ phase] + hap[marker][hi][haplotype];
						if (nf > f) {
							f2 = f;
							f = nf;
						} else
						 f2 = Math.max(f2, nf);
					}
					//retPhase[phase] = f + Misc.getLookup(f - f2);//Misc.LogSumFast(f,f2)
					*/
					retPhase[phase] += f;
				}
			}
			//System.err.println(retPhase1[0] + "\t" + retPhase1[1]  + "\t" + retPhase1[2]  + "\t" + retPhase1[3]);
			if (phasedData)
				ret += retPhase[0];
			else if (selfingPhase)
				ret += Math.max(retPhase[0], retPhase[3]);
			else
				ret += Math.max(Math.max(Math.max(retPhase[0], retPhase[1]), retPhase[2]), retPhase[3]); 
		}
		//System.err.println("score = " + ret);
		if (pfs != null)
			for (PhysicalFamily pf : pfs)
				ret += pf.score(marker, position);
		return ret;	
	}

	
	private static int perm4[][] = Misc.permutation(4);
	private static int perm4_s[][] = Misc.prunePermutation1(perm4);
	private static int perm5[][] = Misc.permutation(5);
	private static int perm5_s[][] = Misc.prunePermutation1(perm5);

	
	final double IMPROVE_TOLERANCE = 0.001;

	//check if distance between two positions is higher than DISTANCE_LIMIT, handles parameter identicalLimit in OrderMarkers2
	private boolean distanceHigher(int pos1, int pos2)
	{
		double d = 0.0;

		//distance based on physical family 
		if (pfs != null) {
			for (PhysicalFamily pf : pfs) {
				double pd = pf.scoreDistance(pos1, pos2); 
				if (pd != 0)
					d += (1.0 - Math.exp(pd)); // is pd is between 0 and 1?
			}
		}
		if (d > DISTANCE_LIMIT)
			return true;
		
		for (float ind[][] : qTables) {
			double maxD = 0.0;
			for (int j = 0; j < 4; ++j) {
				double i1 = ind[pos1][j];
				double i2 = ind[pos2][j];
				if (i1 != i2)
					maxD = Math.max(maxD,  Math.abs(Math.exp(i1) - Math.exp(i2)));
			}
			d += maxD;
			if (d > DISTANCE_LIMIT)
				return true;
		}
		return false;
	}

	private ArrayList<Integer> getUniquePositions(ArrayList<Integer> markers) {
		int numMarkers = markers.size();
		ArrayList<int[][]> phasedHaplotypes = getPhasedHaplotypes(markers, true, false);

		ArrayList<Integer> positions = new ArrayList<Integer>();
		
		int lastPos = 0;
		while (lastPos < numMarkers) { 
			int nextPos = lastPos + 1;
			int delta = 0;
			while (nextPos < numMarkers) {
				for (int ph[][] : phasedHaplotypes) {
					int phh[] = ph[nextPos]; 
					delta += phh[phh.length - 1] + phh[phh.length - 2];
				}
				if (delta > 0)
					break;
				++nextPos;
			}
			positions.add((nextPos + lastPos + 1) / 2);
			lastPos = nextPos;
		}
		return	positions; 
	}
	
	
	private ArrayList<Integer> getPositions(ArrayList<Integer> markers) {
		return getPositions(markers, null);
	}

	private ArrayList<Integer> getPositions(ArrayList<Integer> markers, ArrayList<ArrayList<Integer>> markersInNearestBin) {
		int n = markers.size();
		ArrayList<Integer> ret = new ArrayList<Integer>();
		
		//int numPos = 1;
		int prev = 0;
		ret.add(0);
		for (int i = 1; i <= n; ++i)
			if (i == 1 || i >= n - 1 || distanceHigher(prev, i)) {
				if (markersInNearestBin != null)
					for (int mi = prev; mi < i; ++mi)
						markersInNearestBin.get(prev + 1).add(markers.get(mi));
				ret.add(i);
				prev = i;
				//++numPos;
			}
		//System.err.println("numPos = " + ret.size() + "/" + (n + 1));
		
		return ret;
	}
	
	
	private class polishThread implements Runnable {
		ArrayList<Integer> positions;
		ArrayList<Integer> markers;
		int m1;
		int m2;
		int ret[];
		
		public polishThread(ArrayList<Integer> positions, ArrayList<Integer> markers, int m1, int m2, int ret[]) {
			this.ret = ret;
			this.m1 = m1;
			this.m2 = m2;
			this.markers = markers;
			this.positions = positions;
		}
		@Override
		public void run() {
			for (int mi = m1; mi < m2; ++mi) {
				int m = markers.get(mi);
				int maxPos = 0;
				double maxScore = Double.NEGATIVE_INFINITY;
				
				for (int p : positions) {
					double s = score(m, p); 
					if (s > maxScore) {
						maxScore = s;
						maxPos = p;
					}
				}
				ret[mi] = maxPos;  
			}
			
		}
	
	}

	private ArrayList<Integer> polishTest(ArrayList<Integer> markers)
	{
		double ll = likelihood(markers, true); // to store tables, used in merge...
				 
		//System.err.println("score = " + ll);
		
		ArrayList<Integer> positions = getPositions(markers);

		int numMarkers = markers.size();
		
		ArrayList<ArrayList<Integer>> order = new ArrayList<ArrayList<Integer>>(); 
		for (int i = 0; i <= numMarkers; ++i)
			order.add(new ArrayList<Integer>());

		int numPositions = positions.size();
		
		int nt = (int) calcNumThreadsPolish(numMarkers, numPositions, qTables.length);

		if (nt <= 1) {  
			for (int mi = 0; mi < numMarkers; ++mi) {
				int m = markers.get(mi);
				int maxPos = 0;
				double maxScore = Double.NEGATIVE_INFINITY;
				for (int p: positions) {
					double s = score(m, p); 
					if (s > maxScore) {
						maxScore = s;
						maxPos = p;
					}
					//System.err.print(s + "\t");
				}
				order.get(maxPos).add(m); 
				//System.err.println(maxPos);
				//System.err.println();
			}
		} else { //parallel implementation
			//System.err.println(nt + " threads and " + numMarkers + " markers and " + numPositions + " positions");
			int ret[] = new int[numMarkers];
			
			Thread threads[] = new Thread[nt]; 
			
			for (int t = 0; t < nt; ++t) {
				threads[t] = new Thread(new polishThread(positions, markers, t * markers.size() / nt, (t + 1) * markers.size() / nt, ret));
				threads[t].start();
			}
			try {
				for (int t = 0; t < nt; ++t)
					threads[t].join();
				
			} catch (Exception e) {
				e.printStackTrace();
				Error.error(-99999999);
				System.exit(-1);
			}
			for (int mi = 0; mi < numMarkers; ++mi) {
				int m = markers.get(mi);
				order.get(ret[mi]).add(m); 
			}
		}
		
		markers.clear();
		for (int i = 0; i <= numMarkers; ++i)
			markers.addAll(order.get(i));
		
		return markers;
	}

	private static int perm_all[][][] = 
		{Misc.prunePermutation2(Misc.permutation(3)), Misc.prunePermutation2(Misc.permutation(4)), Misc.prunePermutation2(Misc.permutation(5))};

	private ArrayList<Integer> bestOrder(ArrayList<Integer> markers) {
		int n = markers.size();
		int perm[][] = null;
		if (n > 5 || n <= 2)
			return markers;
		else {
			perm = perm_all[n - 3];
		}
		
		ArrayList<Integer> markers2 = new ArrayList<Integer>();
		markers2.addAll(markers);

		double maxL = Double.NEGATIVE_INFINITY;
		int maxpi = 0;
		int numSolutions = 0;
		for (int pi = 0; pi < perm.length; ++pi) {
			int p[] = perm[pi];
			for (int i = 0; i < n; ++i)
				markers2.set(i, markers.get(p[i]));
			double ll = score(markers2);
			//System.err.println(ll + "\t" + markers2);
			if (ll > maxL - IMPROVE_TOLERANCE) {
				if (ll > maxL + IMPROVE_TOLERANCE)
					numSolutions = 1;
				else 
					++numSolutions;
				if (Misc.random() <= 1.0 / numSolutions)
					maxpi = pi;
				maxL = Math.max(ll, maxL);
			}
		}
		for (int i = 0; i < n; ++i)
			markers2.set(i, markers.get(perm[maxpi][i]));
		
		return markers2;
	}
	
	private ArrayList<Integer> polishFast(ArrayList<Integer> markers, int permFirst[][], int permSecond[][])
	{
		return polishFast(markers, permFirst, permSecond, 1.0);
	}

	private ArrayList<Integer> polishFast(ArrayList<Integer> markers, int permFirst[][], int permSecond[][], double rate)
	{
		int numMarkers = markers.size();

		int k = permFirst[0].length;
		
		if (k > numMarkers)
			return markers;
		
		int tmp[] = new int[k];
		int tmp2[] = new int[k];

		double logL = score(markers);
		//if (k != 4)
		if (rate == 1.0)
			System.err.println("score = " + logL);
		
		for (SingleFamily sf : sfs) // calculate forward tables...
			sf.calcForward(markers);

		for (SingleFamily sf : sfs) // init backward tables...
			sf.initBackward(numMarkers & 1);

		for (int i = numMarkers - k; i >= 0; --i) {

			
			double maxL = Double.NEGATIVE_INFINITY;
			int numSolutions = 0;
			int perm[][] = (i == numMarkers - k) ? permFirst :  permSecond;
			for (int pi = 0; pi < perm.length; ++pi) {
				int p[] = perm[pi]; 
				
				for (int j = 0; j < k; ++j) { 
					tmp[j] = markers.get(i + p[j]);
				}
				
				double ll = 0.0;
				for (SingleFamily sf : sfs) // calculate polish tables...
					ll += sf.polishScore(tmp, i);
				
				if (pfs != null)
					for (PhysicalFamily pf : pfs)
						ll += pf.polishScore(tmp, i, markers);
	
				if (ll >= maxL - IMPROVE_TOLERANCE) {
						if (ll >= maxL + IMPROVE_TOLERANCE)
							numSolutions = 1;
						else
							++numSolutions;
						if (Misc.random() <= 1.0 / numSolutions) { 
							for (int j = 0; j < k; ++j) 
								tmp2[j] = tmp[j];
						}
						maxL = Math.max(ll, maxL);
				}
				while (rate < 1.0 && Misc.random() > rate) // skip with prob 1-rate
					++pi;
			}
			for (int j = 0; j < k; ++j) {
				markers.set(i + j, tmp2[j]);
			}

			for (SingleFamily sf : sfs) {// calculate polish tables...
				sf.recScale.update(tmp2, i); // update recombination scaler
				sf.updateBackward(i + k - 1, 1, markers);
			}
			
		}
		
		ArrayList<Integer> markersRet = new ArrayList<Integer>(); 
		for (int m = 0; m < numMarkers; ++m)
			markersRet.add((markers.get(m)));
		
		return markersRet;
	}
	
	
	private class MergeClass{
		private ArrayList<Integer> positions = null; 
		private ArrayList<Integer> markers = null;
		
		private int numPositions = 0;
		private int numMarkers = 0;

		double f1[] = null;
		double f2[] = null;

		double b1[] = null;
		double b2[] = null;
		
		private ArrayList<ArrayList<Integer>> ret = null;

		ParallelScoreCalculator psc;			
		
		private class ParallelScoreCalculator{
			double scores[] = null;
			int p1, p2, m1, m2, pos;
			int currentP, currentM;
			int nt;
			long numCalculated;
			boolean backward;
			
			private class UpdateScoreThread implements Runnable{
				int currentP, currentM, start, end;

				UpdateScoreThread(int start, int end, int currentP, int currentM) {
					this.currentM = currentM;
					this.currentP = currentP;
					this.start = start;
					this.end = end;
				}
				
				
				@Override
				public void run()
				{
					if (backward) {
					    out:while (currentP >= p1) {
								while (currentM >= m1) {
									scores[start++] = score(markers.get(currentM), positions.get(currentP));
									--currentM;
									if (start >= end)
										break out;
								}
								currentM = m2 - 1;
								--currentP;
							}
					} else {
					    out:while (currentP < p2) {
								while (currentM < m2) {
									scores[start++] = score(markers.get(currentM), positions.get(currentP));
									++currentM;
									if (start >= end)
										break out;
								}
								currentM = m1;
								++currentP;
							}
						}
				}
			}
			
			public ParallelScoreCalculator() {
				scores = scoresM;
			}

			private void updateScores()
			{
				pos = 0;
				int numElements = (int) Math.min(scores.length, ((long) (m2 - m1)) * (p2 - p1) - numCalculated);
				numCalculated += numElements;
				this.nt = calcNumThreadsScore(m2 - m1, p2 - p1, qTables.length, numElements);
				
				if (nt <= 1) { // nt <= 1
					int i = 0;
					if (backward) {
					    out:while (currentP >= p1) {
								while (currentM >= m1) {
									scores[i++] = score(markers.get(currentM), positions.get(currentP));
									--currentM;
									if (i >= scores.length)
										break out;
								}
								currentM = m2 - 1;
								--currentP;
							}
					} else {
				    out:while (currentP < p2) {
							while (currentM < m2) {
								scores[i++] = score(markers.get(currentM), positions.get(currentP));
								++currentM;
								if (i >= scores.length)
									break out;
							}
							currentM = m1;
							++currentP;
						}
					}
				} else {
					//System.err.println(nt + " threads in merge " + numMarkers + " markers and " + numPositions + " positions ");
					//if ((m2 - m1) * (p2 - p1) > scores.length)
					//	System.err.println("split " + currentP + " " + currentM + " (" + p1 + "," + p2 + ") " + "(" + m1 + "," + m2 + ") " + scores.length);
					
					Thread threads[] = new Thread[nt]; 
					
					for (int t = 0; t < nt; ++t) {
						int start = t * numElements / nt;
						int end = (t + 1) * numElements / nt;
						threads[t] = new Thread(new UpdateScoreThread(start, end, currentP, currentM));
						threads[t].start();
						if (backward) {
							currentM -= (end - start) % (m2 - m1);
							currentP -= (end - start) / (m2 - m1);
							if (currentM < m1) {
								currentM += m2 - m1;
								--currentP;
							}
						} else {
							currentM += (end - start) % (m2 - m1); 
							currentP += (end - start) / (m2 - m1);
							if (currentM >= m2) {
								currentM -= m2 - m1;
								++currentP;
							}
						}
					}
					try { 
						for (int t = 0; t < nt; ++t)
							threads[t].join();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
			
			
			public void setArea(int p1, int p2, int m1, int m2, boolean backward)
			{
				this.p1 = p1;
				this.p2 = p2;
				this.m1 = m1;
				this.m2 = m2;
				this.backward = backward;
				this.numCalculated = 0;
				
				if (backward) {
					currentP = p2 - 1; 
					currentM = m2 - 1;
				} else { 
					currentP = p1; 
					currentM = m1;
				}
				updateScores();
				//nt = scores.length / (m2 - m1);  
			}
			
			public double getNextScore(){
				if (pos >= scores.length) {
					updateScores();
				}
				return scores[pos++];
			}
		}
		
		
		MergeClass(ArrayList<Integer> markers1, ArrayList<Integer> markers2) {
			markers = markers2;
			
			ret = new ArrayList<ArrayList<Integer>>(); // there might be an overhead of creating this...
			
			ret.add(new ArrayList<Integer>()); // first place
			for (int m : markers1) {
				ret.add(new ArrayList<Integer>());
				//ret.get(ret.size() - 1).add(m);
			}
			positions = getPositions(markers1, ret);
			//System.err.println("Hiphei");

			numPositions = positions.size();
			numMarkers = markers.size();

			f1 = forwardM1;//new double[numMarkers];
			f2 = forwardM2;//new double[numMarkers];
			b1 = backwardM1;//new double[numMarkers];
			b2 = backwardM2;//new double[numMarkers];

			psc = new ParallelScoreCalculator(); 
			
		}
		
		//boolean print = true;
		//double tmpscore = 0.0;
		
		public ArrayList<Integer> merge()
		{
//			print = true;
//			tmpscore = 0.0;
			
			merge(0, numPositions, 0, numMarkers);
//			System.err.println(tmpscore);

			ArrayList<Integer> ret2 = new ArrayList<Integer>();
			for (ArrayList<Integer> l : ret)
				ret2.addAll(l);
			
			//System.err.println(numPositions + "\t" + ret2);
			
			
			return ret2;
		}
	
		public void merge(int p1, int p2, int m1, int m2)
		{
			if (m2 <= m1) // no markers
				return;
			if (p2 - p1 <= 1) { // only one position
				assert(p2 - p1 == 1);
				int pos = positions.get(p1);
				for (int j = m1; j < m2; ++j) {  
					ret.get(pos).add(markers.get(j));
	//				tmpscore += score(markers.get(j), pos);
				}
				return;
			}
			
			int mid = (p1 + p2) / 2;

			double f[] = calcMergeForward_old(p1, mid, m1, m2);
			double b[] = calcMergeBackward_old(mid, p2, m1, m2);
			double max = Double.NEGATIVE_INFINITY;
			int midMarker = 0;
			
			int numSolutions = 0;
			
			for (int j = m1 - 1; j < m2; ++j) {
				double r = 0.0;
				if (j >= m1) 
					r = f[j];
				if (j < m2 - 1)
					r += b[j + 1];
				if (r > max - IMPROVE_TOLERANCE) {
					if (r > max + IMPROVE_TOLERANCE) {
						numSolutions = 1;
						midMarker = j;
					} else {
						++numSolutions;
						if (Misc.random() <= 1.0 / numSolutions)
							midMarker = j;
					}
					max = Math.max(r, max);
				}
			}
//			if (print) {
//				System.err.println("mid = " + mid + "," + midMarker + " score = " + max);
//				print = false;
//			}
			//System.err.println("mid = " + mid + "," + midMarker + " score = " + max);
			//calcMergeForward_old(p1, p2, m1, m2);
			
			
			merge(p1, mid, m1, midMarker + 1);
			merge(mid, p2, midMarker + 1, m2);
		}
		
		private double[] calcMergeForward_old(int p1, int p2, int m1, int m2)
		{
			psc.setArea(p1, p2, m1, m2, false);
			
			double prevRow[] = f1;
			double row[] = f2;
			
			for (int j = m1; j < m2; ++j)
				prevRow[j] = Double.NEGATIVE_INFINITY;

			for (int i = p1; i < p2; ++i) {
				//row[m1] = Math.max(prevRow[m1], score(markers.get(m1), positions.get(i)));
				row[m1] = Math.max(prevRow[m1], psc.getNextScore());
				
				for (int j = m1 + 1; j < m2; ++j) {
					double v1 = prevRow[j];
					//double v2 = row[j - 1] + score(markers.get(j), positions.get(i));
					double v2 = row[j - 1] + psc.getNextScore();
					if (v1 > v2) {
						row[j] = v1;
					}
					else {
						row[j] = v2;
					}
						
				}
				double tmp[] = prevRow;
				prevRow = row;
				row = tmp;
			}
			//System.err.println("Merge score = " + prevRow[m2 - 1]);
			return prevRow;
		}

		private double[] calcMergeBackward_old(int p1, int p2, int m1, int m2)
		{
			psc.setArea(p1, p2, m1, m2, true);

			double prevRow[] = b1;
			double row[] = b2;
			
			for (int j = m1; j < m2; ++j)
				prevRow[j] = Double.NEGATIVE_INFINITY;

			for (int i = p2 - 1; i >= p1; --i) {
				//row[m2 - 1] = Math.max(prevRow[m2 - 1], score(markers.get(m2 - 1), positions.get(i)));
				row[m2 - 1] = Math.max(prevRow[m2 - 1], psc.getNextScore());
				
				for (int j = m2 - 2; j >= m1; --j) {
					double v1 = prevRow[j];
					//double v2 = row[j + 1] + score(markers.get(j), positions.get(i));
					double v2 = row[j + 1] + psc.getNextScore();
					if (v1 > v2) {
						row[j] = v1;
					}
					else {
						row[j] = v2;
					}
				}
				double tmp[] = prevRow;
				prevRow = row;
				row = tmp;
			}
			//System.err.println("Merge score = " + prevRow[m1]);
			return prevRow;
		}
		
		

		
/*		private double[] calcMergeForward_wrong(int p1, int p2, int m1, int m2)
		{
			int nt = calcNumThreads(p2 - p1, m2 - m1);
			if (nt <= 1)
				return calcMergeForward_old(p1, p2, m1, m2);
			
			for (int i = p1; i < p2; ++i)
				currentPosition[i] = -1;

			if (p1 > 0)
				currentPosition[p1 - 1] = Integer.MAX_VALUE;
			
			for (int j = m1; j < m2; ++j)
				f2[j] = Double.NEGATIVE_INFINITY;
			
			Thread threads[] = new Thread[nt]; 
			
			MergeRow ret = null;
			
			int rett = (p2 - p1 - 1) % nt;
			
			for (int t = 0; t < nt; ++t) {
				MergeRow mt = new MergeRow(p1 + t, p2, m1, m2, nt, ((t & 1) == 0) ? f1 : f2, ((t & 1) == 0) ? f2 : f1, false);
				threads[t] = new Thread(mt);
				threads[t].start();
				if (t == rett)
					ret = mt;
			}
			try { 
				for (int t = 0; t < nt; ++t)
					threads[t].join();
				return ret.getLastLine();
				
			} catch (Exception e) {
				e.printStackTrace();
				return null;
			}
		}*/
		
	}
	
	
	
	private ArrayList<Integer> merge(ArrayList<Integer> markers1, ArrayList<Integer> markers2) {
		MergeClass mc = new MergeClass(markers1, markers2);
		//System.err.println(markers1 + "|" + markers2);
		
		//System.err.println(mc.merge());
		//System.err.println(merge_old(markers1, markers2));
		
		//return merge_old(markers1, markers2);
		return mc.merge();
	}
	

	//TODO: backtrack in linear space..., now path stores |M1|*(|M2| + 1) booleans (bytes), done
	private ArrayList<Integer> merge_old(ArrayList<Integer> markers1, ArrayList<Integer> markers2)
	{
		ArrayList<Integer> positions = getPositions(markers1);
		//int numM1 = markers1.size();
		int numM2 = markers2.size();
		boolean path[][] = new boolean[numM2][positions.size()];

		double prevRow[] = new double[positions.size()];
		double row[] = new double[positions.size()];
		
		for (int j = 0; j < numM2; ++j) {
			double v = score(markers2.get(j), positions.get(0));
			path[j][0] = true;
			row[0] = prevRow[0] + v;
			
			for (int i = 1; i < positions.size(); ++i) {
				double v1 = prevRow[i] + score(markers2.get(j), positions.get(i));
				double v2 = row[i - 1];  
				if (v1 > v2 || (v1 == v2 && Misc.random() < 0.5)) {
					row[i] = v1;
					path[j][i] = true;
				}
				else {
					row[i] = v2;
					path[j][i] = false;
				}
					
			}
			double tmp[] = prevRow;
			prevRow = row;
			row = tmp;
		}
		//System.err.println("Merge score = " + prevRow[positions.size() - 1]);
		ArrayList<Integer> ret = new ArrayList<Integer>(); 		
		
		int prevPosI = markers1.size() - 1;
		int i = positions.size() - 1;
		int j = numM2 - 1;
		
		while (i > 0 || j >= 0) {
			if (j >= 0 && path[j][i]) {
				ret.add(markers2.get(j));
				--j;
			}
			else {
				--i;
				for (int k = prevPosI; k >= positions.get(i); --k)
					ret.add(markers1.get(k));
				prevPosI = positions.get(i) - 1;
	
			}
		}
		return ret;
	}

	//private int totalMarkers = 0;
	
	private ArrayList<Integer> fastOrder(ArrayList<Integer> markers, boolean permutate)
	{
		int numM = markers.size();		
		int perm[] = Misc.randomPermutation(numM);
		ArrayList<Integer> m = new ArrayList<Integer>();
		for (int i = 0; i < numM; ++i)
			if (permutate)
				m.add(markers.get(perm[i]));
			else
				m.add(markers.get(i));

		//totalMarkers = numM;
		m = fastOrder_(m, 0, numM);
		System.err.println(" done");
		
		//printOrder(m, -1, System.err); needs marker names
		return m;
	}
	
	private double computeE(ArrayList<Integer> m)
	{
		int n = m.size();
		int ret1 = 0;
		int ret2 = 0;
		
		for (int i = 0; i < n; ++i)
			System.err.print(markerNames.get(m.get(i)) +"\t");

		for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                            if (markerNames.get(m.get(i)) > markerNames.get(m.get(j)))
                                    ++ret1;
                            else if (markerNames.get(m.get(i)) < markerNames.get(m.get(j)))
                            		++ret2;
            }

		return ((double)Math.max(ret1, ret2)) / (ret1 + ret2);
	}
	
	private ArrayList<Integer> fastOrder_(ArrayList<Integer> markers, int level, int totalMarkers)
	{
		int numM = markers.size();
		if (numM > 2) {
			if (numM <= 5) {
				dataScale = Math.min(1.0, fullDataScale * totalMarkers / numM);
				markers = bestOrder(markers);
			} else {
				ArrayList<Integer> m1 = new ArrayList<Integer>();
				ArrayList<Integer> m2 = new ArrayList<Integer>();
				int r = 0;
				for (int i = 0; i < numM; ++i) {
					if  ((i & 1) == 0) {
						r = (int) (Misc.random() * 2);
						if (r == 0)
							m1.add(markers.get(i));
						else
							m2.add(markers.get(i));
					} else {
						if (r == 1)
							m1.add(markers.get(i));
						else
							m2.add(markers.get(i));
					}
					
//					if  ((i & 1) == 0) 
//						m1.add(markers.get(i));
//					else
//						m2.add(markers.get(i));
				}

				m1 = fastOrder_(m1, level + 1, totalMarkers);
				m2 = fastOrder_(m2, level + 1, totalMarkers);
				dataScale = Math.min(1.0, fullDataScale * totalMarkers / numM);
				
				//initTmpRec(((double) totalMarkers) / numM); // take into account the fact that subproblems have less markers...
				
				double ll = likelihood(m1, true); // to store tables, used in merge...
				//System.err.println("logL (half) = " + ll);
				ArrayList<Integer> merge1 = merge(m1, m2);
				Collections.reverse(m2);
				ArrayList<Integer> merge2 = merge(m1, m2);
		
				merge2 = polishTest(merge2);
				double l2 = score(merge2);
				
				merge1 = polishTest(merge1);
				double l1 = score(merge1);
				
				if (l1 >= l2 || Double.isNaN(l2)) { // NaN check to be safe...
					markers = merge1;
					//System.err.println("logL = " + l1 + " merge1 (" + l2 + ") n=" + numM);
					//System.err.print(l1 + "\t" + l2 + "\t");
				}
				else {
					markers = merge2;
				}
				if (numM <= 32)
					markers = polishFast(markers, perm4, perm4_s, 1.1); //with this extra polish we get more accurate maps but it takes a long time...
			}
		}
		if (level <= 3)
			System.err.print((1 << level));
		//System.err.println("Test\t" + computeE(markers));
		return markers;
	}
	
	
	private double bestScore = 0;
	private boolean refineParentOrder;
	public double getBestScore()
	{
		return bestScore;
	}

	public ArrayList<Integer> findOrder(boolean findInitialOrder)
	{
		ArrayList<Integer> ret = new ArrayList<Integer>();
		ArrayList<Integer> markers = new ArrayList<Integer>();
		for (int i = 0; i < numOriginalMarkers; ++i) {
			markers.add(i);
			ret.add(i);
		}
		bestScore = score(markers);
		double s = bestScore;
		System.err.println("Initial score = " + bestScore);
		
		for (int it = 1; it <= numMergeIterations; ++it) {

			System.err.println("Iteration " + it + " of " + numMergeIterations);
			markers = fastOrder(markers, (it == 1) && findInitialOrder);
			System.err.println("Polishing map...");

			markers = polishFast(markers, perm4, perm4_s); // no need if polishFast is called in fastOrder
			s = score(markers);
			if (s > bestScore) {
				System.err.print("*");
				bestScore = s;
				for (int i = 0; i < numOriginalMarkers; ++i)
					ret.set(i, markers.get(i));
			}			
			System.err.println("score = " + s);
			markers = polishTest(markers);
			markers = polishFast(markers, perm5, perm5_s);
			s = score(markers);

			//System.err.println("Test\t" + markers + "\t" + computeE(markers) + "\t" + score(markers));
			 
			if (s > bestScore) {
				System.err.print("*");
				bestScore = s;
				for (int i = 0; i < numOriginalMarkers; ++i)
					ret.set(i, markers.get(i));
			} else if (Misc.random() < 0.75) {
				for (int i = 0; i < numOriginalMarkers; ++i)
					markers.set(i, ret.get(i));
			}
			System.err.println("score = " + s);
			
		}

		//System.err.println("Final polishing...");
		//ret = polish(ret, perm4first, perm4);
		//bestScore = score(ret);
		
		//improveOrder(ret, 5 * markers.size(), true);
		System.err.println("Final score = " + bestScore);
		
		return ret;
	}

	private double computeLOD(float q1[], float q2[], int phase)
	{
		double f = -23.0;
		double sum1 = 1e-10;
		double sum2 = 1e-10;
		for (int haplotype = 0; haplotype < 4; ++haplotype) {
			double q12 = q1[haplotype] + q2[haplotype ^ phase];
			f = Math.max(f, q12);
			sum1 += Math.exp(q1[haplotype]);
			sum2 += Math.exp(q2[haplotype]);
		}
		return 4.0 * Math.exp(f) / (sum1 * sum2);
	}

	//TODO: Could be more efficient by using StringBuilder
	public void printLOD(ArrayList<Integer> markers, PrintStream stream)
	{
		//ArrayList<Integer> markers = new ArrayList<Integer>();
		//for (int i = 0; i < numMarkers; ++i)
		//	markers.add(i);
		
		double ll = likelihood(markers, true); // to store tables, used in merge...
		//System.err.println("score = " + ll);
		ArrayList<Integer> positions = this.getUniquePositions(markers);//getPositions(markers);
		
		int numF = originalHaplotypes.size();

		stream.print("position");
		for (int i = 0; i < positions.size(); ++i)
			stream.print("\tpos_" + i);
		stream.println();
		
		for (int i = 0; i < markers.size(); ++i) {
			stream.print("pos_" + i);
			for (int j = 0; j < positions.size(); ++j) {
				int fi = 0;
				double loglod2 = 0.0;
				for (int fam = 0; fam < numF; ++fam) {
					int numIndividuals = sfs[fam].numIndividuals;
					double maxLogLod = Double.NEGATIVE_INFINITY;
					for (int phase = 0; phase < 4; ++phase){
						double lod = 1.0;
						double loglod = 0.0;
						for (int hi = 0; hi < numIndividuals; ++hi) {
							float Q1[] = originalHaplotypes.get(fam)[i][hi];//qTables[fi][positions.get(i)];
							float Q2[] = qTables[fi + hi][positions.get(j)];
							lod *= computeLOD(Q1, Q2, phase);
							if (hi == numIndividuals - 1 || (hi & 15) == 15) {
								loglod += Math.log10(lod);
								lod = 1.0;
							}
						}
						maxLogLod = Math.max(maxLogLod, loglod);
					}
					loglod2 += maxLogLod;
					fi += numIndividuals;
				}
				stream.print("\t" + loglod2);
			}
			stream.println();
		}
	}

	//TODO: Could be more efficient by using StringBuilder
	public void printIntervals(ArrayList<Integer> markers, PrintStream stream, double limit)
	{
		//ArrayList<Integer> markers = new ArrayList<Integer>();
		//for (int i = 0; i < numMarkers; ++i)
		//	markers.add(i);
		int numMarkers = markers.size();
		
		double ll = likelihood(markers, true); // to store tables, used in merge...
		ArrayList<int[][]> phasedHaplotypes = getPhasedHaplotypes(markers, true, false);

		ArrayList<Integer> positions = new ArrayList<Integer>();
		ArrayList<Integer>  recombinations = new ArrayList<Integer>();
		
		int lastPos = 0;
		int recombinationPos = 0;
		while (lastPos < numMarkers) { 
			int nextPos = lastPos + 1;
			int delta = 0;
			while (nextPos < numMarkers) {
				for (int ph[][] : phasedHaplotypes) {
					int phh[] = ph[nextPos]; 
					delta += phh[phh.length - 1] + phh[phh.length - 2];
				}
				if (delta > 0)
					break;
				++nextPos;
			}
			positions.add((nextPos + lastPos + 1) / 2);
			recombinations.add(recombinationPos);
			recombinationPos += delta;
			lastPos = nextPos;
		}
		
		int numPositions = positions.size();
		for (int marker = 0; marker < numMarkers; ++marker) {
			stream.print((markerNames.get(marker) + 1));
			double maxLod = Double.NEGATIVE_INFINITY;
			for (int p = 0; p < numPositions; ++p)
				maxLod = Math.max(maxLod, score(marker, positions.get(p)));
			
			int p = 0;
			while (p < numPositions) {
				while (p < numPositions && score(marker, positions.get(p)) < maxLod - limit) {
					++p;
				}
				if (p < numPositions) {
					stream.print("\t" + recombinations.get(p));
					while (p < numPositions && score(marker, positions.get(p)) >= maxLod - limit) {
						++p;
					}
					stream.print("\t" + recombinations.get(p - 1));
				}
			}
			stream.println();
			
		}
		
		//TODO: continue...
		
		
		//System.err.println("score = " + ll);
		
		//stream.println();
	}
	
	

/*	public void refineParents(ArrayList<Integer> markers, boolean unknownOrder[][]) {

		//for (SingleFamily sf : sfs)
		//	sf.refineParents(markers, unknownOrder[sf.getFamilyIndex()]);
			
					
		//00,01,10,11 => 00,10,01,11
		//XX XX
		//XX    XX
		int numF = sfs.length;
		int numM = markers.size();
		for (int f = 0; f < numF; ++f) {
            float[][][] hap = originalHaplotypes.get(fam);
			double ll = sfs[f].phase(markers);
			for (int mi = 0; mi < numM; ++mi) {
				if (unknownOrder[f][mi]) {
					swap(mi); // or swapParents(markers.get(mi));
					double llnew = sfs[f].phase(markers);
					if (llnew <= ll)
						swapParents(mi); // or swapParents(markers.get(mi));
				}
			}
			System.err.println(sfs[f].numIndividuals);
		}
	//TODO: add code to figure parental genotypes...
	private void refineParents(ArrayList<Integer> order)
	{
		ArrayList<Integer> markers = new ArrayList<Integer>();
		for (int m1 : order) {
			boolean informative = false;
			for (Family2 f : data.getFamilies())
				if (f.isFatherInformative(m1) || f.isMotherInformative(m1))
					informative = true;
			if (informative)
				markers.add(m1);
		}
		
		System.err.println("Refining parents: number of markers = " + markers.size());
	
		ArrayList<Integer> markers2 = new ArrayList<Integer>();
		boolean unknownOrder[][] = new boolean[data.getNumFamilies()][data.getNumMarkers()];
		for (int mi = 0; mi < markers.size(); ++mi) {
			int m = markers.get(mi);
			char ipo[] = data.getIgnoreParentOrder(m);
			for (int f = 0; f < data.getNumFamilies(); ++f) 
				unknownOrder[f][mi] = (ipo[f] == '+');
			markers2.add(mi);
		}
		
		of.setHaplotypes(processHaplotypes(markers), informativeMarkers(markers), data, markers);
		of.refineParents(markers2, unknownOrder);
	}
	}*/

	public void setRefineParentOrder(boolean value) {
		refineParentOrder = value;
	}
	public boolean getRefineParentOrder() {
		return refineParentOrder;
	}
	
	
/*	public static void main(String args[])
	{
		OrderFinder of = new OrderFinder(10);
		of.testThreads();
		System.err.println("finished");
	}*/
}
