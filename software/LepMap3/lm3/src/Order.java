
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

//Marker ordering functions
import java.util.*;
import java.text.*;
import java.io.PrintStream;

public class Order {
	private Data2 data;
	private LGMap2 chromosomeMap;
	
	private int mappingFunction = 1;
	
	private boolean sexAveraged = false;
	
	private int numMarkers;
	
	private int printPhased;
	
	private OrderFinder of = new OrderFinder();

	ArrayList<Integer> lastOrder = new ArrayList<Integer>();

	public Order() {
	}

	
	public void setPrintPhased(int value)
	{
		printPhased = value;
	}

	public void setUsePhysical(int value, double prob)
	{
		of.setUsePhysical(value, prob);
	}
	
	public void setKosambi()
	{
		mappingFunction = 2;
	}

	public void setMorgan()
	{
		mappingFunction = 0;
	}

	public void setHaldane()
	{
		mappingFunction = 1;
	}
	public void init(Data2 data, LGMap2 cm) {
		this.data = data;
		this.chromosomeMap = cm;
		numMarkers = data.getNumMarkers();
	}
	public void setDataScale(String s, double maxS, double capS, int mode)
	{
		of.setDataScale(s, maxS, capS, mode);
	}

	public void setSelfingPhase(boolean value)
	{
		of.setSelfingPhase(value);
	}
	
	public void setRecombination(double rec1, double rec2, double inter1, double inter2) {
		of.setRecombination(rec1, rec2, inter1, inter2);
	}
	public void setIdenticalLimit(double limit) {
		of.setIdenticalLimit(limit);
	}
	public void setFilterWindow(int value)
	{
	}
	public void setNumThreads(int numThreads)
	{
		of.setNumThreads(numThreads);
	}
	public void setNumMergeIterations(int numIterations)
	{
		of.setNumMergeIterations(numIterations);
	}
	
	public void computeLODScores(String filename)
	{
		try {
			PrintStream stream = new PrintStream(filename);
			of.printLOD(lastOrder, stream);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void calculateIntervals(String filename, double limit)
	{
		try {
			PrintStream stream = new PrintStream(filename);
			of.printIntervals(lastOrder, stream, limit);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void evaluateOrder(String filename, boolean improveOrder) {
		evaluateOrder(filename, 0, improveOrder);
	}
	

	public void evaluateOrder(String filename, int chr, boolean improveOrder) {
		ArrayList<Integer> order = new ArrayList<Integer>();
		ArrayList<ArrayList<String>> table = Input.loadTable(filename, " \t");
		if (table == null)
			return;
		for (ArrayList<String> line : table) {
//			if (line.size() != 1)
//				Error.error(104);
			if (line.get(0).equals("***"))
				continue;
			int i = Integer.parseInt(line.get(0));
			//System.err.println(i);
			order.add(i - 1);
		}

		if (!improveOrder) {
			ArrayList<Integer> order2 = orderChromosome(order, false, false);
			printOrder(order2, 0, false, System.out, false);
		}	

		if (improveOrder) {
			ArrayList<Integer> order2 = orderChromosome(order, false, true);
			printOrder(order2, 0, false, System.out);
		}
		

	}
	
	private double mappingFunction(double r)
	{
		if (mappingFunction == 0)
			return r;
		else if (mappingFunction == 2)
			return kosambi(r);
		else
			return haldane(r);
	}
	
	private double haldane(double r) {
		return -0.5 * Math.log(1 - 2 * r);
	}
	private double kosambi(double r) {
		return 0.25 * Math.log ((1 + 2 * r) / (1 - 2 * r));
	}

	
	private String getHaplotypeLine(ArrayList<int [][]> hap, int marker)
	{
		char mapping[] = {'0', '1', '-'};
		StringBuilder ret = new StringBuilder();
		for (int h1[][] : hap) {
			ret.append("\t");
			int[] tmp = h1[marker];
			for (int i = 0; i < tmp.length - 2; ++i) // individual 
				ret.append(mapping[tmp[i]]);
			ret.append(" ");
			ret.append(tmp[tmp.length - 2]);
			ret.append(" ");
			ret.append(tmp[tmp.length - 1]);
		}
		return ret.toString();
	}

	private String getQualityLine(ArrayList<float [][]> q, int marker)
	{
		StringBuilder ret = new StringBuilder();
		for (float qf[][] : q) { // family
			float[] tmp = qf[marker];
			ret.append('\t');
			for (int i = 0; i < tmp.length; ++i) {// individual
				ret.append(tmp[i]);
				if (i != tmp.length - 1)
					ret.append(' ');
			}
			
		}
		return ret.toString();
	}
	
	private double[] getRecombinationRates(ArrayList<int [][]> hap, int marker) {
		double ret[] = new double[2];
		int n = 0;
		for (int h1[][] : hap) {
			int[] tmp = h1[marker];
			ret[0] += tmp[tmp.length - 2];
			ret[1] += tmp[tmp.length - 1];
			n += (tmp.length - 2) / 2;
		}
		ret[0] /= n;
		ret[1] /= n;
		return ret;
	}

	private double[] getRecombinationRates(ArrayList<int [][]> hap, int marker, int firstInf1[], int firstInf2[], int lastInf1[], int lastInf2[]) {
		double ret[] = new double[2];
		int n1 = 0;
		int n2 = 0;
		int fam = 0;
		for (int h1[][] : hap) {
			int[] tmp = h1[marker];
			ret[0] += tmp[tmp.length - 2];
			ret[1] += tmp[tmp.length - 1];
			if (marker > firstInf1[fam] && marker <= lastInf1[fam])
				n1 += (tmp.length - 2) / 2;
			if (marker > firstInf2[fam] && marker <= lastInf2[fam])
				n2 += (tmp.length - 2) / 2;
			++fam;
		}
		if (n1 > 0)
			ret[0] /= n1;
		if (n2 > 0)
			ret[1] /= n2;
		return ret;
	}
	
	public void setHyperPhaser(boolean value)
	{
		of.setHyperPhasing(value);
	}

	public void setPhasingIterations(int value)
	{
		of.setPhasingIterations(value);
	}

	public void printOrder(ArrayList<Integer> order, int chr, boolean estimateIndividualErrors, PrintStream stream) {
		printOrder(order, chr, estimateIndividualErrors, stream, true);
	}
	
	
	public void printOrder(ArrayList<Integer> order, int chr, boolean estimateIndividualErrors, PrintStream stream, boolean setHaplotypes) {
		
		if (setHaplotypes) {
			of.setHaplotypes(processHaplotypes(order), informativeMarkers(order), parentOrder(order), data, order);
		}
		
		ArrayList<Integer> order2 = new ArrayList<Integer>();
		for (int i = 0; i < order.size(); ++i)
			order2.add(i);

		lastOrder.clear();
		lastOrder.addAll(order2); 
		
		double ll = of.likelihood(order2);
		
		if (!Double.isNaN(bestScore) && ll != bestScore) {
			System.err.println("Warning: Inconsistencies in phasing, please try using randomPhase=1 and/or hyperPhaser=1 in OrderMarkers2");

			if (ll < bestScore - 1e-6) {
				if (!of.isHyperPhasing()) {
					System.err.println("Enabling hyperPhaser...");
					of.setHyperPhasing(true);
				}
				of.setPhasingIterations(3);
				ll = of.likelihood(order2);
				System.err.println("Phasing likelihoods = " + bestScore + " " + ll);
			}
		}
		
		
		//ArrayList<Double> iee = null; 
		//if (estimateIndividualErrors)
		//	iee = of.getIndividualErrorEstimates(order2);
		
		//of.setHaplotypes(data)
		ArrayList<int [][]> phasedHaplotypes = null;

		phasedHaplotypes = of.getPhasedHaplotypes(order2, false, (printPhased == 2 || printPhased == 4));

		ArrayList<float [][]> quality = null;
		if (printPhased > 2)
			quality = of.getPhasedQualityData(order2);
		
		DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols();
		otherSymbols.setDecimalSeparator('.');
		DecimalFormat df3 = new DecimalFormat("0.000", otherSymbols);
		DecimalFormat df4 = new DecimalFormat("#.####", otherSymbols);
		
		stream.println("#*** LG = " + chr + " likelihood = " + df4.format(ll)) ;
		stream.println("#marker_number\tmale_position\tfemale_position\t( 0.0 )[\tphased data]") ;

		double pos1 = 0.0;
		double pos2 = 0.0;
		int numChr = 1;
		int j = 0;

		int firstInformative1[] = new int[data.getNumFamilies()];
		int firstInformative2[] = new int[data.getNumFamilies()];
		int lastInformative1[] = new int[data.getNumFamilies()];
		int lastInformative2[] = new int[data.getNumFamilies()];
		Arrays.fill(firstInformative1, Integer.MAX_VALUE);
		Arrays.fill(firstInformative2, Integer.MAX_VALUE);
		
		for (int i = 0; i < order.size(); ++i) {
			int m = order.get(i);
			int fi = 0;
			for (Family2 f : data.getFamilies()) {
				if (f.isFatherInformative(m)) {
					if (i < firstInformative1[fi])
						firstInformative1[fi] = i;
					if (i > lastInformative1[fi])
						lastInformative1[fi] = i;
				}
				if (f.isMotherInformative(m)) {
					if (i < firstInformative2[fi])
						firstInformative2[fi] = i;
					if (i > lastInformative2[fi])
						lastInformative2[fi] = i;
				}
				++fi;
			}
		}
		
		for (int i = 0; i < order.size(); ++i) {
			double distance1 = 0.0;
			double distance2 = 0.0;
			double r[] = getRecombinationRates(phasedHaplotypes, i, firstInformative1, firstInformative2, lastInformative1, lastInformative2);//getRecombinationRates(phasedHaplotypes, i);
			if (sexAveraged) {
				distance1 = distance2 = mappingFunction(0.5 * (r[0] + r[1]));
			} else {
					distance1 = mappingFunction(r[0]);
					distance2 = mappingFunction(r[1]);
			}
			pos1 += distance1;
			pos2 += distance2;
			
			if (Double.isInfinite(pos1) || Double.isNaN(pos1) || Double.isInfinite(pos2) || Double.isNaN(pos2)) {
				stream.println("#*** LG = " + chr + "." + (numChr++));
				pos1 = 0.0;
				pos2 = 0.0;
			}
			String qstr = "";
			if (printPhased > 2)
				qstr = "#" + getQualityLine(quality, i);
			
			stream.println((order.get(i) + 1) + "\t" + df3.format(100 * pos1) + "\t" + df3.format(100 * pos2) + "\t( "
					+ df4.format(0.0) + " )" + ((printPhased != 0) ? getHaplotypeLine(phasedHaplotypes, i) : "") + qstr);
			
			
		}
		//stream.println("#COUNT = " + of.getCount());
		//if (estimateIndividualErrors) {
		//	stream.println("#Error estimate for individuals:");
		//	stream.println("#ID\terror_rate:");
		//	int i = 0;
		//	for (Family2 f : data.getFamilies())
		//		for (int child = 0; child < f.getNumChildren(); ++child) 
		//			stream.println("#" + f.getChildID(child) + "\t" + (iee.get(i++)));
		//	assert(i == iee.size());
		//}
		
		
	}

	public void setSexAveraged(boolean sexAveraged)
	{
		this.sexAveraged = sexAveraged;
	}

	public void orderChromosomes(boolean findIntial, boolean improveOrder, boolean estimateIndividualErrors) {
		orderChromosomes(-1, findIntial, improveOrder, estimateIndividualErrors);
	}

	public void orderChromosomes(int chr, boolean findIntial, boolean improveOrder, boolean estimateIndividualErrors) {
		ArrayList<ArrayList<Integer>> chromosomesPI = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> chromosomes = new ArrayList<ArrayList<Integer>>();

		int numInformative[] = new int[data.getNumFamilies()];

		
		for (int i = 0; i < numMarkers; ++i) {
			int fi = 0;
			for (Family2 f : data.getFamilies()) {
				if (f.isFatherInformative(i) || f.isMotherInformative(i))
					++numInformative[fi];
				++fi;
			}
		}
		int perm[] = Misc.randomPermutation(numMarkers);
		
		for (int ii = 0; ii < numMarkers; ++ii) {
			int i = perm[ii];
			boolean informative = false;
			int fi = 0;
			for (Family2 f : data.getFamilies())
				if ((numInformative[fi++] > 1 || !improveOrder) && (f.isFatherInformative(i) || f.isMotherInformative(i)))
					informative = true;
			int nameI = chromosomeMap.getLGName(i);
			for (int j = chromosomes.size(); j <= nameI; ++j) {
				chromosomesPI.add(new ArrayList<Integer>());
				chromosomes.add(new ArrayList<Integer>());
			}
			if (informative)
				chromosomesPI.get(nameI).add(i);
			chromosomes.get(nameI).add(i);
		}

		int numChromosomes = chromosomes.size() - 1;
		System.err.println("Number of LGs = " + numChromosomes);
		
		for (int c = 1; c <= numChromosomes; ++c) {
			if (chr <= 0 || c == chr) {
				
				ArrayList<Integer> markers = new ArrayList<Integer>();
				for (int i : chromosomesPI.get(c)) {
					markers.add(i);
				}
				
				ArrayList<Integer> order = orderChromosome(markers, findIntial, improveOrder);
				printOrder(order, c, estimateIndividualErrors, System.out);
			}
		}
	}
	
	public ArrayList<int[]> informativeMarkers(ArrayList<Integer> markers)
	{
		ArrayList<int[]> ret = new ArrayList<int[]>();
		for (Family2 f : data.getFamilies()) {
			int inf[] = new int[markers.size()];
			int k = 0;
			for (int m : markers) {
				int i = 0;
				if (f.isFatherInformative(m))
					i+=1;
				if (f.isMotherInformative(m))
					i+=2;
				inf[k++] = i;
			}
			ret.add(inf);
		}
		return ret;
	}
	
	public ArrayList<char[]> parentOrder(ArrayList<Integer> markers)
	{
		ArrayList<char[]> ret = new ArrayList<char[]>();
		if (of.getRefineParentOrder()) {
			int numM = markers.size();
			int numF = data.getNumFamilies();
			for (int f = 0; f < numF; ++f) {
				ret.add(new char[numM]);
			}
			for (int m = 0; m < numM; ++m) {
				char po[] = data.getIgnoreParentOrder(markers.get(m));
				for (int f = 0; f < numF; ++f)
					ret.get(f)[m] = po[f];
			}
		}
		return ret;
	}

	
	public ArrayList<float[][][]> processHaplotypes(ArrayList<Integer> markers)
	{
		ArrayList<float[][][]> ret = new ArrayList<float[][][]>();
		for (Family2 f : data.getFamilies()) {
			float h[][][] = new float[markers.size()][f.getNumChildren()][4];
			int k = 0;
			for (int m : markers) {
				float prob[] = f.getProb(m);
				if (prob == null) {
					prob = new float[4 * f.getNumChildren()];
					Arrays.fill(prob, 1.0f);
				}
				for (int i = 0; i < f.getNumChildren(); ++i) {
					float max = Float.MIN_NORMAL;
					for (int j = 0; j < 4; ++j)
						max = Math.max(prob[4 * i + j], max);
					for (int j = 0; j < 4; ++j) {
						h[k][i][j] = prob[4 * i + j] / max;
					}
				}
				++k;
			}
			ret.add(h);
		}
		//System.err.println(ret.get(0).length);
		
		
		/*for (Family2 f : data.getFamilies()) {
		
			byte h1[][] = f.getPaternalHaplotypes(markers);
			byte h2[][] = f.getMaternalHaplotypes(markers);

			for (int mi = 0; mi < markers.size(); ++mi) {
				int m = markers.get(mi);
				if (!f.isFatherInformative(m))
					for (int i = 0; i < h1.length; ++i)
						h1[i][mi] = 2;					
				if (!f.isMotherInformative(m))
					for (int i = 0; i < h2.length; ++i)
						h2[i][mi] = 2;			
				if (f.isFatherInformative(m) && f.isMotherInformative(m))
					for (int i = 0; i < h2.length; ++i)
						if (h1[i][mi] == 2 && h2[i][mi] == 2 && f.getGenotype(m, i) == 2) { 
							h1[i][mi] = 3;
							h2[i][mi] = 3;
						}
			}
			ret.add(Misc.transpose(h1));
			ret.add(Misc.transpose(h2));
		}*/
		return ret;
	}

	private double bestScore = 0;
	
	private ArrayList<Integer> orderChromosome(ArrayList<Integer> markersAll, boolean findInitial, boolean improveOrder) {

		ArrayList<Integer> markers = new ArrayList<Integer>();
		for (int m1 : markersAll) {
			boolean informative = false;
			for (Family2 f : data.getFamilies())
				if (f.isFatherInformative(m1) || f.isMotherInformative(m1))
					informative = true;
			if (informative)
				markers.add(m1);
		}
		
		System.err.println("Number of markers = " + markers.size());

		of.setHaplotypes(processHaplotypes(markers), informativeMarkers(markers), parentOrder(markers), data, markers);

		if (improveOrder) {
			ArrayList<Integer> markers2 = of.findOrder(findInitial);
			bestScore = of.getBestScore();
			for (int i = 0; i < markers.size(); ++i)
				markers2.set(i, markers.get(markers2.get(i)));
			markers = markers2;
		} else
			bestScore = Double.NaN;

		ArrayList<Integer> ret = new ArrayList<Integer>(); 
		for (int m : markers)
			ret.add(m);
		return ret;
	}


	public void setRefineParentOrder(boolean value) {
		// TODO Auto-generated method stub
		of.setRefineParentOrder(value);
	}


	
	
}
