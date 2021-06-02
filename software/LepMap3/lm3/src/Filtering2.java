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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;

import java.util.StringTokenizer;
// Posterior aware (subset of) Filtering
//TODO: removeIndividuals, removeGrandParents, parse pedigree better and detect errors
public class Filtering2 {
	
	private final double SMALL = 1e-20; // to avoid problems with zero probabilities
	private int numIndividuals = 0;
	
	private boolean removeNonInformative = false;
	
	private ArrayList<ArrayList<Integer>> families = new ArrayList<ArrayList<Integer>>(); 
	private ArrayList<String> familyNames = new ArrayList<String>(); 
	private ArrayList<Integer> sex = new ArrayList<Integer>(); 

	HashMap<String, Integer> familyHash = new HashMap<String, Integer>();
	HashMap<Integer, Integer> familyHash2 = new HashMap<Integer, Integer>();
	
	private double dataTolerance = 0.0;
	private double X2Limits[] = new double[4];
	
	private double mafLimit = 0.0;
	
	private boolean outputHWE = false;

	private boolean noSexFiltering = false;
	
	private double missingLimit = 0.0;
	private double missingTolerance = 0.0;

	//rivate double homozygoteLimit = 0.0;
	//private double heterozygoteLimit = 0.0;
	private double familyInformativeLimit = 0.0;
	
	private boolean convert2Biallelic = false;
	
	private final double LIKELIHOOD_TOLERANCE = 0.001; // used in EM and chiX2 tests
	
	private double initDistribution[] = new double[]{0.25,0.25,0.25,0.25};

	
	//private void setHeterozygoteLimit(double value)
	//{
	//	heterozygoteLimit = value;
	//}
	//private void setHomozygoteLimit(double value)
	//{
	//	homozygoteLimit = value;
	//}
	private void setNoSexFiltering(boolean value)
	{
		noSexFiltering = value;
	}

	private void setFamilyInformativeLimit(double value)
	{
		familyInformativeLimit = value;
	}
	
	
	private double tailX2(double x, int df)
	{
		return GammaFunction.incompleteGammaQ(0.5 * df, 0.5 * x);
	}

	
	private	double findLimit(int df, double tolerance)
	{
		if (tolerance <= 0.0)
			return Double.MAX_VALUE;
		//System.err.println();
		
		double x = 0.0;
		double adder = 2.0;
	
		while (adder >= 0.5 * LIKELIHOOD_TOLERANCE) {
			while (tailX2(x, df) >= tolerance)
				x += adder;
			x -= adder;
			adder *= 0.5;
			//System.err.println(x);
		}
		return x;
	}
	
	public void setDataTolerance(double value)
	{
		dataTolerance = value;
		X2Limits[0] = Double.MAX_VALUE;
		X2Limits[1] = findLimit(1, value);
		X2Limits[2] = findLimit(2, value);
		X2Limits[3] = findLimit(3, value);
		System.err.println("chi^2 limits are " +  X2Limits[1] + ", " + X2Limits[2] + ", " + X2Limits[3]);
		
				
	}

	public void setRemoveNonInformative(boolean value)
	{
		removeNonInformative = value;
	}
	
	public void setOutputHWE(boolean value){
		outputHWE = value;
	}
	public void setMAFLimit(double value){
		mafLimit = value;
	}
	public void setMissingLimit(double value){
		missingLimit = value;
	}
	public void setMissingTolerance(double value){
		missingTolerance = value;
	}
	
	private void parsePedigree(String names[], String fathers[], String mothers[]) {
		//TODO: Move the pedigree parsing logic to here
		//TODO: Add pedigree checking here...
		for (int i = 2; i < names.length; ++i) {
			if (fathers[i].equals("0")) { //move parents to the end...
				ArrayList<Integer> findex = families.get(familyHash2.get(i - 2));
				if (!findex.remove(new Integer(i - 2)))
					Error.error(1002);
				findex.add(i - 2);
			}
		}
	}
	
	//also prints...
	public void filter(String filename) {
		DataParser dp = new DataParser();
		dp.loadFile(filename, null, null);
		
		ArrayList<Double> dl = dp.getNextLine(true);

		families = dp.getFamilies();
		//parentIndex = dp.getParentIndex();
		sex = dp.getSex();
		numIndividuals = dp.getNumIndividuals();
		int fam = 0;
		for (ArrayList<Integer> f : families) {
			for (int i : f)
				familyHash2.put(i, fam);
			++fam;
		}
		int numMarkers = 0;
		while (dl != null) {
			++numMarkers;
			//System.err.println(dl);
			filter(dp.getMarkerName(), dl);
			dl = dp.getNextLine(true);
		}
		System.err.println("Number of markers = " + numMarkers);
		//System.err.println("Number of individuals = " + numIndividuals);
		//System.err.println("Number of families = " + families.size());
	}
	
	
	private void printResult(String suffix, ArrayList<Double> posteriors)
	{
		StringBuilder sb = new 	StringBuilder(suffix);
		//System.out.print(suffix);
		
		for (int i = 0; i < numIndividuals; ++i) {
			//System.out.print("\t");
			sb.append('\t');
			double max = SMALL;
			for (int j = 0; j < 10; ++j)
				max = Math.max(max, posteriors.get(10 * i + j));
			
			for (int j = 0; j < 10; ++j) {
				if (max == SMALL)
					//System.out.print(1);
					sb.append('1');
				else {
					double r = posteriors.get(10 * i + j) / max;
					if (r == 0.0)
						//System.out.print(0); // save some space
						sb.append('0');
					else
						//System.out.print(r);
						sb.append(r);
				}
				if (j != 9)
					//System.out.print(" ");
					sb.append(' ');
			}
		}
		System.out.println(sb);
	}
	
	public void convertBiAllelic(ArrayList<Integer> familyIndex, ArrayList<Double> posteriors) {
		
		int numInd = familyIndex.size();
		int fatherIndex = familyIndex.get(numInd - 1);
		int motherIndex = familyIndex.get(numInd - 2);
		if (sex.get(motherIndex) == 1|| sex.get(motherIndex) == 2) {
			int tmp = motherIndex;
			motherIndex = fatherIndex;
			fatherIndex = tmp;
		}

		int motherG = -1;
		int fatherG = -1;
		
		for (int i = 0; i < 10; ++i)
			if (posteriors.get(10 * fatherIndex + i) == 1) {
				if (fatherG >= 0)
					return;
				fatherG = i;
			}
		
		for (int i = 0; i < 10; ++i)
			if (posteriors.get(10 * motherIndex + i) == 1) {
				if (motherG >= 0) 
					return;
				motherG = i;
			}
		
		int fatherA1 = allele1(fatherG);
		int fatherA2 = allele2(fatherG);
		int motherA1 = allele1(motherG);
		int motherA2 = allele2(motherG);
		if (fatherA1 != fatherA2 && motherA1 != motherA2) {
			if (fatherA1 != motherA1 || fatherA2 != motherA2) { // more than two alleles
				for (int i = 0; i < 10; ++i) {
					posteriors.set(10 * motherIndex + i, (i == 1) ? 1.0 : 0.0); // AC  
					posteriors.set(10 * fatherIndex + i, (i == 1) ? 1.0 : 0.0); // AC
				}
				for (int index : familyIndex)
					if (index != fatherIndex && index != motherIndex) {
						double p0 = posteriors.get(10 * index + mapAlleles(fatherA1, motherA1)); 
						double p1 = posteriors.get(10 * index + mapAlleles(fatherA1, motherA2)) + posteriors.get(10 * index + mapAlleles(fatherA2, motherA1)); 
						double p2 = posteriors.get(10 * index + mapAlleles(fatherA2, motherA2));
						for (int i = 0; i < 10; ++i)
							posteriors.set(10 * index + i, 0.0);
						double maxp = Math.max(Math.max(p0, p1), p2);
						posteriors.set(10 * index, p0 / maxp); //AA
						posteriors.set(10 * index + 1, p1 / maxp); //AC
						posteriors.set(10 * index + 4, p2 / maxp); //CC
					}
			}
		}
		//TODO: other cases (do not make difference...)
	} 
	
	public void filter(String suffix, ArrayList<Double> posteriors)
	{
		int numFamilies = families.size();
		double maf[] = new double[4];
		
		boolean skipDistortion = false;
		if (noSexFiltering && suffix.indexOf('*') >= 0)
			skipDistortion = true;
		
		//int family = 0;
		if (outputHWE)
			System.err.print(suffix);
		
		for (ArrayList<Integer> familyIndex : families) {
			if (convert2Biallelic) {
				convertBiAllelic(familyIndex, posteriors);
			}
			boolean filterMarker = isSignificantlyDistorted(posteriors, familyIndex, maf);
			if (skipDistortion)
				filterMarker = false;
			
			double scale = 1.0;
			if (mafLimit >= 1.0)
				scale = familyIndex.size() - 2; 
				
			for (int h = 0; h < 4; ++h)
				if (scale * maf[h] < mafLimit)
					filterMarker = true;

			scale = 1.0;
			if (missingLimit < 1.0)
				scale = 1.0 / (familyIndex.size() - 2);
			if (numMissingGenotypes(posteriors, familyIndex) * scale > missingLimit) 
				filterMarker = true;
			
			
			if (filterMarker)
			{
				for (int i : familyIndex)
					for (int j = 0; j < 10; ++j)
						posteriors.set(10 * i + j, 1.0);
			}
			//++family;
		}
		if (outputHWE)
			System.err.println();


		boolean informative = false;
		int informativeFamilies = 0;
		
		for (ArrayList<Integer> familyIndex : families) {
			int numInd = familyIndex.size();
			int parents[] = {familyIndex.get(numInd - 1), familyIndex.get(numInd - 2)};
			for (int i : parents) {
				int maxj = 0;
				double sum = 0.0;
				for (int j = 0; j < 10; ++j) {
					double p = posteriors.get(10 * i + j);
					sum += p;
					if (p >= 1.0 - SMALL)
						maxj = j;
				}
				if (sum <= 1.0 + SMALL && ("0479".indexOf("" + maxj) < 0)) {
					informative = true;
					++informativeFamilies;
				}
			}
		}
		if ((informative || !removeNonInformative) && informativeFamilies >= familyInformativeLimit)
			printResult(suffix, posteriors);
		
	}
//Returns frequencies for haplotypes 00, 01, 10 and 11
//Used to calculate distortion aware LOD scores...
	public static double[] getFrequencies(double posteriors[][])
	{
		Filtering2 f = new Filtering2();
		double ret[] = new double[4];
		f.em(posteriors, ret, Double.POSITIVE_INFINITY);
		return ret;
	}
	
	public double em(double posteriors[][])
	{
		return em(posteriors, new double[4], Double.POSITIVE_INFINITY);
	}
	
	//allow different distribution for selfing data...
	public void setInitDistribution(double het)
	{
		if (het >= 1.0 || het <= 0.0) {
			System.err.println("Error: Expected distribution must be < 1 and > 0");
			System.exit(-1);
		}
		initDistribution[0] = (1.0 - het) * 0.5;		
		initDistribution[1] = het * 0.5;
		initDistribution[2] = het * 0.5;		
		initDistribution[3] = (1.0 - het) * 0.5;		
	}
	
	public double em(double posteriors[][], double finalP[], double XLimit)
	{
		double p[] = new double[4];
		for (int i = 0; i < 4;++i)
			p[i] = initDistribution[i];
		
		double t[] = new double[]{0, 0, 0, 0};
		double q[] = new double[]{0, 0, 0, 0};
		double initL = Double.NEGATIVE_INFINITY;
		double finalL = Double.NEGATIVE_INFINITY;
		
		for (int iteration = 0; iteration < 100; ++iteration) {
			Arrays.fill(q, 0.0);
			double logL = 0.0;
			double ll = 1.0;
			for (double post[] : posteriors) {
				for (int i = 0; i < 4; ++i)
					t[i] = p[i] * post[i];
				double l = (t[0] + t[1] + t[2] + t[3]);
				if (l < 1e-100) // set missing if very small 
					l = 1e-100;
				double il = 1.0 / l;
				for (int i = 0; i < 4; ++i) 
					q[i] += t[i] * il;
				ll *= l;
				if (ll < 1e-200) { // log when needed
					logL += Math.log(ll);
					ll = 1.0;
				}
			}
			logL += Math.log(ll);
			//System.err.println(logL);
			
			double sum = q[0] + q[1] + q[2] + q[3]; 
			double isum = 1.0 / sum;
			for (int i = 0; i < 4; ++i)
				p[i] = q[i] * isum;
			//System.err.println("ll = " + ll);
			if (iteration == 0)
				initL = logL;
			else {
				double oldL = finalL;  
				finalL = logL;
				if (finalL - oldL < LIKELIHOOD_TOLERANCE || 2 * (finalL - initL) >= XLimit)
					break;
			}
		}
		//System.err.println("p = " + p[0] + "," + p[1] + "," + p[2] + "," + p[3]);
		for (int i = 0; i < 4; ++i)
			finalP[i] = p[i];
		
		return finalL - initL; // likelihood ratio...
	}
	
	public int numMissingGenotypes(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex)
	{
		int ret = 0;
		for (int i = 0; i < familyIndex.size() - 2; ++i) {
			int index = 10 * familyIndex.get(i);		

			double sum = 0.0;
			for (int j = 0; j < 10; ++j)
				sum += posteriors.get(index + j);

			if (sum >= 1.0 + missingTolerance || sum == 0)
				++ret;
		}
		return ret;
	}

	public boolean isSignificantlyDistorted(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex, double maf[]) {
		int numInd = familyIndex.size();
		int fatherIndex = familyIndex.get(numInd - 1);
		int motherIndex = familyIndex.get(numInd - 2);
		if (sex.get(motherIndex) == 1 || sex.get(fatherIndex) == 2) {
			int tmp = motherIndex;
			motherIndex = fatherIndex;
			fatherIndex = tmp;
		}
		fatherIndex *= 10;
		motherIndex *= 10;		

		int motherG = -1;
		int fatherG = -1;
		
		for (int i = 0; i < 10; ++i)
			if (posteriors.get(fatherIndex + i) == 1) {
				if (fatherG >= 0) {
					//clear all
					if (outputHWE)
						System.err.print("\tX2 = " + 0 + "\tdf = " + 0);
					return true;
				}
				fatherG = i;
			}
		
		for (int i = 0; i < 10; ++i)
			if (posteriors.get(motherIndex + i) == 1) {
				if (motherG >= 0) {
					//clear all
					if (outputHWE)
						System.err.print("\tX2 = " + 0 + "\tdf = " + 0);
					return true;
				}
				motherG = i;
			}
		
		int fatherA1 = allele1(fatherG);
		int fatherA2 = allele2(fatherG);
		int motherA1 = allele1(motherG);
		int motherA2 = allele2(motherG);

		
		int df = 0;
		
		if (fatherA1 != fatherA2) {
			if (motherA1 != motherA2) {
				df = 2;
				if (fatherA1 != motherA1 || fatherA2 != motherA2)
					df = 3;
			}
			else
				df = 1;
		} else if (motherA1 != motherA2) {
			df = 1;
		}
		double post[][] = new double[familyIndex.size() - 2][4];
		
		for (int i = 0; i < familyIndex.size() - 2; ++i) {
			int index = 10 * familyIndex.get(i);
			post[i][0] = posteriors.get(index + mapAlleles(fatherA1, motherA1));
			post[i][1] = posteriors.get(index + mapAlleles(fatherA1, motherA2));
			post[i][2] = posteriors.get(index + mapAlleles(fatherA2, motherA1));
			post[i][3] = posteriors.get(index + mapAlleles(fatherA2, motherA2));
		}
		
		double X2 = 2 * em(post, maf, X2Limits[df]);
		if (outputHWE) {
			X2 = 2 * em(post, maf, Double.POSITIVE_INFINITY);
			System.err.print("\tX2 = " + X2 + "\tdf = " + df);
		} 
		return (X2 >= X2Limits[df]);
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
	
	public void setConvert2Biallelic(boolean convert) {
		convert2Biallelic = convert;
	}

	//TODO: "informativeMask", "families", "removeMarkers" 
	//TODO: dataTolerance separately for each family
	//TODO: Give X2 limits (separately for each family)
	
	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(1001);
		String extraParameters = "";
		System.out.print("#java Filtering2");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(1001);
		pp.warning(new String[]{"data", "dataTolerance", "removeNonInformative", "outputHWE", "MAFLimit", "missingLimit", "convert2Biallelic", "familyInformativeLimit", "noSexFiltering", "heterozygoteRate"}); // , "homozygoteLimit", "heterozygoteLimit"
		
		String filename = pp.getValueAsString("data", null);
		if (filename == null)
			Error.error(1001);

		Filtering2 filter = new Filtering2();
		//filter.set
		
		//double post[][] = new double[][]{{1, 1, 0, 1}, {1, 0, 1, 1}, {1, 1, 1, 0}, {0, 1, 1, 0}};
		//System.out.println(filter.em(post));
		//System.out.println();

		double dataTolerance = Double.parseDouble(pp.getValueAsString("dataTolerance", "0.001"));

		filter.setDataTolerance(dataTolerance);
		filter.setRemoveNonInformative(pp.getValueAsString("removeNonInformative", "0").equals("1"));
		filter.setOutputHWE(pp.getValueAsString("outputHWE", "0").equals("1"));
		filter.setMAFLimit(Double.parseDouble(pp.getValueAsString("MAFLimit", "0.0")));

		//filter.setHeterozygoteLimit(Double.parseDouble(pp.getValueAsString("heterozygoteLimit", "0.0")));
		//filter.setHomozygoteLimit(Double.parseDouble(pp.getValueAsString("homozygoteLimit", "0.0")));
		filter.setFamilyInformativeLimit(Integer.parseInt(pp.getValueAsString("familyInformativeLimit", "0")));

		filter.setMissingLimit(Double.parseDouble(pp.getValueAsString("missingLimit", 0, "" + Double.MAX_VALUE)));
		filter.setMissingTolerance(Double.parseDouble(pp.getValueAsString("missingLimit", 1, "0.1")));
		filter.setConvert2Biallelic(pp.getValueAsString("convert2Biallelic", "0").equals("1"));
		
		boolean noSexFiltering = pp.getValueAsString("noSexFiltering", "0").equals("1");
		filter.setNoSexFiltering(noSexFiltering);
		double het = Double.parseDouble(pp.getValueAsString("heterozygoteRate", "0.5"));
		filter.setInitDistribution(het);
		
		filter.filter(filename);
	}
}
