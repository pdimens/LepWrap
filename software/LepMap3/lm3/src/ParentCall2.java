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
import java.util.HashMap;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;

// Posterior aware (subset of) Counts2Genotypes 
public class ParentCall2 {
	
	private final double SMALL = 1e-40; // to avoid problems with zero probabilities
	private int numIndividuals = 0;
	
	private	double ZLimit = Double.NEGATIVE_INFINITY; 
	private	double XLimit = Double.NEGATIVE_INFINITY; 
 	private	double familyLimit = 2;
 	
	private String names[] = null;
	
	private boolean ignoreParentOrder = false;
	private boolean outputParentPosterior = false;
	private boolean halfSibs = false;

	private boolean removeNonInformative = false;
	
	private ArrayList<ArrayList<Integer>> families = null;
	private ArrayList<int[]> parentIndex = null;
	private ArrayList<Integer> sex = null; 
	
	private double gpp[][] = new double[6][4];
	
	private HashMap<Integer, Integer> familyHash2 = new HashMap<Integer, Integer>();

	private HashMap<String, ArrayList<Integer>> halfSibHash = new HashMap<String, ArrayList<Integer>>();
	
	private int maxFather;
	private int maxMother;
	private int maxFather2;
	private int maxMother2;
	private double maxP;
	private double maxP2;
	private double maxP3;
	
	private int numCalledMarkers = 0;
	private int numCalledSexMarkers = 0;
	private int numInformativeMarkers = 0;

	private DataParser dp = null;
	
	private boolean outputRaw = false;
	

	public void setRaw(boolean value)
	{
		outputRaw = value;
	}

	public void setFamilyLimit(double limit)
	{
		familyLimit = limit;
	}
	public void setZLimit(double limit)
	{
		ZLimit = limit;
	}
	public void setXLimit(double limit)
	{
		XLimit = limit;
	}
	public void setIngnoreParentOrder(boolean value)
	{
		ignoreParentOrder = value;
		//if (ignoreParentOrder)
		//	System.err.println("ignoreParentOrder is not currently implemented");
	}
	public void setRemoveNonInformative(boolean value)
	{
		removeNonInformative = value;
	}

	public void setOutputParentPosterior(boolean value)
	{
		outputParentPosterior = value;
	}
	public void setHalfSibs(boolean value)
	{
		halfSibs = value;
		//if (ignoreParentOrder)
		//	System.err.println("ignoreParentOrder is not currently implemented");
	}
	
	
	//prints also...
	public void callParents(String filename, String vcfFile, String posteriorFile, boolean outputLikelihood) {
		dp = new DataParser();
		dp.loadFile(filename, vcfFile, posteriorFile);
		
		ArrayList<Double> dl = dp.getNextLine(true);

		families = dp.getFamilies();
		parentIndex = dp.getParentIndex();
		sex = dp.getSex();
		numIndividuals = dp.getNumIndividuals();
		int fam = 0;
		for (ArrayList<Integer> f : families) {
			for (int i : f)
				familyHash2.put(i, fam);
			int[] gp = parentIndex.get(fam);
			
			// add information about parents to detect half-sibs...
			String father = dp.getIndividualName(gp[0]);
			String mother = dp.getIndividualName(gp[1]);
			if (!halfSibHash.containsKey(father))
				halfSibHash.put(father, new ArrayList<Integer>());
			if (!halfSibHash.containsKey(mother))
				halfSibHash.put(mother, new ArrayList<Integer>());
			halfSibHash.get(father).add(gp[0]);
			halfSibHash.get(mother).add(gp[1]);
			
			for (int i : gp)
				if (i >= 0)
					familyHash2.put(i, fam);
			++fam;
		}
		for (String parent : halfSibHash.keySet())
			if (halfSibHash.get(parent).size() > 1) { 
				System.err.print("parent " + parent + " occurs " + halfSibHash.get(parent).size() + " times in the pedigree");
				System.err.println((halfSibs) ? "": " use halfSibs=1 to utilize this information");
			}
		while (dl != null) {
			//System.err.println(dl);
			if (outputRaw)
				printRaw(dp.getMarkerName(), dl);
			else
				callParents(dp.getMarkerName(), dl, outputLikelihood);
			dl = dp.getNextLine(true);
		}

		System.err.println("Number of called markers = " + numCalledMarkers + " (" + numInformativeMarkers + " informative)");
		System.err.println("Number of called Z/X markers = " + numCalledSexMarkers);
	}
	private void printResult(String suffix, double output[], int bestMother[], int bestFather[])
	{
		//System.out.print(suffix);
		StringBuilder sb = new StringBuilder(suffix);
		
		for (int i = 0; i < numIndividuals; ++i) {
			int f = familyHash2.get(i);
			int mother = bestMother[f];
			int father = bestFather[f];
			int index = 10 * (numIndividuals * (10 * father + mother) + i);
			//System.out.print("\t");
			sb.append('\t');
			double max = SMALL;
			for (int j = 0; j < 10; ++j)
				max = Math.max(max, output[index + j]);
			
/*			int fa1 = allele1(father);
			int fa2 = allele2(father);
			int ma1 = allele1(mother);
			int ma2 = allele2(mother);
			
			System.out.print(output[f][index + mapAlleles(fa1, ma1)] / max + " "); //haplotype 00
			System.out.print(output[f][index + mapAlleles(fa1, ma2)] / max + " "); //haplotype 01
			System.out.print(output[f][index + mapAlleles(fa2, ma1)] / max + " "); //haplotype 10
			System.out.print(output[f][index + mapAlleles(fa2, ma2)] / max);       //haplotype 11
*/
			for (int j = 0; j < 10; ++j) {
				if (max == SMALL)
					//System.out.print(1);
					sb.append('1');
				else {
					double r = output[index + j] / max;
					if (r == 0.0)
						//System.out.print(0); // save some space
						sb.append('0');
					else
						sb.append(r);
						//System.out.print(r);
				}
				if (j != 9)
					sb.append(' ');
					//System.out.print(" ");
			}
		}
		System.out.println(sb);
	}

	private void printRaw(String suffix, ArrayList<Double> dl)
	{
		StringBuilder sb = new StringBuilder(suffix);
		
		for (int i = 0; i < numIndividuals; ++i) {
			sb.append('\t');
			double max = SMALL;
			int index = 10 * i;
			for (int j = 0; j < 10; ++j)
				max = Math.max(max, dl.get(index + j));

			for (int j = 0; j < 10; ++j) {
				if (max == SMALL)
					sb.append('1');
				else {
					double r = dl.get(index + j) / max;
					if (r == 0.0)
						sb.append('0');
					else
						sb.append(r);
				}
				if (j != 9)
					sb.append(' ');
			}
		}
		System.out.println(sb);
	}
	
	
	private void printResult(String suffix, double output[], double likelihoods[][])
	{
		System.out.print(suffix);
		for (double ll[] : likelihoods) { // calculate scales for each family
			double max = Double.NEGATIVE_INFINITY;
			for (int k = 0; k < 100; ++k)
				max = Math.max(max, ll[k]);
			double sum = 0.0;
			for (int k = 0; k < 100; ++k) {
				double e = Misc.exp10(ll[k] - max); 
				ll[k] = e;
				sum += e;
			}
			for (int k = 0; k < 100; ++k)
				ll[k] /= sum;
		} 		
		
		for (int i = 0; i < numIndividuals; ++i) {
			int f = familyHash2.get(i);

			double max = SMALL;
			for (int j = 0; j < 10; ++j) {
				double r = 0.0;
				for (int father = 0; father < 10; ++father)
					for (int mother = 0; mother < 10; ++mother) {
						int index = 10 * (numIndividuals * (10 * father + mother) + i);
						double scale = likelihoods[f][10 * father + mother];
						r += scale * output[index + j];
					}
				max = Math.max(max, r);
			}
			System.out.print("\t");

			for (int j = 0; j < 10; ++j) {
				if (max == SMALL)
					System.out.print(1);
				else {
					double r = 0.0;
					for (int father = 0; father < 10; ++father)
						for (int mother = 0; mother < 10; ++mother) {
							int index = 10 * (numIndividuals * (10 * father + mother) + i);
							double scale = likelihoods[f][10 * father + mother];
							r += scale * output[index + j];
						}
					r = r / max;
					if (r == 0.0)
						System.out.print(0); // save some space
					else
						System.out.print(r);
				}
				if (j != 9)
					System.out.print(" ");
			}
		}
		System.out.println();
	}

	private void update_halfSibs(ArrayList<Double> posteriors, double likelihoods[][], boolean swapSex)
	{
		//marginalize parents out from likelihoods
		
		int numFamilies = likelihoods.length;
		
		double marginal1[][] = new double[numFamilies][10];//father 
		double marginal2[][] = new double[numFamilies][10];//mother
		for (int fam = 0; fam < numFamilies; ++fam){
			Arrays.fill(marginal1[fam], Double.NEGATIVE_INFINITY);
			Arrays.fill(marginal2[fam], Double.NEGATIVE_INFINITY);
		}

		for (int iterations = 1; iterations <= 2; ++iterations) 
			for (int fam = 0; fam < numFamilies; ++fam) { // for each family
				double ll[] = likelihoods[fam];
				
				ArrayList<Integer> familyIndex = families.get(fam); 
				int numInd = familyIndex.size();
				int fatherIndex = familyIndex.get(numInd - 1);
				int motherIndex = familyIndex.get(numInd - 2);
				if ((sex.get(motherIndex) == 1 && !swapSex) || (sex.get(motherIndex) == 2 && swapSex)) {
					int tmp = motherIndex;
					motherIndex = fatherIndex;
					fatherIndex = tmp;
				}
				if (iterations == 1) { // calculate marginals...
					for (int gFather = 0; gFather < 10; ++gFather)
						for (int gMother = 0; gMother < 10; ++gMother) {
							marginal1[fam][gFather] = Misc.logSum(marginal1[fam][gFather], ll[10 * gFather + gMother] * Constants.LN10);
							marginal2[fam][gMother] = Misc.logSum(marginal2[fam][gMother], ll[10 * gFather + gMother] * Constants.LN10);
						}
					double max1 = Double.NEGATIVE_INFINITY;
					double max2 = Double.NEGATIVE_INFINITY;
					for (int parent = 0; parent < 10; ++parent) {
						marginal1[fam][parent] -= Math.log(SMALL + posteriors.get(10 * fatherIndex + parent));					
						marginal2[fam][parent] -= Math.log(SMALL + posteriors.get(10 * motherIndex + parent));				
						max1 = Math.max(max1, marginal1[fam][parent]);
						max2 = Math.max(max2, marginal2[fam][parent]);
					}
					for (int parent = 0; parent < 10; ++parent) {
						marginal1[fam][parent] -= max1;
						marginal2[fam][parent] -= max2;
					}
				} else { // add marginals to parental genotype likelihoods...
					String father = dp.getIndividualName(fatherIndex);
					String mother = dp.getIndividualName(motherIndex);
					ArrayList<Integer> parents1 = halfSibHash.get(father); 
					ArrayList<Integer> parents2 = halfSibHash.get(mother);

					double logP[] = new double[10];
					for (int gFather = 0; gFather < 10; ++gFather)
						logP[gFather] = Math.log(SMALL * 1e-6 + posteriors.get(10 * fatherIndex + gFather));
					double logM[] = new double[10];
					for (int gMother = 0; gMother < 10; ++gMother)
						logM[gMother] = Math.log(SMALL * 1e-6 + posteriors.get(10 * motherIndex + gMother));

					for (int otherFather : parents1)
						if (otherFather != fatherIndex) {
							int fam2 = familyHash2.get(otherFather);
							for (int gFather = 0; gFather < 10; ++gFather)
								logP[gFather] += marginal1[fam2][gFather];
						}
					for (int otherMother : parents2)
						if (otherMother != motherIndex) {
							int fam2 = familyHash2.get(otherMother);
							for (int gMother = 0; gMother < 10; ++gMother)
								logM[gMother] += marginal2[fam2][gMother];								
						}
					double maxLogP = Double.NEGATIVE_INFINITY;
					for (int gFather = 0; gFather < 10; ++gFather)
						maxLogP = Math.max(maxLogP, logP[gFather]);
					double maxlogM = Double.NEGATIVE_INFINITY;
					for (int gMother = 0; gMother < 10; ++gMother)
						maxlogM = Math.max(maxlogM, logM[gMother]);

					for (int gFather = 0; gFather < 10; ++gFather)
						posteriors.set(10 * fatherIndex + gFather, Math.exp(logP[gFather] - maxLogP));
					for (int gMother = 0; gMother < 10; ++gMother)
						posteriors.set(10 * motherIndex + gMother, Math.exp(logM[gMother] - maxlogM));
				}
			}
		//System.out.println(posteriors);
	}
	
	
	public void callParents(String suffix, ArrayList<Double> posteriors, boolean outputLikelihood)
	{
		double lz = Double.NEGATIVE_INFINITY;
		double lx = Double.NEGATIVE_INFINITY;
		double ll = 0.0;
		int numFamilies = families.size();
		
		double output1[] = new double[numIndividuals * 10 * 100];
		double output2[] = new double[numIndividuals * 10 * 100];

		double likelihoods1[][] = new double[numFamilies][100];
		double likelihoods2[][] = new double[numFamilies][100];
		int bestFather[] = new int[numFamilies];
		int bestMother[] = new int[numFamilies];
		int secondBestFather[] = new int[numFamilies];
		int secondBestMother[] = new int[numFamilies];

		int bestFather2[] = new int[numFamilies];
		int bestMother2[] = new int[numFamilies];
		int secondBestFather2[] = new int[numFamilies];
		int secondBestMother2[] = new int[numFamilies];
		
		double maxPs[] = new double[numFamilies];
		double maxP2s[] = new double[numFamilies];
		double maxP3s[] = new double[numFamilies];

		double maxPs2[] = new double[numFamilies];
		double maxP2s2[] = new double[numFamilies];
		double maxP3s2[] = new double[numFamilies];
		
		boolean informative = false;
		boolean orderIgnored = false;
		
		char ignoreOrder[] = new char[numFamilies];
		
		
		int family = 0;
		for (ArrayList<Integer> familyIndex : families) {
			ll += likelihood(posteriors, familyIndex, parentIndex.get(family), output1, likelihoods1[family]);
			
			//System.err.println(maxP - maxP2);
			bestFather[family] = maxFather;
			bestMother[family] = maxMother;
			secondBestFather[family] = maxFather2;
			secondBestMother[family] = maxMother2;

			maxPs[family] = maxP;
			maxP2s[family] = maxP2;
			maxP3s[family] = maxP3;
			
			++family;
		}
		
		if (!Double.isInfinite(ZLimit)) {
			lz = 0.0;
			family = 0;
			for (ArrayList<Integer> familyIndex : families) {
				lz += likelihoodZ(posteriors, familyIndex, parentIndex.get(family), output2, likelihoods2[family]);
				bestFather2[family] = maxFather;
				bestMother2[family] = maxMother;
				secondBestFather2[family] = maxFather2;
				secondBestMother2[family] = maxMother2;
			
				maxPs2[family] = maxP;
				maxP2s2[family] = maxP2;
				maxP3s2[family] = maxP3;
				
				++family;
			}
			//if (lz >= ll + ZLimit)
			//System.err.println(lz + " Z  vs. " + ll + " " + maxFather + " x " + maxMother + " (" + maxFather2 + " x " + maxMother2 + ")");
		} else if (!Double.isInfinite(XLimit)) {
			lx = 0.0;
			family = 0;
			for (ArrayList<Integer> familyIndex : families) {
				lx += likelihoodX(posteriors, familyIndex, parentIndex.get(family), output2, likelihoods2[family]);
				bestFather2[family] = maxFather;
				bestMother2[family] = maxMother;
				secondBestFather2[family] = maxFather2;
				secondBestMother2[family] = maxMother2;

				maxPs2[family] = maxP;
				maxP2s2[family] = maxP2;
				maxP3s2[family] = maxP3;
				
				++family;
			}
			//if (lx >= ll + XLimit) {
			//	System.err.println(lx + " X  vs. " + ll + " " + maxFather + " x " + maxMother + " (" + maxFather2 + " x " + maxMother2 + ")" + " " + bestFather[family - 1] + ","+ bestMother[family - 1]); 				
			//}
				
		}
		
		boolean XMarker = lx - XLimit >= ll;
		boolean ZMarker = lz - ZLimit >= ll;
		boolean sexMarker = (XMarker || ZMarker); 

		if (sexMarker) {
			maxPs = maxPs2; 
			maxP2s = maxP2s2; 
			maxP3s = maxP3s2; 

			bestFather = bestFather2;
			bestMother = bestMother2; 
			secondBestFather = secondBestFather2; 
			secondBestMother = secondBestMother2;
			
			likelihoods1 = likelihoods2;
			output1 = output2;
			
			suffix = suffix + "*"; 
		
			if (outputLikelihood)
				System.err.println(suffix + "\t" + Math.max(lx, lz));
				
		} else {
			if (outputLikelihood)
				System.err.println(suffix + "\t" + ll);
		}
		
		if (halfSibs) {

			if (XMarker)
				update_halfSibs(posteriors, likelihoods1, true);
			else
				update_halfSibs(posteriors, likelihoods1, false);
			
			family = 0;
			Arrays.fill(output1, 0.0);
			for (ArrayList<Integer> familyIndex : families) {
				if (XMarker)
					likelihoodX(posteriors, familyIndex, parentIndex.get(family), output1, likelihoods1[family]);
				else if (ZMarker)
					likelihoodZ(posteriors, familyIndex, parentIndex.get(family), output1, likelihoods1[family]);
				else
					likelihood(posteriors, familyIndex, parentIndex.get(family), output1, likelihoods1[family]);
				
				//System.err.println(maxP - maxP2);
				bestFather[family] = maxFather;
				bestMother[family] = maxMother;
				secondBestFather[family] = maxFather2;
				secondBestMother[family] = maxMother2;
	
				maxPs[family] = maxP;
				maxP2s[family] = maxP2;
				maxP3s[family] = maxP3;
				
				++family;
				
			}
		}		

		for (family = 0; family < numFamilies; ++family) { 
			double maxP = maxPs[family];
			double maxP2 = maxP2s[family];
			double maxP3 = maxP3s[family];

			int maxFather = bestFather[family]; 
			int maxMother = bestMother[family]; 

			int maxFather2 = secondBestFather[family]; 
			int maxMother2 = secondBestMother[family]; 
			
			boolean orderClear = maxP - maxP2 >= familyLimit;
			
			if (orderClear || (ignoreParentOrder && maxMother == maxFather2 && maxFather == maxMother2 && maxP - maxP3 >= familyLimit)) {

				if (("0479".indexOf("" + maxFather) < 0 || "0479".indexOf("" + maxMother) < 0)) // informative 
					informative = true;
				if (sexMarker && maxFather != maxMother)
					informative = true;
				
				if (!orderClear) {
					orderIgnored = true;
					ignoreOrder[family] = '+';
				} else
					ignoreOrder[family] = '_';
				
				//in case of either parents could be informative, put father as informative
				if (!orderClear && ignoreParentOrder && maxP - maxP2 <= SMALL) { // order of parents ignored, no difference in likelihood
					if (("0479".indexOf("" + maxMother) < 0 && "0479".indexOf("" + maxFather) >= 0)) {// mother only informative
						bestFather[family] = maxFather2;
						bestMother[family] = maxMother2;
					}
				}
			} else {
				ArrayList<Integer> familyIndex = families.get(family);
				for (int i : familyIndex) // clear posteriors 
					for (int j = 0; j < 10; ++j)
						output1[10 * numIndividuals * (10 * maxFather + maxMother) + j + 10 * i] = 1.0;
	
				if (outputParentPosterior) {
					for (int k = 0; k < 100; ++k)
						for (int i : familyIndex) // clear posteriors 
							for (int j = 0; j < 10; ++j)
								output1[10 * numIndividuals * k + j + 10 * i] = 1.0;
				}
			
				ignoreOrder[family] = '_';
			}
			//System.err.println(ll + " " + maxFather + " x " + maxMother + " (" + maxFather2 + " x " + maxMother2 + ")");
		}
		
		if (informative)
			++numInformativeMarkers;
		if (informative || !removeNonInformative) {
			if (orderIgnored) {
				if (!outputParentPosterior)
					printResult(suffix + new String(ignoreOrder), output1, bestMother, bestFather);
				else
					printResult(suffix + new String(ignoreOrder), output1, likelihoods1);
			}
			else
				if (!outputParentPosterior)
					printResult(suffix, output1, bestMother, bestFather);
				else
					printResult(suffix, output1, likelihoods1);
			++numCalledMarkers;
			if (sexMarker)
				++numCalledSexMarkers;
		}
	}

	private int mapAlleles(int a1, int a2)
	{
		if (a1 > a2)
			return mapAlleles(a2, a1);	// a1 <= a2
		if (a2 > 3)
			Error.error(505);
		if (a1 == 0)
			return a2;
		if (a1 == 1)
			return 3 + a2;
		if (a1 == 2)
			return 5 + a2;
		if (a1 == 3)
			return 9;
		Error.error(505);
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
		Error.error(505);
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
		Error.error(505);
		return -1;
	}		
	public double likelihoodZ(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex, int grandParentIndex[], double output[], double likelihoods[])
	{
		return likelihoodZX(posteriors, familyIndex, grandParentIndex, output, likelihoods, true);
	}
	public double likelihoodX(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex, int grandParentIndex[], double output[], double likelihoods[])
	{
		return likelihoodZX(posteriors, familyIndex, grandParentIndex, output, likelihoods, false);
	}
	public double likelihoodZX(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex, int grandParentIndex[], double output[], double likelihoods[], boolean isZ) {
		int numInd = familyIndex.size();
		int fatherIndex = familyIndex.get(numInd - 1);
		int motherIndex = familyIndex.get(numInd - 2);
		if ((sex.get(motherIndex) == 1 && isZ) || (sex.get(motherIndex) == 2 && !isZ)) {
			int tmp = motherIndex;
			motherIndex = fatherIndex;
			fatherIndex = tmp;
		}
		fatherIndex *= 10;
		motherIndex *= 10;		
		maxP = Double.NEGATIVE_INFINITY; 
		maxP2 = Double.NEGATIVE_INFINITY; 
		maxP3 = Double.NEGATIVE_INFINITY;
		Arrays.fill(likelihoods, Double.NEGATIVE_INFINITY);
		
		int outputIndex = 0;
		for (int fatherAllele1 = 0; fatherAllele1 < 4; ++fatherAllele1)  
			for (int fatherAllele2 = fatherAllele1; fatherAllele2 < 4; ++fatherAllele2)
				for (int motherAllele1 = 0; motherAllele1 < 4; ++motherAllele1)
					for (int motherAllele2 = motherAllele1; motherAllele2 < 4; ++motherAllele2) {
						if (motherAllele1 == motherAllele2) { 
							int gFather = mapAlleles(fatherAllele1, fatherAllele2);
							int gMother = mapAlleles(motherAllele1, motherAllele2);
							
							double pF = Math.log10(SMALL + posteriors.get(fatherIndex + gFather));						
							double pM = Math.log10(SMALL + posteriors.get(motherIndex + gMother));
							double p = pF + pM;

							//grandparents...
							for (int gpii = 2; gpii < 6; ++gpii) {
								int gpi = 10 * grandParentIndex[gpii];
								if (!isZ)
									gpi = 10 * grandParentIndex[7 - gpii]; // put parentIndex in correct order fox XY system...
								if (gpi >= 0) {
									Arrays.fill(gpp[gpii], 0.0);
									for (int g1 = 0; g1 < 10; ++g1) {
										double post = posteriors.get(gpi + g1);
										if (gpi == fatherIndex) //backcross support
											post = (g1 == gFather) ? 1 : 0;
										if (gpi == motherIndex) //backcross support
											post = (g1 == gMother) ? 1 : 0;

										int a1 = allele1(g1);
										int a2 = allele2(g1);	
										boolean valid = false;
										if (gpii == 2) {
											output[outputIndex + gpi + g1] = post;
											valid = true;
										}
										if (gpii == 3 && a1 == a2 && (a1 == fatherAllele1 || a1 == fatherAllele2)) {
											output[outputIndex + gpi + mapAlleles(a1, motherAllele1)] = post;
											valid = true;
										}
										if (gpii == 4 && (a1 == motherAllele1 || a2 == motherAllele1)) {
											output[outputIndex + gpi + g1] = post;
											valid = true;
										}
										if (gpii == 5 && a1 == a2) {
											output[outputIndex + gpi + mapAlleles(a1, motherAllele1)] = post;
											valid = true;
										}
										if (valid) {
											gpp[gpii][a1] = Math.max(gpp[gpii][a1], post);
											gpp[gpii][a2] = Math.max(gpp[gpii][a2], post);
										}

									}
									
									//normalize gpp
									double max = 0.0;
									for (int a1 = 0; a1 < 4; ++a1) 
										if (gpp[gpii][a1] > max)
											max = gpp[gpii][a1];
									double imax = 1.0;
									if (max > 0.0)
										imax = 1.0 / max;
									
									for (int a1 = 0; a1 < 4; ++a1) 
										gpp[gpii][a1] *= imax;
									
									//if (gpii == 2); // fathers grandfather
									//if (gpii == 3); // fathers grandmother
									//if (gpii == 4); // mothers grandfather
									//if (gpii == 5); // mothers grandmother
									
								} else
									Arrays.fill(gpp[gpii], 1.0);
							}
							
							//System.out.println("**");
							//System.out.println((fatherAllele1+1) + "x" + (fatherAllele2 + 1));
							//Misc.printArray(gpp);
							//System.out.println("**");

							//autosome:
							//p += Math.log10(SMALL + Math.max(gpp[2][fatherAllele1] * gpp[3][fatherAllele2], gpp[2][fatherAllele2] * gpp[3][fatherAllele1]));
							//p += Math.log10(SMALL + Math.max(gpp[4][motherAllele1] * gpp[5][motherAllele2], gpp[4][motherAllele2] * gpp[5][motherAllele1]));
							
							p += Math.log10(SMALL + Math.max(gpp[2][fatherAllele1] * gpp[3][fatherAllele2], gpp[2][fatherAllele2] * gpp[3][fatherAllele1]));
							p += Math.log10(SMALL + gpp[4][motherAllele1]);
							
							//grandparents...
							
							for (int i : familyIndex) {
								int j = 10 * i;
								if (j != motherIndex && j != fatherIndex) {
									int s = sex.get(i);
									if (!isZ && (s == 1 || s == 2))
										s = 3 - s;
									
								if (s == 0) {
	
									int g1 = mapAlleles(fatherAllele1, motherAllele1);
									int g2 = mapAlleles(fatherAllele2, motherAllele1);
									
									int g3 = mapAlleles(fatherAllele1, fatherAllele1);
									int g4 = mapAlleles(fatherAllele2, fatherAllele2);
									
									double p1 = posteriors.get(j + g1);
									double p2 = posteriors.get(j + g2);
									double p3 = posteriors.get(j + g3);
									double p4 = posteriors.get(j + g4);
									
									output[outputIndex + j + g1] += p1;
									output[outputIndex + j + g2] += p2;
									
									if (fatherAllele1 != fatherAllele2) { //note: identical code below...
										if (fatherAllele1 != motherAllele1)
											g3 = mapAlleles(fatherAllele1, motherAllele1);
										if (fatherAllele2 != motherAllele1)
											g4 = mapAlleles(fatherAllele2, motherAllele1);
									}
									
									output[outputIndex + j + g3] += p3;
									output[outputIndex + j + g4] += p4;
									
									p += Math.log10(SMALL + p1 + p2 + p3 + p4) - Constants.LOGFOUR;
									
								} else if (s == 1) {
										int g1 = mapAlleles(fatherAllele1, motherAllele1);
										int g2 = mapAlleles(fatherAllele2, motherAllele1);
										double p1 = posteriors.get(j + g1);
										double p2 = posteriors.get(j + g2);
										
										output[outputIndex + j + g1] += p1;
										output[outputIndex + j + g2] += p2;
										
										p += Math.log10(0.5 * SMALL + p1 + p2) - Constants.LOGTWO;
										
									} else if (s == 2) {
										int g1 = mapAlleles(fatherAllele1, fatherAllele1);
										int g2 = mapAlleles(fatherAllele2, fatherAllele2);
										double p1 = posteriors.get(j + g1);
										double p2 = posteriors.get(j + g2);
		
										if (fatherAllele1 != fatherAllele2) { //note: identical code above...
											if (fatherAllele1 != motherAllele1)
												g1 = mapAlleles(fatherAllele1, motherAllele1);
											if (fatherAllele2 != motherAllele1)
												g2 = mapAlleles(fatherAllele2, motherAllele1);
										}
										output[outputIndex + j + g1] += p1;
										output[outputIndex + j + g2] += p2;
										p += Math.log10(0.5 * SMALL + p1 + p2) - Constants.LOGTWO;

										
									} else {
										Error.error(506);
									}
								}
							}
							if (p > maxP) {
								maxP3 = maxP2;
								
								maxP2 = maxP;
								maxFather2 = maxFather;
								maxMother2 = maxMother;
								
								maxP = p;
								maxFather = gFather;
								maxMother = gMother;
							} else if (p > maxP2) {
								maxP3 = maxP2;
								
								maxP2 = p;
								maxFather2 = gFather;
								maxMother2 = gMother;
							} else if (p > maxP3)
								maxP3 = p;

							likelihoods[10 * gFather + gMother] = p;
						
							if (fatherAllele1 == fatherAllele2 && fatherAllele1 != motherAllele1) {
								gMother = mapAlleles(fatherAllele1, motherAllele1);
							}
	
							output[outputIndex + fatherIndex + gFather] += 1;
							output[outputIndex + motherIndex + gMother] += 1;
						} // end of if (motherAllele1 == motherAllele2)
						outputIndex += 10 * numIndividuals;											
					}
		return maxP;
		//System.err.println(maxMother + " x " + maxFather);
	}	

	public double likelihood(ArrayList<Double> posteriors, ArrayList<Integer> familyIndex, int grandParentIndex[], double output[], double likelihoods[]) {
		int numInd = familyIndex.size();
		int fatherIndex = familyIndex.get(numInd - 1);
		int motherIndex = familyIndex.get(numInd - 2);
		
		if (sex.get(motherIndex) == 1) { // is this needed?
			int tmp = motherIndex;
			motherIndex = fatherIndex;
			fatherIndex = tmp;
		}
		//System.err.println("motherIndex = " + motherIndex);
		//System.err.println("fatherIndex = " + fatherIndex);
		fatherIndex *= 10;
		motherIndex *= 10;
		maxP = Double.NEGATIVE_INFINITY; 
		maxP2 = Double.NEGATIVE_INFINITY; 
		maxP3 = Double.NEGATIVE_INFINITY;
		int outputIndex = 0;
		for (int fatherAllele1 = 0; fatherAllele1 < 4; ++fatherAllele1)  
			for (int fatherAllele2 = fatherAllele1; fatherAllele2 < 4; ++fatherAllele2)
				for (int motherAllele1 = 0; motherAllele1 < 4; ++motherAllele1)  
					for (int motherAllele2 = motherAllele1; motherAllele2 < 4; ++motherAllele2) {
						int gFather = mapAlleles(fatherAllele1, fatherAllele2);
						int gMother = mapAlleles(motherAllele1, motherAllele2);
						double pF = Math.log10(SMALL + posteriors.get(fatherIndex + gFather));						
						double pM = Math.log10(SMALL + posteriors.get(motherIndex + gMother));
						double p = pF + pM;

						//grandparents...
						for (int gpii = 2; gpii < 6; ++gpii) {
							int gpi = 10 * grandParentIndex[gpii];
							if (gpi >= 0) {
								Arrays.fill(gpp[gpii], 0.0);
								
								for (int g1 = 0; g1 < 10; ++g1) {
									double post = posteriors.get(gpi + g1);

									//System.err.println(gpi + "\t" + fatherIndex);
									if (gpi == fatherIndex) { //backcross support 
										post = (g1 == gFather) ? 1 : 0;
										//System.err.println("Hiphei");
									}
									if (gpi == motherIndex) { //backcross support
										post = (g1 == gMother) ? 1 : 0;
									}
									
									output[outputIndex + gpi + g1] = post; //keep the posteriors for grandparents

									gpp[gpii][allele1(g1)] = Math.max(gpp[gpii][allele1(g1)], post);
									gpp[gpii][allele2(g1)] = Math.max(gpp[gpii][allele2(g1)], post);
								}
								//if (gpii == 2); // fathers grandfather
								//if (gpii == 3); // fathers grandmother
								//if (gpii == 4); // mothers grandfather
								//if (gpii == 5); // mothers grandmother
							} else
								Arrays.fill(gpp[gpii], 1.0);
						}
								
						p += Math.log10(SMALL + Math.max(gpp[2][fatherAllele1] * gpp[3][fatherAllele2], gpp[2][fatherAllele2] * gpp[3][fatherAllele1]));
						p += Math.log10(SMALL + Math.max(gpp[4][motherAllele1] * gpp[5][motherAllele2], gpp[4][motherAllele2] * gpp[5][motherAllele1]));
						
						//grandparents...					
						
						for (int i : familyIndex) {
							int j = 10 * i;
							if (j != motherIndex && j != fatherIndex) {
								int g1 = mapAlleles(fatherAllele1, motherAllele1);
								int g2 = mapAlleles(fatherAllele1, motherAllele2);
								int g3 = mapAlleles(fatherAllele2, motherAllele1);
								int g4 = mapAlleles(fatherAllele2, motherAllele2);
								double p1 = posteriors.get(j + g1);
								double p2 = posteriors.get(j + g2);
								double p3 = posteriors.get(j + g3);
								double p4 = posteriors.get(j + g4);
								output[outputIndex + j + g1] += p1;
								output[outputIndex + j + g2] += p2;
								output[outputIndex + j + g3] += p3;
								output[outputIndex + j + g4] += p4;
								p += Math.log10(SMALL + p1 + p2 + p3 + p4) - Constants.LOGFOUR;

							}
						}
						output[outputIndex + fatherIndex + gFather] += 1;
						output[outputIndex + motherIndex + gMother] += 1;
						outputIndex += 10 * numIndividuals;
						if (p > maxP) {
							maxP3 = maxP2;
							
							maxP2 = maxP;
							maxFather2 = maxFather;
							maxMother2 = maxMother;
							
							maxP = p;
							maxFather = gFather;
							maxMother = gMother;
						} else if (p > maxP2) {
							maxP3 = maxP2;
							
							maxP2 = p;
							maxFather2 = gFather;
							maxMother2 = gMother;
						} else if (p > maxP3)
							maxP3 = p;
						likelihoods[10 * gFather + gMother] = p;
					}
		return maxP;
		//System.err.println(maxMother + " x " + maxFather);
	}

	public static void main(String args[])
	{
		if (args.length == 0)
			Error.error(501);
		String extraParameters = "";
		System.out.print("#java ParentCall2");
		for (int i = 0; i < args.length; ++i) {
			extraParameters += args[i] + " ";
			System.out.print(" " + args[i]);
		}
		System.out.println();
		
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters))
			Error.error(501);
		pp.warning(new String[]{"data", "familyLimit", "ZLimit", "XLimit", "ignoreParentOrder", "removeNonInformative", "vcfFile", "posteriorFile", "outputParentPosterior", "halfSibs", "outputSNPLikelihood", "outputRaw"});
		
		String filename = pp.getValueAsString("data", null);
		if (filename == null)
			Error.error(501);
		
		String vcfFile = pp.getValueAsString("vcfFile", null);
		String posteriorFile = pp.getValueAsString("posteriorFile", null);
		if (vcfFile != null && posteriorFile != null)
			Error.error(514);

		ParentCall2 pg = new ParentCall2();

		double familyLimit = Double.parseDouble(pp.getValueAsString("familyLimit", "2")); 
		pg.setFamilyLimit(familyLimit);
		pg.setZLimit(Double.parseDouble(pp.getValueAsString("ZLimit", "" + Double.POSITIVE_INFINITY)));
		pg.setXLimit(Double.parseDouble(pp.getValueAsString("XLimit", "" + Double.POSITIVE_INFINITY)));
		if (pp.getValueAsString("ignoreParentOrder", "0").equals("1"))
			pg.setIngnoreParentOrder(true);
		if (pp.getValueAsString("removeNonInformative", "0").equals("1"))
			pg.setRemoveNonInformative(true);

		if (pp.getValueAsString("outputParentPosterior", "0").equals("1"))
			pg.setOutputParentPosterior(true);

		if (pp.getValueAsString("halfSibs", "0").equals("1"))
			pg.setHalfSibs(true);

		if (pp.getValueAsString("outputRaw", "0").equals("1"))
			pg.setRaw(true);
		
		
		pg.callParents(filename, vcfFile, posteriorFile, pp.getValueAsString("outputSNPLikelihood", "0").equals("1"));
	}
}
