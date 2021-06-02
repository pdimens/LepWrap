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
    along with Lep-MAP.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/
import java.util.ArrayList;
import java.util.Arrays;

public class Family2 {
	
	private String name;
	
	private int numChildren;
	private int numMarkers;
	
	private int numInfFather = 0;
	private int numInfMother = 0;
	private int numThreeAlleles = 0;
	private int numFourAlleles = 0;	

	//use these instead of every 64 in LOD calculation
	private final double MAX_VALUE = 1e200;
	private final double MIN_VALUE = -1e200;
	
	private boolean phasedData = false;
	

	// how double informative markers are treated (typically 2 alleles ==> haplotypes not known for heterozygotes)
	private int lod3Mode = 1; // 1 normal, 2 identical, 3 identical or complementary

	// how double informative markers are treated (typically 2 alleles ==> haplotypes not known for heterozygotes)
	private boolean distortionLodMode = false; // segregation distortion aware LOD score 
	
	private ArrayList<float[]> haplotypeProb = new ArrayList<float[]>();
	
	private ArrayList<Integer> motherInf = new ArrayList<Integer>(); 
	private ArrayList<Integer> fatherInf = new ArrayList<Integer>(); 
	//private ArrayList<Integer> numAlleles = new ArrayList<Integer>(); 

	private double[][] freqs = null; 


	//double genotypeProbabilities[][];		
	
	public String getChildID(int c)
	{
		return "" + c;
	}
	
	public float[] getProb(int marker)
	{
		return haplotypeProb.get(marker);
	}
	
	
	public int getNumChildren()
	{
		return numChildren;
	}
	
	public String getName()
	{
		return name;
	}
	
	public void setLod3Mode(int mode)
	{
			lod3Mode = mode;
	}
	public void setPhasedData(boolean value)
	{
		phasedData = value;
	}
	
	
	public void setDistortionLodMode()
	{
		distortionLodMode = true;
		
		System.err.print("Calculating allele frequencies...");
		//calculate haplotype frequencies...
		freqs = new double[numMarkers][4];
		double post[][] = new double[numChildren][4];
		for (int i = 0; i < numMarkers; ++i)
			if (isFatherInformative(i) || isMotherInformative(i)) {
				float hp[] = haplotypeProb.get(i);
				for (int j = 0; j < numChildren; ++j) {
					post[j][0] = hp[4 * j + 0];
					post[j][1] = hp[4 * j + 1];
					post[j][2] = hp[4 * j + 2];
					post[j][3] = hp[4 * j + 3];
				}
				
				freqs[i] = Filtering2.getFrequencies(post);
				//System.err.println(freqs[i][0] + " " + freqs[i][1] + " " + freqs[i][2] + " "  + freqs[i][3]); 
			}
		System.err.println("done");
		
	}
	
	
	public int getNumMarkers()
	{
		return numMarkers;
	}
	
	public boolean isMotherInformative(int marker) {
		return (motherInf.get(marker) == 1);
	}
	public boolean isFatherInformative(int marker) {
		return (fatherInf.get(marker) == 1);
	}

	public void maskInformative(String informativeMask)
	{
		for (int m = 0; m < numMarkers; ++m) {
			int info = (isFatherInformative(m)?1:0) + (isMotherInformative(m)?2:0);
			if (informativeMask.indexOf("" + info) < 0)
				removeMarker(m);
		}
	}
	
	private int firstAllele(int genotype)
	{
		if (genotype < 0)
			Error.error(-995);
		if (genotype < 4)
			return 0;
		if (genotype < 7)
			return 1;
		if (genotype < 9)
			return 2;
		if (genotype == 9)
			return 3;
		Error.error(-995);
		return -1;
	}

	private int secondAllele(int genotype)
	{
		if (genotype < 0)
			Error.error(-995);
		if (genotype < 4)
			return genotype;
		if (genotype < 7)
			return genotype - 3;
		if (genotype < 9)
			return genotype - 5;
		if (genotype == 9)
			return 3;
		Error.error(-995);
		return -1;
	}

	private int mapAlleles(int a1, int a2)
	{
		if (a1 > a2)
			return mapAlleles(a2, a1);	// a1 <= a2
		if (a2 > 3)
			Error.error(-995);
		if (a1 == 0)
			return a2;
		if (a1 == 1)
			return 3 + a2;
		if (a1 == 2)
			return 5 + a2;
		if (a1 == 3)
			return 9;
		Error.error(-995);
		return -1;
	}
	
	
	private boolean isInformative(int genotype) 
	{
		return !(genotype == 0 || genotype == 4 || genotype == 7 || genotype == 9); 
	}

	public Family2(String name)
	{
		this.name = name;
	}
	public void updateNumChildren(int newNumChildren)
	{
		if (numChildren == 0)
			numChildren = newNumChildren;
		if (numChildren != newNumChildren)
			Error.error(4404);
	}

	public void addPedigree(ArrayList<Integer> familyIndex, ArrayList<String> names, ArrayList<Integer> sex, ArrayList<String> traits)
	{
		int numInd = familyIndex.size();
		int fatherIndex = familyIndex.get(numInd - 2);
		int motherIndex = familyIndex.get(numInd - 1);

		if (sex.get(motherIndex) == 1) {
			int tmp = fatherIndex;
			fatherIndex = motherIndex;
			motherIndex = tmp;
		}
		assert((sex.get(motherIndex) == 2) && (sex.get(fatherIndex) == 1));
		for (int i : familyIndex) {
			if (i != fatherIndex && i != motherIndex) {
			}
		}
	}
	public void addMarker() // add a dummy marker 
	{
		++numMarkers;
		fatherInf.add(0);
		motherInf.add(0);
		//numAlleles.add(2);
		haplotypeProb.add(null); // to save memory...
	}

	//calculate the probability of allele a from genotype posteriors
	private double pAllele(int a, ArrayList<Double> prob, int index)
	{
		if (index < 0)
			return 1.0;
		double ret = 0;
		index *= 10;
		for (int a1 = 0; a1 < 4; ++a1)
			for (int a2 = a1; a2 < 4; ++a2)
				if (a1 == a || a2 == a)
					ret += prob.get(index + mapAlleles(a1, a2));
		return ret;
	}
	
	public void addMarker(ArrayList<Double> alleleProb, ArrayList<Integer> familyIndex, ArrayList<Integer> sex, String informativeMask, boolean grandparentPhase, int parentIndex[])
	{

		int numInd = familyIndex.size();
		
		updateNumChildren(numInd - 2);
	
		int fatherIndex = familyIndex.get(numInd - 2);
		int motherIndex = familyIndex.get(numInd - 1);

		if (sex.get(motherIndex) == 1) {
			int tmp = fatherIndex;
			fatherIndex = motherIndex;
			motherIndex = tmp;
		}
		//System.err.println();
		
		assert((sex.get(motherIndex) == 2) && (sex.get(fatherIndex) == 1));
		
		int fatherA = 0;
		int motherA = 0;
		double sumMother = 0;
		double sumFather = 0;
		for (int i = 0; i < 10; ++i) {
			if (alleleProb.get(10 * fatherIndex + fatherA) < alleleProb.get(10 * fatherIndex + i))
				fatherA = i;
			if (alleleProb.get(10 * motherIndex + motherA) < alleleProb.get(10 * motherIndex + i))
				motherA = i;
			sumFather += alleleProb.get(10 * motherIndex + i);
			sumMother += alleleProb.get(10 * motherIndex + i);
		}
		final double FAMILY_LIMIT = 0.1; // how certain parental genotypes must be...
		
		int info = (sumFather <= 1 + FAMILY_LIMIT && isInformative(fatherA)?1:0) + (sumMother <= 1 + FAMILY_LIMIT && isInformative(motherA)?2:0);
		
		if (informativeMask.indexOf("" + info) < 0) { // to save memory...
			addMarker();
			return;
		}		

		int fa1 = firstAllele(fatherA);
		int fa2 = secondAllele(fatherA);
		int ma1 = firstAllele(motherA);
		int ma2 = secondAllele(motherA);
		
		if (grandparentPhase && info > 0) { //phasing logic with grandparents
			int p1 = parentIndex[0]; //parent1 
			int p2 = parentIndex[1]; //parent2
			int gp11 = parentIndex[2]; //grandparent1 of parent1
			int gp12 = parentIndex[3]; //grandparent2 of parent1
			int gp21 = parentIndex[4]; //grandparent1 of parent2
			int gp22 = parentIndex[5]; //grandparent2 of parent2
			
			final double PHASE_LIMIT = 0.01; // how certain phase must be...
			
			assert(p1 == fatherIndex && p2 == motherIndex);
			
			if (info == 1 || info == 3) { //father informative
				double at1 = pAllele(fa1, alleleProb, gp11) * pAllele(fa2, alleleProb, gp12);
				double ta1 = pAllele(fa2, alleleProb, gp11) * pAllele(fa1, alleleProb, gp12);
				if (at1 < ta1 * PHASE_LIMIT) { // second, first
					int tmp = fa1;
					fa1 = fa2;
					fa2 = tmp;
				} else if (ta1 < at1 * PHASE_LIMIT) { // first, second
					//default way will do
				} else {
					addMarker(); // remove marker
					return;
				}
			}
			if (info == 2 || info == 3) { // mother informative
				double at2 = pAllele(ma1, alleleProb, gp21) * pAllele(ma2, alleleProb, gp22);
				double ta2 = pAllele(ma2, alleleProb, gp21) * pAllele(ma1, alleleProb, gp22);

				if (at2 < ta2 * PHASE_LIMIT) { // second, first
					int tmp = ma1;
					ma1 = ma2;
					ma2 = tmp;
				} else if (ta2 < at2 * PHASE_LIMIT) { // first, second
					//default way will do
				} else {
					addMarker(); // remove marker
					return;
				}
				
			}
		} 

		fatherInf.add((isInformative(fatherA)) ? 1:0);
		motherInf.add((isInformative(motherA)) ? 1:0);
		
		int c[] = {-1, -1, -1, -1};

		int alleles = 2;
		
		if (isInformative(fatherA)) {
			++numInfFather;
			if (isInformative(motherA)) {
				++numInfMother;
				if (fatherA == motherA) // two alleles
					;
				else if (ma1 == fa1 || ma1 == fa2 || ma2 == fa1 || ma2 == fa2) {// three alleles
					++numThreeAlleles;
					alleles = 3;
				}
				else {// four alleles
					++numFourAlleles;
					alleles = 4;
				}
					
			}
		} else if (isInformative(motherA))
			++numInfMother;

		//numAlleles.add(alleles);

		c[0] = mapAlleles(fa1, ma1); 
		c[1] = mapAlleles(fa1, ma2); 
		c[2] = mapAlleles(fa2, ma1); 
		c[3] = mapAlleles(fa2, ma2);
		
		
		int numTable = 4 * (numInd - 2);
		float f[] = new float[numTable];
		

		for (int i = 0; i < numTable; i+=4)
			for (int j = 0; j < 4; ++j)
				f[i + j] = new Float(alleleProb.get(10 * familyIndex.get(i / 4) + c[j]));
				//if (c[j] >= 0) 
				//	f[i + j] = new Float(alleleProb.get(10 * familyIndex.get(i / 4) + c[j]));
				//else
				//	f[i + j] = 1.0f;
		normalize(f);
		
		
		haplotypeProb.add(f);
		/*
		for (int i = 0; i < numTable; i+=4) { 
			for (int j = 0; j < 4; ++j)
				System.err.print(f[i + j] + " ");
			System.err.print("\t");
		}
		System.err.println();
		*/
		++numMarkers;

	}
	
	public void printStatistics()
	{
		System.err.println("Family " + name + ":");		
		
		System.err.println("Number of paternally informative markers = " + numInfFather);
		System.err.println("Number of maternally informative markers = " + numInfMother);

		System.err.println("Number of double informative markers with three alleles = " + numThreeAlleles);
		System.err.println("Number of double informative markers with four alleles = " + numFourAlleles);
	}
	
	public void callGenotypes(int marker, double limit)
	{
		float f[] = getProb(marker);
		int tableSize = 4 * numChildren;

		if ((lod3Mode==2 || lod3Mode==3) && isFatherInformative(marker) && isMotherInformative(marker)) {
			//TODO : check lod3Mode==3

			for (int ci = 0; ci < tableSize; ci+=4) {
				double g0 = f[ci + 0];
				double g1 = f[ci + 1] + f[ci + 2];
				double g2 = f[ci + 3];
				if (g0 > 1 - limit) {
					f[ci + 0] = 1;					
					f[ci + 1] = 0;					
					f[ci + 2] = 0;					
					f[ci + 3] = 0;					
				} else if (g1 > 1 - limit) {
					f[ci + 0] = 0;					
					f[ci + 1] /= g1;					
					f[ci + 2] /= g1;					
					f[ci + 3] = 0;					
				} else if (g2 > 1 - limit) {
					f[ci + 0] = 0;					
					f[ci + 1] = 0;					
					f[ci + 2] = 0;					
					f[ci + 3] = 1;					
				}
			}
			return;
		}
		
		for (int ci = 0; ci < tableSize; ci+=4) {
			double pat0 = f[ci] + f[ci + 1];
			if (pat0 <= limit) {
				f[ci + 0] = 0;
				f[ci + 1] = 0;
				f[ci + 2] /= (1 - pat0);
				f[ci + 3] /= (1 - pat0);	
			}
			else if (pat0 >= 1 - limit) {
				f[ci + 0] /= pat0;
				f[ci + 1] /= pat0;		
				f[ci + 2] = 0;
				f[ci + 3] = 0;
			}
			double mat0 = f[ci] + f[ci + 2];

			if (mat0 <= limit) {
				f[ci + 0] = 0;
				f[ci + 2] = 0;
				f[ci + 1] /= (1 - mat0);
				f[ci + 3] /= (1 - mat0);	
			}
			else if (mat0 >= 1 - limit) {
				f[ci + 0] /= mat0;
				f[ci + 2] /= mat0;		
				f[ci + 1] = 0;
				f[ci + 3] = 0;
			}
		}
	}

	public void normalize(int m1)
	{
		normalize(haplotypeProb.get(m1));
	}

	private void normalize(float prob[])
	{
		for (int i = 0; i < prob.length; i+=4) {
			double sum = prob[i] + prob[i + 1] + prob[i + 2] + prob[i + 3] + 1e-20;
			for (int j = 0; j < 4; ++j) 
				prob[i + j] = (float) ((prob[i + j] + 0.25e-20) / sum);
		}
	}
	
	private final int mapping[][] = {{0,1,2,3}, {2,3,0,1}, {1,0,3,2}, {3,2,1,0}};
	
	public void scale(double value) {
		if (value != 1.0)
			for (float h[] : haplotypeProb) {
				if (h != null)
					for (int i = 0; i < h.length; ++i)
						h[i] = (float) Math.pow(h[i], value);
			}
	}
	
	private void swap(int i, int j, float t[])
    {
            float tmp = t[i];
            t[i] = t[j];
            t[j] = tmp;                     
    }

	public void randomPhase() {
		for (float h[] : haplotypeProb) {
			if (h != null) {
				if (Math.random() < 0.5)
					for (int i = 0; i < h.length; i+=4) {
						swap(i + 0, i + 2, h);
						swap(i + 1, i + 3, h);
					}
				if (Math.random() < 0.5)
					for (int i = 0; i < h.length; i+=4) {
						swap(i + 0, i + 1, h);
						swap(i + 2, i + 3, h);
					}
			}
		}
	}


	public void minError(float value) {
		if (value != 0.0)
			for (float h[] : haplotypeProb) {
				if (h != null)
					for (int i = 0; i < h.length; ++i)
						h[i] = (h[i] < value) ? value : h[i];
			}
	}
	
	double computeLOD(double theta1, double theta2, int m1, int m2)
	{
		if (isMotherInformative(m1) && isMotherInformative(m2)) { 
			if (isFatherInformative(m1) && isFatherInformative(m2)) { 
				//TODO: LOD3,LOD3'
				if (lod3Mode == 1) {
					//	return computeLOD(theta1, theta2, m1, m2, 0);
					double lod1 = computeLOD(theta1, theta2, m1, m2, 0);
					if (phasedData)
						return lod1;
					double lod2 = computeLOD(theta1, theta2, m1, m2, 1);
					double lod3 = computeLOD(theta1, theta2, m1, m2, 2);
					double lod4 = computeLOD(theta1, theta2, m1, m2, 3);
					return Math.max(Math.max(Math.max(lod1, lod2), lod3), lod4);
				} else if (lod3Mode == 2 && theta1 == theta2) {
					double lod1 = computeLOD2(theta1, m1, m2, 0);
					if (phasedData)
						return lod1;
					double lod2 = computeLOD2(theta1, m1, m2, 1);
					return Math.max(lod1, lod2);
				} else if (lod3Mode == 3 && theta1 == theta2) {
					double lod1 = computeLOD3(theta1, m1, m2, 0);
					if (phasedData)
						return lod1;
					double lod2 = computeLOD3(theta1, m1, m2, 1);
					return Math.max(lod1, lod2);
				} else
					Error.error(1014);
			} else {
				double lod1 = computeLOD_2(theta2, m1, m2, 0);
				if (phasedData)
					return lod1;
				double lod2 = computeLOD_2(theta2, m1, m2, 2);
				return Math.max(lod1, lod2);
			}
		} else if (isFatherInformative(m1) && isFatherInformative(m2)) {
			double lod1 = computeLOD_1(theta1, m1, m2, 0);
			if (phasedData)
				return lod1;
			double lod2 = computeLOD_1(theta1, m1, m2, 1);
			return Math.max(lod1, lod2);
		}
		return 0.0;
	}
	
	private void computeLNoLinkage12(int m1)
	{
		//freqs[];
	}
	private void computeLNoLinkage3(int phaseMapping[])
	{
				
	}
	private void computeLNoLinkage4(int phaseMapping[])
	{
				
	}
	
	//TODO: consider lod3Mode
	public double maxLodScore(double theta1, double theta2, int m1)
	{
		float f1[] = getProb(m1);
		if (f1 == null)
			return 0.0;
		
		int tableSize = 4 * numChildren;
		
		double lNoLinkage1 = 2.0;
		double lNoLinkage2 = 2.0;
		double lNoLinkage3 = 3.0;
		double lNoLinkage4 = 2.0;
		if (distortionLodMode) { // TODO: figure out how to implement here...
			lNoLinkage1 = 2.0;
			lNoLinkage2 = 2.0;
			lNoLinkage3 = 3.0;
			lNoLinkage4 = 2.0;
		}
		double logRet = 0.0;
		double ret = 1.0;
		for (int ci = 0; ci < tableSize; ci+=4)
		{
			double sum1 = f1[ci + 0] + f1[ci + 1];
			double pat = Math.max(sum1, (1 - sum1)); 
 
			sum1 = f1[ci + 0] + f1[ci + 2];
			double mat = Math.max(sum1, (1 - sum1));
						 
			ret *= lNoLinkage1 * (theta1 * (1 - pat) + (1 - theta1) * pat);
			ret *= lNoLinkage2 * (theta2 * (1 - mat) + (1 - theta2) * mat);
			if ((ci + 4 >= tableSize) || ret > MAX_VALUE || ret < MIN_VALUE ) {
				logRet += Math.log10(ret);
				ret = 1.0;
			}
		}
		if (lod3Mode == 2 && theta1 == theta2 && isMotherInformative(m1) && isFatherInformative(m1)) {
			double ret2 = 1.0;			

			double logRet2 = 0.0;
			for (int ci = 0; ci < tableSize; ci+=4)
			{
				double sum1 = 2 * f1[ci + 0];
				double sum2 = f1[ci + 1] + f1[ci + 2]; // TODO:Check
				double sum3 = 2 * f1[ci + 3];
				double sum = sum1 + sum2 + sum3; 
				
				double max = Math.max(Math.max(sum1, sum2), sum3) / sum;
				
				ret2 *= lNoLinkage3 * (theta1 * (1 - max) + (1 - theta1) * max);

				if ((ci + 4 >= tableSize) || ret2 > MAX_VALUE || ret2 < MIN_VALUE ) {
					logRet2 += Math.log10(ret2);
					ret2 = 1.0;
				}
			
			}
			return Math.max(logRet, logRet2); 
		}
		if (lod3Mode == 3 && theta1 == theta2 && isMotherInformative(m1) && isFatherInformative(m1)) {
			double ret2 = 1.0;
			double logRet2 = 0.0;
			for (int ci = 0; ci < tableSize; ci+=4)
			{
				double sum1 = f1[ci + 0] + f1[ci + 3];
				double max = Math.max(sum1, 1 - sum1);
				ret2 *= lNoLinkage4 * (theta1 * (1 - max) + (1 - theta1) * max);
				if ((ci + 4 >= tableSize) ||  ret2 > MAX_VALUE || ret2 < MIN_VALUE) {
					logRet2 += Math.log10(ret2);
					ret2 = 1.0;
				}
			}
			return Math.max(logRet, logRet2); 
		}
		return logRet;
	}

	//TODO: faster by evaluating maximum LOD score obtained 
	private double computeLOD(double theta1, double theta2, int m1, int m2, int phase)
	{
		float f1[] = getProb(m1);
		float f2[] = getProb(m2);
		
		int tableSize = 4 * numChildren;
		int map[] = mapping[phase];
	
		double ret = 1.0;
		double logRet = 0.0;
		
		double lNoLinkage1 = 2.0;
		double lNoLinkage2 = 2.0;
		if (distortionLodMode) { // distortion aware lod score 
			double p = freqs[m1][0]+freqs[m1][1];
			double q = freqs[m2][map[0]]+freqs[m2][map[1]];
			lNoLinkage1 = p * q + (1 - p) * (1 - q);
			if (lNoLinkage1 <= 0.5)
				lNoLinkage1 = 2.0;
			else
				lNoLinkage1 = 1.0 / lNoLinkage1;
			
			p = freqs[m1][0]      + freqs[m1][2];
			q = freqs[m2][map[0]] + freqs[m2][map[2]];
			lNoLinkage2 = p * q + (1 - p) * (1 - q);
			if (lNoLinkage2 <= 0.5)
				lNoLinkage2 = 2.0;
			else
				lNoLinkage2 = 1.0 / lNoLinkage2; 
		}
		
		for (int ci = 0; ci < tableSize; ci+=4)
		{
			double sum1 = f1[ci + 0] + f1[ci + 1];
			double sum2 = f2[ci + map[0]] + f2[ci + map[1]];
			double pat = sum1 * sum2 + (1 - sum1) * (1 - sum2);
			
			sum1 = f1[ci + 0] + f1[ci + 2];
			sum2 = f2[ci + map[0]] + f2[ci + map[2]];
			double mat = sum1 * sum2 + (1 - sum1) * (1 - sum2);
						 
			ret *= lNoLinkage1 * (theta1 * (1 - pat) + (1 - theta1) * pat);
			ret *= lNoLinkage2 * (theta2 * (1 - mat) + (1 - theta2) * mat);
			if ((ci + 4 >= tableSize) ||  ret > MAX_VALUE || ret < MIN_VALUE) {
				logRet += Math.log10(ret);
				ret = 1.0;
			}
		}
		//TODO: support for larger and smaller LOD scores... (maybe by taking log at every K:th individual)
		return logRet;
	}
	private double computeLOD_1(double theta1, int m1, int m2, int phase)
	{
		float f1[] = getProb(m1);
		float f2[] = getProb(m2);
		
		int tableSize = 4 * numChildren;
		int map[] = mapping[phase];
		
		double lNoLinkage1 = 2.0;
		if (distortionLodMode) { // distortion aware lod score
			double p = freqs[m1][0]+freqs[m1][1];
			double q = freqs[m2][map[0]]+freqs[m2][map[1]];
			lNoLinkage1 = p * q + (1 - p) * (1 - q);
			if (lNoLinkage1 <= 0.5)
				lNoLinkage1 = 2.0;
			else
				lNoLinkage1 = 1.0 / lNoLinkage1;
		}
	
		double ret = 1.0;
		double logRet = 0.0;
		for (int ci = 0; ci < tableSize; ci+=4)
		{
			double sum1 = f1[ci + 0] + f1[ci + 1];
			double sum2 = f2[ci + map[0]] + f2[ci + map[1]];
			double pat = sum1 * sum2 + (1 - sum1) * (1 - sum2);
			
			ret *= lNoLinkage1 * (theta1 * (1 - pat) + (1 - theta1) * pat);
			if ((ci + 4 >= tableSize) ||  ret > MAX_VALUE || ret < MIN_VALUE ) {
				logRet += Math.log10(ret);
				ret = 1.0;
			}
		}
		//TODO: support for larger and smaller LOD scores... (maybe by taking log at every K:th individual)
		return logRet;
	}	
	
	private double computeLOD_2(double theta2, int m1, int m2, int phase)
	{
		float f1[] = getProb(m1);
		float f2[] = getProb(m2);
		
		int tableSize = 4 * numChildren;
		int map[] = mapping[phase];
		
		double lNoLinkage2 = 2.0;
		if (distortionLodMode) { // distortion aware lod score
			double p = freqs[m1][0]+freqs[m1][2];
			double q = freqs[m2][map[0]]+freqs[m2][map[2]];
			lNoLinkage2 = p * q + (1 - p) * (1 - q);
			if (lNoLinkage2 <= 0.5)
				lNoLinkage2 = 2.0;
			else
				lNoLinkage2 = 1.0 / lNoLinkage2; 
		}
	
		double ret = 1.0;
		double logRet = 0.0;
		for (int ci = 0; ci < tableSize; ci+=4)
		{
			double sum1 = f1[ci + 0] + f1[ci + 2];
			double sum2 = f2[ci + map[0]] + f2[ci + map[2]];
			double mat = sum1 * sum2 + (1 - sum1) * (1 - sum2);
						 
			ret *= lNoLinkage2 * (theta2 * (1 - mat) + (1 - theta2) * mat);
			if ((ci + 4 >= tableSize) ||  ret > MAX_VALUE || ret < MIN_VALUE  ) {
				logRet += Math.log10(ret);
				ret = 1.0;
			}
		}
		//TODO: support for larger and smaller LOD scores... (maybe by taking log at every K:th individual)
		return logRet;
	}

	// LOD score for double informative markers, 
	// mode 2
	
	private final int mapping3[][] = {{0,3}, {3,0}};
	
	//+log(3) bit for identical genotype
	private double computeLOD2(double theta, int m1, int m2, int phase)
	{
		float f1[] = getProb(m1);
		float f2[] = getProb(m2);
		
		int tableSize = 4 * numChildren;
		int map[] = mapping3[phase];

		double lNoLinkage3 = 3.0;
		if (distortionLodMode) {
			double p1 = 2.0 * freqs[m1][0];
			double p2 = 2.0 * freqs[m1][3];
			double p3 = freqs[m1][1] + freqs[m1][2];
			double sump = p1 + p2 + p3;

			double q1 = 2.0 * freqs[m1][map[0]];
			double q2 = 2.0 * freqs[m1][map[1]];
			double q3 = freqs[m1][1] + freqs[m1][2];
			double sumq = q1 + q2 + q3;
			
			lNoLinkage3 = (p1 * q1 + p2 * q2 + p3 * q3) / (sump * sumq);
			if (lNoLinkage3 <= Constants.ONETHIRD)
				lNoLinkage3 = 3.0;
			else 
				lNoLinkage3 = 1.0 / lNoLinkage3; 
		}
		
	
		double ret = 1.0;
		double logRet = 0.0;
		for (int ci = 0; ci < tableSize; ci+=4)
		{
			double sum11 = 2 * f1[ci + 0];
			double sum21 = 2 * f2[ci + map[0]];

			double sum12 = 2 * f1[ci + 3];
			double sum22 = 2 * f2[ci + map[1]];
			
			double sum13 = f1[ci + 1] + f1[ci + 2]; 
			double sum23 = f2[ci + 1] + f2[ci + 2];
			
			double sum1 = sum11 + sum12 + sum13; 
			double sum2 = sum21 + sum22 + sum23;

			double s = 1.0 / (sum1 * sum2);			
			
			double t0 = (sum11 * sum21) * s; 
			double t1 = (sum12 * sum22) * s;
			double t2 = (sum13 * sum23) * s;
			
			double t = t0 + t1 + t2;

			ret *= lNoLinkage3 * (theta * (1 - t) + (1 - theta) * t);
			if ((ci + 4 >= tableSize) ||  ret > MAX_VALUE || ret < MIN_VALUE  ) {
				logRet += Math.log10(ret);
				ret = 1.0;
			}
		}
		//TODO: support for larger and smaller LOD scores... (maybe by taking log at every K:th individual)
		return logRet;
	}	

	// LOD score for double informative markers, 
	// mode 3
	//+1 bit for homozygote-heterozygote... (phase=1)
	//+1 bit for homozygote-homozygote... (phase=0)
	private double computeLOD3(double theta, int m1, int m2, int phase)
	{
		float f1[] = getProb(m1);
		float f2[] = getProb(m2);
		
		int tableSize = 4 * numChildren;
	
		double lNoLinkage4 = 2.0;
		if (distortionLodMode) {
			double p = freqs[m1][0]+freqs[m1][3];
			double q = (phase == 1) ? freqs[m2][1] + freqs[m2][2] : freqs[m2][0] + freqs[m2][3];
			lNoLinkage4 = p * q + (1 - p) * (1 - q);
			if (lNoLinkage4 <= 0.5)
				lNoLinkage4 = 2.0;
			else
				lNoLinkage4 = 1.0 / lNoLinkage4;
		}
		
		double ret = 1.0;
		double logRet = 0.0;
		for (int ci = 0; ci < tableSize; ci+=4)
		{
			double sum11 = f1[ci + 0] + f1[ci + 3];
			double sum21 = (phase == 1) ? f2[ci + 1] + f2[ci + 2] : f2[ci + 0] + f2[ci + 3];
			double sum12 = 1 - sum11; 
			double sum22 = 1 - sum21;

			double t0 = sum11 * sum21;
			double t1 = sum12 * sum22;
			double t = t0 + t1;
			ret *= lNoLinkage4 * (theta * (1 - t) + (1 - theta) * t);
			if ((ci + 4 >= tableSize) ||  ret > MAX_VALUE || ret < MIN_VALUE ) {
				logRet += Math.log10(ret);
				ret = 1.0;
			}
		}
		//TODO: support for larger and smaller LOD scores... (maybe by taking log at every K:th individual)
		return logRet;
	}
	
	//private final double MINPROB_COMBINE = 0.001;
	
	// from m2 to m1 
	public void setProb(int m1, int m2)
	{
/*		float f1[] = haplotypeProb.get(m1);
		float f2[] = haplotypeProb.get(m2);
		int tableSize = 4 * numChildren;
		for (int i = 0; i < tableSize; ++i)
			f1[i] = f2[i];*/
		haplotypeProb.set(m1, haplotypeProb.get(m2));
	}
	public void setProb(int m1, float f[])
	{
/*		float f1[] = haplotypeProb.get(m1);
		float f2[] = haplotypeProb.get(m2);
		int tableSize = 4 * numChildren;
		for (int i = 0; i < tableSize; ++i)
			f1[i] = f2[i];*/
		haplotypeProb.set(m1, f);
	}
	
	

	private final double NORMALIZE_SMALL = 0.25e-20; 
	public void combineAndNormalize(int m1, int m2, double theta1, double theta2)
	{
		if (theta1 != 0 && theta2 != 0)
			return;

		int tableSize = 4 * numChildren;
		float f1[] = getProb(m1);
		if (f1 == null) {
			f1 = new float[4 * numChildren];
			Arrays.fill(f1, 0.25f);
			haplotypeProb.add(m1, f1);
		}
			
		float f2[] = getProb(m2);
		if (f2 == null) {
			f2 = new float[4 * numChildren];
			Arrays.fill(f2, 0.25f);
			haplotypeProb.add(m2, f2);
		}

		
		if (theta1 == theta2 && lod3Mode == 2 && isMotherInformative(m1) && isMotherInformative(m2) && isFatherInformative(m1) && isFatherInformative(m2)) {
			double table[] = {
					computeLOD2(theta1, m1, m2, 0),
					computeLOD2(theta1, m1, m2, 1)};
			int max = table[0] > table[1] ? 0 : 1;
			if (max == 0)
				for (int i = 0; i < tableSize; i+=4) {
					double sum = 4 * NORMALIZE_SMALL;
					for (int j = 0; j < 4; ++j) {
						sum += f1[i + j] * f2[i + j];
				}
				for (int j = 0; j < 4; ++j)
					f1[i + j] = (float)((f1[i + j] * (double) f2[i + j] + NORMALIZE_SMALL) / sum);
			} else {
				for (int i = 0; i < tableSize; i+=4) {
					double sum = 4 * NORMALIZE_SMALL;
					for (int j = 0; j < 4; ++j)
						sum += f1[i + j] * f2[i + ((j == 0 || j == 3) ? 3 - j : j)];
					for (int j = 0; j < 4; ++j)
						f1[i + j] = (float)((f1[i + j] * (double) f2[i + ((j == 0 || j == 3) ? 3 - j : j)] + NORMALIZE_SMALL) / sum);
				}
			}
			return;
		}
		
		if (theta1 == theta2 && lod3Mode == 3 && isMotherInformative(m1) && isMotherInformative(m2) && isFatherInformative(m1) && isFatherInformative(m2)) {
			double table[] = {
					computeLOD3(theta1, m1, m2, 0),
					computeLOD3(theta1, m1, m2, 1)};
			
			int max = table[0] >= table[1] ? 0 : 1;

			if (max == 0)
				for (int i = 0; i < tableSize; i+=4) {
					double sum = 4 * NORMALIZE_SMALL;
					for (int j = 0; j < 4; ++j) {
						sum += f1[i + j] * f2[i + j];
				}
				for (int j = 0; j < 4; ++j)
					f1[i + j] = (float)((f1[i + j] * (double) f2[i + j] + NORMALIZE_SMALL) / sum);
			} else {
				for (int i = 0; i < tableSize; i+=4) {
					double hom1 =  f1[i + 0] + f1[i + 3];
					double hom2 =  f2[i + 0] + f2[i + 3];
					
					double sum = hom1 * (1 - hom2) + (1 - hom1) * hom2 + 4 * NORMALIZE_SMALL; 
					
					f1[i + 0] = (float)((f1[i + 0] * (1 - hom2) + NORMALIZE_SMALL) / sum);
					f1[i + 1] = (float)((f1[i + 1] * (hom2) + NORMALIZE_SMALL) / sum);
					f1[i + 2] = (float)((f1[i + 2] * (hom2) + NORMALIZE_SMALL) / sum);
					f1[i + 3] = (float)((f1[i + 3] * (1 - hom2) + NORMALIZE_SMALL) / sum);
				}
				
			}
			
			return;
		}
		
		
		double table[] = {
				computeLOD(theta1, theta2, m1, m2, 0),
				computeLOD(theta1, theta2, m1, m2, 1),
				computeLOD(theta1, theta2, m1, m2, 2),
				computeLOD(theta1, theta2, m1, m2, 3)};
		
		int maxi = 0;
		for (int i = 1; i < 4; ++i)
			if (table[i] > table[maxi])
				maxi = i;
		
		int map[] = mapping[maxi];
		
		
		
		if (theta1 == 0 && theta2 == 0) 
			for (int i = 0; i < tableSize; i+=4) {
				
				double sum = 4 * NORMALIZE_SMALL;
				for (int j = 0; j < 4; ++j)
					sum += f1[i + j] * (double) f2[i + map[j]];
				
				for (int j = 0; j < 4; ++j) {
					//if (p < MINPROB_COMBINE)
					//	p = MINPROB_COMBINE;
					f1[i + j] = (float)((f1[i + j] * (double) f2[i + map[j]] + NORMALIZE_SMALL) / sum);
				}
			}
		
		else if (theta1 == 0) {
			for (int i = 0; i < tableSize; i+=4) {
				double p = f2[i + map[0]];
				p += f2[i + map[1]];
				double ret[] = new double[4]; 
				//if (p < MINPROB_COMBINE)
				//	p = MINPROB_COMBINE;
				for (int j = 0; j < 4; ++j)
					ret[j] = f1[i + j] * ((j < 2) ? p : 1 - p);

				double sum = 4 * NORMALIZE_SMALL;
				for (int j = 0; j < 4; ++j)
					sum += ret[j];
				for (int j = 0; j < 4; ++j)
					f1[i + j] = (float)((ret[j] + NORMALIZE_SMALL) / sum); 
			}
		}
		else if (theta2 == 0) {
			for (int i = 0; i < tableSize; i+=4) {
				double p = f2[i + map[0]];
				p += f2[i + map[2]];
				double ret[] = new double[4]; 
				//if (p < MINPROB_COMBINE)
				//	p = MINPROB_COMBINE;
				for (int j = 0; j < 4; ++j)
					ret[j] = f1[i + j] * ((j == 0 || j == 2) ? p : 1 - p);

				double sum = 4 * NORMALIZE_SMALL;
				for (int j = 0; j < 4; ++j)
					sum += ret[j];
				for (int j = 0; j < 4; ++j)
					f1[i + j] = (float)((ret[j] + NORMALIZE_SMALL) / sum); 
			}
		}
		
		//normalize(f1);
		
	}

	
	public String pattern(int marker)
	{
		int tableSize = 4 * numChildren;
		
		final double LIMIT = 0.1;
		
		float f[] = getProb(marker);

		String s1 = "";
		String s2 = "";
	
		if (f == null)
			for (int i = 0; i < tableSize; i+=4) {
				s1 = s1 + "?";
				s2 = s2 + "?";
			}
		else if ((lod3Mode == 2 || lod3Mode == 3) && (isFatherInformative(marker) && isMotherInformative(marker))) {
			for (int i = 0; i < tableSize; i+=4) {
				double p0 = f[i + 0]; 
				double p1 = f[i + 1] + f[i + 2]; 
				double p2 = f[i + 3];
				if (p0 >= 1 - LIMIT)
					s1 = s1 + "0";
				else if (p1 >= 1 - LIMIT)
					s1 = s1 + "1";
				else if (p2 >= 1 - LIMIT)
					s1 = s1 + "2";
				else 
					s1 = s1 + "?";
			}
			s2 = s1;
			
		} else
			for (int i = 0; i < tableSize; i+=4) {
				double pat = f[i + 0] + f[i + 1]; 
				double mat = f[i + 0] + f[i + 2];
				if (pat >= 1 - LIMIT)
					s1 = s1 + "0";
				else if (pat <= LIMIT)
					s1 = s1 + "1";
				else
					s1 = s1 + "?";
		
				if (mat >= 1 - LIMIT)
					s2 = s2 + "0";
				else if (mat <= LIMIT)
					s2 =  s2 + "1";
				else
					s2 = s2 + "?";
			}
		
		return s1 + "\t" + s2;
	}

	private void printPattern(int marker)
	{
		System.err.println(pattern(marker));
	}

	public void removeMarker(int marker)
	{
		motherInf.set(marker, motherInf.get(marker) + 2);
		fatherInf.set(marker, fatherInf.get(marker) + 2);
		haplotypeProb.set(marker, null);
		
	}
	
	

}
