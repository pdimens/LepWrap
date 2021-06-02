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

	Copyright (C) 2020 Pasi Rastas, pasi.rastas@helsinki.fi, University of Helsinki	
*/
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

//TODO: vcf output, indels, output depth, ploidy

public class Pileup2Likelihoods {
	
	int ploidy = 2;
	double minQ = 0.001;
	boolean callIndels = false;
	
	float logP2[][][] = new float[3][94][94]; //33-126
	ArrayList<int []> numAlleles = null;
	ArrayList<int []> numAllelesT = null; // store array in transpose as well
	
	
	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
		logP2 = new float[ploidy + 1][94][94];
	}

	public void setCallIndels(boolean value) {
		callIndels = value;
	}
	
	public void setMinQ(double value)
	{
		minQ = value;
	}
	
	private final double logTerm = -0.1 * Constants.LN10;
	
	private double quality(int q) {
		if (q == 0)
			return 0.75; // all 4 bases equal prob
		if (q == 1)
			return 0.725; // between q=1 and q=2
		return Math.exp(logTerm * q);
	}

	private double combineQ(int q1, int q2) {
        double p1 = quality(q1);
        double p2 = quality(q2);
        double ret = p1 + p2 - p1 * p2;
        if (ret > 0.75)
                ret = 0.75;
        return ret;
	}

	public void initTables()
	{
		double iploidy = 1.0 / ploidy; 
		for (int i = 0; i < 94; ++i) {
			boolean fill = false;
			for (int j = i; j < 94; ++j) {
				if (fill) {
					for (int pl = 0; pl <= ploidy; ++pl)
						logP2[pl][j][i] = logP2[pl][i][j] = logP2[pl][i][j - 1]; 
				} else {
	                double q = combineQ(i, j);
	                if (q < minQ) {
	                        q = minQ;
	                        fill = true;
	                }
	                double q3 = q / 3.0;
	                
					for (int pl = 0; pl <= ploidy; ++pl)
						logP2[pl][j][i] = logP2[pl][i][j] = (float) Math.log((pl * (1 - q) + (ploidy - pl) * q3) * iploidy);

					//logP2[0]
					//logP2[1]
					//logP2[2]

					//logP[j][i] = logP[i][j] = (float) Math.log(q / 3);
	                //logNP2[j][i] = logNP2[i][j] = (float) (Math.log(q / 3 + (1 - q)) + Constants.LNHALF);
	                //logNP[j][i] = logNP[i][j] = (float) Math.log(1 - q); 
				}
			}

        }
		int allele[] = new int[ploidy]; 
		//for (int pl = 0; pl < ploidy; ++pl)
        //    allele[pl] = 1;
		numAlleles = new ArrayList<int []>();
		numAlleles.add(new int[]{ploidy,0,0,0}); //A..A
		int k = ploidy - 1;
        while (k >= 0) {
            if (allele[k] < 3) {
                ++allele[k];
                for (int l = k + 1; l < ploidy; ++l)
                        allele[l] = allele[k];
                int na[] = new int[4];
                for (int l = 0; l < ploidy; ++l)
                        ++na[allele[l]];
                numAlleles.add(na);
                k = ploidy - 1;
            }
            else
                --k;
        }
        //print genotypes
        int g = 0;
        for (int[] list : numAlleles) {
        	String alleles = "ACGT";
        	for (int i = 0; i < 4; ++i)
        		for (int j = 0; j < list[i]; ++j)
        			System.err.print(alleles.charAt(i));
        	System.err.println("\t" + (g + 1));
        	++g;
        }
        numAllelesT = new ArrayList<int []>(); // transpose
        for (int i = 0; i < 4; ++i) {
        	int na[] = new int[numAlleles.size()];
        	int j = 0;
            for (int[] list : numAlleles)
            	 na[j++] = list[i];
        	numAllelesT.add(na);
        }
        

        System.err.println("number of different genotypes is " + numAlleles.size());
	}
	
	public Pileup2Likelihoods() {
	}
	int numIndividuals = 0;

	private int count[];
	private int count2[] = new int[4];

	private int intMapping[];
	
	private int coveragePerIndividual = 0;
	private int numLowerCoverage = 0;
	private int minimumAlleleFreq = 0;
	private int minCoverageSum = 0;
	
	private void setMinCoverageSum(int value) {
		if (value < 0)
			minCoverageSum = (numIndividuals - numLowerCoverage) * coveragePerIndividual;   
		else
			minCoverageSum = value;
	}
	
	private void setCoveragePerIndividual(int value) {
		coveragePerIndividual = value;
	}

	private void setMinimumAlleleFreq(int value) {
		minimumAlleleFreq = value;
	}

	private void setMinimumAlleleFreqRatio(double value) {
		minimumAlleleFreq = (int)(0.5 + value * numIndividuals);
	}

	private void setNumLowerCoverage(int value) {
		numLowerCoverage = value;
	}

	private void setNumLowerCoverageRatio(double value) {
		numLowerCoverage = (int) (0.5 + numIndividuals * value);
	}

	ArrayList<String> individualNames = new ArrayList<String>();

    private void processMapping(String mappingFile)
    {
    	ArrayList<ArrayList<String>> it = Input.loadTable(mappingFile, "\t ");
    	ArrayList<String> values = new ArrayList<String>();
    	
    	if (it == null || it.size() == 0)
    		Error.error(3002);

    	int n = 0;
		for (ArrayList<String> row : it) {
			if (row.size() == 0 || (n > 0 && row.size() != n))
				Error.error(3002);
			values.addAll(row);
			n = row.size();
		}
		if (it.size() > 1 && n > 1)
			Error.error(3002);
    	
    	
    	HashMap<String, Integer> hm = new HashMap<String, Integer>();
    	    	
    	numIndividuals = 0;
    	for (String m : values)
    		if (!hm.containsKey(m)) {
    			hm.put(m, numIndividuals++);
    			individualNames.add(m);
    		}
    	
    	intMapping = new int[values.size()];
    	int i = 0;
    	for (String m : values)
    		intMapping[i++] = hm.get(m);
    	
		count = new int[numIndividuals];
    	
    	System.err.println("Number of bams is " + values.size());
    	System.err.println("Number of individuals is " + numIndividuals);
    }
    private void processPileup(String filename)
    {
		System.out.print("CHR\tPOS");
		for (int i = 0; i < numIndividuals; ++i) {
			System.out.print("\t" + individualNames.get(i));
		}
		System.out.println();
    	
		try {
			BufferedReader br = null;
			if (filename.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(filename));
			
			ArrayList<String> line = Input.loadTableRow(br, "\t");
			while (line != null) {
				processOneRow(line);
				line = Input.loadTableRow(br, "\t");
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
    }

    // save storage space (small=>0, 1.0=>1) and slightly faster
    private String myExp(double value){
    	if (value == 0.0)
    		return "1";
    	else if (value < -46.0)
    		return "0";
    	else
    		return "" + (float) Math.exp(value);
    }
    
    
    private void processOneRow(ArrayList<String> row)
    {
    	//ArrayList<String> row2 = new ArrayList<String>();
    	//row2.addAll(row);
    	ArrayList<ArrayList<Integer>> indels = null;
    	if (callIndels) {
    		indels = new ArrayList<ArrayList<Integer>>();
    		for (int i = 0; i < numIndividuals; ++i)
    			indels.add(new ArrayList<Integer>());
		}
    	
     	int sum = 0;
		Arrays.fill(count, 0);
		
		int mIndex = 0;
		for (int i = 3; i < row.size(); i+=4) {
			int c = Integer.parseInt(row.get(i));
			sum += c;
			count[intMapping[mIndex++]] += c;
			if (c == 0){ // remove "*" if no reads
				while (row.size() <= i + 3)
					row.add("");
				row.set(i + 1, "");
				row.set(i + 2, "");
				row.set(i + 3, "");
			}
		}
		if (sum < minCoverageSum)
			return;
		
		int missing = 0;
		for (int i = 0; i < numIndividuals; ++i)
			if (count[i] < coveragePerIndividual)
				++missing;
		
		if (missing > numLowerCoverage)
			return;
		//filtering...
		
		Arrays.fill(count2, 0);
		for (int i = 4; i < row.size(); i+=4) {
			char ca[] = row.get(i).toCharArray();
			int caIndex = 0;
			
			int n = ca.length;
			for (int j = 0; j < n; ++j) {
				char c = ca[j];
				switch (c) {
				case '*': 
						ca[caIndex++]=4; // missing
						break;
				case '+': 
				case '-': 
					int l = 0; 
					++j; 
					while (j < n && ca[j] >= '0' && ca[j] <= '9') {
						l = 10 * l + (ca[j] - '0');
						++j;
					}
					j += (l - 1); // correct for last j++
					if (callIndels) {
						if (c == '-')
							indels.get(intMapping[(i >> 2) - 1]).add(-l);
						else
							indels.get(intMapping[(i >> 2) - 1]).add(l);
					}
					break;
				case '^': ++j; 
				case '$':
					break;

				case 'a':
				case 'A': ++count2[0];ca[caIndex++]=0;
					break;
				case 'c':
				case 'C': ++count2[1];ca[caIndex++]=1;
					break;
				case 'g':
				case 'G': ++count2[2];ca[caIndex++]=2;
					break;
				case 't':
				case 'T': ++count2[3];ca[caIndex++]=3;
					break;
				}
			}
			row.set(i, new String(ca, 0, caIndex));
		}
		int alleles = 0;
		for (int i = 0; i < 4; ++i)
			if (count2[i] >= minimumAlleleFreq)
				++alleles;
		
		//System.err.println(rowNumber + " " + row[0]);
		
		if (alleles >= 2) {
			int numGenotypes = numAlleles.size();
			double prob[][] = new double[numIndividuals][numGenotypes];
			
			mIndex = 0;

			for (int i = 4; i < row.size(); i+=4) {
				String s = row.get(i);
				String q1 = (row.size() > i + 1) ? row.get(i + 1) : "";
				String q2 = (row.size() > i + 2) ? row.get(i + 2) : "";

				//System.err.println(s);
				//System.err.println(q1);
				//System.err.println(q2);

				//assert(s.length() == q1.length() && q1.length() == q2.length());
				if(s.length() != q1.length() || q1.length() != q2.length()) { //quality and aligment have different number of values...
					//System.err.println(row2.get(i));
					
					for (int j = 0; j < s.length(); ++j) 
							System.err.print("acgtN".charAt(s.charAt(j)));
					System.err.println();
					System.err.println(q1);
					System.err.println(q2);
					Error.warning(3003);
					int l = Math.min(Math.min(s.length(), q1.length()), q2.length());
					s = s.substring(0, l);
					q1 = q1.substring(0, l);
					q2 = q2.substring(0, l);
				}
				int n = s.length();
				int ind = intMapping[mIndex++];
				for (int j = 0; j < n; ++j) {
					char c = s.charAt(j);
					if (c < 4) {
						char qc1 = q1.charAt(j);
						char qc2 = q2.charAt(j);
						int nat[] = numAllelesT.get(c);
	                    for (int g = 0; g < numGenotypes; ++g)
                            prob[ind][g] += logP2[nat[g]][qc1-33][qc2-33];					
                    }
				}
				// get qualities and output genotype likelihoods...
			}
			StringBuilder sb = new StringBuilder();
			sb.append(row.get(0));
			sb.append('\t');
			sb.append(row.get(1));
			for (int ind = 0; ind < numIndividuals; ++ind) {
				double maxp = Double.NEGATIVE_INFINITY;
                for (int g = 0; g < numGenotypes; ++g)
                    maxp = Math.max(maxp, prob[ind][g]);
                for (int g = 0; g < numGenotypes; ++g) {
    				if (g == 0) 
    					sb.append('\t');
    				else
    					sb.append(' ');
    				sb.append(myExp(prob[ind][g] - maxp));
            	}
			}
			//print result
           System.out.println(sb);
		} else if (alleles == 1 && callIndels) {
			ArrayList<Integer> lengths = new ArrayList<Integer>(); 
			for (ArrayList<Integer> list : indels)
				lengths.addAll(list);
			if (lengths.size() > minimumAlleleFreq && lengths.size() >= 1) {
				Collections.sort(lengths);
				int numL = lengths.size();
				int cid = 1;
				int id = 0;
				for (int i = 1; i < numL; ++i) {
					if (lengths.get(i - 1) == lengths.get(i)) {
						++cid;
					} else {
						if (cid >= minimumAlleleFreq) {
							++alleles;
							id = lengths.get(i - 1);
						}
						cid = 1;
					}
				}
				if (cid >= minimumAlleleFreq) {
					++alleles;
					id = lengths.get(lengths.size() - 1);
				}
				if (alleles == 2) { // two values
					int numGenotypes = numAlleles.size();

					StringBuilder sb = new StringBuilder();
					sb.append(row.get(0));
					sb.append('\t');
					sb.append(row.get(1));
					
					for (int ind = 0; ind < numIndividuals; ++ind) {
						ArrayList<Integer> list = indels.get(ind);
						int c = 0;
						for (int l : list)
							if (l == id)
								++c;

							double prob[] = count2likelihood(count[ind] - c, c); 
							for (int g = 0; g < numGenotypes; ++g) {
			    				if (g == 0) 
			    					sb.append('\t');
			    				else
			    					sb.append(' ');
			    				if (ind == 0)
			    					sb.append((float) Math.exp(prob[g]));
			    				else
			    					sb.append(myExp(prob[g]));
			            	}
					}
					System.out.println(sb);
				}
			}
		}
    }
    private double[] count2likelihood(int c1, int c2){
    	int numGenotypes = numAlleles.size();

    	double prob[] = new double[numGenotypes];
    	
    	int nat[] = numAllelesT.get(0);
        for (int g = 0; g < numGenotypes; ++g)
            prob[g] += c1 * logP2[nat[g]][20][20];
        
    	nat = numAllelesT.get(1);
        for (int g = 0; g < numGenotypes; ++g)
            prob[g] += c2 * logP2[nat[g]][20][20];
        
		double maxp = Double.NEGATIVE_INFINITY;
        for (int g = 0; g < numGenotypes; ++g)
            maxp = Math.max(maxp, prob[g]);

        for (int g = 0; g < numGenotypes; ++g)
            prob[g] = prob[g] - maxp;
    	
    	return prob;
    }
	
	public static void main(String args[]){ 	
		ParameterParser pp = new ParameterParser();
		StringBuilder extraParameters = new StringBuilder();
		for (int i = 0; i < args.length; ++i) {
			extraParameters.append(args[i] + " ");
		}
		if (!pp.init(extraParameters.toString()))
			Error.error(3001);
		pp.warning(new String[]{"minCoverage", "numLowerCoverage", "minAlleleFreq", "minQuality", "minCoverageSum", "mappingFile", "pileup", "ploidy", "callIndels"});
		
		Pileup2Likelihoods p2l = new Pileup2Likelihoods();
		int ploidy = Integer.parseInt(pp.getValueAsString("ploidy", "2"));
		if (ploidy != 2)
			p2l.setPloidy(ploidy);
		
		int cov = Integer.parseInt(pp.getValueAsString("minCoverage", "3"));
		
		double numLow = Double.parseDouble(pp.getValueAsString("numLowerCoverage", "0.3"));
		
		double minAlleleFreq = Double.parseDouble(pp.getValueAsString("minAlleleFreq", "0.1"));
		//double minSum = Double.parseDouble(pp.getValueAsString("minCoverageSum", "0"));

		double minQ = Double.parseDouble(pp.getValueAsString("minQuality", "0.001"));

		p2l.processMapping(pp.getValueAsString("mappingFile", "mapping.txt"));
		
		p2l.setCoveragePerIndividual(cov);
		if (numLow >= 1)
			p2l.setNumLowerCoverage((int) (0.5 + numLow));
		else
			p2l.setNumLowerCoverageRatio(numLow);
		if (minAlleleFreq >= 1)
			p2l.setMinimumAlleleFreq((int)(0.5 + minAlleleFreq));
		else
			p2l.setMinimumAlleleFreqRatio(minAlleleFreq);
		
		p2l.setMinQ(minQ);
	

		p2l.setMinCoverageSum(Integer.parseInt(pp.getValueAsString("minCoverageSum", "-1")));

		if (pp.getValueAsString("callIndels", "0").equals("1"))
			p2l.setCallIndels(true);
		
		p2l.initTables();
		p2l.processPileup(pp.getValueAsString("pileup", "-"));
		
	}

}
