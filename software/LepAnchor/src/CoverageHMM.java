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
//User interface and implementation of ScaffoldHMM  
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

public class CoverageHMM {

	private int column = 2;
	private int numMixture = 0;
	private int MAX_ITERATIONS = 1000;
	private int numScaffolds;
	
	private double minProb = 0.001;
	private double scale = 0.001;
	private double sample = 0.1;
	
	private ArrayList<String> scaffolds = new ArrayList<String>();
	
	public void setColumn(int value) 
	{
		column = value;
	}

	public void setMinProb(double value) 
	{
		minProb = value;
	}

	public void setSample(double value) 
	{
		sample = value;
	}

	public void setScale(double value) 
	{
		scale = value;
	}

	private ArrayList<double[]> mixture = new ArrayList<double[]>();
	private HashMap<String, ArrayList<Integer>> scaffoldMap = new HashMap<String, ArrayList<Integer>>();
	
	public void loadMixture(String mixtureFile)
	{
		numMixture = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(mixtureFile));
			
			ArrayList<String> line = Input.loadTableRow(br, " \t");
			while (line != null) {
				if (numMixture == 0)
					numMixture = line.size() - 3;
				else
					if (numMixture != line.size() - 3) {
						System.err.println("Error: Different number of columns in the input");
						System.exit(-1);
					}
				int value = Integer.parseInt(line.get(0));
				if (mixture.size() != value) {
					System.err.println("Error: Incorrect mixture file");
					System.exit(-1);
				}
				double m[] = new double[numMixture];
				for (int i = 0; i < numMixture; ++i)
					m[i] = Double.parseDouble(line.get(i + 3));
				mixture.add(m);
				line = Input.loadTableRow(br, " \t");
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
	public void loadDepth(String depthFile)
	{
		try {
			BufferedReader br = null;
			if (depthFile.equals("-"))	
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(depthFile));
			
			ArrayList<String> line = Input.loadTableRow(br, " \t");
			double pick = 0.0;
			while (line != null) {
				if (line.size() <= column) {
					System.err.println("Warning: skipping " + line); 
					continue; 
				}
				pick += sample;
				if (pick >= 1.0) {
					pick -= 1.0;
					String scaffold = line.get(0);
					long position = Long.parseLong(line.get(1));
					if (position > 0xffffffffl) {
						System.err.println("Error: Too long scaffold");
						System.exit(-1);
					}
					int cov = Integer.parseInt(line.get(column));
	
					if (!scaffoldMap.containsKey(scaffold))
						scaffoldMap.put(scaffold, new ArrayList<Integer>());
					ArrayList<Integer> list = scaffoldMap.get(scaffold);
					list.add((int) position);
					list.add(cov);
				}
				line = Input.loadTableRow(br, " \t");
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		numScaffolds = scaffoldMap.keySet().size();
//		double markerDensityScales[][] = new double[numScaffolds][]; 
		scaffolds = new ArrayList<String>();
		for (String scaffold : scaffoldMap.keySet())
			scaffolds.add(scaffold);

		for (String scaffold : scaffolds) {
			ArrayList<Integer> positionsAndCov = scaffoldMap.get(scaffold); 
			int numM = positionsAndCov.size();
			for (int i = 2; i < numM; i+=2) {
				if (positionsAndCov.get(i - 2) > positionsAndCov.get(i)) {
					System.err.println("physical positions must be increasing...");
					System.exit(0);
				}
			}
		}
		System.err.println("Number of scaffolds = " + numScaffolds);
	}
	
	public void hmm()
	{
		DecimalFormat df = new DecimalFormat("#0.000"); 
		
		for (double m[] : mixture)
			for (int i = 0; i < numMixture; ++i) {
				double v = m[i];
				if (v < minProb)
					m[i] = minProb;
				m[i] = Math.pow(m[i], scale);
			}
		
		//Init transition
		double t[][] = new double[numMixture][numMixture];
	
		for (double tmp[] : t)
			Arrays.fill(tmp, 1.0 / numMixture);
		
		for (int i = 0; i < numMixture; ++i) {
			double sum = 0.0;
			for (int j = 0; j < numMixture; ++j) {
				double r = (i == j) ? 1.0 : 0.1 * scale * Math.random();
				sum += r;
				t[i][j] = r;
			}
			double iSum = 1.0 / sum;
			for (int j = 0; j < numMixture; ++j)
				t[i][j] *= iSum;
		}

		boolean converge = false;
		double oldLogL = Double.NEGATIVE_INFINITY;
		for (int iterations = 0; iterations <= MAX_ITERATIONS; ++iterations) {			
			int maxCount = mixture.size() - 1;
			double QT[][] = new double[numMixture][numMixture];
			double logL = 0.0;
			double l = 1.0;
			for (int scaffold=0; scaffold < numScaffolds; ++scaffold) {
				String scaffoldS = scaffolds.get(scaffold);
				ArrayList<Integer> positionsAndCov = scaffoldMap.get(scaffoldS);

				int numM = positionsAndCov.size() / 2;
				double forward[][] = new double[numM + 1][numMixture];
				double scale[] = new double[numM + 1];
				scale[0] = 1.0;
				
				double iM = 1.0 / numMixture;
				for (int i = 0; i < numMixture; ++i)
					forward[0][i] = iM;

				for (int p = 0; p < numM; ++p) {
					int cov = positionsAndCov.get(p + p + 1);
					if (cov >= maxCount)
						cov = maxCount;
					
					double m[] = mixture.get(cov); 
					
					for (int j = 0; j < numMixture; ++j) {
						double fp1j = 0.0;
						double f_[] = forward[p];
						for (int i = 0; i < numMixture; ++i)
							fp1j += f_[i] * t[i][j];
						forward[p + 1][j] = fp1j * m[j]; 
					}
					
					
					double sum = 0.0;
					for (int j = 0; j < numMixture; ++j)
						sum += forward[p + 1][j];
					double isum = 1.0 / sum;
					
					for (int j = 0; j < numMixture; ++j)
						forward[p + 1][j] *= isum;
					scale[p + 1] = isum;
					
					l *= sum;
					if (l < 1e-100) {
						logL += Math.log10(l); 
						l = 1.0;
					}
				}
				double backward[][] = new double[numM + 1][numMixture];
				for (int j = 0; j < numMixture; ++j)
					backward[numM][j] = 1.0;

				for (int p = numM - 1; p >= 0; --p) {
					int cov = positionsAndCov.get(p + p + 1);
					if (cov >= maxCount)
						cov = maxCount;
					
					double m[] = mixture.get(cov);

					
					for (int j = 0; j < numMixture; ++j) {
						double tmp = backward[p + 1][j] * m[j] * scale[p + 1];
						for (int i = 0; i < numMixture; ++i)
							backward[p][i] += t[i][j] * tmp;
					}
				}
				//double sum = 0.0;
				//for (int j = 0; j < numMixture; ++j)
				//	sum += backward[1][j] * forward[1][j];
				//System.err.println(sum);
				
				if (iterations == MAX_ITERATIONS || converge) {
					for (int p = 0; p < numM; ++p) {
						long pos = positionsAndCov.get(p + p) & 0xffffffffl;
						int cov = positionsAndCov.get(p + p + 1);
						double max = 0.0;
						double max2 = 0.0;
						int maxj = 0;
						for (int j = 0; j < numMixture; ++j) {
							double pj = forward[p][j] * backward[p][j];
							if (pj > max) {
								max2 = max;
								max = pj;
								maxj = j;
							} else if (pj > max2)
								max2 = pj;
						}
						StringBuilder sb = new StringBuilder();
						sb.append(scaffoldS);
						sb.append('\t');
						sb.append(pos);
						sb.append('\t');
						sb.append(cov);
						sb.append('\t');
						sb.append(maxj);
						sb.append('\t');
						sb.append(df.format(Math.log10(max) - Math.log10(max2)));
						System.out.println(sb.toString());
					}
					
				} else
				for (int p = 1; p < numM; ++p) {
					int cov = positionsAndCov.get(p + p + 1);
					if (cov >= maxCount)
						cov = maxCount;
					double m[] = mixture.get(cov);
					for (int j = 0; j < numMixture; ++j) {
						double tmp = scale[p + 1] * backward[p + 1][j] * m[j];
						for (int i = 0; i < numMixture; ++i)
							QT[i][j] += forward[p][i] * t[i][j] * tmp;
					}
				}
				
			}
			if (iterations == MAX_ITERATIONS || converge)
				break;
			
			for (int i = 0; i < numMixture; ++i) {
				double sum = 0.0;
				for (double q : QT[i])
					sum += q;
				double iSum = 1.0 / sum;
				for (int j = 0; j < numMixture; ++j)
					t[i][j] = QT[i][j] * iSum;
			}
			
			//for (int i = 0; i < numMixture; ++i) {
			//	for (int j = 0; j < numMixture; ++j)
			//		System.err.println(i + "\t" + j + "\t" + t[i][j]);
			//}

			logL = logL + Math.log10(l);
			System.err.println("logL = " + logL);
			if (oldLogL > logL - 0.01)
				converge = true;
					
			oldLogL = logL;
		}

				
/*
				//backward
				double backward[][] = new double[numM + 1][numChromosomes];
				for (int j = 0; j < numChromosomes; ++j)
					backward[numM][j] = 1.0;

				for (int i = numM - 1; i >= 0; --i) {
					long pos = positionsAndChr.get(3 * i);
					int chr = positionsAndChr.get(3 * i + 1).intValue();
					int density = positionsAndChr.get(3 * i + 2).intValue();
					double p1 = 1.0;
					if (i > 0) {
							prevPos = positionsAndChr.get(3 * i - 3);
							p1 = transition(pos - prevPos);
					}
					double p2 = (1.0 - p1) / numChromosomes;

					double sum = 0.0;
					for (int j = 0; j < numChromosomes; ++j)
						sum += backward[i + 1][j];
					
					double sump2 = sum * p2;
					
					for (int j = 0; j < numChromosomes; ++j)
						backward[i][j] = (sump2 + p1 * backward[i + 1][j]) * scale[i + 1] * emission(chr, j, density);
				}
				//for (int j = 0; j < numChromosomes; ++j) {
				//	System.err.println(backward[0][j]);
				//}

				//QE
				for (int i = 0; i < numM; ++i) {
					//int pos = positionsAndChr.get(3 * i);
					int chr = positionsAndChr.get(3 * i + 1).intValue();
					int density = positionsAndChr.get(3 * i + 2).intValue();
					for (int j = 0; j < numChromosomes; ++j) {
						if (chr < 0) {
							//QE[0] += ?
							//QE[1] += ?
						} else {
							if (chr == j)
								QE[0] += backward[i + 1][j] * forward[i + 1][j] / density;
							else
								QE[1] += backward[i + 1][j] * forward[i + 1][j] / density;
						}
					}
				}

				//QT
				for (int i = 0; i < numM; ++i) {
					long pos = positionsAndChr.get(3 * i);
					if (i > 0) {
						int chr = positionsAndChr.get(3 * i + 1).intValue();
						int density = positionsAndChr.get(3 * i + 2).intValue();

						double sum  = 0.0;
						for (int j = 0; j < numChromosomes; ++j)
							sum += backward[i + 1][j];

						double p1 = transition(pos - prevPos);
						double p2 = (1.0 - p1) / numChromosomes;
						//if (ts != 1.0)
						//	System.err.println(ts);
						for (int j = 0; j < numChromosomes; ++j) {
								
							double se = scale[i + 1] * emission(chr, j, density);
							QT[0] += p1 * forward[i][j] * backward[i + 1][j] * se;
							QT[1] += p2 * forward[i][j] * sum * se; 
						}
						
					}
					prevPos = pos;
				}
				
				if (iterations == MAX_ITERATIONS) {
					// Viterbi
					double viterbi[][] = new double[numM + 1][numChromosomes];
					int path[][] = new int[numM + 1][numChromosomes];

					for (int j = 0; j < numChromosomes; ++j)
						viterbi[0][j] = 1.0 / numChromosomes;
					
					prevPos = 0;

					for (int i = 0; i < numM; ++i) {
						long pos = positionsAndChr.get(3 * i);
						int chr = positionsAndChr.get(3 * i + 1).intValue();
						int density = positionsAndChr.get(3 * i + 2).intValue();

						double p1 = 1.0;
						if (i > 0)
							p1 = transition((int)(pos - prevPos));

						double p2 = (1.0 - p1) / numChromosomes;
						
						int maxj = 0;
						for (int j = 1; j < numChromosomes; ++j)
							if (viterbi[i][j] > viterbi[i][maxj])
								maxj = j;

						double sump2 = p2 * viterbi[i][maxj];
						
						for (int j = 0; j < numChromosomes; ++j) {
							double sump1 = p1 * viterbi[i][j]; 
							if (sump2 > sump1) {
								viterbi[i + 1][j] = sump2 * emission(chr, j, density);
								path[i + 1][j] = maxj;
							}
							else {
								viterbi[i + 1][j] = sump1 * emission(chr, j, density);;
								path[i + 1][j] = j;
							}
						}

						double max = 0.0;
						for (int j = 0; j < numChromosomes; ++j)
							max = Math.max(max, viterbi[i + 1][j]);
						
						double imax = 1.0 / max;
						for (int j = 0; j < numChromosomes; ++j)
							viterbi[i + 1][j] *= imax;

						prevPos = pos;
					}
					//backtrack path
					int finalPath[] = new int[numM];

					int maxj = 0;
					for (int j = 1; j < numChromosomes; ++j)
						if (viterbi[numM][j] > viterbi[numM][maxj])
							maxj = j;
					
					for (int i = numM; i >= 1; --i) {
						finalPath[i - 1] = maxj;
						maxj = path[i][maxj];
					}
					
					for (int i = 0; i < numM; ++i) {
						int maxJ = finalPath[i]; // Viterbi path
						double max2 = Double.NEGATIVE_INFINITY;
						for (int j = 0; j < numChromosomes; ++j) {
							if (j != maxJ) {
								max2 = Math.max(max2, backward[i + 1][j] * forward[i + 1][j]);
							}
						}
						max2 = Math.log10(backward[i + 1][maxJ] * forward[i + 1][maxJ]) - Math.log10(max2);
						long pos = positionsAndChr.get(3 * i);
						int chr = positionsAndChr.get(3 * i + 1).intValue();
						if (chr >= 0) // do not print padding markers
							System.out.println(scaffolds.get(scaffold) + "\t" + pos + "\t" + (maxJ + 1) + "\t" + (chr + 1) + "\t" + max2);
					}
				}
	
				//System.err.println(scaffolds.get(scaffold) + "\t" + logScale + "\t" + QE[0] + "\t" + QE[1]);
			}
			logL = logL + Math.log10(l);
			System.err.println("logL = " + logL);
			System.err.println(P_ERROR + "\t" + P_CHIMERA);
			P_ERROR = QE[1] / (QE[0] + QE[1]);
			P_CHIMERA = QT[1] / (QT[0] + QT[1]);*/
	}

	private double emissionTable1[] = new double[1000]; 
	private double emissionTable2[] = new double[1000]; 
	
	private void initParameters()
	{
		//double e = P_ERROR / (double) (numChromosomes - 1); // P_ERROR/0 does not matter!
		for (int i = 0; i < emissionTable1.length; ++i) {
			//emissionTable1[i] = Math.pow(e, 1.0 / i);
			//emissionTable2[i] = Math.pow(1.0 - P_ERROR, 1.0 / i);
		}
		
	}
	
	private double emission(int chr, int state, int density){
		if (chr < 0)
			return 1;

		assert(density >= 1);
		
		if (density < emissionTable1.length) {
			if (chr != state)
				return emissionTable1[density];
			else
				return emissionTable2[density];
			
		}
		return	Math.pow(emission(chr, state, 1), 1.0 / density); 
	}

	private static void usageInfo()
	{
		System.err.println("Usage: java CoverageHMM depth=depth.txt mixture=mixture.txt [options]");
		System.err.println("options:");
		System.err.println("       depth=file            output from samtools -a depth");
        System.err.println("       mixture=file          output from CoverageAnalyser");
        System.err.println("       column=NUM            Take column NUM from depth file [3]");
        System.err.println("       scale=NUM             One position counts this much [0.01]");
        System.err.println("       sample=NUM            Take only NUM fraction of data [0.1]");
        System.err.println("       (default scale*sample=0.001, 1 per 1kb)");
        System.err.println("       minProb=NUM           One position can be this sure [0.001] (1:1000)");
	}
	
	
	public static void main(String args[])
	{
		ParameterParser pp = new ParameterParser();

		String extraParameters = "";
	    for (int i = 0; i < args.length; ++i) {
	            extraParameters +=  " " + args[i];
	    }
	    
		if (args.length == 0 || !pp.init(extraParameters)) {
			usageInfo();
            System.exit(0);
		}
		pp.warning(new String[]{"depth", "mixture", "column", "scale", "minProb", "sample"});
		
		
		System.out.println("#java CoverageHMM" + extraParameters);
			
		int column = Integer.parseInt(pp.getValueAsString("column", "3"));
		CoverageHMM ch = new CoverageHMM();
		ch.setColumn(column - 1);
		String df = pp.getValueAsString("depth", null);
		String mf = pp.getValueAsString("mixture", null);
		if (df == null || mf == null) {
			usageInfo();
            System.exit(0);
		}
		ch.setScale(Double.parseDouble(pp.getValueAsString("scale", "0.001")));
		ch.setMinProb(Double.parseDouble(pp.getValueAsString("minProb", "0.001")));
		ch.setSample(Double.parseDouble(pp.getValueAsString("sample", "0.1")));

		ch.loadMixture(mf);
		ch.loadDepth(df);
		ch.hmm();
	}
}
