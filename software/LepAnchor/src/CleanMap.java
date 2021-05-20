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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class CleanMap {

	int markerDistance = 50;
	int chimericDistance = 1000;
	int numChromosomes;
	
	int numMarkerTypes = 1; // how many different marker types, used for having different error parameters for repeats or different maps... 
	
	double P_CHIMERA = 0.01;
	//double P_ERROR = 0.01;//new double[]{0.01};
	double P_ERROR[] = new double[]{0.01};
	
	int MAX_ITERATIONS = 20;
	
	public void setMarkerDistance(int value)
	{
		markerDistance = value;
	}

	public void setChimericDistance(int value)
	{
		chimericDistance = value;
	}

	public void setInitialValues(double p_c, double p_e)
	{
		for (int type = 0; type < numMarkerTypes; ++type)
			P_ERROR[type] = p_e;
		//P_ERROR = p_e;
		P_CHIMERA = p_c;
	}
	private class CleanMarker {
		public CleanMarker(long position, int chr, int density, int type){
			this.position = position;
			this.chr = chr;
			this.density = density;
			this.type = type;
		}
		
		public long getPosition()
		{
			return position;
		}
		public int getChr()
		{
			return chr;
		}
		public int getDensity()
		{
			return density;
		}
		public int getType()
		{
			return type;
		}
		public void setDensity(int value)
		{
			density = value;
		}
		public void setType(int value)
		{
			type = value;
		}
		
		long position;
		int chr;
		int density;
		int type;
	}
	
	
	public void cleanMap(String mapFile, int annotation)
	{
		
		HashMap<String, ArrayList<CleanMarker>> scaffoldMap = new HashMap<String, ArrayList<CleanMarker>>();
		
		HashMap<Integer,Integer> numAnnotations = new HashMap<Integer,Integer>(); 

		numChromosomes = 0;
		
		int noAnnotation = 0;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(mapFile));
			
			ArrayList<String> line = Input.loadTableRow(br, " \t");
			while (line != null) {
				if (!line.get(0).equals("CHR") && !line.get(0).equals("CHROM")) { 
					String scaffold = line.get(0);
					long position = InputData.myParseLong(line.get(1));
					
					int chr = Integer.parseInt(line.get(2));
					if (chr > 0) {
						if (!scaffoldMap.containsKey(scaffold))
							scaffoldMap.put(scaffold, new ArrayList<CleanMarker>());
						ArrayList<CleanMarker> list = scaffoldMap.get(scaffold);

						//load possible annotation...
						int annot = 0;
						if (annotation > 0 && line.size() >= annotation) {
							annot = Integer.parseInt(line.get(annotation - 1));
							if (annot >= 0) { // no negative annotations... 
								if (!numAnnotations.containsKey(annot))
									numAnnotations.put(annot, 0);
								numAnnotations.put(annot, numAnnotations.get(annot) + 1);
							} else
								++noAnnotation;
						} else
							++noAnnotation;
						
						list.add(new CleanMarker(position, chr - 1, 1, annot));
						// code chr as 0..(c-1), negative missing
						// density = 1, for calculating marker density...
						numChromosomes = Math.max(numChromosomes, chr);
					}
				}
				line = Input.loadTableRow(br, " \t");
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		
		int numScaffolds = scaffoldMap.keySet().size();
//		double markerDensityScales[][] = new double[numScaffolds][]; 
		ArrayList<ArrayList<CleanMarker>> positionsAndChrs = new ArrayList<ArrayList<CleanMarker>>();
		ArrayList<String> scaffolds = new ArrayList<String>();
		for (String scaffold : scaffoldMap.keySet())
			scaffolds.add(scaffold);

		int totalPads = 0;
		for (String scaffold : scaffolds) {
			ArrayList<CleanMarker> positionsAndChr = scaffoldMap.get(scaffold); 
			positionsAndChrs.add(positionsAndChr);
			int numM = positionsAndChr.size();
//			markerDensityScales[si] = new double[numM];
			for (int i = 1; i < numM; ++i) {
				if (positionsAndChr.get(i - 1).getPosition() > positionsAndChr.get(i).getPosition()) {
					System.err.println("physical positions must be increasing...");
					System.exit(0);
				}
			}
			int j = 0;
			int k = 1;
			for (int i = 0; i < numM; ++i) {
				while (positionsAndChr.get(i).getPosition() - positionsAndChr.get(j).getPosition() > markerDistance)
					++j;
				while (k < numM && positionsAndChr.get(k).getPosition() - positionsAndChr.get(i).getPosition() <= markerDistance)
					++k;
				positionsAndChr.get(i).setDensity(k - j);
				//markerDensityScales[si][i] = 1.0 / (k - j);
			}

			ArrayList<CleanMarker> positionsAndChr_padded = new ArrayList<CleanMarker>();
			
			for (int i = 0; i < numM; ++i) {
				if (i > 0) {
					long prevP = positionsAndChr.get(i - 1).getPosition();
					long pos = positionsAndChr.get(i).getPosition();
							
					long numPads = (pos - prevP) / chimericDistance;
					if (numPads >= 2)
						totalPads += numPads - 1;
					for (j = 1; j < numPads; ++j) {
						double p = ((double) j) / numPads;
						// chr missing
						// density is 1
						positionsAndChr_padded.add(new CleanMarker((long)(pos + p * (prevP - pos)), -1, 1, 0));
					}
				}
				positionsAndChr_padded.add(positionsAndChr.get(i));
			}
			positionsAndChr.clear();
			positionsAndChr.addAll(positionsAndChr_padded);
			
		}
		System.err.println("Number of scaffolds = " + numScaffolds);
		System.err.println("Number of chromosomes = " + numChromosomes);
		
		System.err.println("Number of padding markers = " + totalPads);
		
		System.err.println("Number of markers without annotation = " + noAnnotation);
		if (numAnnotations.size() > 0) {
			for (int key : numAnnotations.keySet()) {
				numMarkerTypes = Math.max(numMarkerTypes, key + 1);
				System.err.println("Number of markers with annotation " + key + " = " + numAnnotations.get(key));
			}
		}
		
		DecimalFormat df = new DecimalFormat("#0.000"); 
		
		for (int iterations = 0; iterations <= MAX_ITERATIONS; ++iterations) {
			initParameters();
			
			double QE[][] = new double[numMarkerTypes][2];
//			double QE[] = new double[2];
			double QT[] = new double[2];
			double logL = 0.0;
			double l = 1.0;
			for (int scaffold=0; scaffold < numScaffolds; ++scaffold) {
				ArrayList<CleanMarker> positionsAndChr = positionsAndChrs.get(scaffold);
//				double markerDensityScale[] = markerDensityScales[scaffold];
				int numM = positionsAndChr.size();
				double forward[][] = new double[numM + 1][numChromosomes];

				for (int j = 0; j < numChromosomes; ++j)
					forward[0][j] = 1.0 / numChromosomes;
				
				long prevPos = 0;
				double scale[] = new double[numM + 1];
				scale[0] = 1.0;

				for (int i = 0; i < numM; ++i) {
					CleanMarker m = positionsAndChr.get(i);
					long pos = m.getPosition();
					int chr = m.getChr();
					int density = m.getDensity();
					int type = m.getType();

					double p1 = 1.0;
					if (i > 0)
						p1 = transition((int)(pos - prevPos));

					double p2 = (1.0 - p1) / numChromosomes;
					
					double sum = 0.0;
					for (int j = 0; j < numChromosomes; ++j)
						sum += forward[i][j];
					double sump2 = p2 * sum;
					
					for (int j = 0; j < numChromosomes; ++j)
						forward[i + 1][j] = (sump2 + p1 * forward[i][j]) * emission(chr, j, density, type);

					sum = 0;
					for (int j = 0; j < numChromosomes; ++j)
						sum += forward[i + 1][j];
					double isum = 1.0 / sum;
					
					for (int j = 0; j < numChromosomes; ++j)
						forward[i + 1][j] *= isum;
					scale[i + 1] = isum;
					
					l *= sum;
					if (l < 1e-100) {
						logL += Math.log10(l); 
						l = 1.0;
					}

					prevPos = pos;
				}

				//backward
				double backward[][] = new double[numM + 1][numChromosomes];
				for (int j = 0; j < numChromosomes; ++j)
					backward[numM][j] = 1.0;

				for (int i = numM - 1; i >= 0; --i) {
					CleanMarker m = positionsAndChr.get(i);
					long pos = m.getPosition();
					int chr = m.getChr();
					int density = m.getDensity();
					int type = m.getType();
					double p1 = 1.0;
					if (i > 0) {
							prevPos = positionsAndChr.get(i - 1).getPosition();
							p1 = transition(pos - prevPos);
					}
					double p2 = (1.0 - p1) / numChromosomes;

					double sum = 0.0;
					for (int j = 0; j < numChromosomes; ++j)
						sum += backward[i + 1][j];
					
					double sump2 = sum * p2;
					
					for (int j = 0; j < numChromosomes; ++j)
						backward[i][j] = (sump2 + p1 * backward[i + 1][j]) * scale[i + 1] * emission(chr, j, density, type);
				}
				//for (int j = 0; j < numChromosomes; ++j) {
				//	System.err.println(backward[0][j]);
				//}

				//QE
				for (int i = 0; i < numM; ++i) {
					CleanMarker m = positionsAndChr.get(i);
					int chr = m.getChr();
					int density = m.getDensity();
					int type = m.type;
					for (int j = 0; j < numChromosomes; ++j) {
						if (chr < 0) {
							//QE[0] += ?
							//QE[1] += ?
						} else {
							if (chr == j)
								QE[type][0] += backward[i + 1][j] * forward[i + 1][j] / density;
							else
								QE[type][1] += backward[i + 1][j] * forward[i + 1][j] / density;
						}
					}
				}

				//QT
				for (int i = 0; i < numM; ++i) {
					CleanMarker m = positionsAndChr.get(i);
					long pos = m.getPosition();
					if (i > 0) {
						int chr = m.getChr();
						int density = m.getDensity();
						int type = m.getType();

						double sum  = 0.0;
						for (int j = 0; j < numChromosomes; ++j)
							sum += backward[i + 1][j];

						double p1 = transition(pos - prevPos);
						double p2 = (1.0 - p1) / numChromosomes;
						//if (ts != 1.0)
						//	System.err.println(ts);
						for (int j = 0; j < numChromosomes; ++j) {
								
							double se = scale[i + 1] * emission(chr, j, density, type);
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
						CleanMarker m = positionsAndChr.get(i);
						long pos = m.getPosition();
						int chr = m.getChr();
						int density = m.getDensity();
						int type = m.getType();

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
								viterbi[i + 1][j] = sump2 * emission(chr, j, density, type);
								path[i + 1][j] = maxj;
							}
							else {
								viterbi[i + 1][j] = sump1 * emission(chr, j, density, type);
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
					
					StringBuilder sb = new StringBuilder();
					for (int i = 0; i < numM; ++i) {
						/*int maxJ = 0;
						for (int j = 1; j < numChromosomes; ++j) {
							if (backward[i + 1][j] * forward[i + 1][j] > backward[i + 1][maxJ] * forward[i + 1][maxJ])
								maxJ = j;
						}*/
						int maxJ = finalPath[i]; // Viterbi path
						double max2 = Double.NEGATIVE_INFINITY;
						for (int j = 0; j < numChromosomes; ++j) {
							if (j != maxJ) {
								max2 = Math.max(max2, backward[i + 1][j] * forward[i + 1][j]);
							}
						}
						max2 = Math.log10(backward[i + 1][maxJ] * forward[i + 1][maxJ] / max2);
						CleanMarker m = positionsAndChr.get(i);
						long pos = m.getPosition();
						int chr = m.getChr();
						if (chr >= 0) {// do not print padding markers
							sb.append(scaffolds.get(scaffold));
							sb.append('\t');
							sb.append(pos);
							sb.append('\t');
							sb.append((maxJ + 1));
							sb.append('\t');
							sb.append((chr + 1));
							sb.append('\t');
							sb.append(df.format(max2));
							sb.append('\n');
							if (sb.length() >= 10000) {
								System.out.print(sb.toString());
								sb.setLength(0);
							}
							//System.out.println(scaffolds.get(scaffold) + "\t" + pos + "\t" + (maxJ + 1) + "\t" + (chr + 1) + "\t" + max2);
						}
					}
					if (sb.length() > 0)
						System.out.print(sb.toString());
				}
				//System.err.println(scaffolds.get(scaffold) + "\t" + logScale + "\t" + QE[0] + "\t" + QE[1]);
			}
			logL = logL + Math.log10(l);
			System.err.println("logL = " + logL);
			for (int type = 0; type < numMarkerTypes; ++type)
				System.err.print(P_ERROR[type] + "\t");
			System.err.println(P_CHIMERA);
			for (int type = 0; type < numMarkerTypes; ++type)
				P_ERROR[type] = QE[type][1] / (QE[type][0] + QE[type][1]);
			P_CHIMERA = QT[1] / (QT[0] + QT[1]);
		}
	}	

	private double transition(long distance)
	{
		return (1 - P_CHIMERA);
	}

	private double emissionTable1[][] = null; 
	private double emissionTable2[][] = null; 
	
	private void initParameters()
	{
		if (P_ERROR.length == 1 && numMarkerTypes > 1) {
			double oldValue = P_ERROR[0];
			P_ERROR = new double[numMarkerTypes];
			Arrays.fill(P_ERROR, oldValue);
		}
		if (emissionTable1 == null) {
			emissionTable1 = new double[numMarkerTypes][1000];
			emissionTable2 = new double[numMarkerTypes][1000];
		}
		
		for (int type = 0; type < numMarkerTypes; ++type) {
			double e = P_ERROR[type] / (double) (numChromosomes - 1); // P_ERROR/0 does not matter!
			for (int i = 0; i < emissionTable1[type].length; ++i) {
				emissionTable1[type][i] = Math.pow(e, 1.0 / i);
				emissionTable2[type][i] = Math.pow(1.0 - P_ERROR[type], 1.0 / i);
			}
		}
		
	}
	
	private double emission(int chr, int state, int density, int type){
		if (chr < 0)
			return 1;

		assert(density >= 1);
		
		if (density < emissionTable1[type].length) {
			if (chr != state)
				return emissionTable1[type][density];
			else
				return emissionTable2[type][density];
			
		}
		return	Math.pow(emission(chr, state, 1, type), 1.0 / density); 
	}

	private static void usage()
	{
		System.err.println("Usage: java CleanMap map=contig_pos_map.txt [options]");
        System.err.println("       map=file               a map file containing columns contig, position and chromosome");
        System.err.println("       markerDistance=NUM     Minimum effective distance of markers [50]");
        System.err.println("       chimericDistance=NUM   The default distance of markers for transition (change of chromosome) [1000]");
        System.err.println("       initialValues=NUM NUM  Initial values for theta and epsilon [0.01 0.01]");

        System.err.println("       markerAnnotation=NUM   Load marker annotation from input map column NUM");
	}
	
	
	public static void main(String args[])
	{
		ParameterParser pp = new ParameterParser();

		String extraParameters = "";
	    for (int i = 0; i < args.length; ++i) {
	            extraParameters +=  " " + args[i];
	    }
	    
		if (args.length == 0 || !pp.init(extraParameters)) {
			usage();
            System.exit(0);
		}
		pp.warning(new String[]{"map", "markerDistance", "chimericDistance", "initValues", "markerAnnotation"});
		
		
		System.out.println("#java CleanMap" + extraParameters);
			
		CleanMap cm = new CleanMap();
		cm.setChimericDistance(Integer.parseInt(pp.getValueAsString("chimericDistance", "1000")));
		cm.setMarkerDistance(Integer.parseInt(pp.getValueAsString("markerDistance", "50")));
		cm.setInitialValues(Double.parseDouble(pp.getValueAsString("initialValues",1,"0.01")), Double.parseDouble(pp.getValueAsString("initialValues", 2,"0.01")));
		
		cm.cleanMap(pp.getValueAsString("map", null), Integer.parseInt(pp.getValueAsString("markerAnnotation", "-1")));
	}

}
