import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.locks.ReentrantLock;

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
// User interface for SeparateChromosomes module
//Errorcodes 44xx
public class Separate2 {
	private final ReentrantLock lock = new ReentrantLock();
	
	private double theta1;
	private double theta2;
	
	private double minLod = Double.NEGATIVE_INFINITY;
	
	private Data2 data = new Data2();
	
	private UnionFind chromosomes;
	private UnionFind identicals;

	public void keepFamilies(ArrayList<String> families) {
		data.keepFamilies(families);
	}
	
	public void setPhasedData(boolean value)
	{
		for (Family2 f : data.getFamilies())
			f.setPhasedData(value);
		data.setPhasedData(value);
	}
	
	public int getNumMarkers()
	{
		return data.getNumMarkers();
	}
	
	public void setTheta1(double theta)
	{
		theta1 = theta;
		//logTheta1 = Math.log(theta);		
		//log1mTheta1 = Math.log(1.0 - theta);		
	}
	public void setTheta2(double theta)
	{
		theta2 = theta;
		//logTheta2 = Math.log(theta);		
		//log1mTheta2 = Math.log(1.0 - theta);		
	}
	
	//public void addFamilyFromFile(String filename) {
	//	data.addFamilyFromFile(filename);
	//}
	public void addFamilyFromFile(String filename, String informativeMask) {
		data.addFamilyFromFile(filename, informativeMask, null, false);
	}
	
	public void addFamilyFromFile(String filename, String informativeMask, boolean mask[], boolean gpPhase) {
		data.addFamilyFromFile(filename, informativeMask, mask, gpPhase);
	}
	public void subsample(double rate) {
		data.subSampleData(rate);
	}
	
	public double computeLOD(int m1, int m2) {
		double ret = 0.0;
		for (Family2 f : data.getFamilies()) {
			double l = f.computeLOD(theta1, theta2, m1, m2);
			if (l < minLod)
				l = minLod;
			ret += l;
		}
		return ret;
	}
	
	public void maskInformative(String informativeMask) {
		for (Family2 f : data.getFamilies())
			f.maskInformative(informativeMask);
	}
	public void maskLG(LGMap2 map, int lg) {
		for (Family2 f : data.getFamilies())
			for (int i = 0; i < map.getNumMarkers(); ++i)
				if (map.getLGName(i) != lg)
					f.removeMarker(i);
	}

	public void maskLG(LGMap2 map, boolean lgs[]) {
		for (Family2 f : data.getFamilies())
			for (int i = 0; i < map.getNumMarkers(); ++i)
				if (!lgs[map.getLGName(i)])
					f.removeMarker(i);
	}
	
	public void setLod3Mode(int mode)
	{
		for (Family2 f : data.getFamilies())
			f.setLod3Mode(mode);
	}	
	
	public void setDistortionLodMode()
	{
		for (Family2 f : data.getFamilies())
			f.setDistortionLodMode();
	}	
	
	
	
/*	private void joinIdenticals(LGMap2 map)
	{
		joinIdenticals(map, true, theta1, theta2);
	}
	private void joinIdenticals(LGMap2 map, boolean printPatterns)
	{
		joinIdenticals(map, printPatterns, theta1, theta2);
	}
	private void joinIdenticals(LGMap2 map, double theta1, double theta2)
	{
		joinIdenticals(map, true, theta1, theta2);
	}
	*/
	
	private void joinIdenticals(ArrayList<Integer> markers, double theta1, double theta2, double lodLimit, UnionFind uf)
	{
		ArrayList<Integer> infoMarkers = new ArrayList<Integer>();
		
		ArrayList<Double> maxScores = new ArrayList<Double>(); 

		for (int m1 : markers) {
			boolean informative = false;
			for (Family2 f : data.getFamilies())
				if (f.isFatherInformative(m1) || f.isMotherInformative(m1)) {
					informative = true;
					break; 
				}
			if (informative) {
				double maxLODScore = 0.0; 
				for (Family2 f : data.getFamilies())
					maxLODScore += f.maxLodScore(theta1, theta2, m1);
				
				if (maxLODScore >= lodLimit) {				
					infoMarkers.add(m1);
					maxScores.add(maxLODScore);
				}
			}
		}
		if (infoMarkers.size() == 0) {
			markers.clear();
			return;
		}

		Misc.ArrayIndexComparator<Double> comparator = new Misc.ArrayIndexComparator<Double>(maxScores);
		Integer[] indexes = comparator.createIndexArray();
		Arrays.sort(indexes, comparator);
		
		ArrayList<Integer> infoMarkers2 = new ArrayList<Integer>();
		for (int i : indexes) {
			infoMarkers2.add(infoMarkers.get(i));
			//System.err.println(maxScores.get(i));
		}
		//System.err.println();
		infoMarkers = infoMarkers2;
		Collections.reverse(infoMarkers);
		
		while (true) {
			ArrayList<ArrayList<Integer>> cliques = new ArrayList<ArrayList<Integer>>();
			cliques.add(new ArrayList<Integer>());
			for (int m : infoMarkers) {
				boolean canAdd = true;
				//int perm[] = Misc.randomPermutation(cliques.size());
				//for (int cliquei : perm) {
				//	ArrayList<Integer> clique = cliques.get(cliquei); //Randomize the clique order
				for (ArrayList<Integer> clique : cliques) {
					canAdd = true;
					for (int cm : clique)
						if (computeLOD(m, cm) < lodLimit) {
							canAdd = false;
							break;
						}
					if (canAdd) {
						clique.add(m);
						break;
					}
				}
				if (!canAdd) {
					ArrayList<Integer> al = new ArrayList<Integer>();
					al.add(m);
					cliques.add(al);
				}
			}
			//System.err.println(cliques.size() +"\t" + infoMarkers.size());
			if (infoMarkers.size() == cliques.size())
				break;
			infoMarkers.clear();
			for (ArrayList<Integer> clique : cliques) { 
				int c0 = clique.get(0);
				infoMarkers.add(c0);
				for (int cm : clique)
					if (cm != c0) {
						for (Family2 f : data.getFamilies()) {
							f.combineAndNormalize(c0, cm, theta1, theta2);
						}
						uf.union(cm, c0);
					}
				for (int cm : clique)
					if (cm != c0)
						for (Family2 f : data.getFamilies())
							f.setProb(cm, c0);
			}
		}
		markers.clear();
		markers.addAll(infoMarkers);
	}
	
	private void printPatterns(int name, int size, int marker)
	{
		for (Family2 f : data.getFamilies()) {
			int type = 0;
			if (f.isFatherInformative(marker))
				type += 1;
			if (f.isMotherInformative(marker))
				type += 2;
			
			System.err.print(name + "\t" + type + "\t" + size + "\t" + f.pattern(marker) + "\t");
		}
		System.err.println();
	}
	private void printPatterns(int name, ArrayList<Integer> markers)
	{
		for (Family2 f : data.getFamilies()) {
			int numP = 1;
			for (int m : markers) {
				System.err.print(name + "." + numP + "\t" + f.pattern(m) + "\t");				
				++numP;
			}
		}
		System.err.println();
	}

	
	public void joinSingles2Identicals(double lodLimit_[], double lodDifference, LGMap2 map, boolean betweenSameTypeOnly, int numThreads)
	{
		//joinIdenticals(map, true, theta1, theta2, lodLimit[0]);
		if (betweenSameTypeOnly && (data.getNumFamilies() > 1 || lodLimit_.length != 3))
			Error.error(701);

		if (!betweenSameTypeOnly && lodLimit_.length > 1)
			Error.error(702);
		
		@SuppressWarnings("unchecked")
		ArrayList<Double> smap[] = new ArrayList[getNumMarkers()];
		for (int i = 0; i < getNumMarkers(); ++i) {
			smap[i] = new ArrayList<Double>();
			int chr = map.getLGName(i); 
			if (chr != 0)
				smap[i].add((double) chr); 
		}
		//int joined = 0;
		//int joined2 = 0;
		
		String masks[] = (betweenSameTypeOnly) ? new String[]{"1", "2", "3"}:new String[]{"123"};
		
		joinIdenticalsFast(map);
		
		for (String mask : masks) {
			double lodLimit = (betweenSameTypeOnly) ? lodLimit_[Integer.parseInt(mask) - 1] : lodLimit_[0]; 
		
			ArrayList<Integer> infoMarkers = new ArrayList<Integer>();
			for (int m1 = 0; m1 < getNumMarkers(); ++m1)
				if (map.getLGName(m1) == 0) {
					boolean informative = false;
					for (Family2 f : data.getFamilies()) {
						int imode = 0;
						if (f.isFatherInformative(m1))
							++imode;
						if (f.isMotherInformative(m1))
							imode +=2;
						if (mask.contains("" + imode)) {
							informative = true;
							break; 
						}
					}
					if (informative) {
						double maxLODScore = 0.0; 
						for (Family2 f : data.getFamilies())
							maxLODScore += f.maxLodScore(theta1, theta2, m1);
						if (maxLODScore >= lodLimit) 		
							infoMarkers.add(m1);
					}
				}
			ArrayList<ArrayList<Integer>> markers = map.getMarkersOfLGs(getNumMarkers());
	
			ArrayList<Integer> groupsWithMarkers = new ArrayList<Integer>(); 
			for (int c = 1; c < markers.size(); ++c)
				if (markers.get(c).size() > 0) {
					if (betweenSameTypeOnly) {
						int m1 = markers.get(c).get(0);
						int imode = 0;
						for (Family2 f : data.getFamilies()) {
							if (f.isFatherInformative(m1))
								++imode;
							if (f.isMotherInformative(m1))
								imode +=2;
						}
						if (mask.contains("" + imode))
							groupsWithMarkers.add(c);
					} else
						groupsWithMarkers.add(c);
				}
			
			ArrayList<ArrayList<Integer>> markersMasked = new ArrayList<ArrayList<Integer>>();

			for (ArrayList<Integer> mc : markers) {
				ArrayList<Integer> mcm = new ArrayList<Integer>();
				markersMasked.add(mcm);
				for (int m2 : mc) {
					int imode = 0;
					for (Family2 f : data.getFamilies()) {
						if (f.isFatherInformative(m2))
							++imode;
						if (f.isMotherInformative(m2))
							imode +=2;
					}
					if (mask.contains("" + imode)) {
						mcm.add(m2);
						break;
					}
				}
			}
				
			Thread threads[] = new Thread[numThreads];
			
			for (int t = 0; t < numThreads; ++t) {
				threads[t] = new Thread(new JoinSinglesThread(infoMarkers, markersMasked, lodLimit, lodDifference, t, numThreads, smap, t == 0));
				threads[t].start();
			}
			try {
				for (int t = 0; t < numThreads; ++t)
					if (threads[t].isAlive())
						threads[t].join();
			} catch (Exception e) {
				e.printStackTrace();
				Error.error(-99999999);
			}
		}
		for (int i = 0; i < getNumMarkers(); ++i) {
			StringBuilder str = new StringBuilder();
			if (smap[i].size() == 0)
				str.append('0');
			if (smap[i].size() == 1)
				str.append(smap[i].get(0).intValue());
			else
				for (int j = 0; j < smap[i].size(); j+=2) {
					str.append(smap[i].get(j).intValue());
					str.append('\t');
					str.append(smap[i].get(j+1));
				}
			System.out.println(str);
		}
		//System.err.println("#joined = " + joined + " #joined2 = " + joined2);
	}
	
	private class JoinSinglesThread implements Runnable
	{
		ArrayList<Integer> infoMarkers;
		ArrayList<ArrayList<Integer>> markers;
		double lodLimit, lodDifference;
		int numThreads, start;
		ArrayList<Double> smap[];
		boolean print;
		
		public JoinSinglesThread(ArrayList<Integer> infoMarkers, ArrayList<ArrayList<Integer>> markers, double lodLimit, double lodDifference, int start, int numThreads,ArrayList<Double> smap[], boolean print)
		{
			this.infoMarkers = infoMarkers;
			this.markers = markers;
			this.lodLimit = lodLimit;
			this.lodDifference = lodDifference;
			this.numThreads = numThreads;
			this.start = start;
			this.smap = smap;
			this.print = print;
		}

		@Override
		public void run() {
			int progress = 0;
			int prevProgress = 0;
			
			for (int m2i = start; m2i < infoMarkers.size(); m2i+=numThreads) {
				int m2 = infoMarkers.get(m2i);
				if (print) {
					progress = (10 * m2i / infoMarkers.size());
					if (progress != prevProgress) {
						System.err.print(progress + ".");
						prevProgress = progress;
					}
				}
				
				
				double maxScore = Double.NEGATIVE_INFINITY;
				double maxScore2 = Double.NEGATIVE_INFINITY;
				
				
				for (int c = 1; c < markers.size(); ++c) {
					double score = Double.NEGATIVE_INFINITY;					
					for (int m1 : markers.get(c)) { 
						double l = computeLOD(m1, m2);
						if (l > score)
							score = l;
					}
					if (score > maxScore) {
						maxScore2 = maxScore;
						maxScore = score;
					} else
						maxScore2 = Math.max(maxScore2, score);
					

					if (score >= lodLimit) {

						if (smap[m2].size() == 0) {
							smap[m2].add((double) c);
							smap[m2].add(score);
						}
						else  {
							double bestScore = smap[m2].get(1); 
							if (score > bestScore) {
								smap[m2].add(smap[m2].get(0));
								smap[m2].add(bestScore);
								smap[m2].set(0, (double) c);
								smap[m2].set(1, score);
							}
							else {
								smap[m2].add((double)c);
								smap[m2].add(score);
							}
						}
							
					}
				}
				if (maxScore - maxScore2 < lodDifference) {
					smap[m2].clear();
				}
			}
			if (print)
				System.err.println(" done!");

		}
	}
	
	private class JoinSinglesThreadMaxDistance implements Runnable
	{
		ArrayList<Integer> infoMarkers;
		double lodLimit, lodDifference;
		int numThreads, start;
		ArrayList<Double> smap[];
		boolean print;
		int maxDistance;
		
		public JoinSinglesThreadMaxDistance(ArrayList<Integer> infoMarkers, double lodLimit, double lodDifference, int start, int numThreads,ArrayList<Double> smap[], boolean print, int maxDistance)
		{
			this.infoMarkers = infoMarkers;
			this.lodLimit = lodLimit;
			this.lodDifference = lodDifference;
			this.numThreads = numThreads;
			this.start = start;
			this.smap = smap;
			this.print = print;
			this.maxDistance = maxDistance;
		}

		@Override
		public void run() {
			int progress = 0;
			int prevProgress = 0;
			
			for (int m2i = start; m2i < infoMarkers.size(); m2i+=numThreads) {
				int m2 = infoMarkers.get(m2i);
				if (print) {
					progress = (10 * m2i / infoMarkers.size());
					if (progress != prevProgress) {
						System.err.print(progress + ".");
						prevProgress = progress;
					}
				}
				if (smap[m2].size() == 1)
					continue;
				
				double maxScore = Double.NEGATIVE_INFINITY;
				int maxChr = -1;
				
				for (int m1i = Math.max(0, m2i - maxDistance); m1i < m2i + maxDistance && m1i < infoMarkers.size(); ++m1i) {
					int m1 = infoMarkers.get(m1i);
					if (smap[m1].size() != 1)
						continue;
					double l = computeLOD(m1, m2);
					if (l > maxScore) {
						int chr = smap[m1].get(0).intValue();
						maxScore = l;
						maxChr = chr;
					}
				}
				double maxScore2 = Double.NEGATIVE_INFINITY;
				for (int m1i = Math.max(0, m2i - maxDistance); m1i < m2i + maxDistance && m1i < infoMarkers.size(); ++m1i) {
					int m1 = infoMarkers.get(m1i);
					if (smap[m2].size() != 1)
						continue;
					int chr = smap[m1].get(0).intValue();
					if (chr != maxChr)
						maxScore2 = Math.max(maxScore2, computeLOD(m1, m2));
				}
				if (maxScore >= lodLimit) {
					if (smap[m2].size() == 0) {
						smap[m2].add((double) maxChr);
						smap[m2].add(maxScore);
					}
					else  {
						double bestScore = smap[m2].get(1); 
						if (maxScore > bestScore) {
							smap[m2].add(smap[m2].get(0));
							smap[m2].add(bestScore);
							smap[m2].set(0, (double) maxChr);
							smap[m2].set(1, maxScore);
						}
						else {
							smap[m2].add((double) maxChr);
							smap[m2].add(maxScore);
						}
					}
							
				}
				if (maxScore - maxScore2 < lodDifference) {
					smap[m2].clear();
				}
			}
			if (print)
				System.err.println(" done!");
		}
	}
	
	public void joinSingles(double lodLimit_[], double lodDifference, LGMap2 map, boolean betweenSameTypeOnly, boolean iterate, final int numThreads, int maxDistance)
	{
		//joinIdenticals(map, true, theta1, theta2, lodLimit[0]);
		if (betweenSameTypeOnly && (data.getNumFamilies() > 1 || lodLimit_.length != 3))
			Error.error(701);

		if (!betweenSameTypeOnly && lodLimit_.length > 1)
			Error.error(702);
		
		@SuppressWarnings("unchecked")
		ArrayList<Double> smap[] = new ArrayList[getNumMarkers()];
		for (int i = 0; i < getNumMarkers(); ++i) {
			smap[i] = new ArrayList<Double>();
			int chr = map.getLGName(i); 
			if (chr != 0)
				smap[i].add((double) chr); 
		}
		
		String masks[] = (betweenSameTypeOnly) ? new String[]{"1", "2", "3"}:new String[]{"123"};

		ArrayList<ArrayList<Integer>> markers = map.getMarkersOfLGs(getNumMarkers());
		
		for (String mask : masks) {
			double lodLimit = (betweenSameTypeOnly) ? lodLimit_[Integer.parseInt(mask) - 1] : lodLimit_[0];
		
			final ArrayList<Integer> infoMarkers = new ArrayList<Integer>();
			for (int m1 = 0; m1 < getNumMarkers(); ++m1)
				if (map.getLGName(m1) == 0 || maxDistance < Integer.MAX_VALUE) {
					boolean informative = false;
					for (Family2 f : data.getFamilies()) {
						int imode = 0;
						if (f.isFatherInformative(m1))
							++imode;
						if (f.isMotherInformative(m1))
							imode +=2;
						if (mask.contains("" + imode)) {
							informative = true;
							break; 
						}
					}
					if (informative) {
						double maxLODScore = 0.0; 
						for (Family2 f : data.getFamilies())
							maxLODScore += f.maxLodScore(theta1, theta2, m1);
						if (maxLODScore >= lodLimit) 		
							infoMarkers.add(m1);
					}
				}

			ArrayList<ArrayList<Integer>> markersMasked = null;
			if (betweenSameTypeOnly) {
				markersMasked = new ArrayList<ArrayList<Integer>>();
				for (ArrayList<Integer> mc : markers) {
					ArrayList<Integer> mcm = new ArrayList<Integer>();
					markersMasked.add(mcm);
					for (int m2 : mc) {
						int imode = 0;
						for (Family2 f : data.getFamilies()) {
							if (f.isFatherInformative(m2))
								++imode;
							if (f.isMotherInformative(m2))
								imode +=2;
						}
						if (mask.contains("" + imode))
							mcm.add(m2);
					}
					
				}
			}
			else 
				markersMasked = markers;

			boolean changed = true;
			int iteratation = 1;
			do { 
				Thread threads[] = new Thread[numThreads];
				
				for (int t = 0; t < numThreads; ++t) {
					if (maxDistance < Integer.MAX_VALUE)
						threads[t] = new Thread(new JoinSinglesThreadMaxDistance(infoMarkers, lodLimit, lodDifference, t, numThreads, smap, t == 0, maxDistance));
					else
						threads[t] = new Thread(new JoinSinglesThread(infoMarkers, markersMasked, lodLimit, lodDifference, t, numThreads, smap, t == 0));
					threads[t].start();
				}
				try {
					for (int t = 0; t < numThreads; ++t)
						if (threads[t].isAlive())
							threads[t].join();
				} catch (Exception e) {
					e.printStackTrace();
					Error.error(-99999999);
				}
				if (iterate) {
					changed = false;
					for (ArrayList<Integer> mc : markersMasked)
						mc.clear();
					ArrayList<Integer> infoMarkers2 = new ArrayList<Integer>();
					for (int m : infoMarkers)
						if (smap[m].size() == 0)
							infoMarkers2.add(m);
						else {
							int chr = smap[m].get(0).intValue(); 
							markersMasked.get(chr).add(m);
							changed = true;
						}
					if (changed)
						System.err.println("Starting iteration " + (++iteratation));

					infoMarkers.clear();
					infoMarkers.addAll(infoMarkers2);  
				}
				
			} while (changed && iterate);
		}
				
		for (int i = 0; i < getNumMarkers(); ++i) {
			StringBuilder str = new StringBuilder();
			if (smap[i].size() == 0)
				str.append('0');
			if (smap[i].size() == 1)
				str.append(smap[i].get(0).intValue());
			else
				for (int j = 0; j < smap[i].size(); j+=2) {
					if (j > 0)
						str.append('\t');
					str.append(smap[i].get(j).intValue());
					str.append('\t');
					str.append(smap[i].get(j+1));
				}
			System.out.println(str);
		}
	}
	
	
	public void joinIdenticals(LGMap2 map)
	{
		UnionFind uf = new UnionFind(getNumMarkers());
		//joinIdenticals(....,uf)
		ArrayList<ArrayList<Integer>> markers = map.getMarkersOfLGs(getNumMarkers());
		for (int c1 = 1; c1 < markers.size(); ++c1) {
			ArrayList<Integer> mc = markers.get(c1);
			joinIdenticals(mc, 0.0, 0.0, 0.0, uf); // quite slow with many markers...
			if (mc.size() > 0)
				printPatterns(c1, mc.size(), mc.get(0));
		}
		
	}
	public void joinIdenticalsFast(LGMap2 map) //avoid quadratic time consumption...
	{
		if (theta1 != 0.0 || theta2 != 0.0)
			return;
		
		//UnionFind uf = new UnionFind(getNumMarkers());
		//joinIdenticals(....,uf)
		ArrayList<ArrayList<Integer>> markers = map.getMarkersOfLGs(getNumMarkers());
		for (int c1 = 1; c1 < markers.size(); ++c1) {
			ArrayList<Integer> mc = markers.get(c1);
			if (mc.size() > 0) {
				int mc0 = mc.get(0);
				for (int m : mc)
					if (m != mc0)
						join(mc0, m);
				for (int m : mc)
					if (m != mc0)
						for (Family2 f : data.getFamilies())
							f.setProb(m, mc0);
				
			}
			//joinIdenticals(mc, 0.0, 0.0, 0.0, uf); // quite slow with many markers...
			if (mc.size() > 0)
				printPatterns(c1, mc.size(), mc.get(0));
		}
		
	}
	
	
	
	public void joinLGs(double lodLimit[], LGMap2 map, int sizeLimit)
	{
		//joinIdenticals(map, true, 0.0, 0.0, lodLimit[0]);
	
/*		UnionFind uf = new UnionFind(getNumMarkers());
		//joinIdenticals(....,uf)
		ArrayList<ArrayList<Integer>> markers = map.getMarkersOfLGs(getNumMarkers());
		for (int c1 = 1; c1 < markers.size(); ++c1) {
			ArrayList<Integer> mc = markers.get(c1);
			ArrayList<Integer> tmp = new ArrayList<Integer>();
			tmp.addAll(mc);
			//joinIdenticals(mc, 0.0, 0.0, 0.0, uf);
			if (mc.size() > 0)
				printPatterns(c1, mc.get(0));
			markers.set(c1, tmp);
			//if (mc.size() > 0) {
			//	int m0 = mc.get(0);
			//	for (int m: mc)
			//		if (m != m0)
			//			join(m0, m);
			//	printPatterns(c1, m0);
			//}
		}*/

		joinIdenticalsFast(map);
		
		ArrayList<ArrayList<Integer>> markers = map.getMarkersOfLGs(getNumMarkers());
		
		UnionFind chromosomes = new UnionFind(markers.size());
		
		double maxScore = Double.NEGATIVE_INFINITY;

		for (int c1 = 1; c1 < markers.size(); ++c1) {
			if (markers.get(c1).size() > sizeLimit) {
				int m2 = markers.get(c1).get(0); 
				int f1 = chromosomes.find(c1);
				for (int c2 = c1 + 1; c2 < markers.size(); ++c2)
					if (markers.get(c2).size() > sizeLimit && chromosomes.find(c2) != f1) {
						int m1 = markers.get(c2).get(0);
						double l = computeLOD(m1, m2);
						double maxLODScore1 = 0.0;
						double maxLODScore2 = 0.0;
						for (Family2 f : data.getFamilies()) {
							maxLODScore1 += f.maxLodScore(theta1, theta2, m1);
							maxLODScore2 += f.maxLodScore(theta1, theta2, m1);
						}
						int j = 0;
						while (j < lodLimit.length - 1 && (maxLODScore1 < lodLimit[j] || maxLODScore2 < lodLimit[j]))
							++j;
						
						if (l > lodLimit[j]) {
							f1 = chromosomes.union(c1, c2);
						}
						else 
							maxScore = Math.max(maxScore, l);
					}
				}
		}


		int names2[] = new int[getNumMarkers()];
		
		int names[] = new int[markers.size()];
		int numLGs = 0;

		for (int c1 = 1; c1 < markers.size(); ++c1) {
			if (markers.get(c1).size() > 0) {
				int f1 = chromosomes.find(c1);
				if (names[f1] == 0)
					names[f1] = ++numLGs;
				for (int m: markers.get(c1))
					names2[m] = names[f1];
					
			}
		}
		for (int i = 0; i < getNumMarkers(); ++i) {
			System.out.println(names2[i] + "\t" + map.getLGName(i));
		}
		
		System.err.println((markers.size() - 1) + "->" + numLGs);
		System.err.println("Next join at LOD=" + maxScore);
	}
	
	
	private void join(int m1, int m2)
	{
		for (Family2 f : data.getFamilies()){
			f.combineAndNormalize(m1, m2, theta1, theta2);
			//f.normalize(m1);
		}
	}

	private class SeparateIdenticalThread implements Runnable{
		private ArrayList<Integer> infoMarkers_;
		private double lodLimit;
		private boolean print;
		private double keepRate = 0.5;
		
		public SeparateIdenticalThread(ArrayList<Integer> markers, double lodLimit, boolean print, double keepRate)
		{
			infoMarkers_ = markers;
			this.lodLimit = lodLimit;
			this.print = print;
			this.keepRate = keepRate;
		}
		
		public void run()
		{
			ArrayList<Integer> infoMarkers = new ArrayList<Integer>();
			infoMarkers.addAll(infoMarkers_);
			
			int oldProgress = 0;
			if (print)
				System.err.println("computing pairwise LOD scores");
			for (int m1i = 0; m1i < infoMarkers.size() - 1; ++m1i) {
				int m1 =  infoMarkers.get(m1i); 
				int c1 = chromosomes.find(m1);

				if (print) {
					int newProgress = 10 * m1i / infoMarkers.size();
					for (int i = oldProgress; i < newProgress; ++i)
						System.err.print(i + 1);
					oldProgress = newProgress;
				}
				
				int infoIndex = m1i + 1;
				
				for (int m2i = m1i + 1; m2i < infoMarkers.size(); ++m2i) 
					{
						int m2 = infoMarkers.get(m2i);
						int c2 = chromosomes.find(m2);
						boolean didJoin = false;
						if (c1 != c2) {
							double lod = computeLOD(m1, m2);
							if (lod >= lodLimit) {
								c1 = chromosomes.union(c1, c2);
								didJoin = true;
							}
						}
						if (!didJoin || Math.random() < keepRate) {
							infoMarkers.set(infoIndex, m2);
							++infoIndex;
						}
					}
				for (int i = infoMarkers.size() - 1; i >= infoIndex; --i)
					infoMarkers.remove(i);
					
			}
			
			infoMarkers.clear();
			
			ArrayList<ArrayList<Integer>> markers = new ArrayList<ArrayList<Integer>>();

			int name[] = new int[getNumMarkers()];
			int numC = 0;
			
			for (int m : infoMarkers_) {
				int c = chromosomes.find(m);
				if (name[c] == 0)
					name[c] = ++numC;
				while (markers.size() < name[c]) {
					markers.add(new ArrayList<Integer>());
				}
				markers.get(name[c] - 1).add(m);
			}
			
			for (ArrayList<Integer> al : markers) {
				joinIdenticals(al, theta1, theta2, lodLimit, identicals);
				for (int m : al)
					infoMarkers.add(m);
			}
			
			
			
			infoMarkers_.clear();
			infoMarkers_.addAll(infoMarkers);
			
			if (print)
				System.err.println(" done!");			
		}
	}

	
	public LGMap2 separateIdenticals(double lodLimit_[], int numParts, int numThreads, boolean removeSingles, double keepRate, boolean betweenSameTypeOnly)
	{
		
		if (betweenSameTypeOnly && (data.getNumFamilies() > 1 || lodLimit_.length != 3))
			Error.error(701);

		if (!betweenSameTypeOnly && lodLimit_.length > 1)
			Error.error(702);
		
		String masks[] = (betweenSameTypeOnly) ? new String[]{"1", "2", "3"}:new String[]{"123"}; 
		
		numThreads = Math.min(numThreads, numParts);
		
		chromosomes = new UnionFind(getNumMarkers());
		identicals = new UnionFind(getNumMarkers());
		
		for (String mask : masks) {
			double lodLimit = (betweenSameTypeOnly) ? lodLimit_[Integer.parseInt(mask) - 1] : lodLimit_[0];
		
			ArrayList<Integer> infoMarkers = new ArrayList<Integer>();
			ArrayList<Double> maxScores = new ArrayList<Double>(); 
	
			for (int m1 = 0; m1 < getNumMarkers(); ++m1) {
				boolean informative = false;
				for (Family2 f : data.getFamilies()) {
					int imode = 0;
					if (f.isFatherInformative(m1))
						++imode;
					if (f.isMotherInformative(m1))
						imode +=2;
					if (mask.contains("" + imode)) {
						informative = true;
						break; 
					}
				}
				if (informative) {
					double maxLODScore = 0.0; 
					for (Family2 f : data.getFamilies())
						maxLODScore += f.maxLodScore(theta1, theta2, m1);
					
					if (maxLODScore >= lodLimit) {				
						infoMarkers.add(m1);
						maxScores.add(maxLODScore);
					}
				}
			}
	
			
			int markerIndex[] = new int[getNumMarkers()];
			int index = 0;
			for (int m : infoMarkers) {
				markerIndex[m] = index >> 5; // keep 32 adjacent markers together (physical order)
				if ((++index >> 5) >= numParts)
					index = 0;
			}
			
			Misc.ArrayIndexComparator<Double> comparator = new Misc.ArrayIndexComparator<Double>(maxScores);
			Integer[] indexes = comparator.createIndexArray();
			Arrays.sort(indexes, comparator);
			
			ArrayList<Integer> infoMarkers2 = new ArrayList<Integer>();
			
			
			for (int i : indexes) {
				infoMarkers2.add(infoMarkers.get(i));
				//System.err.println(maxScores.get(i));
			}
			infoMarkers = infoMarkers2;
			Collections.reverse(infoMarkers);
			
			ArrayList<ArrayList<Integer>> infoMarkersForEachThread = new ArrayList<ArrayList<Integer>>();
			
			
	
			for (int i = 0; i < numParts; ++i)
				infoMarkersForEachThread.add(new ArrayList<Integer>());
			
			for (int m : infoMarkers) {
				infoMarkersForEachThread.get(markerIndex[m]).add(m); // keep 32 markers in the physical order together...
			}
			
			Thread threads[] = new Thread[numThreads];  
	
			for (int j = 0; j < numParts; j+= numThreads) {  
				for (int i = 1; i < numThreads; ++i) {
					if (i + j < numParts) { 
						threads[i] = new Thread(new SeparateIdenticalThread(infoMarkersForEachThread.get(i + j), lodLimit, false, keepRate)); 
						threads[i].start();
					}
				}
				new SeparateIdenticalThread(infoMarkersForEachThread.get(j), lodLimit, true, keepRate).run();
				int t = 1;
				while (t < numThreads) {
					while (threads[t].isAlive())
						try {
							threads[t].join();
						} catch (Exception e) {
							e.printStackTrace();
							Error.error(-99999999);
						}
					++t;
				}
			}
			
			if (numParts >= 1) {
				System.err.println("combining all");
	
				System.err.print(infoMarkers.size() + "->");
				
				if (removeSingles) {
					int sizes[] = new int[getNumMarkers()];
					int sizes2[] = new int[getNumMarkers()];
					for (int m = 0; m < getNumMarkers(); ++m) {
						++sizes[identicals.find(m)];
						++sizes2[chromosomes.find(m)];
					}
		
					for (int i = 0; i < numParts; ++i) {
						ArrayList<Integer> ml = infoMarkersForEachThread.get(i);
						ArrayList<Integer> ml2 = new ArrayList<Integer>();
						for (int m : ml)
							if (sizes2[chromosomes.find(m)] > 1 && sizes[identicals.find(m)] >= 1) {
								ml2.add(m);
								sizes[identicals.find(m)] = 0;
							}
						infoMarkersForEachThread.set(i, ml2);					
					}
				}
				
				infoMarkers.clear();
				int maxM = 0;
				for (int i = 0; i < numParts; ++i)
					maxM = Math.max(maxM, infoMarkersForEachThread.get(i).size());
		
				for (int m = 0; m < maxM; ++m)
					for (int i = 0; i < numParts; ++i) {
						ArrayList<Integer> ml = infoMarkersForEachThread.get(i);
						if (ml.size() > m)
							infoMarkers.add(ml.get(m));
					}
				System.err.println(infoMarkers.size());
				
				SeparateIdenticalThread st = new SeparateIdenticalThread(infoMarkers, lodLimit, true, 1.0);
				int oldSize = infoMarkers.size();
				st.run();
				System.err.println(oldSize + "->" + infoMarkers.size());
				st.run();
			}
		}
			
		int sizes[] = new int[getNumMarkers()];
		int names[] = new int[getNumMarkers()];
		int numChromosomes = 0;
		
		for (int m = 0; m < getNumMarkers(); ++m) {
			int c = identicals.find(m);
			++sizes[c];
			if (names[c] == 0 && sizes[c] > 1)
				names[c] = ++numChromosomes;
		}
		//System.err.println("numChromosomes = " + numChromosomes);

		LGMap2 lGMap = new LGMap2(getNumMarkers());
		for (int m = 0; m < getNumMarkers(); ++m)
			lGMap.setLGName(m, names[identicals.find(m)]);

		lGMap.renameLGs();
		for (int m = 0; m < getNumMarkers(); ++m) {
			int c = identicals.find(m);
			if (names[c] > 0) {
				printPatterns(lGMap.getLGName(m), sizes[c], m);				
				names[c] = 0;
			}
		}				
		
		//lGMap.printLGAssignment(numMarkers);
		return lGMap;
	}
	
	public void outputData()
	{
		outputData(0);
	}
	public void outputData(int sizeLimit)
	{

		System.out.print("CHR\tPOS");
		for (Family2 f : data.getFamilies()) {
			for (int i = 0; i < f.getNumChildren() + 2; ++i)
				System.out.print("\t" + f.getName());
		}
		System.out.println();

		System.out.print("CHR\tPOS");
		int fam = 1;
		for (Family2 f : data.getFamilies()) {
			for (int i = 0; i < f.getNumChildren() + 2; ++i)
				System.out.print("\t" + (i + fam));
			fam += f.getNumChildren() + 3;
		}			
		System.out.println();
		
		System.out.print("CHR\tPOS");
		fam = 1;
		for (Family2 f : data.getFamilies()) {
			System.out.print("\t0\t0");
			for (int i = 0; i < f.getNumChildren(); ++i)
				System.out.print("\t" + fam);
			fam += f.getNumChildren() + 3;

		}			
		System.out.println();

		System.out.print("CHR\tPOS");
		fam = 1;
		for (Family2 f : data.getFamilies()) {
			System.out.print("\t0\t0");
			for (int i = 0; i < f.getNumChildren(); ++i)
				System.out.print("\t" + (fam + 1));
			fam += f.getNumChildren() + 3;
		}			
		System.out.println();

		System.out.print("CHR\tPOS");
		for (Family2 f : data.getFamilies()) {
			System.out.print("\t1\t2");
			for (int i = 0; i < f.getNumChildren(); ++i)
				System.out.print("\t0");
		}
		System.out.println();
		
		System.out.print("CHR\tPOS");
		for (Family2 f : data.getFamilies()) {
			for (int i = 0; i < f.getNumChildren() + 2; ++i)
				System.out.print("\t0");
		}
		System.out.println();
			
		HashMap<Object, Integer> hm = new HashMap<Object, Integer>();
		
		for (Family2 f : data.getFamilies())
			for (int m1 = 0; m1 < getNumMarkers(); ++m1) {
				float p[] = f.getProb(m1);
				if (!hm.containsKey(p))
					hm.put(p, 0);
				hm.put(p, hm.get(p) + 1);
			}
		
		for (int m1 = 0; m1 < getNumMarkers(); ++m1) {
			//System.out.print("CHR\tPOS");		
			//System.out.print(data.getMarkerName(m1));
			StringBuilder sb = new StringBuilder(data.getMarkerName(m1));
			for (Family2 f : data.getFamilies()) {
				float p[] = f.getProb(m1);
				if (hm.get(p) < sizeLimit)
					p = null;
				if (p == null || (!f.isFatherInformative(m1) && !f.isMotherInformative(m1))) {
					for (int i = 0; i < f.getNumChildren() + 2; ++i)
						//System.out.print("\t1 1 1 1 1 1 1 1 1 1");
						sb.append("\t1 1 1 1 1 1 1 1 1 1");
				} else {
					int tableSize = p.length;
					if (f.isFatherInformative(m1) && f.isMotherInformative(m1)) {
						//System.out.print("\t0 1 0 0 0 0 0 0 0 0");
						//System.out.print("\t0 1 0 0 0 0 0 0 0 0");
						sb.append("\t0 1 0 0 0 0 0 0 0 0\t0 1 0 0 0 0 0 0 0 0");
						for (int i = 0; i < tableSize; i+=4) {
							double r0 = p[i];
							double r2 = p[i + 3];
							double r1 = 1 - r0 - r2;
							double max = Math.max(r0, Math.max(r1, r2));
							//System.out.print("\t" + r0 / max + " " + r1 / max + " 0 0 " + r2/max + " 0 0 0 0 0");
							sb.append("\t" + r0 / max + " " + r1 / max + " 0 0 " + r2/max + " 0 0 0 0 0");
						}
					}
					else if (f.isFatherInformative(m1)) {
						//System.out.print("\t0 1 0 0 0 0 0 0 0 0");
						//System.out.print("\t1 0 0 0 0 0 0 0 0 0");
						sb.append("\t0 1 0 0 0 0 0 0 0 0\t1 0 0 0 0 0 0 0 0 0");
						for (int i = 0; i < tableSize; i+=4) {
							double r = (p[i] + p[i + 1]);
							double max = Math.max(r, 1 - r);
							sb.append("\t" + r / max + " " + (1 - r) / max + " 0 0 0 0 0 0 0 0");
							//System.out.print("\t" + r / max + " " + (1 - r) / max + " 0 0 0 0 0 0 0 0");
						}

					}
					else if (f.isMotherInformative(m1)) {
						//System.out.print("\t1 0 0 0 0 0 0 0 0 0");
						//System.out.print("\t0 1 0 0 0 0 0 0 0 0");
						sb.append("\t1 0 0 0 0 0 0 0 0 0\t0 1 0 0 0 0 0 0 0 0");
						for (int i = 0; i < tableSize; i+=4) {
							double r = (p[i] + p[i + 2]);
							double max = Math.max(r, 1 - r);
							sb.append("\t" + r / max + " " + (1 - r) / max + " 0 0 0 0 0 0 0 0");
							//System.out.print("\t" + r / max + " " + (1 - r) / max + " 0 0 0 0 0 0 0 0");
						}
					}
				}
				
			}
			//System.out.println();
			System.out.println(sb);
		}
	}
	
	private final int MAX_TABLE_SIZE = 1000000;
	
	private final int RAND_TABLE_SIZE = 1024;
	
	private class SeparateThread implements Runnable{
		private ArrayList<Integer> infoMarkers;
		private double lodLimit;
		private boolean print;

		private int numResults = 0;
		private int resultsM1[] = new int[MAX_TABLE_SIZE];
		private int resultsM2[] = new int[MAX_TABLE_SIZE];
		private int prevR1 = 0;
		private int prevR2 = 0;
		private int rand[] = null;
		private boolean sample = false;
		
		private int m1i, m2i, start, numThreads;
		
		int oldProgress = 0;
		
		
		public SeparateThread(ArrayList<Integer> markers, int start, int numThreads, double lodLimit, boolean print, double samplingRate)
		{
			this.infoMarkers = markers;
			this.lodLimit = lodLimit;
			this.print = print;
			this.start = start;
			this.numThreads = numThreads;
			this.m1i = 0;
			this.m2i = start + 1;
			this.sample = false;

			if (samplingRate < 1.0) {
				sample = true;
				rand = new int[RAND_TABLE_SIZE];
				for (int i = 0; i < RAND_TABLE_SIZE; ++i) {
					double r = 0.0;
					while (r < 1e-20)
						r = Math.random();
					int j = 1;
					while (r >= samplingRate) {
						r *= (1 - samplingRate);
						++j;
					}
					rand[i] = j * numThreads;
				}
			}
			
		}
		
		
		public void run()
		{
			int rindex = 0;
			numResults = 0;
			while (m1i < infoMarkers.size() - 1) {
				int m1 =  infoMarkers.get(m1i); 
				int c1 = chromosomes.find(m1);
				if (print) {
					int newProgress = 10 * m1i / infoMarkers.size();
					for (int i = oldProgress; i < newProgress; ++i)
						System.err.print(i + 1);
					oldProgress = newProgress;
				}
				while (m2i < infoMarkers.size()) {
					int m2 = infoMarkers.get(m2i);
					int c2 = chromosomes.find(m2);
					if (c1 != c2 && ((prevR1 != c1 || prevR2 != c2) && (prevR1 != c2 || prevR2 != c1))) {
						double lod = computeLOD(m1, m2);
						if (lod >= lodLimit) {
							resultsM1[numResults] = c1;
							resultsM2[numResults] = c2;
							prevR1 = c1;
							prevR2 = c2;
							++numResults;
							if (numResults == MAX_TABLE_SIZE) {
								processResults();
								c1 = chromosomes.find(m1);
							}
							//if (numResults == MAX_TABLE_SIZE)
							//	return;
						}
					}
					if (sample)
						m2i += rand[(rindex++) & (RAND_TABLE_SIZE - 1)];
					else
						m2i += numThreads;
				}
				++m1i;
				m2i = start + m1i + 1;
			}

			processResults();
			
			if (print)
				System.err.println(" done!");			
		}

		public void processResults()
		{
			lock.lock();
			try {
				for (int r = 0; r < numResults; ++r)
					chromosomes.union(resultsM1[r], resultsM2[r]);
				numResults = 0;
			} finally {
				lock.unlock();
			}
		}
	}
	
	public LGMap2 separateChromosomes(double lodLimit, int numThreads, double samplePairs)
	{
		chromosomes = new UnionFind(getNumMarkers());
		ArrayList<Integer> infoMarkers = new ArrayList<Integer>();
		
		for (int m1 = 0; m1 < getNumMarkers(); ++m1) {
			boolean informative = false;
			for (Family2 f : data.getFamilies())
				if (f.isFatherInformative(m1) || f.isMotherInformative(m1)) {
					informative = true;
					break; 
				}			
			
			if (informative) {
				double maxLODScore = 0.0; 
				for (Family2 f : data.getFamilies())
					maxLODScore += f.maxLodScore(theta1, theta2, m1);
				if (maxLODScore >= lodLimit)
					infoMarkers.add(m1);
			}
		}
		System.err.println("computing pairwise LOD scores");

		if (samplePairs < 0.0001) {
			samplePairs = 0.0001;
			System.err.println("Warning: samplePairs rounded to " + samplePairs);
		}
		
		SeparateThread sepT[] = new SeparateThread[numThreads];
		for (int t = 0; t < numThreads; ++t)
			sepT[t] = new SeparateThread(infoMarkers, t, numThreads, lodLimit, t == 0, samplePairs);
		
		Thread threads[] = new Thread[numThreads];

		if (numThreads > 1) {
			for (int t = 0; t < numThreads; ++t) {
				threads[t] = new Thread(sepT[t]);
				threads[t].start();
			}
			try {
				for (int t = 0; t < numThreads; ++t)
					if (threads[t].isAlive())
						threads[t].join();
			} catch (Exception e) {
				e.printStackTrace();
				Error.error(-99999999);
			}
		} else
			sepT[0].run();
		System.err.println(" done!");
	
		int sizes[] = new int[getNumMarkers()];
		for (int m = 0; m < getNumMarkers(); ++m)
			++sizes[chromosomes.find(m)];

		int names[] = new int[getNumMarkers()];
		int numChromosomes = 0;
		int numSingles = 0;
		for (int m = 0; m < getNumMarkers(); ++m) {
			if (sizes[m] > 1)
				names[chromosomes.find(m)] = ++numChromosomes;
			else if (sizes[m] == 1)
				++numSingles;
		}
		System.err.println("number of LGs = " + numChromosomes + " singles = " + numSingles);
		LGMap2 lGMap = new LGMap2(getNumMarkers());
		for (int m = 0; m < getNumMarkers(); ++m)
			lGMap.setLGName(m, names[chromosomes.find(m)]);
		lGMap.renameLGs();
	
		return lGMap;
	}

	public void setMinLod(double parseDouble) {
		// TODO Auto-generated method stub
		minLod = parseDouble;
	}
	

}
