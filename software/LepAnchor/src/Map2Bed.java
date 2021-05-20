import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

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
public class Map2Bed {
	private int markerSupport = 2;
	private int splitLength = 10000;
	private double minQuality = 2;
	
	public static void usage()
	{
		System.err.println("Usage: java Map2Bed map=contig_pos_map.txt contigLengths=contig.sizes [options] >map_clean.bed");
        System.err.println("       map=file           a map file containing columns contig, position and chromosome ");
        System.err.println("                          in sorted order (typically from CleanMap)");

        System.err.println("       contigLength=file  a file containing columns contig name and its length");
        
        System.err.println("       markerSupport=NUM  at least NUM markers are needed for splitting a contig [2]");
        System.err.println("       minSplitLength=NUM do not split shorter regions than NUM [10000]");
        System.err.println("       minQuality=NUM     only consider markers with quality >= NUM [2]");
        System.err.println("                          (column 5 from CleanMap)");
	}
	private void makeBed(String mapFile, String lengthFile)
	{
		
		HashMap<String, Long> scaffoldLength = new HashMap<String, Long>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(lengthFile));
			
			ArrayList<String> line = Input.loadTableRow(br, " \t");
			while (line != null) {
				String scaffold = line.get(0);
				long length = Long.parseLong(line.get(1));
				scaffoldLength.put(scaffold, length);
				line = Input.loadTableRow(br, " \t");
			}
		} catch (Exception e) {
			System.err.println("Error loading contig lengths");
			e.printStackTrace();
			System.exit(-1);
		}
		
		int numChromosomes = 0;
		String prevScaffold = "";
		ArrayList<Long> list = new ArrayList<Long>();
		
		HashMap<String, Long> scaffoldMap = new HashMap<String, Long>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(mapFile));
			
			ArrayList<String> line = Input.loadTableRow(br, " \t");
			while (line != null) {
				if (!line.get(0).equals("CHR") && !line.get(0).equals("CHROM")) { 
					String scaffold = line.get(0);
					long position = Long.parseLong(line.get(1));
					int chr = Integer.parseInt(line.get(2));
					double quality = minQuality;
					if (line.size() >= 5)
						quality = Double.parseDouble(line.get(4));
					if (chr > 0 && quality >= minQuality) {
						if (scaffoldMap.containsKey(scaffold) && scaffoldMap.get(scaffold) > position) {
							System.err.println("Error: map file must be sorted");
							System.exit(-1);
						}
						if (!scaffoldLength.containsKey(scaffold)) {
							System.err.println("Error: contig " + scaffold + " does not have length in " + lengthFile);
							System.exit(-1);
						}
							
						if (!prevScaffold.equals("") && !prevScaffold.equals(scaffold)) {
							processScaffold(prevScaffold, scaffoldLength.get(prevScaffold), list);
							list.clear();
						}
						scaffoldMap.put(scaffold, position);
						prevScaffold = scaffold;

						list.add(position);
						list.add((long)chr);
						list.add((long)(100 * quality));
						numChromosomes = Math.max(numChromosomes, chr);
					}
				}
				line = Input.loadTableRow(br, " \t");
			}
			processScaffold(prevScaffold, scaffoldLength.get(prevScaffold), list);
		} catch (Exception e) {
			System.err.println("Error loading map file");
			e.printStackTrace();
			System.exit(-1);
		}
	}
	private class Interval{
		int startIndex;
		int endIndex;
		long quality;
		long position;
		long startPosition;
		long chr;
		int markers;
		public Interval(int startIndex, long position, long quality, long chr){
			this.startIndex = startIndex;
			this.endIndex = startIndex;
			this.quality = quality;
			this.chr = chr;
			this.position = position;
			startPosition = position;
			markers = 1;
		}
		public void addMarker(int index, long position, long quality) {
			endIndex = index;
			this.quality += quality;
			this.position = position;
			++markers;
		}

		public void addMarkerLeft(int index, long position, long quality) {
			startIndex = index;
			this.quality += quality;
			this.startPosition = position;
			++markers;
		}

		public void addInterval(Interval i) {
			endIndex = i.endIndex;
			quality += i.quality;
			position = i.position;
			markers += i.markers;
		}
		public int getMarkers() {
			return markers;
		}
		public long getLength() {
			return position - startPosition + 1;
		}
		public long getQuality() {
			return quality;
		}
		public long getChr() {
			return chr;
		}
		public int getStartIndex() {
			return startIndex;
		}
		public int getEndIndex() {
			return endIndex;
		}
		public long getStartPosition() {
			return startPosition;
		}

		public long getPosition() {
			return position;
		}
		public String toString()
		{
			//return chr + "\t" + markers + "\t" + quality;
			return startPosition + "\t" + position + "\t" + chr;
		}
	}
	
	private void processScaffold(String scaffold, long scaffoldLength, ArrayList<Long> list)
	{
		
		//create intervals of adjacent markers in the same chromosome
		ArrayList<Interval> intervals = new ArrayList<Interval>();
		long prevChr = 0;
		for (int i = 0; i < list.size(); i+=3) {
			long position = list.get(i);
			long chr = list.get(i + 1);
			long quality = list.get(i + 2);
			if (i == 0 || prevChr != chr)
				intervals.add(new Interval(i, position, quality, chr));
			else
				intervals.get(intervals.size() - 1).addMarker(i, position, quality);
			prevChr = chr;
		}
		
		if (intervals.size() > 1) { //more than one interval, prune
			do {
				int removed = -1;
				for (int ii = 0; ii < intervals.size(); ++ii) {
					Interval i = intervals.get(ii);
					if (i.getMarkers() < markerSupport || i.getLength() < splitLength) { // < markerSupport markers and/or length < splitLength
						if (removed < 0)
							removed = ii;
						else {
							long l = i.getLength();
							long lr = intervals.get(removed).getLength();
							if (l < lr || (l == lr && i.getQuality() < intervals.get(removed).getQuality()))
								removed = ii;
						}
					}
				}
				
				if (removed >= 0) { //remove removed (too short and/or too few markers)
					if (removed > 0 && removed + 1 < intervals.size() && intervals.get(removed - 1).getChr() == intervals.get(removed + 1).getChr()) {
						Interval prev = intervals.get(removed - 1);
						Interval next = intervals.get(removed + 1);
						for (int li = prev.getEndIndex() + 3; li < next.getStartIndex(); li+=3) {
							long chr = list.get(li + 1);
							if (chr == prev.getChr()) {
								long position = list.get(li);
								long quality = list.get(li + 2);
								prev.addMarker(li, position, quality);
							}
						}
						prev.addInterval(next);	
						intervals.remove(removed + 1);
						intervals.remove(removed);
					} else
						intervals.remove(removed);
				}
				else
					break;
			} while (true);

			//try to extend intervals
			for (int ii = 0; ii + 1 < intervals.size(); ++ii) {
				Interval prev = intervals.get(ii);
				Interval next = intervals.get(ii + 1);
				
				for (int li = prev.getEndIndex() + 3; li < next.getStartIndex(); li+=3) {
					long chr = list.get(li + 1);
					if (chr == next.getChr()) // conflict
						break;
					if (chr == prev.getChr()) {
						//System.err.println("Extend");
						long position = list.get(li);
						long quality = list.get(li + 2);
						prev.addMarker(li, position, quality);
					}
				}

				for (int li = next.getStartIndex() - 3; li > prev.getEndIndex(); li-=3) {
					long chr = list.get(li + 1);
					if (chr == prev.getChr()) // conflict
						break;
					if (chr == next.getChr()) {
						//System.err.println("Extend");
						long position = list.get(li);
						long quality = list.get(li + 2);
						next.addMarkerLeft(li, position, quality);
					}
				}
			}
			
		}
		
		for (int ii = 0; ii < intervals.size(); ++ii) {
			Interval i = intervals.get(ii);
			long start = 1;
			if (ii > 0)
				start =  intervals.get(ii - 1).getPosition() + 1;

			String end = scaffoldLength + "*";
			if (ii + 1 < intervals.size())
				end = "" + (intervals.get(ii + 1).getStartPosition() - 1);
			System.out.println(scaffold + "\t" + start + "-" + i.getStartPosition() + "\t" + i.getPosition() + "-" + end + "\t?\t" + i.getChr());
		}
		
		//System.err.println(scaffold + "\t" +  intervals);
	}
	
	private void setMarkerSupport(int parseInt) {
		markerSupport = parseInt;
	}
	private void setSplitLength(int parseInt) {
		splitLength = parseInt;
	}
	private void setMinQuality(double parseDouble) {
		minQuality = parseDouble;
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
		
		
		String mapFile = pp.getValueAsString("map", null);
		String lengthFile = pp.getValueAsString("contigLength", null);
		if (mapFile == null) {
			System.err.println("Please specify a map file");
			usage();
            System.exit(0);
			
		}
		if (lengthFile == null) {
			System.err.println("Please specify a contig length file");
			usage();
            System.exit(0);
			
		}

		pp.warning(new String[]{"map", "contigLength", "markerSupport", "minSplitLength", "minQuality"});

		
		Map2Bed m2p = new Map2Bed();
		m2p.setMarkerSupport(Integer.parseInt(pp.getValueAsString("markerSupport", "2")));
		m2p.setSplitLength(Integer.parseInt(pp.getValueAsString("minSplitLength", "10000")));
		m2p.setMinQuality(Double.parseDouble(pp.getValueAsString("minQuality", "2.0")));

		System.out.println("#java Map2Bed" + extraParameters);
		m2p.makeBed(mapFile, lengthFile);
	}
}
