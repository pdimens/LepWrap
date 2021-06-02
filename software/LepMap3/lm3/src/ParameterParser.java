/**
    This file is part of Lep-MAP.

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

	Copyright (C) 2013 Pasi Rastas, pasi.rastas@helsinki.fi, University of Helsinki
*/
//Parses parameters given as "x=a y=b c d z = e"
import java.util.StringTokenizer;
import java.util.ArrayList;
import java.util.HashMap;

public class ParameterParser {
	private HashMap<String, ArrayList<String>> hm;

	public ParameterParser()
	{
		
	}
	public ParameterParser(String parameters)
	{
		if (!init(parameters))
			hm = null;
	}
	
	public boolean init(ArrayList<String> tokens) 
	{
		hm = new HashMap<String, ArrayList<String>>();
		int n = tokens.size();
		int i = 0;
		while (i < n) {
			String key = tokens.get(i);
			if (key.equals("="))
				return false;
			++i;
			if (i < n - 1)  {
				if (!tokens.get(i).equals("="))
					return false;
				++i;
			} else
				return false;
			
			ArrayList<String> value = new ArrayList<String>();
			while (i == n - 1 || (i < n - 1 && !tokens.get(i + 1).equals("="))) {
				String nv = tokens.get(i);
				if (nv.equals("="))
					return false;
				value.add(nv);
				++i;
			}
			if (hm.containsKey(key))
				return false;
			hm.put(key, value);
		}
		return true;
	}
	public boolean init(String parameters)
	{
		String delims = "\t =;";
		StringTokenizer st = new StringTokenizer(parameters, delims, true);
		
		ArrayList <String> tokens = new ArrayList <String>();  
		while (st.hasMoreTokens()) {
			String nt = st.nextToken();
			if (delims.indexOf(nt) >= 0) {
				if (nt.equals("=") || nt.equals(";"))
					tokens.add(nt);
			  } else
				  tokens.add(nt);
			
		}
		return init(tokens);
	}
	public ArrayList<String> getValue(String variable)
	{
		return hm.get(variable);
	}
	public String getValueAsString(String variable, String defaultValue)
	{
		ArrayList<String> ret = getValue(variable);
		if (ret != null) {
			String rets = "";
			for (String s : ret)
				rets += " " + s;
			return rets.substring(Math.min(1, rets.length())); // remove first space
		}
		return defaultValue;
	}
	public int getNumberOfValues(String variable)
	{
		ArrayList<String> ret = getValue(variable);
		if (ret != null)
			return ret.size();
		else
			return 0;
	}
	public String getValueAsString(String variable, int index, String defaultValue)
	{
		ArrayList<String> ret = getValue(variable);
		if (ret != null) {
			//if (ret.size() == 1 && index >= 0)
				//return ret.get(0);
			if (index + 1 > ret.size()) // this makes more sense?
				return defaultValue;
			else
				return ret.get(index);
		}
		return defaultValue;
	}
	public ArrayList<String> getValuesAsList(String variable)
	{
		return getValuesAsList(variable, 1, null);
	}
	
	public ArrayList<String> getValuesAsList(String variable, int size, String defaultValue)
	{
		ArrayList<String> ret = new ArrayList<String>(); 
		if (getNumberOfValues(variable) == 1 && getValueAsString(variable, null).startsWith("file:")) {
			ArrayList<ArrayList<String>> matrix = Input.loadTable(getValueAsString(variable, null).substring(5), "[\t ]");
			if (matrix.size() == 1)
				ret.addAll(matrix.get(0));
			if (matrix.size() > 1) { // transpose  matrix
				ArrayList<String> tp = new ArrayList<String>();
				boolean ok = true;
//				System.err.println("Reading matrix");
				for (ArrayList<String> als : matrix) {
					if (als.size() != 1)
						ok = false;
					else
						tp.add(als.get(0));
				}
				if (ok)
					ret.addAll(tp);
//				System.err.println("Reading matrix ok" + ok);
			}
		} else
			for (int i = 0; i < getNumberOfValues(variable); ++i)
				ret.add(getValueAsString(variable, i, null));
		if (ret.size() == 0)
			ret.add(defaultValue);
		if (ret.size() == 1)
			for (int i = 0; i < size - 1; ++i)
				ret.add(new String(ret.get(0)));
		return ret;
	}
	

	public ArrayList<ArrayList<String>> getValuesAsMatrix(String variable)
	{
		if (getNumberOfValues(variable) == 1 && getValueAsString(variable, null).startsWith("file:")) {
//			System.err.println("loading matrix " + getValueAsString(variable, null).substring(5));
			ArrayList<ArrayList<String>> matrix = Input.loadTable(getValueAsString(variable, null).substring(5), "[\t ]");
			return matrix;
		}
		ArrayList<ArrayList<String>> ret = new ArrayList<ArrayList<String>>();
		if (getNumberOfValues(variable) > 0)
			ret.add(new ArrayList<String>());
		int line = 0;

		for (int i = 0; i < getNumberOfValues(variable); ++i) {
			if (getValueAsString(variable, i, null).equals(";")) {
				++line;
				ret.add(new ArrayList<String>());
			} else
			ret.get(line).add(getValueAsString(variable, i, null));
		}
		return ret;
	}
	public boolean warning(String keyWords[])
	{
		HashMap<String, Integer> hm2 = new HashMap<String, Integer>();
		
		for (String s : hm.keySet())
			hm2.put(s, 1);
		
		for (String s : keyWords)
			hm2.remove(s);
		
		if (hm2.size() > 0) {
			System.err.println("Error! Unknown parameters: "); 
			for (String s : hm2.keySet())
				System.err.print(s + " ");
			System.err.println();
			System.exit(-1);
		}

		return false;
	}
	
	public static void main(String args[])
	{
		ParameterParser pp = new ParameterParser("a=40 10 30 b=10 d = 1 2;3 4 6 e=file:ped1");
		System.err.println("a = " + pp.getValue("a"));
		System.err.println("b = " + pp.getValue("b"));
		System.err.println("c = " + pp.getValue("c"));
		System.err.println("d = " + pp.getValue("d"));
		System.err.println("number of values for a = " + pp.getNumberOfValues("a"));
		System.err.println("number of values for b = " + pp.getNumberOfValues("b"));
		System.err.println("number of values for c = " + pp.getNumberOfValues("c"));
		System.err.println(pp.getValuesAsList("d"));
		System.err.println(pp.getValuesAsMatrix("d"));
		System.err.println(pp.getValuesAsMatrix("e"));
	}
	
}
