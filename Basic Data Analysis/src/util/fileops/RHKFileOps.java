package util.fileops;

import image.ImageEditing;
import image.imageIO;

import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import schrodinger.MovieMaker;
import util.ArrayOps;
import util.FieldOps;
import util.LayerUtil;
import util.Sort;
import util.TopomapUtil;
import util.fourier.FFTOps;
import util.geom.Distance;
import main.SRAW;

public class RHKFileOps {

	//INSTRUCTIONS:
	//To create a map bin file from a "topomap" directory, first put all the txt files in one folder.
	//Then run doRHKDir to create the folder index file and a bunch of useless image files (but important binary files?).
	//Then run writeTopoFile(dir) to create the topomap index file. Then go into the file and fix the bias.
	//Then run Topomap t = getTopoTopo(dir, the appendage of the file (usually "didv_r"))
	// and use Topomap.writeBin(t, fc);
	
	static String DIR;//Topomap.stddir;//RHKFileOpsBinary.rhkdir;	
	public static void main(String[] args)
	{
		Topomap.setStdDir();
//		DIR = RHKFileOpsBinary.rhkdir;
		DIR = Topomap.stddir;
		String dir = "C:\\Users\\Charlie\\Desktop\\spm data\\Practice\\processed data\\Topomap - Run #488 SnTe on Bi2Te3 90 mins Area 2 0T\\topomap";
		JFileChooser fc = new JFileChooser(Topomap.stddir);
//		dir = FileOps.selectDir(fc);
//		doRHKDir(dir);
//		writeTopoFile(dir);
//		Topomap t = getTopoTopo(dir, "didv_r");
//		Topomap.writeBIN(t, fc);
//		t = getTopoTopo(dir, "topo_r");
//		Topomap.writeBIN(t, fc);
//		Topomap t = getTopoTopo(dir, "dIdV");
//		Topomap.writeBIN(t, fc);
//		t = getTopoTopo(dir, "topo");
//		Topomap.writeBIN(t, fc);
//		String[] things = {"didv_r", "didv_f", "curr_r", "curr_f", "topo_r", "topo_f"};
//		String[] things = {/*"didv_r", "didv_f",*/ "topo"};
//		String[] things = {".txtdidv=back", ".txtdidv=for", ".txtcurr=back", ".txtcurr=for", ".txttopo=back", ".txttopo=for"};
//		String[] things = {"didv", "curr", "topo"};
//		killFirstEntriesAndHeader();
//		writeITX(t.data[390], "Map Layer", fc, 210, 175.5, 100);
//		writeDatDirToITX3(dir, "", 150, "topomap_906");
//		writeDatDirToITX3(dir, "", 25, 15, "topomap_513");
		
//		dir = "C:\\data\\analysis\\Cu Data\\Drift Test\\txt\\";
//		doRHKDir(dir);
//		Topomap t = getTopoFromTxtList(dir, "topo_f", false);
//		writeTopoFile(dir);
//		writeTopoFile("C:\\data\\analysis\\Cu Data\\211\\tm_1229\\txt\\");
//		Topomap t = getTopoTopo("C:\\data\\analysis\\Cu Data\\211\\tm_1229\\txt\\", "ch4_r");
//		Topomap.writeBIN(t, fc);
		
//		doOneFile();
//		doOneFileLayer();
//		convertPointSpecTXTtoBIN();
		convertMapTXTtoBIN();
		//		System.out.println(Math.acos(+1));
//		Topomap.setStdDir();
//		write3Files(new JFileChooser(Topomap.stddir), getFromGreyBMP());
//		double[][] data = getFromGreyBMP();
//		double[][] down4 = FieldOps.reduce(4, data);
//		RHKFileOps.write3Files(fc, down4);
		
//		readSM4File("C:\\data\\raw\\2013\\Run #300 PbSnSe0.3\\Mar.5,2012,5K\\didv map\\01-25-2013--5k-0328.SM4");
//
//		Layer t = Layer.openFree(fc);
//		doUserFitting(t);
//		Layer.writeBIN(t, fc);
	}

	//The ASCII file has 13 lines before the table. Also, each table line has some spaces
	//followed by the row number starting from zero, followed by the data.
	public static void makeAsciiTable(String ascii, String table, String dira, String dirt)
	{
		File asc = new File(dira + ascii);
		File out = new File(dirt + table);
		String line;
		PrintStream output = null;
		Scanner in = null;
		try {
			in = new Scanner(asc);
			output = new PrintStream(out);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int i = 0;
		for (i = 0; i < 13; i++)
			in.nextLine();
		
		i = 0;
		int log = 1, pow = 10;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			if (line.trim().equals(""))
				break;
			
			log = 1; pow = 10;
			while(pow <= i){
				pow *= 10;
				log++;
			}
			
			output.append(line.substring(line.indexOf("" + i) + log) + "\r\n");
			i++;
			
		}
	}
	public static double[][] getTable(String ascii, String dir)
	{
		File asc = new File(dir + ascii);
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int i = 0;
		for (i = 0; i < 13; i++)
			in.nextLine();
		
		i = 0;
		int log = 1, pow = 10;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			nlines++;
			line = in.nextLine();
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < 13; i++)
			in.nextLine();
		i = 0;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length - 1;
			if (data == null) data = new double[nperline][nlines];
			for (int j = 1; j < words.length; j++)
			{
				data[j-1][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		return data;
	}
	public static double[][] getTable(File asc)
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int i = 0;
		for (i = 0; i < 13; i++)
			in.nextLine();
		
		i = 0;
		int log = 1, pow = 10;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			nlines++;
			line = in.nextLine();
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < 13; i++)
			in.nextLine();
		i = 0;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length - 1;
			if (data == null) data = new double[nperline][nlines];
			for (int j = 1; j < words.length; j++)
			{
				data[j-1][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		return data;
	}
	public static Layer getLayer(File asc)
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int i = 0;
		for (i = 0; i < 13; i++)
			in.nextLine();
		
		i = 0;
		int log = 1, pow = 10;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			nlines++;
			line = in.nextLine();
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		double v = 0, current = 0;
		int nx = 0, ny = 0;
		double xscale = 1, yscale = 1;
		double[] x = null, y = null;
		String s; String[] tok;
		for (i = 0; i < 13; i++){
			line = in.nextLine();
			//read the header information in the most general possible form.
			line = line.trim();
			if (line.contains("bias"))
			{
				s = line.substring(line.indexOf("bias")+4, line.indexOf(","));
				if (s.endsWith("V"))
					v = Double.parseDouble(s.substring(0, s.length()-2));
				else if (s.endsWith("mV"))
					v = Double.parseDouble(s.substring(0, s.length()-3))/1000;
				else v = 0;
			}
			if (line.contains("current"))
			{
				s = line.substring(line.indexOf("current")+7);
				if (s.endsWith("pA"))
					current = Double.parseDouble(s.substring(0, s.length()-3))/1000;
				else if (s.endsWith("nA"))
					current = Double.parseDouble(s.substring(0, s.length()-3));
				else current = 0;
			}
			if (line.contains("page size"))
			{
				s = line.substring(9, line.indexOf("pixels"));
				s = s.trim();
				tok = s.split(" by ");
				nx = Integer.parseInt(tok[0]);
				ny = Integer.parseInt(tok[1]);
			}
			if (line.contains("page dimensions"))
			{
				s = line.substring("page dimensions ".length());
				s = s.trim();
				tok = s.split(" by ");
				if (tok[0].endsWith("nm"))
					xscale = Double.parseDouble(tok[0].substring(0, tok[0].length()-2).trim())*Math.pow(10, -9);
				else JOptionPane.showMessageDialog(null, "The unit of length is not nm!!! Lengths not set properly.");
				if (tok[1].endsWith("nm"))
					yscale = Double.parseDouble(tok[1].substring(0, tok[1].length()-2).trim())*Math.pow(10, -9);
			}
		}
		//now do some stuff:
		data = new double [nx][ny];
		x = new double [nx];
		y = new double [ny];
		for (int j = 0; j < nx; j++)
			x[j] = j*xscale/nx;
		for (int j = 0; j < ny; j++)
			y[j] = -j*yscale/ny;
		
		i = 0;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length - 1;
			if (data == null) data = new double[nperline][nlines];
			for (int j = 1; j < words.length; j++)
			{
				data[j-1][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		ArrayOps.flipX(data);
		return new Layer(data, x, y, v, current);
	}
	public static double[][] getTableDIDV(File asc)
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int nskiplines = 20;
		int i = 0;
		for (i = 0; i < nskiplines; i++)
			in.nextLine();
		
		i = 0;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			if(!line.trim().equals(""))
				nlines++;
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < nskiplines; i++)
			in.nextLine();
		i = 0;
		while(i < nlines)
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length-1;
			if (data == null) data = new double[nperline][nlines];
			for (int j = 1; j < words.length; j++)
			{
				data[j-1][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		return data;
	}
	
	public static void topoDIDVTextToBin()
	{
		JOptionPane.showMessageDialog(null, "Open the text file.");
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		Topomap t = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			t = getTopoDIDV(fc.getSelectedFile());
		}
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			Topomap.writeBIN(t, fc.getSelectedFile().toString() + ".bin");
		}
	}
	public static void topoDIDVTextToBins()
	{
		JOptionPane.showMessageDialog(null, "Open the text file.");
		JFileChooser fc = new JFileChooser("C:\\data\\analysis\\");
		Topomap t = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			t = getTopoDIDV(fc.getSelectedFile());
		}
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			Topomap.writeBINsImages(t, fc.getCurrentDirectory().toString(), false);
		}
	}
	/**
	 * Note: If the angle of the 
	 */
	public static void convertMapTXTtoBIN()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		try{
			Topomap t = getTopoDIDVWindow(FileOps.selectOpen(fc));
			Topomap.writeBIN(t, fco);
			JOptionPane.showMessageDialog(null, "Done.");
			System.exit(0);
		}
		catch (Exception e)
		{
			JOptionPane.showMessageDialog(null, e.toString());
			e.printStackTrace();
			System.exit(0);
		}
	}
	public static void convertPointSpecTXTtoBIN()
	{
		JFileChooser fc = new JFileChooser(DIR);
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		try{
			PointSpectra t = getPointSpectra(FileOps.selectOpen(fc));
			PointSpectra.writeBIN(t, fco);
			JOptionPane.showMessageDialog(null, "Done.");
			System.exit(0);
		}
		catch (Exception e)
		{
			JOptionPane.showMessageDialog(null, e.toString());
			System.exit(0);
		}
	}
	/**
	 * Deprecated
	 * @param asc
	 * @return
	 */
	public static Topomap getTopoDIDV(File asc)
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int nskiplines = 20;
		int i = 0;
		for (i = 0; i < nskiplines; i++)
			in.nextLine();
		
		i = 0;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			if(!line.trim().equals(""))
				nlines++;
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < nskiplines-2; i++)
			in.nextLine();
		//read the x and y lines
		line = in.nextLine();
		line = line.trim();
		words = line.split("\t");
		nperline = words.length-1;
		int size = (int)Math.sqrt(nperline);
		double[] x = new double [size], y = new double [size];
		for (int j = 0; j < nperline; j++)
		{
			x[j%size] = Double.parseDouble(words[j%size + 1]);
		}
		line = in.nextLine();
		line = line.trim();
		words = line.split("\t");
		nperline = words.length-1;
		for (int j = 0; j < nperline; j++)
		{
			y[j/size] = Double.parseDouble(words[j/size + 1]);
		}
		double[] v = new double [nlines];
		                         
		                         
		i = 0;
		while(i < nlines)
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length-1;
			if (data == null) data = new double[nperline][nlines];
			v[i] = Double.parseDouble(words[0]);
			for (int j = 1; j < words.length; j++)
			{
				data[j-1][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		size = (int)Math.sqrt(nperline);
		double[][][] topomap = new double[nlines][size][size];
		for (i = 0; i < nlines; i++)
			for (int j = 0; j < nperline; j++)
			{
				topomap[i][j%size][j/size]=data[j][i];
			}
		for (i = 0; i < nlines; i++)
			ArrayOps.flipX(topomap[i]);
		ArrayOps.flip(x);
		String[] names = new String[v.length];
		
		for (i = 0; i < v.length; i++)
			names[i] = MovieMaker.fromInt(i);
		
		return new Topomap(topomap, v, x, y, names);
//		return data;
	}
	/**
	 * This method will scale the results incorrectly if the angle of the map is not Zero degrees. This ought to be pointed out. 
	 * @param asc
	 * @return
	 */
	public static Topomap getTopoDIDVWindow(File asc)
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int nskiplines = 20;
		int i = 0;
		for (i = 0; i < nskiplines; i++)
			in.nextLine();
		
		i = 0;
		int nlines = 0;
		String[] words;
		JFrame jf = new JFrame();
		jf.setSize(350, 50);
		jf.setLocation(600, 500);
		jf.setTitle("Determining the number of lines...");
		jf.show();
		while(in.hasNextLine())
		{
			line = in.nextLine();
			if(!line.trim().equals(""))
				nlines++;
		}
		//done.
		int nperline;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < nskiplines-2; i++)
			in.nextLine();
		//read the x and y lines
		line = in.nextLine();
		line = line.trim();
		words = line.split("\t");
		nperline = words.length-1;
//		int size = (int)Math.sqrt(nperline);
		int size = (int)Math.sqrt(nperline);
		int sizex = size, sizey = size;
		boolean weirdDoubleMap = true;
		if ((double)size != Math.sqrt(nperline)) //assume it is wider rather than round
		{
//			if (Double.parseDouble(words[0]) == Double.parseDouble(words[1])){
//				sizex = (int)Math.sqrt(nperline*2);
//				sizey = (int)Math.sqrt(nperline/2);
//			}
			if (weirdDoubleMap) //This basically repeats the entire method 
			{
				i = 0;
				sizex = (int)Math.sqrt(nperline/2);
				sizey = (int)Math.sqrt(nperline/2);
				double[] x = new double [sizex], y = new double [sizey];
				for (int j = 0; j < nperline/2; j++)
				{
					x[j%sizex] = Double.parseDouble(words[(j*2)%sizex + 1]);
				}
				line = in.nextLine();
				line = line.trim();
				words = line.split("\t");
				nperline = words.length-1;
				for (int j = 0; j < nperline/2; j++)
				{
					y[j/sizex] = Double.parseDouble(words[j*2 + 1]);
				}
				double[] v = new double [nlines];
				double[][][] topomap = new double[nlines][sizex][sizey];
				double temp;
				while(i < nlines)
				{
					jf.setTitle("Reading " + i + " out of " + nlines);
					line = in.nextLine();
					line = line.trim();
					words = line.split("\t");
					nperline = words.length-1;
					v[i] = Double.parseDouble(words[0]);
					for (int j = 1; j < words.length-1; j+=2)
					{
						temp = (Double.parseDouble(words[j]) + Double.parseDouble(words[j+1]))/2;
						topomap[i][((j-1)/2)%sizex][((j-1)/2)/sizex] = temp;
					}
					i++;
				}
				jf.setTitle("Flipping images...");
				for (i = 0; i < nlines; i++)
					ArrayOps.flipX(topomap[i]);
				ArrayOps.flip(x);
				String[] names = new String[v.length];
				
				jf.setTitle("Writing strings...");
				for (i = 0; i < v.length; i++)
					names[i] = MovieMaker.fromInt(i);
				jf.hide();
				jf.dispose();
				
				//if the voltage is reversed, flip it:
				if (v[v.length-1] < v[0])
				{
					double[] tv = new double [v.length];
					for (i = 0; i < v.length; i++)
						tv[i] = v[v.length-1-i];
					v = tv;
					double[][][] layers = new double[v.length][][];
					for (i = 0; i < v.length; i++)
					{
						layers[i] = topomap[v.length-1-i];
					}
					topomap = layers;
				}
				
				return new Topomap(topomap, v, x, y, names);
				
			}
		}
		double[] x = new double [sizex], y = new double [sizey];
		for (int j = 0; j < nperline; j++)
		{
			x[j%sizex] = Double.parseDouble(words[j%sizex + 1]);
		}
		line = in.nextLine();
		line = line.trim();
		words = line.split("\t");
		nperline = words.length-1;
		for (int j = 0; j < nperline; j++)
		{
			y[j/sizex] = Double.parseDouble(words[j + 1]);
		}
		double[] v = new double [nlines];
		                         
		                         
		i = 0;
//		size = (int)Math.sqrt(nperline);
		double[][][] topomap = new double[nlines][sizex][sizey];
		
		while(i < nlines)
		{
			jf.setTitle("Reading " + i + " out of " + nlines);
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length-1;
			v[i] = Double.parseDouble(words[0]);
			for (int j = 1; j < words.length; j++)
			{
				topomap[i][(j-1)%sizex][(j-1)/sizex] = Double.parseDouble(words[j]);
			}
			i++;
		}
		jf.setTitle("Flipping images...");
		for (i = 0; i < nlines; i++)
			ArrayOps.flipX(topomap[i]);
		ArrayOps.flip(x);
		String[] names = new String[v.length];
		
		jf.setTitle("Writing strings...");
		for (i = 0; i < v.length; i++)
			names[i] = MovieMaker.fromInt(i);
		jf.hide();
		jf.dispose();
		
		//if the voltage is reversed, flip it:
		if (v[v.length-1] < v[0])
		{
			double[] tv = new double [v.length];
			for (i = 0; i < v.length; i++)
				tv[i] = v[v.length-1-i];
			v = tv;
			double[][][] layers = new double[v.length][][];
			for (i = 0; i < v.length; i++)
			{
				layers[i] = topomap[v.length-1-i];
			}
			topomap = layers;
		}
		
		return new Topomap(topomap, v, x, y, names);
//		return data;
	}
	public static PointSpectra getPointSpectra(File asc)
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int nskiplines = 20;
		int i = 0;
		for (i = 0; i < nskiplines; i++)
			in.nextLine();
		
		i = 0;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			if(!line.trim().equals(""))
				nlines++;
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < nskiplines-2; i++)
			in.nextLine();
		//read the x and y lines
		line = in.nextLine();
		line = line.trim();
		words = line.split("\t");
		nperline = words.length-1;
		double[] x = new double [nperline], y = new double [nperline];
		for (int j = 0; j < nperline; j++)
		{
			x[j] = Double.parseDouble(words[j + 1]);
		}
		line = in.nextLine();
		line = line.trim();
		words = line.split("\t");
		nperline = words.length-1;
		for (int j = 0; j < nperline; j++)
		{
			y[j] = Double.parseDouble(words[j + 1]);
		}
		double[] v = new double [nlines];
		                         
		                         
		i = 0;
		while(i < nlines)
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length-1;
			if (data == null) data = new double[nperline][nlines];
			v[i] = Double.parseDouble(words[0]);
			for (int j = 1; j < words.length; j++)
			{
				data[j-1][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		
		return new PointSpectra(data, v, x, y);
//		return data;
	}
	public static void doOneFile()
	{
//		JOptionPane.showMessageDialog(null, "Open the text file.");
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		double[][] data = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			data = RHKFileOps.getTable(fc.getSelectedFile());
		}
		int o = JOptionPane.showConfirmDialog(null, "To subtract the average value of each horizontal line, click Yes.\r\n" +
				"To subtract the average value of each vertical line, click No.\r\n" +
				"To not subtract anything, click Cancel.", "Line Delta Z Subtract", JOptionPane.YES_NO_CANCEL_OPTION);
		if (o == JOptionPane.YES_OPTION)
			FieldOps.subtactLineAvg(data, false);
		else if (o == JOptionPane.NO_OPTION)
			FieldOps.subtactLineAvg(data, true);
		else;
		
		o = JOptionPane.showConfirmDialog(null, "To flip the image horizontally, click Yes.\r\n" +
				"To flip vertically, click No.\r\n" +
				"To not flip at all, click Cancel.", "Flip the image.", JOptionPane.YES_NO_CANCEL_OPTION);
		if (o == JOptionPane.YES_OPTION)
			ArrayOps.flipX(data);
		else if (o == JOptionPane.NO_OPTION)
			ArrayOps.flipY(data);
		else;
			
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			write3Files(fc.getSelectedFile().toString(), data);
		}
	}
	public static void doOneFileLayer()
	{
//		JOptionPane.showMessageDialog(null, "Open the text file.");
		JFileChooser fc = new JFileChooser(DIR);
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		Layer l = Layer.openFree(fc);
//		int o = JOptionPane.showConfirmDialog(null, "To subtract the average value of each horizontal line, click Yes.\r\n" +
//				"To subtract the average value of each vertical line, click No.\r\n" +
//				"To not subtract anything, click Cancel.", "Line Delta Z Subtract", JOptionPane.YES_NO_CANCEL_OPTION);
//		if (o == JOptionPane.YES_OPTION)
//			FieldOps.subtactLineAvg(l.data, false);
//		else if (o == JOptionPane.NO_OPTION)
//			FieldOps.subtactLineAvg(l.data, true);
//		else;
		doUserFitting(l);
		doUserFlip(l);
		askUserOrigin(l);
		
		if (fco.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			Layer.writeBIN(l, fco.getSelectedFile().toString() + ".bin");
			write2Images(fco.getSelectedFile().toString(), l.data);
		}
	}
	public static void doMultiFileLayersMap()
	{
//		JOptionPane.showMessageDialog(null, "Open the text file.");
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		Layer l = Layer.openFree(fc);
		File dir = fc.getCurrentDirectory();
		
		String name;
		File[] files = dir.listFiles();
		File outdir = new File(dir.toString() + "\\output\\");
		if (!outdir.exists()) outdir.mkdir();
		Layer[] data = new Layer[files.length];
//		int o = JOptionPane.showConfirmDialog(null, "To subtract the average value of each horizontal line, click Yes.\r\n" +
//				"To subtract the average value of each vertical line, click No.\r\n" +
//				"To not subtract anything, click Cancel.", "Line Delta Z Subtract", JOptionPane.YES_NO_CANCEL_OPTION);
//		if (o == JOptionPane.YES_OPTION)
//			FieldOps.subtactLineAvg(l.data, false);
//		else if (o == JOptionPane.NO_OPTION)
//			FieldOps.subtactLineAvg(l.data, true);
//		else;
		int fit = getUserFittingChoice();
		double[] origin = getOriginFromUser();
		
		for (int i = 0; i < files.length; i++)
		{
			data[i] = Layer.openFree(files[i]);
			doFitting(data[i], fit);
			data[i].makeOrigin(origin);
//			if (i % 2 != 0)
//				data[i].flipX(); //This is assumptive
			Layer.writeBIN(data[i], outdir + "\\" + files[i].getName().substring(0, files[i].getName().length() - 4) + ".bin");
			RHKFileOps.write2Images(outdir + "\\" + files[i].getName().substring(0, files[i].getName().length() - 4), data[i].data);
		}
		
		Topomap t = Topomap.newTopomap(data);
		Topomap.writeBIN(t, fc);
		
//		doUserFitting(l);
//		doUserFlip(l);
	}
	
	public static void askUserOrigin(Layer l)
	{
		int o;
		o = JOptionPane.showConfirmDialog(null, "Would you like to specify an origin at this time?\r\n" +
				"(This allows you to set the location of the image)", "Set the origin?", JOptionPane.YES_NO_OPTION);
		if (o == JOptionPane.YES_OPTION)
			addOrigin(l);
	}

	public static void doUserFitting(Layer l)
	{
		int o =getUserFittingChoice();
		doFitting(l, o);
	}
	public static void doUserFitting(Topomap  t)
	{
		int o =getUserFittingChoice();
		doFittingTopo(t, o);
	}
	public static int getUserFittingChoice()
	{
		int o = Integer.parseInt(JOptionPane.showInputDialog(null,
				"Select a line correction method (enter the number):\r\n"+
						"0) Do nothing.\r\n" +
						"1) Subract average from horizontal lines.\r\n" +
						"2) Subtract average from vertical lines.\r\n" +
						"3) Subtract best-fit line from horizontal lines. \r\n" +
						"4) Subtract best-fit line from vertical lines. \r\n" +
						"5) Do 3 and 4, (3 first).\r\n" +
						"6) Fourier-filter the horizontal and vertical lines (like 5 but less drastic in FFT)\r\n" +
						"7) Fourier-filter the horizontal and vertical lines, but leave the overall average.\r\n" +
						"8) Fourier-filter vertical line only\r\n" +
						"9) Foutier-filter horizontal line only\r\n" + 
						"10) Subtract plane fit\r\n" +
						"11) Subtract parabolic sheet fit\r\n" +
						"12) Subtract overall average value\r\n" + 
						"13) Replace with the gradient magnitude\r\n" +
						"14) Replace with -1*laplacian\r\n" + 
						"15) Smooth with a gaussian\r\n" + 
						"16) Replace with autocorrelation\r\n" +
						"17) Replace each line with its FFT\r\n" +
						"18) Split by step edge and fit each shelf to a parabola\r\n" +
						"19) Replace with the cross-correlation with another layer to be opened (after taking Laplacian of both)\r\n" +
						"20) Subtract an N-polynomial sheet\r\n" + 
						"21) Subtract an N-polynomial line (for eliminating z-drift)\r\n" +  
						"22) Crop out everything outside a certain box\r\n" + 
						"23) Take the natural logarithm\r\n" + 
						"24) Subtract another layer\r\n" + 
						"25) Subtract gaussian smoothing\r\n" +
						"26) \"Suppress bad pixels\"\r\n" + 
						"27) Take inverse laplacian\r\n" + 
						"28) Do curvature thing (like laplacian but allegedly better)\r\n" +
						"29) Do some kind of Fourier-filter\r\n" +
						"30) Normalize to the range [0, 1]\r\n" + 
						"31) Take the radial derivative\r\n" + 
						"32) Get the local value of a Q-vector (specify length scale)\r\n" + 
						"33) Shift by a certain amount leaving blank space at the edge\r\n"
						, "Scan Correction", JOptionPane.OK_CANCEL_OPTION));
		return o;
	}
	public static void doFitting(Layer l, int o)
	{
		switch (o)
		{
		case 0: break;
		case 1:
			FieldOps.subtactLineAvg(l.data, false); break;
		case 2:
			FieldOps.subtactLineAvg(l.data, true); break;
		case 3:
			FieldOps.subtractLineFit(l.data, true); break;
		case 4:
			FieldOps.subtractLineFit(l.data, false); break;
		case 5:
			FieldOps.subtractLineFit(l.data, true); 
			FieldOps.subtractLineFit(l.data, false); break;
		case 6:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldCross(l.nx, l.ny, true)); break;
		case 7:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldCrossNotO(l.nx, l.ny, true)); break;
		case 8:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldLine(l.nx, l.ny, true, true)); break;
		case 9:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldLine(l.nx, l.ny, false, true)); break;
		case 10:
			FieldOps.subtractPlaneFit(l.data);
			break;
		case 11:
			FieldOps.subtractParabolicFit(l.data); break;
		case 12:
			FieldOps.subtractAvgFactor(l.data, 0.999); break;
		case 13: 
			FieldOps.replaceWithGradMag(l.data); break;
		case 14:
			l.data = FieldOps.laplace(l.data); FieldOps.negate(l.data); break;
		case 15:
			l.data = FieldOps.gaussSmooth(l.data, Double.parseDouble(JOptionPane.showInputDialog("Smoothing length in pixels?"))); break;
		case 16:
			l.data = FieldOps.getAutocorrelationFourier(l.data); break;
//			l.data = FieldOps.getJDOS_sameSize(l.data); break;
		case 17:
			l.data = FieldOps.getReplacedRowFFT(l.data, true, true);
			try{
				double lineTime = Double.parseDouble(JOptionPane.showInputDialog("Line time?"));
				double[] freq = new double [l.nx];
				for (int i = 0; i < l.nx; i++)
				{
					freq[i] = i/(lineTime);
				}
				l.x = freq;
				l.resetCoords();
			}
			catch(Exception e)
			{
				
			}
			break;
		case 18:
			LayerUtil.fitEachPlateauToParabolaAndZeroTheEdges(l, ImageEditing.loadFromImage(FileOps.selectOpen(null).toString(), java.awt.Color.WHITE));
			break;
//		case 19:
//			FieldOps.abs(l.data);
//			FieldOps.log(l.data);
//			break;
		case 19:
			l.data = FieldOps.getCrosscorrelationFourier(FieldOps.gradMag(l.data), FieldOps.gradMag(Layer.openFree().data));
			break;
		case 20:
			FieldOps.subtractNSheetFit(l.data, Integer.parseInt(JOptionPane.showInputDialog("Enter the power of the sheet.")));
//			l.data = FieldOps.getVorticity(l.data, Math.PI);
//			BufferedImage grad = ImageEditing.getBufferedImage(FieldOps.curl(FieldOps.gradientNM2(l.data)), null);
//			SimpleImageDrawer drawer = new SimpleImageDrawer(grad);
			break;
		case 21:
			FieldOps.subtractLinePolynomialFit(l.data, Integer.parseInt(JOptionPane.showInputDialog("Enter the power.")), JOptionPane.showConfirmDialog(null, "Insert blank space for the other layer?") == JOptionPane.YES_OPTION);
			break;
		case 22:
			int xi, xf, yi, yf;
			String i = JOptionPane.showInputDialog("Enter the coordinates of the upper-left corner separated by commas, in pixels.");
			String[] tok = i.split(",");
			xi = Integer.parseInt(tok[0].trim());
			yi = Integer.parseInt(tok[1].trim());
			i = JOptionPane.showInputDialog("Enter the coordinates of the lower-right corner separated by commas, in pixels.");
			tok = i.split(",");
			xf = Integer.parseInt(tok[0].trim());
			yf = Integer.parseInt(tok[1].trim());
			FieldOps.changeToAverageOutsideBox(l.data, xi, xf, yi, yf);
			break;
		case 23:
			FieldOps.log(l.data);
			break;
		case 24:
			double[][] aux = Layer.openFree().data;
			double factor = Double.parseDouble(JOptionPane.showInputDialog("Multiply by?", 1));
			FieldOps.timesEquals(aux, factor);
			FieldOps.minusEquals(l.data, aux);
			break;
		case 25:
			l.data = FieldOps.subtractGaussSmooth(l.data, Double.parseDouble(JOptionPane.showInputDialog("Smoothing length in pixels?"))); break;
		case 26:
			double d = Double.parseDouble(JOptionPane.showInputDialog("Smoothing length in pixels?"));
			double minperc = Double.parseDouble(JOptionPane.showInputDialog("Minimum distrubition cutoff (0 to 1)?"));
			double maxperc = Double.parseDouble(JOptionPane.showInputDialog("Minimum distrubition cutoff (0 to 1)?"));
			double[][][] smoothStuff = FieldOps.subtractGaussSmooth_ReturnBoth(l.data, d);
			double[][] flattened = FieldOps.copy(smoothStuff[0]);
			FieldOps.cutOffExtremes(flattened, minperc, maxperc);
			boolean save = JOptionPane.showConfirmDialog(null, "Save the differences?") == JOptionPane.YES_OPTION;
			if (save){
				double[][] diff = FieldOps.minus(smoothStuff[0], flattened);
				Layer.writeBIN(Layer.newLayer(l, diff), FileOps.selectSave(null).toString());
			}
			FieldOps.add(flattened, smoothStuff[1], l.data);
			break;
		case 27: //inverse laplace.
			l.data = FieldOps.inverseLaplaceFourier(l.data); FieldOps.negate(l.data); break;
		case 28:
			l.data = FieldOps.getCurvatureThing(l.data, Double.parseDouble(JOptionPane.showInputDialog("Enter the parameter (usual range 0 to a lot)" +
					"\r\nSee Rev. Sci. Instrum. 82, 043712 (2011)\r\n" +
					"Infinity is basically laplacian\r\n" +
					"Zero is some weird curvature thing\r\n", "" + 0.5)));
			FieldOps.negate(l.data);
			break;
		case 29:
		{
			int choice = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice:\r\n"
					+ "1 - Use a circle centered at the origin\r\n"
					+ "2 - Exclude a strip directed along a certain vector"));
			boolean[][] suppress = null;
			if (choice == 1)
				suppress = TopomapUtil.FourierFilterMethods.getSuppressionFieldCircle(Double.parseDouble(JOptionPane.showInputDialog("Enter the radius")), l.nx, l.ny);
			else if (choice == 2)
			{
				double[] vector = getTwoDoubles();
				double thickness = Double.parseDouble(JOptionPane.showInputDialog("Enter the thickness."));
				suppress = TopomapUtil.FourierFilterMethods.getSupressionFieldStripThroughOrigin(l.nx, l.ny, thickness, Distance.unitVector(vector), true);
			}
			else if (choice == 3)
			{
				int ncircles = Integer.parseInt(JOptionPane.showInputDialog(null, "How many circles?"));
				int[][] pts = new int [2*ncircles][];
				int[] radius = new int [2*ncircles];
				for (int j = 0; j < ncircles; j++)
				{
					radius[j] = Integer.parseInt(JOptionPane.showInputDialog("Enter the radius."));
					radius[j+ncircles] = radius[j];
					String inp = JOptionPane.showInputDialog("Enter the center of the circle, w.r.t. the center.");
					pts[j] = new int[] {Integer.parseInt(inp.split(",")[0]), Integer.parseInt(inp.split(",")[1])};
					pts[j+ncircles] = new int [] {-pts[j][0], -pts[j][1]};
				}
				suppress = TopomapUtil.FourierFilterMethods.getSuppressionFieldCircles(radius, pts, l.nx, l.ny);
			}
			if (JOptionPane.showConfirmDialog(null, "To include, click Yes. To exclude, click No.") == JOptionPane.YES_OPTION)
				FieldOps.negate(suppress);
			l.data = FFTOps.getFourierFilteredIFFT(l.data, suppress);
			break;
		}
		case 30:
		{
			double min = FieldOps.min(l.data);
			double max = FieldOps.max(l.data);
			FieldOps.minusEquals(l.data, min);
			FieldOps.timesEquals(l.data, 1/(max-min));
			break;
		}
		case 31:
		{
			double[][][] gradNM2 = FieldOps.gradientNM2(l.data);
			for (int k = 0; k < l.nx; k++)
				for (int j = 0; j < l.ny; j++)
				{
					double[] unitV = Distance.unitVector(new double[] {k-l.nx/2, j-l.ny/2});
					
					l.data[k][j] = gradNM2[k][j][0]*unitV[0] + gradNM2[k][j][1]*unitV[1];
				}
			break;
		}
		case 32:
		{
			String[] token = JOptionPane.showInputDialog("Enter the Q-vector peak w.r.t. center, in pixels, comma separated.").split(",");
			double[] q = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
			double length = Double.parseDouble(JOptionPane.showInputDialog("Enter the length scale in pixles"));
			l.data = FieldOps.getSignalAmplitude(l.data, q, length);
			break;
		}
		case 33:
		{
			double[][] temp = FieldOps.copy(l.data);
			String[] stuff = JOptionPane.showInputDialog("Enter the comma-separated integer vector.").split(",");
			int dx = Integer.parseInt(stuff[0]);
			int dy = Integer.parseInt(stuff[1]);
			FieldOps.shiftAvgExcluded(l.data, temp, dx, dy);
			l.data = temp;
			break;
		}
		}
		l.setMean();
	}
	
	public static double[] getTwoDoubles(){
		String i = JOptionPane.showInputDialog("Enter the vector separated by commas.");
		String[] tok = i.split(",");
		double xi = Double.parseDouble(tok[0].trim());
		double yi = Double.parseDouble(tok[1].trim());
		return new double[] {xi, yi};
	}
	public static void doFittingTopo(Topomap t, int o)
	{
		Layer l;

		double lineTime = 0;
		double[] freq = null;
		String chosenFile = "";
		switch(o)
		{
		case 15:
			double d = Double.parseDouble(JOptionPane.showInputDialog("Smoothing length in pixels?"));
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
				l.data = FieldOps.gaussSmooth(l.data, d);
				t.data[i] = l.data;
			}
			break;
		case 17:
			lineTime = Double.parseDouble(JOptionPane.showInputDialog("Line time?"));
			freq = new double [t.nx];
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
				l.data = FieldOps.getReplacedRowFFT(l.data, true, true);
				for (int j = 0; j < l.nx; j++)
				{
					freq[i] = i/(lineTime);
				}
				l.x = freq;
				l.resetCoords();
				t.data[i] = l.data;
			}
			break;
		case 18:
			boolean[][] isStep = ImageEditing.loadFromImage(FileOps.selectOpen(null).toString(), java.awt.Color.WHITE);
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
				LayerUtil.fitEachPlateauToParabolaAndZeroTheEdges(l, isStep);
				t.data[i] = l.data;
			}
			break;
		case 19:
			Layer l2 =Layer.openFree();	
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
		
				l.data = FieldOps.getCrosscorrelationFourier(FieldOps.gradMag(l.data), FieldOps.gradMag(l2.data));
				t.data[i] = l.data;
			}
			break;
		case 20:{
			int power = Integer.parseInt(JOptionPane.showInputDialog("Enter the power of the sheet."));
			double[][] fit  = new double [t.nlayers][];
			boolean savesheet = JOptionPane.showConfirmDialog(null, "Save the sheet?") == JOptionPane.YES_OPTION;
			File sv = savesheet ? FileOps.selectSave(null) : null;
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
				fit[i] = FieldOps.fitToNSheet(l.data, null, power);
				FieldOps.subtractNSheetFit(l.data, fit[i], power);
			}
			if (savesheet)
			{
				double[][][] sheet = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
					FieldOps.putNSheetFit(sheet[i], fit[i], power);
				
				Topomap.writeBIN(Topomap.newTopomap(t, sheet), sv.toString());
			}
			
			break;}
		case 21:{
			int power = Integer.parseInt(JOptionPane.showInputDialog("Enter the power of the fit."));
			boolean blankSpace = JOptionPane.showConfirmDialog(null, "Insert blank space for the other layer?") == JOptionPane.YES_OPTION;
				for (int i = 0; i < t.nlayers; i++)
					FieldOps.subtractLinePolynomialFit(t.data[i], power, blankSpace);
			break;}
		case 24:
			double[][] aux = Layer.openFree().data;
			double factor = Double.parseDouble(JOptionPane.showInputDialog("Multiply by?", 1));
			FieldOps.timesEquals(aux, factor);
			for (int i = 0; i < t.nlayers; i++)
				FieldOps.minusEquals(t.data[i], aux);
			break;

		case 25:
			double d1 = Double.parseDouble(JOptionPane.showInputDialog("Smoothing length in pixels?"));
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
				l.data = FieldOps.subtractGaussSmooth(l.data, d1);
				t.data[i] = l.data;
			}
			break;
		case 26:
			double d11 = Double.parseDouble(JOptionPane.showInputDialog("Smoothing length in pixels?"));
			double minperc = Double.parseDouble(JOptionPane.showInputDialog("Minimum distrubition cutoff (0 to 1)?"));
			double maxperc = Double.parseDouble(JOptionPane.showInputDialog("Maximum distrubition cutoff (0 to 1)?"));
			boolean save = JOptionPane.showConfirmDialog(null, "Save the differences?") == JOptionPane.YES_OPTION;
			double[][][] diff = null;
			if (save) diff = new double [t.nlayers][t.nx][t.ny];
			for (int i = 0; i < t.nlayers; i++){
				double[][][] smoothStuff = FieldOps.subtractGaussSmooth_ReturnBoth(t.data[i], d11);
				double[][] flattened = FieldOps.copy(smoothStuff[0]);
				FieldOps.cutOffExtremes(flattened, minperc, maxperc);
				FieldOps.add(flattened, smoothStuff[1], t.data[i]);
				if (save)
					diff[i] = FieldOps.minus(smoothStuff[0], flattened);
			}
			
			if (save){
				Topomap.writeBIN(Topomap.newTopomap(t, diff), FileOps.selectSave(null).toString());
			}
			break;
		case 28:
			double curvparam = Double.parseDouble(JOptionPane.showInputDialog("Enter the parameter (usual range 0 to a lot)" +
					"\r\nSee Rev. Sci. Instrum. 82, 043712 (2011)\r\n"			+		"Infinity is basically laplacian\r\n" +
					"Zero is some weird curvature thing\r\n", "" + 0.5));
			for (int i = 0; i < t.nlayers; i++){
				t.data[i] = FieldOps.getCurvatureThing(t.data[i], curvparam);
				FieldOps.negate(t.data[i]);
			}
			break;
		case 29:
		{
			int choice = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice:\r\n"
					+ "1 - Use a circle centered at the origin\r\n"
					+ "2 - Exclude a strip directed along a certain vector\r\n"
					+ "3 - Input filter circles by hand\r\n"));
			boolean[][] filter = null;
			if (choice == 1){ double rad = Double.parseDouble(JOptionPane.showInputDialog("Enter the radius"));
				filter = TopomapUtil.FourierFilterMethods.getSuppressionFieldCircle(rad, t.nx, t.ny);
				if (JOptionPane.showConfirmDialog(null, "To include, click Yes. To exclude, click No.") == JOptionPane.YES_OPTION)
					FieldOps.negate(filter);
				
			}
			else if (choice == 2)
			{
				double[] vector = getTwoDoubles();
				double thickness = Double.parseDouble(JOptionPane.showInputDialog("Enter the thickness."));
				filter = TopomapUtil.FourierFilterMethods.getSupressionFieldStripThroughOrigin(t.nx, t.ny, thickness, Distance.unitVector(vector), true);
			}
			else if (choice == 3)
			{
				int ncircles = Integer.parseInt(JOptionPane.showInputDialog(null, "How many circles?"));
				int[][] pts = new int [ncircles][];
				int[] radius = new int [ncircles];
				for (int i = 0; i < ncircles; i++)
				{
					radius[i] = Integer.parseInt(JOptionPane.showInputDialog("Enter the radius."));
					String inp = JOptionPane.showInputDialog("Enter the center of the circle, w.r.t. the center.");
					pts[i] = new int[] {Integer.parseInt(inp.split(",")[0]), Integer.parseInt(inp.split(",")[1])};
				}
				filter = TopomapUtil.FourierFilterMethods.getSuppressionFieldCircles(radius, pts, t.nx, t.ny);
			}
			
			if (filter != null) {
				for (int i = 0; i < t.nlayers; i++)
					t.data[i] = FFTOps.getFourierFilteredIFFT(t.data[i], filter);
				ImageEditing.copyToClipboard(ImageEditing.getBufferedImage(filter));
			}
			
			break;
		}
		case 32:
		{
			String[] token = JOptionPane.showInputDialog("Enter the Q-vector peak w.r.t. center, in pixels, comma separated.").split(",");
			double[] q = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
			double length = Double.parseDouble(JOptionPane.showInputDialog("Enter the length scale in pixles"));
			for (int i = 0; i < t.nlayers; i++){ t.data[i] = FieldOps.getSignalAmplitude(t.data[i], q, length); System.out.print(" " + i);}
			System.out.println();
			break;
		}
		case 33:
		{
			double[][] temp = FieldOps.copy(t.data[0]);
			String[] stuff = JOptionPane.showInputDialog("Enter the comma-separated integer vector.").split(",");
			int dx = Integer.parseInt(stuff[0]);
			int dy = Integer.parseInt(stuff[1]);
			for (int i = 0; i < t.nlayers; i++)
			{
				FieldOps.shiftAvgExcluded(t.data[i], temp, dx, dy);
				t.data[i] = FieldOps.copy(temp);
			}
			break;
		}

		default:
			for (int i = 0; i < t.nlayers; i++)
			{
				l = t.getLayer(i);
				doFitting(l, o);
				t.data[i] = l.data;
			}
		}
		
	}
	public static void doUserFlip(Layer l)
	{
		int o = JOptionPane.showConfirmDialog(null, "To flip the image horizontally, click Yes.\r\n" +
				"To flip vertically, click No.\r\n" +
				"To not flip at all, click Cancel.", "Flip the image.", JOptionPane.YES_NO_CANCEL_OPTION);
		doFlip(l, o);
	}
	public static void doFlip(Layer l, int o)
	{
		if (o == JOptionPane.YES_OPTION)
			l.flipX();
		else if (o == JOptionPane.NO_OPTION)
			l.flipY();
		else;
	}
	public static void doOperationsToLayer()
	{
//		JOptionPane.showMessageDialog(null, "Open the text file.");
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		Layer l = Layer.openFree(fc);
//		int o = JOptionPane.showConfirmDialog(null, "To subtract the average value of each horizontal line, click Yes.\r\n" +
//				"To subtract the average value of each vertical line, click No.\r\n" +
//				"To not subtract anything, click Cancel.", "Line Delta Z Subtract", JOptionPane.YES_NO_CANCEL_OPTION);
//		if (o == JOptionPane.YES_OPTION)
//			FieldOps.subtactLineAvg(l.data, false);
//		else if (o == JOptionPane.NO_OPTION)
//			FieldOps.subtactLineAvg(l.data, true);
//		else;
		int o = Integer.parseInt(JOptionPane.showInputDialog(null,
				"Select a line correction method (enter the number):\r\n"+
						"0) Do nothing.\r\n" +
						"1) Subract average from horizontal lines.\r\n" +
						"2) Subtract average from vertical lines.\r\n" +
						"3) Subtract best-fit line from horizontal lines. \r\n" +
						"4) Subtract best-fit line from vertical lines. \r\n" +
						"5) Do 3 and 4, (3 first).\r\n" +
						"6) Fourier-filter the horizontal and vertical lines (like 5 but less drastic in FFT)\r\n" +
						"7) Fourier-filter the horizontal and vertical lines, but leave the overall average.\r\n" +
						"8) Fourier-filter vertical line only\r\n" +
						"9) Foutier-filter horizontal line only\r\n" + 
						"10) Subtract plane fit"
						, "Scan line Correction", JOptionPane.OK_CANCEL_OPTION));
		switch (o)
		{
		case 0: break;
		case 1:
			FieldOps.subtactLineAvg(l.data, false); break;
		case 2:
			FieldOps.subtactLineAvg(l.data, true); break;
		case 3:
			FieldOps.subtractLineFit(l.data, true); break;
		case 4:
			FieldOps.subtractLineFit(l.data, false); break;
		case 5:
			FieldOps.subtractLineFit(l.data, true); 
			FieldOps.subtractLineFit(l.data, false); break;
		case 6:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldCross(l.nx, l.ny, true)); break;
		case 7:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldCrossNotO(l.nx, l.ny, true)); break;
		case 8:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldLine(l.nx, l.ny, true, true)); break;
		case 9:
			LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldLine(l.nx, l.ny, false, true)); break;
		case 10:
			FieldOps.subtractPlaneFit(l.data); break;
		}
		o = JOptionPane.showConfirmDialog(null, "To flip the image horizontally, click Yes.\r\n" +
				"To flip vertically, click No.\r\n" +
				"To not flip at all, click Cancel.", "Flip the image.", JOptionPane.YES_NO_CANCEL_OPTION);
		if (o == JOptionPane.YES_OPTION)
			l.flipX();
		else if (o == JOptionPane.NO_OPTION)
			l.flipY();
		else;
			
		o = JOptionPane.showConfirmDialog(null, "Would you like to specify an origin at this time?\r\n" +
				"(This allows you to set the location of the image)", "Set the origin?", JOptionPane.YES_NO_OPTION);
		if (o == JOptionPane.YES_OPTION)
			addOrigin(l);

		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			Layer.writeBIN(l, fc.getSelectedFile().toString() + ".bin");
			write2Images(fc.getSelectedFile().toString(), l.data);
		}
	}
	public static void write3Files(String dir, double[][] layer, String name)
	{
		ColumnIO.writeBin(layer, dir + name + ".dat");
		FFTOps.writeFFTBMPCent(dir + name + "fft", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(dir + name + ".bmp", layer);
	}
	public static void write3Files(String dir, String fftdir, double[][] layer, String name)
	{
		File f = new File (dir);
		if (!f.exists()) f.mkdir();
		f = new File(fftdir);
		if (!f.exists()) f.mkdir();
		
		ColumnIO.writeBin(layer, dir + name + ".dat");
		FFTOps.writeFFTBMPCent(fftdir + name + "fft", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(dir + name + ".bmp", layer);
	}
	public static void write3Files(String path, double[][] layer)
	{
		ColumnIO.writeBin(layer, path + ".dat");
		FFTOps.writeFFTBMPCent(path + "fft", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(path + ".bmp", layer);
	}
	public static void write2Images(String path, double[][] layer)
	{
//		ColumnIO.writeBin(layer, path + ".dat");
		FFTOps.writeFFTBMPCent(path + "fft", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(path + ".bmp", layer);
	}
	public static void write3Files(JFileChooser fc, double[][] layer)
	{
		String path = FileOps.selectSave(fc).toString();
		ColumnIO.writeBin(layer, path + ".dat");
		FFTOps.writeFFTBMPCent(path + "fft", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(path + ".bmp", layer);
	}
	public static void write3FilesBothFFT(JFileChooser fc, double[][] layer)
	{
		String path = FileOps.selectSave(fc).toString();
		ColumnIO.writeBin(layer, path + ".dat");
		FFTOps.writeFFTBMPCent(path + "fftlin", FFTOps.obtainFFT(layer), false, true);
		FFTOps.writeFFTBMPCent(path + "fftlog", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(path + ".bmp", layer);
	}
	public static void write3Files(JFileChooser fc, double[][] layer, boolean ftlog)
	{
		String path = FileOps.selectSave(fc).toString();
		ColumnIO.writeBin(layer, path + ".dat");
		FFTOps.writeFFTBMPCent(path + "fft", FFTOps.obtainFFT(layer), ftlog, true);
		SRAW.writeImage(path + ".bmp", layer);
	}
	public static void writeFilesFFTSeparate(String path, double[][] layer)
	{
		ColumnIO.writeBin(layer, path + ".dat");
//		FFTOps.writeFFTBMPCent(path + "fft", FFTOps.obtainFFT(layer), true, true);
		SRAW.writeImage(path + ".bmp", layer);
	}
	public static Topomap getTopoTopo(String dir, String appendage)
	{
		Scanner in = null;
		try {
			in = new Scanner(new File(dir + "Topomap.txt"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ArrayList<String> names = new ArrayList<String>();
		ArrayList<Double> v = new ArrayList<Double>();
		String line;
		String[] words;
		String[] namesArray;
		in.nextLine();
		double dx = 0, dy = 0;
		while (in.hasNextLine())
		{
			line = in.nextLine();
			words = line.split("\t");
			names.add(words[0].toString());
			v.add(Double.parseDouble(words[1]));
			dx = Double.parseDouble(words[2]);
			dy = Double.parseDouble(words[3]);
		}
		//transfer names to array
		namesArray = new String[names.size()];
		for (int i = 0; i < names.size(); i++)
			namesArray[i] = names.get(i);
		
		//Now read the files:
		double[][] temp = null;
		double[][][] data = null;
		double[] bias = null;
		double[] x = null, y = null;
		System.out.print("Files read (out of " + (names.size()) + "):");
		for (int i = 0; i < names.size(); i++)
		{
//			temp = ColumnIO.readSquareTable(dir + appendage + "\\" + names.get(i) + ".dat");	//if files are in split folders
//			temp = ColumnIO.readSquareTable(dir + names.get(i) + appendage + ".dat");
			temp = RHKFileOps.getLayer(new File(dir + names.get(i) + appendage + ".txt")).data;
//			temp = ColumnIO.readAllColumns(new File(dir + names.get(i) + appendage + ".TXT"), null);
			System.out.print("" + (i+1) + " ");
			if (data == null)
				data = new double[names.size()][temp.length][temp[0].length];
			if (bias == null)
				bias = new double[names.size()];
			if (x == null)
				x = new double [temp.length];
			if (y == null)
				y = new double [temp[0].length];
			FieldOps.copy(temp, data, i);
			bias[i] = v.get(i);
			if (Math.abs(bias[i]) > 10) bias[i]/=1000;
		}
		System.out.println();
		for (int i = 0; i < x.length; i++)
			x[i] = i*dx;
		for (int i = 0; i < y.length; i++)
			y[i] = i*dy;
		return new Topomap(data, bias, x, y, namesArray);
	}
	public static double[][][] getTopoMap(double[][] data)
	{
		int sizesq = data.length;
		int nlayers = data[0].length;
		int size = (int)Math.sqrt(sizesq);
		double[][][] topomap = new double[nlayers][size][size];
		for (int i = 0; i < nlayers; i++)
			for (int j = 0; j < sizesq; j++)
			{
				topomap[i][j%size][j/size]=data[j][i];
			}
		for (int i = 0; i < nlayers; i++)
			ArrayOps.flipX(topomap[i]);
		return topomap;
		
	}
	
	public static double[][] getTableSwitched(File asc) //this flips the x-axis
	{
		String line;
		Scanner in = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int i = 0;
		for (i = 0; i < 13; i++)
			in.nextLine();
		
		i = 0;
		int log = 1, pow = 10;
		int nlines = 0;
		String[] words;
		while(in.hasNextLine())
		{
			nlines++;
			line = in.nextLine();
		}
		//done.
		int nperline;
		double[][] data = null;
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (i = 0; i < 13; i++)
			in.nextLine();
		i = 0;
		while(in.hasNextLine())
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			nperline = words.length - 1;
			if (data == null) data = new double[nperline][nlines];
			for (int j = 1; j < words.length; j++)
			{
				data[(words.length-1) - j][i] = Double.parseDouble(words[j]);
			}
			i++;
		}
		return data;
	}
	public static void saveWithRHKRowStarts()
	{
		
	}
	
	public static void doIt(String file, String outdir, String name)
	{
		double[][] data = getTableSwitched(new File(file));
		ColumnIO.writeBin(data, outdir + name + ".dat");
		SRAW.writeImage(outdir + name + ".bmp", data);
		FFTOps.writeFFTBMPs(outdir + name + "fft", FFTOps.obtainFFT(data), true, true);
	}
	//Assumes that this is the only file in the directory.
	public static void doItDir1(String dir, String outdir, String name)
	{
		double[][] data = getTable(new File(dir).listFiles()[0]);
		ColumnIO.writeBin(data, outdir + name + ".dat");
		SRAW.writeImage(outdir + name + ".bmp", data);
		FFTOps.writeFFTBMPs(outdir + name + "fft", FFTOps.obtainFFT(data), true, true);
	}
	
	public static void doTxtFilesIn(String dir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name = "";
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".txt"))
			{
				data = getTable(files[i]);
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				ColumnIO.writeBin(data, name + ".dat");
				SRAW.writeImage(dir + name, data);
				FFTOps.writeFFTBMPCent(dir + name + "fft", FFTOps.obtainFFT(data), true, true);
			}
		}
	}
	public static void doRHKDir(String dir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name = "";
		String filetable = "File Table.txt";
		String content = "Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		int stupidSpaces = 0;
		int o = getUserFittingChoice();
		
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".txt"))
			{
				data = getTable(files[i]);
				Layer l = Layer.getFreeLayer(data);
				name = files[i].getName().substring(0, files[i].getName().length() - 4);
				if (name.contains("topo")){
					doFitting(l, o);
				}
				data = l.data;
				if (name.contains("ch4") || name.contains("ch5") || name.contains("didv") || name.contains("x"))
					stupidSpaces = 1;
				else if (name.contains("bias")  || name.contains("curr"))
					stupidSpaces = 2;
				else
					stupidSpaces = 0;
				content += name + "\t" + TopoInfo.getFromFile(files[i].toString(), stupidSpaces).toString() + "\r\n";
//				Layer.writeBIN(l, name + ".bin");
//				ColumnIO.writeBin(data, name + ".dat");
//				SRAW.writeImage(name, data);
//				FFTOps.writeFFTBMPCent(name + "fft", FFTOps.obtainFFT(data), true, true);
			}
		}
		ColumnIO.writeString(content, dir + filetable);
	}
	public static void doRHKDir(String dir, int[] toBeFlipped)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name = "";
		String filetable = "File Table.txt";
		String content = "Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		boolean flip = false;
		int stupidSpaces = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".txt") && !files[i].toString().contains(filetable))
			{
				flip = false;
				data = getTable(files[i]);
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				for (int j = 0; j < toBeFlipped.length; j++)
					if (name.endsWith("" + toBeFlipped[j])) flip = true;
				if(flip) ArrayOps.flipX(data);
				if (name.contains("didv") || name.contains("curr"))
					stupidSpaces = 1;
				else if (name.contains("bias"))
					stupidSpaces = 2;
				else
					stupidSpaces = 0;
				content += name.substring(name.lastIndexOf("\\")) + "\t" + TopoInfo.getFromFile(files[i].toString(), stupidSpaces).toString() + "\r\n";
				ColumnIO.writeBin(data, name + ".dat");
				SRAW.writeImage(name, data);
				FFTOps.writeFFTBMPCent(name + "fft", FFTOps.obtainFFT(data), true, true);
			}
		}
		ColumnIO.writeString(content, dir + filetable);
	}
	public static void doRHKDirSubDirs(String dir, String[] appendages)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name = "";
		String subname;
		String filetable = "File Table.txt";
		String content = "Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		boolean flip = false;
		for (int i = 0; i < appendages.length; i++)
			if (!(new File(dir + appendages[i] + "\\").exists()))
					new File(dir + appendages[i] + "\\").mkdir();
		int stupidSpaces = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".txt") && !files[i].toString().contains(filetable))
			{
				flip = false;
				data = getTable(files[i]);
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				for (int j = 0; j < appendages.length; j++)
					if (name.endsWith("" + appendages[j])){
						subname = name.substring(name.lastIndexOf("\\") + 1, name.length() - appendages[j].length());
						ColumnIO.writeBin(data, dir + appendages[j] + "\\" + subname + ".dat");
						SRAW.writeImage(dir + appendages[j] + "\\" + subname, data);
						FFTOps.writeFFTBMPCent(dir + appendages[j] + "\\" + subname + "fft", FFTOps.obtainFFT(data), true, true);
					}
				if (name.contains("didv") || name.contains("curr"))
					stupidSpaces = 1;
				else if (name.contains("bias"))
					stupidSpaces = 2;
				else
					stupidSpaces = 0;
				content += name.substring(name.lastIndexOf("\\")) + "\t" + TopoInfo.getFromFile(files[i].toString(), stupidSpaces).toString() + "\r\n";
			}
		}
		ColumnIO.writeString(content, dir + filetable);
	}
	public static void doRHKDirSubDirs(String dir, String[] appendages, String[] toBeFlipped, boolean flipxtrue)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name = "";
		String subname;
		String filetable = "File Table.txt";
		String filetable2 = "Topomap.txt";
		String content = "Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		ArrayList<Double> biases;
		ArrayList<String> tablenames;
		
		boolean flip = false;
		for (int i = 0; i < appendages.length; i++)
			if (!(new File(dir + appendages[i] + "\\").exists()))
					new File(dir + appendages[i] + "\\").mkdir();
		
		int stupidSpaces = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".txt") && !files[i].toString().contains(filetable)&& !files[i].toString().contains(filetable2))
			{
				flip = false;
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				for (int j = 0; j < appendages.length; j++)
					if (name.endsWith("" + appendages[j])){
						data = getTable(files[i]);
						flip = false;
						subname = name.substring(name.lastIndexOf("\\") + 1, name.length() - appendages[j].length());
						for (int k = 0; k < toBeFlipped.length; k++)
							if (subname.endsWith("" + toBeFlipped[k])) flip = true;
						if (flip && flipxtrue) ArrayOps.flipX(data);
						else if (flip && !flipxtrue) ArrayOps.flipY(data);
						ColumnIO.writeBin(data, dir + appendages[j] + "\\" + subname + ".dat");
						SRAW.writeImage(dir + appendages[j] + "\\" + subname, data);
						FFTOps.writeFFTBMPCent(dir + appendages[j] + "\\" + subname + "fft", FFTOps.obtainFFT(data), true, true);
					}
				if (name.contains("didv") || name.contains("curr") || name.contains("bias"))
					stupidSpaces = 1;
				else if (name.contains("phse"))
					stupidSpaces = 2;
				else
					stupidSpaces = 0;
				content += name.substring(name.lastIndexOf("\\")+1) + "\t" + TopoInfo.getFromFile(files[i].toString(), stupidSpaces).toString() + "\r\n";
			}
		}
		ColumnIO.writeString(content, dir + filetable);
		writeTopoFile(dir);
	}
	public static void writeTopoFile(String dir)
	{
		File table = new File(dir + "File Table.txt");
		Scanner s = null;
		try {
			s = new Scanner(table);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		s.nextLine();
		ArrayList<Double> v = new ArrayList<Double>();
		ArrayList<Object> name = new ArrayList<Object>();
		String line;
		String[] words;
		double d;
		double dx= 0, dy = 0;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			words = line.split("\t");
			d = Double.parseDouble(words[5]);
			if (v.contains(d));
			else
			{
				v.add(d);
				name.add(words[0].substring(0, words[0].length()-6));
			}
			dx = Double.parseDouble(words[3])/Double.parseDouble(words[1]);
			dy = Double.parseDouble(words[4])/Double.parseDouble(words[2]);
		}
		Sort.bubbleSort(v, name);
		String content = "Name \tBias \tnm/pixel x \t nm/pixel y \r\n";
		for (int i = 0; i< v.size(); i++)
			content += name.get(i) + "\t" + v.get(i) + "\t" + dx + "\t" + dy + "\r\n";
		ColumnIO.writeString(content, dir + "Topomap.txt");
	}
	
	//This is to change the "order" of the files
	public static void doRHKDirSubDirs(String dir, String[] appendages, String[] toBeFlipped, boolean flipxtrue, String[] newnames)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name = "";
		String subname = null, altsubname = "";
		String filetable = "File Table.txt";
		String content = "Name\tNew Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		int n = 0, total = 0;
		boolean flip = false;
		for (int i = 0; i < appendages.length; i++)
			if (!(new File(dir + appendages[i] + "\\").exists()))
					new File(dir + appendages[i] + "\\").mkdir();

		int stupidSpaces = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".txt") && !files[i].toString().contains(filetable))
			{
				flip = false;
				data = getTable(files[i]);
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				for (int j = 0; j < appendages.length; j++)
					if (name.endsWith("" + appendages[j])){
						flip = false;
						subname = newnames[n];
						altsubname = name.substring(name.lastIndexOf("\\") + 1, name.length() - appendages[j].length());
						for (int k = 0; k < toBeFlipped.length; k++)
							if (altsubname.endsWith("" + toBeFlipped[k])) flip = true;
						if (flip && flipxtrue) ArrayOps.flipX(data);
						else if (flip && !flipxtrue) ArrayOps.flipY(data);
						ColumnIO.writeBin(data, dir + appendages[j] + "\\" + subname + ".dat");
						SRAW.writeImage(dir + appendages[j] + "\\" + subname, data);
						FFTOps.writeFFTBMPCent(dir + appendages[j] + "\\" + subname + "fft", FFTOps.obtainFFT(data), true, true);
					}
				if (name.contains("didv") || name.contains("curr"))
					stupidSpaces = 1;
				else if (name.contains("bias"))
					stupidSpaces = 2;
				else
					stupidSpaces = 0;
				content += name.substring(name.lastIndexOf("\\")) + "\t" + subname + "\t" + TopoInfo.getFromFile(files[i].toString(), stupidSpaces).toString() + "\r\n";
				total++;
				if (total == appendages.length) {total = 0; n++;}
			}
		}
		ColumnIO.writeString(content, dir + filetable);
	}
	public static void killButTxt(String dir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		for (int i = 0; i < files.length; i++)
		{
			if (!files[i].toString().endsWith(".txt"))
				files[i].delete();
		}
	}
	public static class TopoInfo
	{
		double lx, ly; //in nm
		double bias; //mv
		double setpoint; //pA.
		int px, py;
		
		public TopoInfo(double lx, double ly, double bias, double setpoint,
				int px, int py) {
			super();
			this.lx = lx;
			this.ly = ly;
			this.bias = bias;
			this.setpoint = setpoint;
			this.px = px;
			this.py = py;
		}
		
		public String toString()
		{
			return "" + px + "\t" + py + "\t" + lx + "\t" + ly + "\t" + bias + "\t" + setpoint;
		}
		public double[] toDarray()
		{
			return new double[] {lx, ly, bias, setpoint};
		}
		public static TopoInfo getFromFile(String filename, int nspaces)
		{
			Scanner in = null;
			File asc = new File(filename);
			try {
				in = new Scanner(asc);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			for (int i = 0; i < 4; i++)
				in.nextLine();
			
			String[] line5 = in.nextLine().split(" ");
			int i = 2, j = 5;
			for (int k = 0; k < nspaces; k++)
				{i++; j++;}
			double bias = Double.parseDouble(line5[i]);
			double setpoint = Double.parseDouble(line5[j]);
			
			in.nextLine();
			String[] line7 = in.nextLine().split(" ");
		
			int px = Integer.parseInt(line7[2]), py = Integer.parseInt(line7[4]);
			System.out.println(bias + "\t" + setpoint + "\t" + px + "\t" + py);
			String[] line8 = in.nextLine().split(" ");
			//assumy always nm:
			double lx = Double.parseDouble(line8[2]), ly = Double.parseDouble(line8[5]);
			
			return new TopoInfo(lx, ly, bias, setpoint, px, py);
		}
	}
	public static double[] readTopoInfo(String topoInfo)
	{
		double[] result = new double [4];
		Scanner in = null;
		File asc = new File(topoInfo);
		try {
			in = new Scanner(asc);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (int i = 0; i < 4; i++)
			result[i] = in.nextDouble();
		return result;
	}
	
	public static String getITX(double[][] data, double topolengthX, double topolengthY, String name)
	{
		int N = data.length, M = data[0].length;
		String content = "";
		content += "IGOR" + "\r\n";
		content += "WAVES/D/N=(" + N + "," + M + ") " + name + "\r\n";
		content += "BEGIN/r/n";
		for (int j = 0; j < M; j++){System.out.print(j);
			for (int i = 0; i < N-1; i++)
				content += "" + data[i][j] + "\t";
			content += data[N-1][j] + "\r\n";
		}
		content += "END\r\n";
		content += "X SetScale/P x 0," + topolengthX + ",\"\", "+name + "; SetScale/P y 0," + topolengthY + ",\"\", "+name + "; SetScale d 0,0,\"\", " + name; 
		return content;
	}
	public static String getITX3D(double[][][] data, double topolengthX, double topolengthY, String name)
	{
		int ntps = data.length, N = data[0].length, M = data[0][0].length;
		String content = "";
		content += "IGOR" + "\r\n";
		content += "WAVES/D/N=(" + N + "," + M + "," + ntps + ") " + name + "\r\n";
		content += "BEGIN/r/n";
		for (int k = 0; k < ntps; k++)
			for (int j = 0; j < M; j++){System.out.print(j);
				for (int i = 0; i < N-1; i++)
					content += "" + data[k][i][j] + "\t";
			content += data[k][N-1][j] + "\r\n";
		}
		content += "END\r\n";
		content += "X SetScale/P x 0," + topolengthX + ",\"\", "+name + "; SetScale/P y 0," + topolengthY + ",\"\", "+name + "; Setscale z/P 0,1,\"\"," + name + "; SetScale d 0,0,\"\", " + name; 
		return content;
	}
	public static String getIgorNote(String name, double bias, double setpoint, String biasU, String setU)
	{
		return "X Note w" + name + "," + "\"\\rSet V(" + biasU + ")="+bias +"\\rSet I("+setU+")="+setpoint+"\\r\"";
	}
	public static void writeITX(double[][] data, String dir, String name, double topolength, double bias, double setpoint)
	{
		int N = data.length, M = data[0].length;
//		String ans = getITX(data, topolength, topolength, name) + "\r\n";
		
		String header = "";
		header += "IGOR" + "\r\n";
		header += "WAVES/D/N=(" + N + "," + M + ") w" + name + "\r\n";
		header += "BEGIN\r\n";
//		ans += getIgorNote(name, bias, setpoint, "mV", "pA");
//		ColumnIO.writeString(ans, dir + name + ".itx");
		PrintStream firstStream = null;
		String tail = "";
		tail += "END\r\n";
		tail += "X SetScale/P x 0," + topolength + ",\"\", w"+name + "; SetScale/P y 0," + topolength + ",\"\", w"+name + "; SetScale d 0,0,\"\", w" + name; 
		tail += "\r\n" + getIgorNote(name, bias, setpoint, "mV", "pA");
		
		try{
			firstStream = new PrintStream(new File(dir + name + ".itx"));
		}
		catch (FileNotFoundException e)
		{
			System.out.println("No new file.");
		}
		firstStream.append(header);
		for (int i = 0; i < M; i++){
			for (int j = 0; j < N-1; j++)
				firstStream.append(data[j][i] + "\t");
			firstStream.append(data[N-1][i]+"\r\n");
		}
		firstStream.append(tail);
		firstStream.close();
		
	}
	public static void writeITX(double[][] data, String name, JFileChooser fc, double topolength, double bias, double setpoint)
	{
		int N = data.length, M = data[0].length;
//		String ans = getITX(data, topolength, topolength, name) + "\r\n";
		
		String header = "";
		header += "IGOR" + "\r\n";
		header += "WAVES/D/N=(" + N + "," + M + ") w" + name + "\r\n";
		header += "BEGIN\r\n";
//		ans += getIgorNote(name, bias, setpoint, "mV", "pA");
//		ColumnIO.writeString(ans, dir + name + ".itx");
		PrintStream firstStream = null;
		String tail = "";
		tail += "END\r\n";
		tail += "X SetScale/P x 0," + topolength + ",\"\", w"+name + "; SetScale/P y 0," + topolength + ",\"\", w"+name + "; SetScale d 0,0,\"\", w" + name; 
		tail += "\r\n" + getIgorNote(name, bias, setpoint, "mV", "pA");
		
		File f = null;
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			f = fc.getSelectedFile();
		try{
			firstStream = new PrintStream(f);
		}
		catch (FileNotFoundException e)
		{
			System.out.println("No new file.");
		}
		firstStream.append(header);
		for (int i = 0; i < M; i++){
			for (int j = 0; j < N-1; j++)
				firstStream.append(data[j][i] + "\t");
			firstStream.append(data[N-1][i]+"\r\n");
		}
		firstStream.append(tail);
		firstStream.close();
	}
	public static void writeITX3D(double[][][] data, String dir, String name, double topolength)
	{
		writeITX3D(data, dir, name, topolength, topolength);
	}
	public static void writeITX3D(double[][][] data, String dir, String name, double topolengthX, double topolengthY)
	{
		int N = data[0].length, M = data[0][0].length, ntps = data.length;
//		String ans = getITX(data, topolength, topolength, name) + "\r\n";
		
		String header = "";
		header += "IGOR" + "\r\n";
		header += "WAVES/D/N=(" + N + "," + M + "," + ntps + ") w" + name + "\r\n";
		header += "BEGIN\r\n";
//		ans += getIgorNote(name, bias, setpoint, "mV", "pA");
//		ColumnIO.writeString(ans, dir + name + ".itx");
		PrintStream firstStream = null;
		String tail = "";
		tail += "END\r\n";
		tail += "X SetScale/P x 0," + (topolengthX/N) + ",\"\", w"+name + "; SetScale/P y 0," + (topolengthY/M) + ",\"\", w"+name + "; Setscale/P z 0,1,\"\", w" + name + "; SetScale/P d 0,0,\"\", w" + name; 
		try{
			firstStream = new PrintStream(new File(dir + name + ".itx"));
		}
		catch (FileNotFoundException e)
		{
			System.out.println("No new file.");
		}
		firstStream.append(header);
		for (int k = 0; k < ntps; k++)
			for (int i = 0; i < M; i++){
				for (int j = 0; j < N-1; j++)
					firstStream.append(data[k][i][j] + "\t");
				firstStream.append(data[k][i][N-1]+"\r\n");
			}
		firstStream.append(tail);
		firstStream.close();
		
	}
	public static void writeDatDirToITX3(String dir, String keyword, double topolength, String finalname)
	{
		File[] files = new File(dir).listFiles();
		double[][] data = null;
		double[][][] alltopos = null;
		int N;
		String name = "";
//		String filetable = "File Table.txt";
//		String content = "Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		int filecount = 0;
		for (int i = 0; i < files.length; i++)
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
				filecount++;
		
		int f = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
			{
				if (data == null) {
					data = ColumnIO.readSquareTable(files[i].toString());
					N = data.length; 
					alltopos = new double[filecount][N][N];
					FieldOps.copy(data, alltopos[f]);
				}
				else
					alltopos[f] = ColumnIO.readSquareTable(files[i].toString());
				name = files[i].toString().substring(files[i].toString().lastIndexOf("\\") + 1, files[i].toString().length() - 4);
				System.out.println(name);
				f++;
			}
		}
		writeITX3D(alltopos, dir, finalname, topolength);
//		ColumnIO.writeString(content, dir + filetable);
	}
	public static void writeDatDirToITX3(String dir, String keyword, double topolengthX, double topolengthY, String finalname)
	{
		File[] files = new File(dir).listFiles();
		double[][] data = null;
		double[][][] alltopos = null;
		int N;
		String name = "";
//		String filetable = "File Table.txt";
//		String content = "Name\tx Pix\ty Pix\tx length(nm)\ty length(nm)\tBias\tSetpoint" + "\r\n";
		int filecount = 0;
		for (int i = 0; i < files.length; i++)
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
				filecount++;
		
		int f = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
			{
				if (data == null) {
					data = ColumnIO.readSquareTable(files[i].toString());
					N = data.length; 
					alltopos = new double[filecount][N][N];
					FieldOps.copy(data, alltopos[f]);
				}
				else
					alltopos[f] = ColumnIO.readSquareTable(files[i].toString());
				name = files[i].toString().substring(files[i].toString().lastIndexOf("\\") + 1, files[i].toString().length() - 4);
				System.out.println(name);
				f++;
			}
		}
		writeITX3D(alltopos, dir, finalname, topolengthX, topolengthY);
//		ColumnIO.writeString(content, dir + filetable);
	}
	
	public static Topomap getTopoFromTxtList(String dir, String suffix, boolean useNormalV)
	{
		File[] files = new File(dir).listFiles();
		int o = Integer.parseInt(JOptionPane.showInputDialog(null,
				"Select a line correction method (enter the number):\r\n"+
						"0) Do nothing.\r\n" +
						"1) Subract average from horizontal lines.\r\n" +
						"2) Subtract average from vertical lines.\r\n" +
						"3) Subtract best-fit line from horizontal lines. \r\n" +
						"4) Subtract best-fit line from vertical lines. \r\n" +
						"5) Do 3 and 4, (3 first).\r\n" +
						"6) Fourier-filter the horizontal and vertical lines (like 5 but less drastic in FFT)\r\n" +
						"7) Fourier-filter the horizontal and vertical lines, but leave the overall average."
						, "Scan line Correction", JOptionPane.OK_CANCEL_OPTION));
		Layer l;
		ArrayList<Layer> ls = new ArrayList <Layer>();
		String name;
		for (int i = 0; i < files.length; i++)
		{
			name = files[i].toString().substring(0, files[i].toString().length() - 4);
			if (files[i].toString().endsWith(".txt") && name.endsWith(suffix))
			{
				l = RHKFileOps.getLayer(files[i]);
				switch (o)
				{
				case 0: break;
				case 1:
					FieldOps.subtactLineAvg(l.data, false); break;
				case 2:
					FieldOps.subtactLineAvg(l.data, true); break;
				case 3:
					FieldOps.subtractLineFit(l.data, true); break;
				case 4:
					FieldOps.subtractLineFit(l.data, false); break;
				case 5:
					FieldOps.subtractLineFit(l.data, true); 
					FieldOps.subtractLineFit(l.data, false); break;
				case 6:
					LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldCross(l.nx, l.ny, true)); break;
				case 7:
					LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldCrossNotO(l.nx, l.ny, true)); break;
				case 8:
					LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldLine(l.nx, l.ny, true, true)); break;
				case 9:
					LayerUtil.fourierFilter(l, TopomapUtil.FourierFilterMethods.getSupressionFieldLine(l.nx, l.ny, false, true)); break;
				}
				
				if (!useNormalV)
					l.v = Double.parseDouble(name.substring(name.indexOf(suffix)-4, name.indexOf(suffix)));

				ls.add(l);
				System.out.println(name);
			}
		}
		Layer[] layers = new Layer[ls.size()];
		for (int i = 0; i < layers.length; i++)
			layers[i] = ls.get(i);
		return Topomap.newTopomap(layers);
	}
	public static void killFirstEntriesAndHeader()
	{
		JFileChooser fc = new JFileChooser();
		File f;
		double[][] data;
		JOptionPane.showMessageDialog(null, "Select the first RHK file. To stop, click cancel on the file chooser.");
		while(fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			f = fc.getSelectedFile();
			data = getTable(f);
			ColumnIO.writeTable(data, f.toString());
		}
	}
	
	public static void addOrigin(Layer t)
	{
		double xc = Double.parseDouble(JOptionPane.showInputDialog("Enter the x-coordinate of the center in nanometers (as in RHK tooltip)."));
		double yc = Double.parseDouble(JOptionPane.showInputDialog("Enter the y-coordinate of the center in nanometers (as in RHK tooltip)."));
		t.makeOrigin(new double[] {xc, yc});
	}
	public static double[] getOriginFromUser()
	{
		double xc = Double.parseDouble(JOptionPane.showInputDialog("Enter the x-coordinate of the center in nanometers (as in RHK tooltip)."));
		double yc = Double.parseDouble(JOptionPane.showInputDialog("Enter the y-coordinate of the center in nanometers (as in RHK tooltip)."));
		return new double[] {xc, yc};
	}
	public static double[][] getFromGreyBMP()
	{
		int[][][] pixels = imageIO.loadImage();
		double[][] data = new double [pixels.length][pixels[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				data[j][i] = (ArrayOps.sum(pixels[i][j])-256)/765.0;
		return data;
	}

//	public static void readSM4File(String path)
//	{
//		File file = new File(path);
//		FileInputStream inf = null;
//		BufferedInputStream inbuff = null;
//		DataInputStream ind = null;
//		byte[] data = new byte[(int)file.length()];
//		try {
//			inf = new FileInputStream(file);
//			inbuff = new BufferedInputStream(inf, (int)file.length());
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		ind = new DataInputStream(inbuff);
//
//		try{
//			short s = ind.readShort();
//			String name = ind.
////			ind.read(data);
////			for (int i = 0; i < 56; i++)
////				System.out.println(data[i]);
//			byte[] sig = new byte[36];
//			for (int i = 0; i < 36; i++)
//				sig[i] = data[i+2];
//			String s = new String(sig, "UTF-16LE");
//			System.out.println(s);
//			
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//
//	}
}
