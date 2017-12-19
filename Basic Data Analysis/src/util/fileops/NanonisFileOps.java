package util.fileops;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import drawing.GraphDrawerCart;
import drawing.LayerViewer;
import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.Printer;
import util.TopomapUtil;
import util.fourier.FFTOps;

public class NanonisFileOps {
	
	public static String nanonisDir = "H:\\Work computer\\data\\raw\\nanonis unisoku\\Run #467 Mn-doped PbSnSe 0.3% T-0037\\Area 2\\";
//	public static boolean changeNaN = false;
	public static int nanThreshold = 200*200; //number of square pixels. Above this value we do not check for NaN as it would take too long.
	
	public static String[] defaultNames = {"Z", "ZBack", "LockX", "LockXback", "Current", "Currentback"};
	
	public static void main(String[] args)
	{
		int choice = JOptionPane.showConfirmDialog(null, "To open a scan, click Yes.\r\n" +
				"To open grid spectroscopy, click No.\r\n" +
				"To do something with the time signal, click cancel.");
		if (choice == JOptionPane.YES_OPTION)
			convertScanToBins(false);
		else if (choice == JOptionPane.NO_OPTION)
			convert3dsToBins();
		else if (choice == JOptionPane.CANCEL_OPTION)
		{
			String path = FileOps.selectOpen(new JFileChooser(nanonisDir)).toString();
			doTimeSignalStuff(path);
		}
//		System.exit(0);
//		do3dsConversions();
//		convertScanToBins(false);
//		convertScanToBins("C:\\data\\raw\\nanonis\\Run #18 Bi4Se3\\scan015.sxm", "C:\\data\\analysis\\Bi4X3\\Se\\big topo\\", true);
		
//		loadLayerFromScan("C:\\data\\raw\\nanonis\\Run #18 Bi4Se3\\scan015.sxm");
//		String dir = "C:\\data\\Run #18 Bi4Se3\\nice topos and point spectra\\";
//		Topomap[] t = loadTopomapFrom3ds("C:\\data\\Run #18 Bi4Se3\\nice topos and point spectra\\Grid Spectroscopy001.3ds");
//		String[] names = getLayerNamesFrom3ds("C:\\data\\Run #18 Bi4Se3\\nice topos and point spectra\\Grid Spectroscopy001.3ds");
//
//		for (int i = 0; i < names.length; i++)
//		{
//			if (names[i].contains("\\")) names[i] = names[i].split("\\")[0];
//			if (names[i].contains("/")) names[i] = names[i].split("/")[0];
//			Topomap.writeBIN(t[i], dir + "A_" + names[i] + ".bin");
//		}
	}
	
	public static void doTimeSignalStuff(String path)
	{
		double[][] filter = {/*{14.7, 14.85},*/ 
//				{8.955, 8.985}, 
				{2.13, 2.15},
//				{5.986, 5.988},
				{59.95, 60.05},
				{119.9, 121.1},
				{179.9, 181.1},
				{239.85, 240.1},
				{299.9, 300.1}};
//				};
		int o = Integer.parseInt(JOptionPane.showInputDialog("Enter the code:\r\n"+
				"0 - Convert and filter the layer\r\n"+
				"1 - Look at FFT\r\n"+
				"2 - Look at real space signal\r\n" +
				"3 - Get the amplitude of a frequency as a function of time\r\n"+
				"4 - Get the Fourier transform as a function of time at a certain frequency\r\n"+
				"5 - Look at the time derivative of a real-space signal\r\n"));
		TimeSignal ts;
		int layerCode; 		
		
		
		switch(o)
		{
		case 0:
			convertScanToBinsTimeFilter(false, getFilterRangesFromUser(path), path);
			break;
		case 1:{
			layerCode = getLayerCodeFromUser();
			boolean scanDown = NanonisFileOps.scanWasDown(path);
			ts = getTimeSignal(path, layerCode, scanDown);
			GraphDrawerCart.plotGraph(ts.frequencies, ts.fftmag);
			break;}
		case 2:{
			layerCode =  getLayerCodeFromUser();
			boolean scanDown = NanonisFileOps.scanWasDown(path);
			ts = getTimeSignal(path, layerCode, scanDown);
			GraphDrawerCart.plotGraph(ts.time, ts.values);
			break;}
		case 3:
			layerCode =  getLayerCodeFromUser();
			ts = getTimeSignal(path, layerCode);
//			ts = TimeSignal.getSineWave(2000, 15, 512*512);
			double freq = Double.parseDouble(JOptionPane.showInputDialog("Enter the frequency in Hertz."));
			double length = Double.parseDouble(JOptionPane.showInputDialog("Enter decay length in seconds."));
			double[] result = new double [ts.values.length];
			double[] copy = ArrayOps.subtractPolynomialFit(ts.time, ts.values, 2)[0];
			double[] cosineTimes = new double [result.length];
			double k = ts.getNatrualFrequency(freq);
			for (int i = 0; i < cosineTimes.length; i++)
				cosineTimes[i] = copy[i]*Math.cos(k*i);
			double[] gaussDecay = ArrayOps.getGaussDecayingSine(0, length/ts.pixelTime, 4);
//			GraphDrawerCart.plotGraph(ArrayOps.generateArrayNotInclUpper(0, gaussDecay.length, gaussDecay.length), gaussDecay);
			for (int i = 0; i < result.length; i++)
				result[i] = ArrayOps.getProductWithEvenArray(cosineTimes, gaussDecay, i);
			GraphDrawerCart.plotGraph(ts.time, result);
			GraphDrawerCart.plotGraph(ts.frequencies, FFTOps.get1DFFTMag(result));
			break;
		case 4:
			layerCode =  getLayerCodeFromUser();
			ts = getTimeSignal(path, layerCode);
			double minFreq = Double.parseDouble(JOptionPane.showInputDialog("Enter the minimum of the band frequency in Hertz."));
			double maxFreq = Double.parseDouble(JOptionPane.showInputDialog("Enter the minimum of the band frequency in Hertz."));
			int minI = ArrayOps.indexOf(ts.frequencies, minFreq, true);
			int maxI = ArrayOps.indexOf(ts.frequencies, maxFreq, true);
			int nf  = maxI - minI;
			JOptionPane.showMessageDialog(null, "There will be " + nf + " frequency pixels. (" + minI + " to " + maxI + " out of " + ts.values.length + ")");
			length = Double.parseDouble(JOptionPane.showInputDialog("Enter decay length in seconds."));
			gaussDecay = ArrayOps.getGaussDecayingSine(0, length/ts.pixelTime, 4);
			copy = ArrayOps.subtractPolynomialFit(ts.time, ts.values, 2)[0];
			
			int nt = Integer.parseInt(JOptionPane.showInputDialog("The map contains " + (int)(ts.totalTime/length) + " decay lengths and (probably) " + (int)Math.sqrt(ts.values.length/2) + " lines.\r\n"
					+ "How many time points do you want to take?"));
			
			int[] is = new int [nt];
			for (int i = 0; i < nt; i++)
				is[i] = FieldOps.round(((double)(ts.values.length-1)/nt)*i);
			
			double[][] ans = new double [nf][nt];
			double[] freqUsed = new double [nf];
			for (int j = 0; j < nf; j++)
				freqUsed[j] = ts.frequencies[minI+j];
			double[] timeUsed = ArrayOps.generateArrayInclBoth(0, ts.totalTime, nt);
			double[] tempFFT;
			for (int i = 0; i < nt; i++)
			{
				result = ArrayOps.getFilteredWithEvenArray(copy, gaussDecay, is[i]);
				if (i == nt/2){
					GraphDrawerCart.plotGraph(ArrayOps.generateArrayNotInclUpper(0, ts.values.length, ts.values.length), result);
					GraphDrawerCart.plotGraph(null, gaussDecay);
				}
				tempFFT = FFTOps.get1DFFTMag(result);
				System.out.print(i + " ");
				for (int j = 0; j < nf; j++)
					ans[j][i] = tempFFT[j+minI];
			}
			System.out.println();
			Layer layer = new Layer(ans, freqUsed, timeUsed, 1, 1);
			LayerViewer lv = new LayerViewer(layer, Topomap.stddir, 512);
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			break;
		case 5:
			layerCode =  getLayerCodeFromUser();
			double[] tempArray = getTimeSignalValues(path, layerCode);
			double[] derivative = ArrayOps.getDerivative(tempArray);
			
			ts = getTimeSignalFromValues(path, layerCode, derivative);
			GraphDrawerCart.plotGraph(ts.time, ts.values);
			break;
		}
		
//		if (o == JOptionPane.YES_OPTION){
////			System.exit(0);
//		}
//		else if (o == JOptionPane.NO_OPTION){
////			JFileChooser fci = new JFileChooser(nanonisDir);
//			TimeSignal ts = getTimeSignal(path, 0);
//			GraphDrawerCart.plotGraph(ts.frequencies, ts.fftmag);
//		}
//		else if (o == JOptionPane.CANCEL_OPTION){
////			JFileChooser fci = new JFileChooser(nanonisDir);
//			TimeSignal ts = getTimeSignal(path, 0);
//			GraphDrawerCart.plotGraph(ts.time, ts.values);
//		}

	}
	public static int getLayerCodeFromUser()
	{
		int layerCode = Integer.parseInt(JOptionPane.showInputDialog("Enter the code:\r\n"+
				"0 - Z\r\n"+
				"1 - Current\r\n"+
				"2 - Lockin\r\n"));
		if (layerCode == 0) return TimeSignal.Z_CODE;
		else if (layerCode == 1 ) return TimeSignal.CURRENT_CODE;
		else if (layerCode == 2) return TimeSignal.LOCKIN_CODE;
		return 0;

	}
	//This extracts the "maps" of each channel as well as the topography, which is in the final "map."
	//It also saves them to a user-chosen file.
	
	public static void do3dsConversions()
	{
		String base = "G:\\Run #18 Bi4Se3\\nice topos and point spectra\\";
		String[] names = new String [2];
		for (int i = 1; i <= names.length; i++)
			names[i-1] = base + "Grid Spectroscopy" + "00" + i + ".3ds";
		
		String[] dests = new String[names.length];
		dests[0] = "C:\\data\\analysis\\Bi4X3\\Se\\map 1 -300 to 100\\";
		dests[1] = "C:\\data\\analysis\\Bi4X3\\Se\\map 2 -100 to 400\\";
		for (int i = 0; i < names.length; i++)
			convert3dsToBins(names[i], dests[i]);
	}
	public static void convert3dsToBins()
	{
		Topomap.setStdDir();
		JFileChooser fci = new JFileChooser(nanonisDir);
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		String path = FileOps.selectOpen(fci).toString();
		String dest = FileOps.selectSave(fco).toString();
		boolean cropIt = JOptionPane.showConfirmDialog(null, "Crop any pixels?") == JOptionPane.YES_OPTION;
		if (!cropIt)
			convert3dsToBins(path, dest);
		else
			convert3dsToBinsCropped(path, dest);
	}	
	/**
	 * This assumes that the directory contains NOTHING except the files in the topomap.
	 */
	public static void convertDirToTopomap(String dir)
	{
		File[] f  = new File(dir).listFiles();
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		String dest = FileOps.selectSave(fco).toString();
		boolean useFilter = JOptionPane.showConfirmDialog(null, "Attempt to fourier-filter the layers?") == JOptionPane.YES_OPTION;
		if (!useFilter) convertScansToTopomaps(f, dest);
		else convertScansToTopomapsWithTimeFilter(f, dest);
	}
	public static void convertDirToTopomap2(String dir)
	{
		File[] f  = new File(dir).listFiles();
		FileWithVoltage[] fv = new FileWithVoltage[f.length];
		for (int i = 0; i < f.length; i++)
			fv[i] = new FileWithVoltage(f[i]);
		ArrayOps.quicksort(fv);
		f = FileWithVoltage.getFiles(fv);
		
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		String dest = FileOps.selectSave(fco).toString();
		boolean useFilter = JOptionPane.showConfirmDialog(null, "Attempt to fourier-filter the layers?") == JOptionPane.YES_OPTION;
		if (!useFilter) convertScansToTopomaps(f, dest);
		else convertScansToTopomapsWithTimeFilter(f, dest);
	}
	public static class FileWithVoltage implements Comparable
	{
		public File f;
		double v;
		
		public FileWithVoltage(File f)
		{
			this.f = f;
			v = NanonisFileOps.getVoltageFromScan(f.toString());
		}
		@Override
		public int compareTo(Object arg0) {
			// TODO Auto-generated method stub
			if (((FileWithVoltage)arg0).v < v) return -1;
			if (((FileWithVoltage)arg0).v > v) return 1;
			return 0;
		}
		public static File[] getFiles(FileWithVoltage[] fv)
		{
			File[] f = new File[fv.length];
			for (int i = 0; i < f.length; i++)
				f[i] = fv[i].f;
			return f;
		}
	}
	public static void convert3dsToBins(String sourceFilePath, String dest)
	{
		Topomap[] t = loadTopomapFrom3ds(sourceFilePath);
		String[] names = getLayerNamesFrom3ds(sourceFilePath);
		for (int i = 0; i < names.length; i++)
		{
//			names[i].replaceAll("/", "S");
//			names[i].replaceAll("\\", "");
			if (names[i].contains("\\")) names[i] = names[i].split("\\")[0] + names[i].split("\\")[1];
			if (names[i].contains("/")) names[i] = names[i].split("/")[0] + names[i].split("/")[1];
			Topomap.writeBIN(t[i], dest + names[i] + ".bin");
		}
		Layer.writeBIN(t[names.length].getLayer(0), dest + "topo.bin");
		ColumnIO.writeString(sourceFilePath, dest + "_source file.txt");
	}
	public static void convert3dsToBinsCropped(String sourceFilePath, String dest)
	{
		Topomap[] t = loadTopomapFrom3ds(sourceFilePath);
		int xi; int xf; int yi; int yf;
		String input = JOptionPane.showInputDialog("Enter the minimum and maximum of x, comma separated.", "" + 0 + "," + t[0].nx);
		String[] token = input.split(",");
		xi = Integer.parseInt(token[0]); xf = Integer.parseInt(token[1]);
		input = JOptionPane.showInputDialog("Enter the minimum and maximum of y, comma separated.", "" + 0 + "," + t[0].ny);
		token  = input.split(",");
		yi = Integer.parseInt(token[0]); yf = Integer.parseInt(token[1]);
		
		String[] names = getLayerNamesFrom3ds(sourceFilePath);
		for (int i = 0; i < names.length; i++)
		{
//			names[i].replaceAll("/", "S");
//			names[i].replaceAll("\\", "");
			if (names[i].contains("\\")) names[i] = names[i].split("\\")[0] + names[i].split("\\")[1];
			if (names[i].contains("/")) names[i] = names[i].split("/")[0] + names[i].split("/")[1];
			Topomap.writeBIN(Topomap.getCropped(t[i], xi, xf, yi, yf), dest + names[i] + ".bin");
		}
		Layer.writeBIN(Layer.getCropped(t[names.length].getLayer(0), xi, xf, yi, yf), dest + "topo.bin");
		ColumnIO.writeString(sourceFilePath
				+ "\r\n" + xi + "\t" + xf + "\r\n" + yi + "\t" + yf, dest + "_source file.txt");
	}
	public static void convert3dsToBinsLineCut_PointSpec(String sourceFilePath, String dest)
	{
		Topomap[] t = loadTopomapFrom3ds(sourceFilePath);
		int xi; int xf; int yi; int yf;
		String input = JOptionPane.showInputDialog("Enter the minimum and maximum of x, comma separated.", "" + 0 + "," + t[0].nx);
		String[] token = input.split(",");
		xi = Integer.parseInt(token[0]); xf = Integer.parseInt(token[1]);
		input = JOptionPane.showInputDialog("Enter the minimum and maximum of y, comma separated.", "" + 0 + "," + t[0].ny);
		token  = input.split(",");
		yi = Integer.parseInt(token[0]); yf = Integer.parseInt(token[1]);
		
		String[] names = getLayerNamesFrom3ds(sourceFilePath);
		ArrayList<Topomap> wanted = new ArrayList<Topomap>();
		
		for (int i = 0; i < names.length; i++)
		{
//			names[i].replaceAll("/", "S");
//			names[i].replaceAll("\\", "");
			t[i] = Topomap.getCropped(t[i], xi, xf, yi, yf);
			if (names[i].contains("\\")) names[i] = names[i].split("\\")[0] + names[i].split("\\")[1];
			if (names[i].contains("/")) names[i] = names[i].split("/")[0] + names[i].split("/")[1];
			//This condition signifies that the data is dI/dV, but not the average:
			if (names[i].contains("Lockin") && names[i].contains("X") && !names[i].contains("[AVG]"))
			{	
				wanted.add(t[i]);
				System.out.println(i + "\t" + names[i] + "\t" + t[i]);
			}
		}
		
		int nPerPoint = wanted.size();
		int nPoints = t[0].nx*t[0].ny;
		int nspec = nPerPoint*nPoints;
		double [][] data = new double [nspec][t[0].nlayers];
		double[] v = t[0].v;
		double[] x = new double [nspec], y = new double [nspec];
		int p = 0;
		for (int i = 0; i < t[0].nx; i++)
			for (int j = 0; j < t[0].ny; j++)
			{
				for (int k = 0; k < wanted.size(); k++){
					data[p*nPerPoint+k] = wanted.get(k).getSpectrum(i, j);
					x[p] = wanted.get(k).x[i];
					y[p] = wanted.get(k).y[j];
				}
				p++;
			}
		PointSpectra.writeBIN(new PointSpectra(data, v, x, y), dest + ".bin");
//		Layer.writeBIN(Layer.getCropped(t[names.length].getLayer(0), xi, xf, yi, yf), dest + "topo.bin");
		ColumnIO.writeString(sourceFilePath
				+ "\r\n" + xi + "\t" + xf + "\r\n" + yi + "\t" + yf, dest + "_source file.txt");
		
	}
	public static void convertScansToTopomaps(File[] sources, String dest)
	{
		Topomap t;
		String [] lines = new String [sources.length];
		
		for (int i = 0; i < defaultNames.length; i++)
		{
			t = getMapFromScanSeries(sources, i);
			if (t.v[0] > t.v[t.nlayers-1])
				TopomapUtil.flipEnergy(t);
			Topomap.writeBIN(t, dest + defaultNames[i] + ".bin");
		}
		for (int i = 0; i < sources.length; i++)
			lines[i] = sources[i].toString();
		
		ColumnIO.writeLines(lines, dest + "_source file.txt");
	}
	public static void convertScansToTopomapsWithTimeFilter(File[] sources, String dest)
	{
		Topomap[] all = new Topomap[6];
		double[][][][] data = new double[6][sources.length][][];
		String[] lines = new String [sources.length];
		String[] names = {"Z", "ZBack", "Curr", "CurrBack", "X", "XBack"};
//		double[][][] filterRanges = getFilterRanges3FromUser(sources[0].toString());
		double[][] filterRangebasic = getFilterRangesFromUser(sources[0].toString());
		double[][][] filterRanges = new double [][][] {filterRangebasic,filterRangebasic,filterRangebasic};
		Layer[][] set = null;
		double[] v = new double [sources.length];
		for (int i = 0; i < sources.length; i++)
		{
			set = getWithTimeFilter(sources[i].toString(), filterRanges);
			data[0][i] = set[0][0].data;
			data[1][i] = set[0][1].data;
			data[2][i] = set[1][0].data;
			data[3][i] = set[1][1].data;
			data[4][i] = set[2][0].data;
			data[5][i] = set[2][1].data;
			v[i] = set[0][0].v;
		}
		for (int i = 0; i < 6; i++)
		{
			all[i] = new Topomap(data[i], v.clone(), set[0][0].x, set[0][0].y, null);
			if (all[i].v[0] > all[i].v[all[i].nlayers-1])
				TopomapUtil.flipEnergy(all[i]);
			Topomap.writeBIN(all[i], dest + names[i] + ".bin");
		}
		for (int i = 0; i < sources.length; i++)
			lines[i] = sources[i].toString();
		
		ColumnIO.writeLines(lines, dest + "_source file.txt");
	}
	public static Topomap[] loadTopomapFrom3ds(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Topomap[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/4);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (!s.contains("HEADER_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//read the header
		//First line: nx, ny.
		int nx, ny;
		String thing = headerLines.get(0).substring(headerLines.get(0).indexOf("=")+2, headerLines.get(0).length() - 1);
		String[] dim = thing.split(" x ");
		nx = Integer.parseInt(dim[0]);
		ny = Integer.parseInt(dim[1]);
		
		String signal = headerLines.get(3).substring(headerLines.get(3).indexOf("=")+1, headerLines.get(3).length());
		
		//0 if bias, 1 if z
		int signal_type = 0;
		if (signal.contains("Z (m)"))
			signal_type = 1;
		
		//second line: grid parameters: (this assumes an angle of zero degrees).
		thing = headerLines.get(1).substring(headerLines.get(1).indexOf("=")+1, headerLines.get(1).length());
		String[] grid = thing.split(";");
		double cx = Double.parseDouble(grid[0]);
		double cy = Double.parseDouble(grid[1]);
		double lx = Double.parseDouble(grid[2]);
		double ly = Double.parseDouble(grid[3]);
		double angle = Double.parseDouble(grid[4]);
		
		double[] x = ArrayOps.generateArrayInclBoth(cx-lx/2, cx+lx/2, nx);
		double[] y = ArrayOps.generateArrayInclBoth(cy-ly/2, cy+ly/2, ny);
		
		//Extract the voltage.
		thing = headerLines.get(8).substring(headerLines.get(8).indexOf("=")+1, headerLines.get(8).length());
		int nlayers = Integer.parseInt(thing);
		double[] v = null;
		double[][] topo = new double [nx][ny]; 
		//the sweep parameters may be in the header, if it was included.
		double vstart = 0, vstop = 0;
		if (signal_type == 0){
			for (int i = 0; i < headerLines.size(); i++)
			{
				if (headerLines.get(i).contains("Bias Spectroscopy>Sweep Start (V)"))
				{
					vstart = Double.parseDouble(headerLines.get(i).substring(headerLines.get(i).indexOf("=")+1));
				}
				if (headerLines.get(i).contains("Bias Spectroscopy>Sweep End (V)"))
				{
					vstop = Double.parseDouble(headerLines.get(i).substring(headerLines.get(i).indexOf("=")+1));
				}
			}
			if (vstart != 0 || vstop != 0){
				v = ArrayOps.generateArrayInclBoth(vstart, vstop, nlayers);
			}
		}
		if (signal_type == 1)
		{
			for (int i = 0; i < headerLines.size(); i++)
			{
				if (headerLines.get(i).contains("Z Spectroscopy>Initial Z-offset (m)"))
				{
					vstart = Double.parseDouble(headerLines.get(i).substring(headerLines.get(i).indexOf("=")+1));
				}
				if (headerLines.get(i).contains("Z Spectroscopy>Sweep distance (m)"))
				{
					vstop = Double.parseDouble(headerLines.get(i).substring(headerLines.get(i).indexOf("=")+1));
				}
			}
			if (vstart != 0 || vstop != 0){
				v = ArrayOps.generateArrayInclBoth(vstart, vstop, nlayers);
			}
		}
		//save the number of parameter bytes before each spectrum:
		thing = headerLines.get(6).substring(headerLines.get(6).indexOf("=")+1, headerLines.get(6).length());
		int nPtParam = Integer.parseInt(thing);
		
		//Print the names of the channels:
		thing = headerLines.get(9).substring(headerLines.get(9).indexOf("=")+2, headerLines.get(9).length()-1);
		int nchannels = thing.split(";").length;
		
		
		//If the voltage array is not in the header, extract it from the first data point:
		if (v == null)
		{
			try {
				vstart = (double)ind.readFloat();
				vstop = (double)ind.readFloat();
				v = ArrayOps.generateArrayInclBoth(vstart, vstop, nlayers);
				System.out.println(vstart + "\t" + vstop);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else //clear those two floats anyway
		{
			try {
				ind.readFloat();
				ind.readFloat();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		double[][][][] dataset = new double [nchannels][nlayers][nx][ny];
		int count = 0;
		
//		while (count < 100)
//		{
//			count++;
//			try {
//				System.out.println(ind.readFloat());
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		boolean fileNotFull = false;
		float[] flt = new float [nPtParam-2];
		try {
			for (int n = 0; n < ny; n++)
				for (int m = 0; m < nx; m++)
				{
					if (count != 0)
					{
						ind.readFloat();
						ind.readFloat();
					}
					//discard the parameter values except z and load only the data
					for (int i = 0; i < nPtParam-2; i++)
						flt[i] = ind.readFloat();
					topo[m][n] = flt[2];
					for (int i = 0; i < nchannels; i++)
						for (int j = 0; j < nlayers; j++)
						{
							dataset[i][j][m][n] = (double)ind.readFloat();
						}
					count++;
				}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.out.println("Exception with " + filepath);
			e.printStackTrace();
			//Hopefully this is the end of file error:
			fileNotFull = true;
		}

		for (int i = 0; i < nchannels; i++)
			for (int j = 0; j < nlayers; j++)
			{
//				if (nx == ny)
//					dataset[i][j] = FieldOps.transpose(dataset[i][j]);
				if (angle != 180) ArrayOps.flipY(dataset[i][j]);
				
				if (fileNotFull)
					FieldOps.changeZeroToAverage(dataset[i][j]);
				
				if (angle == 90 && nx == ny){
					FieldOps.rotatePlus90(dataset[i][j]);
					FieldOps.putTranspose(dataset[i][j]);
				}
			}
//		if (nx == ny)
//			topo = FieldOps.transpose(topo);
		if (angle == 90 && nx == ny)
		{
			FieldOps.rotatePlus90(topo);
			FieldOps.putTranspose(topo);
		}
//		if (angle != 180){
//			ArrayOps.flipY(topo);
//			ArrayOps.flip(y);
//		}
		if (fileNotFull) FieldOps.changeZeroToAverage(topo);
		//if the voltage is flipped, flip it
		if (v[v.length-1] < v[0])
		{
			double[] tv = new double [v.length];
			for (int i = 0; i < v.length; i++)
				tv[i] = v[v.length-1-i];
			v = tv;
			double[][][] layers = new double[v.length][][];
			for (int j = 0; j < nchannels; j++){
				for (int i = 0; i < v.length; i++)
				{
					layers[i] = dataset[j][v.length-1-i];
				}
				dataset[j] = layers.clone();
			}
		}

		
		//now, make the topomaps.
		t = new Topomap[nchannels+1];
		for (int i = 0; i < nchannels; i++){
			t[i] = new Topomap(dataset[i], v, x, y, null);
			t[i].fileIsComplete = !fileNotFull;
		}
		
		
		t[nchannels] = new Topomap(new double[][][] {topo}, new double[] {1}, x, y, null);
		t[nchannels].fileIsComplete = !fileNotFull;

		//check for NaN
//		for (int i = 0; i < t.length; i++)
//			for (int j = 0; j < t[i].nlayers; j++)
//			{
//				FieldOps.changeNaNToAverage(t[i].data[j]);
//			}

		
		return t;
	}
	public static String[] getLayerNamesFrom3ds(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Topomap t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/16);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("HEADER_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//Print the names of the channels:
		String thing = headerLines.get(9).substring(headerLines.get(9).indexOf("=")+2, headerLines.get(9).length()-1);
		return thing.split(";");
	}
	public static Layer[] loadLayerFromScan(String filepath, boolean evaluateNan)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/4);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("SCANIT_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		String[] names = null;
		int nchannels = 0;
		int nx=0, ny=0;
		double lx=0, ly=0, cx=0, cy=0;
		double bias=0;
		double current=0;
		String thing;
		String[] tokens;
		String direction = null;
		double angle = 0;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("SCAN_PIXELS"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				nx = Integer.parseInt(tokens[0]);
				ny = Integer.parseInt(tokens[tokens.length-1]);
			}
			if (headerLines.get(i).contains("SCAN_RANGE"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				lx = Double.parseDouble(tokens[0]);
				ly = Double.parseDouble(tokens[tokens.length-1]);
			}
			if (headerLines.get(i).contains("SCAN_OFFSET"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				cx = Double.parseDouble(tokens[0]);
				cy = Double.parseDouble(tokens[tokens.length-1]);
			}
			if (headerLines.get(i).contains("SCAN_ANGLE"))
			{
				thing = headerLines.get(i+1).trim();
				angle = Double.parseDouble(thing);
			}
			if (headerLines.get(i).contains("SCAN_DIR"))
			{
				direction = headerLines.get(i+1).trim();
			}
			if (headerLines.get(i).contains(":BIAS:"))
			{
				thing = headerLines.get(i+1).trim();
				bias = Double.parseDouble(thing);
			}
			if (headerLines.get(i).contains("Current>Current (A)"))
			{
				thing = headerLines.get(i+1).trim();
				current = Double.parseDouble(thing);
			}
			if (headerLines.get(i).contains("Scan>channels"))
			{
				thing = headerLines.get(i+1).trim();
				nchannels = thing.split(";").length;
				names = new String [nchannels*2];
			}
			if (headerLines.get(i).contains("DATA_INFO"))
			{
				for (int j = i+2; j < i+2+nchannels; j++)
				{
					tokens = headerLines.get(j).trim().split("\t");
					names[j-(i+2)] = tokens[1];
					names[j-(i+2) + nchannels] = tokens[1] + "back";
				}
			}
		}
		
		//if we forgot to save the Scan parameters (we should NEVER do this)
		//The nchannels will be the number of lines after DATA_INFO
		if (nchannels == 0)
		{
			
		}
		double[] x = ArrayOps.generateArrayInclBoth(cx-lx/2, cx+lx/2, nx);
		double[] y = ArrayOps.generateArrayInclBoth(cy+ly/2, cy-ly/2, ny);
		double[][][] dataset = new double [nchannels*2][nx][ny];
		boolean full = true;
		try {
			ind.readFloat();
			for (int k = 0; k < nchannels; k++){
				for (int j = 0; j < ny; j++)
					for (int i = 0; i < nx; i++)
					{
						dataset[k][i][j] = ind.readFloat();
					}
				for (int j = 0; j < ny; j++)
					for (int i = 0; i < nx; i++)
					{
						dataset[k+nchannels][i][j] = ind.readFloat();
					}

			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			full = false;
		}
		
		try {
			inbuff.close();
			ind.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		double[] time = getLineTimeFromScan(filepath);
		//check for NaN;
		String str;
		String nan = "NaN";
		if (full && (nx*ny < nanThreshold || evaluateNan))
//			for (int i = 0; i < nx && full; i++)
//				for (int j = 0; j < ny && full; j++)
//				{
//					str = "" + dataset[0][i][j];
//					if (str.equals(nan))
//						full = false;
//				}
			{//just evaluate the last pixel
				str = "" + dataset[0][nx-1][ny-1];
				if (str.equals(nan))
					full = false;
			}
		//All pixels after nNan are considered as being NaN.
		
//		System.out.println(i + "\t" + j + "\t" + dataset[0][i][j]);
		if (!full && nx*ny < nanThreshold  || evaluateNan)
			for (int i = 0; i < 2*nchannels; i++)
				FieldOps.changeNaNToAverage(dataset[i]);
		
		for (int i = nchannels; i < 2*nchannels; i++)
			ArrayOps.flipX(dataset[i]);
		
		for (int i = 0; i < 2*nchannels; i++)
			if (angle == 90)
				dataset[i] = FieldOps.rotateMinus90(dataset[i]);
		//now, make the topomaps.
		t = new Layer[nchannels*2];
		for (int i = 0; i < nchannels*2; i++)
		{
			t[i] = new Layer(dataset[i], x, y, bias, current);
			t[i].fileIsComplete = full;
			t[i].lineTime = time[0];
			if (direction.contains("up"))
				t[i].flipY();
		}//		t[nchannels] = new Topomap(new double[][][] {topo}, new double[] {1}, x, y, null);
		return t;
	}
	public static boolean scanWasDown(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/4);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("SCANIT_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		String direction = null;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("SCAN_DIR"))
			{
				direction = headerLines.get(i+1).trim();
				break;
			}
		}
		return (!direction.contains("up"));
	}
	/**
	 * This loads only the layer specified by index, to save memory. 
	 * @param filepath
	 * @param index 0=Z, 1 = Zback, etc.
	 * @return
	 */
	public static Layer loadSpecificLayerFromScan(String filepath, int index)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/4);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("SCANIT_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		String[] names = null;
		int nchannels = 0;
		int nx=0, ny=0;
		double lx=0, ly=0, cx=0, cy=0;
		double bias=0;
		double current=0;
		String thing;
		String[] tokens;
		String direction = null;
		double angle = 0;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("SCAN_PIXELS"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				nx = Integer.parseInt(tokens[0]);
				ny = Integer.parseInt(tokens[tokens.length-1]);
			}
			if (headerLines.get(i).contains("SCAN_RANGE"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				lx = Double.parseDouble(tokens[0]);
				ly = Double.parseDouble(tokens[tokens.length-1]);
			}
			if (headerLines.get(i).contains("SCAN_OFFSET"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				cx = Double.parseDouble(tokens[0]);
				cy = Double.parseDouble(tokens[tokens.length-1]);
			}
			if (headerLines.get(i).contains("SCAN_ANGLE"))
			{
				thing = headerLines.get(i+1).trim();
				angle = Double.parseDouble(thing);
			}
			if (headerLines.get(i).contains("SCAN_DIR"))
			{
				direction = headerLines.get(i+1).trim();
			}
			if (headerLines.get(i).contains(":BIAS:"))
			{
				thing = headerLines.get(i+1).trim();
				bias = Double.parseDouble(thing);
			}
			if (headerLines.get(i).contains("Current>Current (A)"))
			{
				thing = headerLines.get(i+1).trim();
				current = Double.parseDouble(thing);
			}
			if (headerLines.get(i).contains("Scan>channels"))
			{
				thing = headerLines.get(i+1).trim();
				nchannels = thing.split(";").length;
				names = new String [nchannels*2];
			}
			if (headerLines.get(i).contains("DATA_INFO"))
			{
				for (int j = i+2; j < i+2+nchannels; j++)
				{
					tokens = headerLines.get(j).trim().split("\t");
					names[j-(i+2)] = tokens[1];
					names[j-(i+2) + nchannels] = tokens[1] + "back";
				}
			}
		}
		
		double[] x = ArrayOps.generateArrayInclBoth(cx-lx/2, cx+lx/2, nx);
		double[] y = ArrayOps.generateArrayInclBoth(cy+ly/2, cy-ly/2, ny);
		double[][] dataset = new double [nx][ny];
		boolean full = true;
		try {
			ind.readFloat();
			for (int k = 0; k < 2*nchannels; k++){
				if (k == index)
					for (int j = 0; j < ny; j++)
						for (int i = 0; i < nx; i++)
							dataset[i][j] = ind.readFloat();
				else
					for (int j = 0; j < ny; j++)
						for (int i = 0; i < nx; i++)
							ind.readFloat();

			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			full = false;
		}
		
		try {
			inbuff.close();
			ind.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		double[] time = getLineTimeFromScan(filepath);
		//check for NaN;
		String str;
		String nan = "NaN";
		if (full && nx*ny < nanThreshold || false)
			for (int i = 0; i < nx && full; i++)
				for (int j = 0; j < ny && full; j++)
				{
					str = "" + dataset[i][j];
					if (str.equals(nan))
						full = false;
				}
		//All pixels after nNan are considered as being NaN.
		
//		System.out.println(i + "\t" + j + "\t" + dataset[0][i][j]);
		if (!full && (nx*ny < nanThreshold || true))
			FieldOps.changeNaNToAverage(dataset);
		
		if (index %2 == 1)
			ArrayOps.flipX(dataset);
		
		if (angle == 90)
			dataset = FieldOps.rotateMinus90(dataset);
		//now, make the topomaps.
		t = new Layer(dataset, x, y, bias, current);
		t.fileIsComplete = full;
		t.lineTime = time[0];
		if (direction.contains("up"))
			t.flipY();
		return t;
	}
	public static String[] loadLayerNamesFromScan(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/64);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("SCANIT_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		String[] names = null;
		int nchannels = 0;
		String thing;
		String[] tokens;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("Scan>channels"))
			{
				thing = headerLines.get(i+1).trim();
				nchannels = thing.split(";").length;
				names = new String [nchannels*2];
			}
			if (headerLines.get(i).contains("DATA_INFO"))
			{
				for (int j = i+2; j < i+2+nchannels; j++)
				{
					tokens = headerLines.get(j).trim().split("\t");
					names[j-(i+2)] = tokens[1];
					names[j-(i+2) + nchannels] = tokens[1] + "back";
				}
			}
		}
		return names;
	}
	public static double[] getScanSize(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/64);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("SCANIT_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		String[] names = null;
		int nchannels = 0;
		String thing;
		String[] tokens;
		//find all important data in the header
		double lx = 0, ly = 0;
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("SCAN_RANGE"))
			{
				thing = headerLines.get(i+1).trim();
				tokens = thing.split(" ");
				lx = Double.parseDouble(tokens[0]);
				ly = Double.parseDouble(tokens[tokens.length-1]);
			}
		}
		return new double [] {lx, ly};
	}
	public static double[] getLineTimeFromScan(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length()/64);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !s.contains("SCANIT_END"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		double[] time = new double [2];
		int nchannels = 0;
		String thing;
		String[] tokens;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("SCAN_TIME"))
			{
				thing = headerLines.get(i+1).trim();
				nchannels = thing.split(";").length;
				tokens = thing.split(" ");
				time[0] = Double.parseDouble(tokens[0]);
				time[1] = Double.parseDouble(tokens[tokens.length-1]);
			}
		}
		return time;
	}
	public static double getAngleFromFile(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, 1024);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !(s.contains("SCANIT_END") || s.contains("HEADER_END")))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		double angle = -1;
		int nchannels = 0;
		String thing;
		String[] tokens;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains("SCAN_ANGLE"))
			{
				thing = headerLines.get(i+1).trim();
				angle = Double.parseDouble(thing);
			}
		}
		return angle;
	}
	public static double getVoltageFromScan(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, 1024);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (s.length() < 1000 && !(s.contains("SCANIT_END") || s.contains("HEADER_END")))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		for (int i  = 0; i < headerLines.size(); i++)
//			System.out.println(headerLines.get(i));
		double bias = -1;
		int nchannels = 0;
		String thing;
		String[] tokens;
		//find all important data in the header
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).contains(":BIAS:"))
			{
				thing = headerLines.get(i+1).trim();
				bias = Double.parseDouble(thing);
			}
		}
		return bias;
	}
	public static void convertScanToBins(boolean writeImages)
	{
		Topomap.setStdDir();
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		JFileChooser fci = new JFileChooser(nanonisDir);
		String path = FileOps.selectOpen(fci).toString();
		String dest = FileOps.selectSave(fco).toString();
		convertScanToBins(path, dest, writeImages);
	}
	public static Layer convertScanToBins(String sourceFilePath, String dest, boolean writeImages)
	{
		Layer[] t = loadLayerFromScan(sourceFilePath, true);
		String[] names = loadLayerNamesFromScan(sourceFilePath);
		for (int i = 0; i < names.length; i++)
		{
			if (names[i].contains("\\")) names[i] = names[i].split("\\")[0];
			if (names[i].contains("/")) names[i] = names[i].split("/")[0];
			if (writeImages)
				RHKFileOps.write2Images(dest + names[i], t[i].data);
			Layer.writeBIN(t[i], dest + names[i] + ".bin");
		}
		ColumnIO.writeString(sourceFilePath, dest + "_source file.txt");
		return t[0];
	}
	public static Layer getBlankLayer(String sourceFilePath)
	{
		Layer[] t = loadLayerFromScan(sourceFilePath, true);
		return t[0];
	}
	
	public static double[][] getFilterRangesFromUser(String path)
	{
		if (JOptionPane.showConfirmDialog(null, "Get the range from a file?") == JOptionPane.YES_OPTION)
		{
			JFileChooser temp = new JFileChooser(path);
			String[] lines  = FileOps.openLines(temp);
			double[][] filter = new double [lines.length][];
			String[] token;
			for (int i = 0; i < lines.length; i++)
			{
				token = lines[i].split(",");
				filter[i] = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
			}
			return filter;
			
		}
		int nRanges = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of frequency ranges you wish to remove."));
		double[][] filter = new double [nRanges][];
		String input;
		String[] token;
		for (int i = 0; i < nRanges; i++)
		{
			input = JOptionPane.showInputDialog("Enter the minimum and maximum of the range, comma separated.");
			token = input.split(",");
			filter[i] = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
		}
		return filter;
	}
	public static double[][][] getFilterRanges3FromUser(String path)
	{
		JFileChooser temp = new JFileChooser(path);
		double[][][] filterRanges = new double [3][][];
		for (int k = 0; k < 3; k++)
		{
			if (JOptionPane.showConfirmDialog(null, "Get the range from a file?") == JOptionPane.YES_OPTION)
			{
				String[] lines  = FileOps.openLines(temp);
				double[][] filter = new double [lines.length][];
				
				String[] token;
				for (int i = 0; i < lines.length; i++)
				{
					token = lines[i].split(",");
					filter[i] = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
				}
				filterRanges[k] = filter;
			}
			else{
				int nRanges = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of frequency ranges you wish to remove."));
				double[][] filter = new double [nRanges][];
				String input;
				String[] token;
				for (int i = 0; i < nRanges; i++)
				{
					input = JOptionPane.showInputDialog("Enter the minimum and maximum of the range, comma separated.");
					token = input.split(",");
					filter[i] = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
				}
				filterRanges[k] = filter;
			}
		}
		return filterRanges;
	}
	public static void convertScanToBinsTimeFilter(boolean writeImages, double[][] filterRanges)
	{
		Topomap.setStdDir();
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		JFileChooser fci = new JFileChooser(nanonisDir);
		String path = FileOps.selectOpen(fci).toString();
		String dest = FileOps.selectSave(fco).toString();
		convertScanToBinsWithTimeFilter(path, dest, writeImages, filterRanges);
	}
	public static void convertScanToBinsTimeFilter(boolean writeImages, double[][] filterRanges, String path)
	{
		Topomap.setStdDir();
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		String dest = FileOps.selectSave(fco).toString();
		convertScanToBinsWithTimeFilter(path, dest, writeImages, filterRanges);
//		System.exit(0);
	}
	public static void convertScanToBinsTimeFilter(boolean writeImages, double[][][] filterRanges, String path)
	{
		Topomap.setStdDir();
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		String dest = FileOps.selectSave(fco).toString();
		convertScanToBinsWithTimeFilter(path, dest, writeImages, filterRanges);
//		System.exit(0);
	}
	public static void convertScanToBinsWithTimeFilter(String sourceFilePath, String dest, boolean writeImages, double[][] filterRanges)
	{
		Layer blank = convertScanToBins(sourceFilePath, dest, false);
		blank.data = null;
		TimeSignal ts;
		File destf = new File(dest);
		String destdir = destf.getParent();
		String destName = destf.getName();
		destdir = destdir + "\\" + destName + " frequency filtered\\";
		File f = new File(destdir);
		if (!f.exists()) f.mkdir();
		String note = "All frequencies in the following ranges were suppressed (Hz):\r\n";
		note += "From\t to\r\n";
		
		boolean scanWasDown = NanonisFileOps.scanWasDown(sourceFilePath);
		
		Layer[] pair;
		String[] suffix = {"", "back"};
		
		for (int i = 0; i < 3; i++)
		{
			ts = getTimeSignal(sourceFilePath, i, scanWasDown);
//			if (i == 0) GraphDrawerCart.plotGraph(ts.frequencies.clone(), ts.fftmag.clone());
			for (int j = 0; j < filterRanges.length; j++)
			{
//				Robo.wait(10000);
				ts.suppressFrequencies(filterRanges[j][0], filterRanges[j][1]);
				if (i == 0) note += Printer.arrayLnHorizontal(filterRanges[j]);
//				if (i == 0) GraphDrawerCart.plotGraph(ts.frequencies.clone(), ts.fftmag.clone());
			}
			pair = ts.getLayersIFFT(blank);
			
			for (int j = 0; j < 2; j++)
			{
				Layer.writeBIN(pair[j], destdir + TimeSignal.codes[i] + suffix[j] + ".bin");
				if (writeImages)
					RHKFileOps.write2Images(destdir + TimeSignal.codes[i] + suffix[j], pair[j].data);
			}
		}
		ColumnIO.writeString(note, destdir + "Filtering note.txt");
		ColumnIO.writeString(sourceFilePath, dest + "_source file.txt");
	}
	public static Layer[][] getWithTimeFilter(String sourceFilePath, double[][][] filterRanges)
	{
		Layer blank = getBlankLayer(sourceFilePath);
//		blank.data = null;
		TimeSignal ts;
//		File destf = new File(dest);
		String note = "All frequencies in the following ranges were suppressed (Hz):\r\n";
		note += "From\t to\r\n";
		
		Layer[] pair;
		Layer[][] ans = new Layer[3][2];
		String[] suffix = {"", "back"};
		boolean scanDown = NanonisFileOps.scanWasDown(sourceFilePath);
		for (int i = 0; i < 3; i++)
		{
			ts = getTimeSignal(sourceFilePath, i, scanDown);
//			if (i == 0) GraphDrawerCart.plotGraph(ts.frequencies.clone(), ts.fftmag.clone());
			for (int j = 0; j < filterRanges[i].length; j++)
			{
//				Robo.wait(10000);
				ts.suppressFrequencies(filterRanges[i][j][0], filterRanges[i][j][1]);
				if (i == 0) note += Printer.arrayLnHorizontal(filterRanges[i][j]);
//				if (i == 0) GraphDrawerCart.plotGraph(ts.frequencies.clone(), ts.fftmag.clone());
			}
			pair = ts.getLayersIFFT(blank);
//			if (!scanDown)
//			{
//				pair[0].flipY();
//				pair[1].flipY();
//			}
			ans[i] = pair;
		}
		return ans;
	}
	
	/**
	 * Does Fourier-filtering of the time signal with different filters for Z, current, and Lockin.
	 * the array is {filterRangesZ, filterRangesCurrent, filterRangesLockin} i.e. the index matches the LAYER_CODE.
	 * @param sourceFilePath
	 * @param dest
	 * @param writeImages
	 * @param filterRangesZ
	 * @param filterRangesCurrent
	 * @param filterRangesX
	 */
	public static void convertScanToBinsWithTimeFilter(String sourceFilePath, String dest, boolean writeImages, double[][][] filterRanges)
	{
		Layer blank = convertScanToBins(sourceFilePath, dest, false);
		blank.data = null;
		TimeSignal ts;
		File destf = new File(dest);
		String destdir = destf.getParent();
		destdir = destdir + "\\frequency filtered\\";
		File f = new File(destdir);
		if (!f.exists()) f.mkdir();
		String note = "All frequencies in the following ranges were suppressed (Hz):\r\n";
		note += "From\t to\r\n";
		
		Layer[] pair;
		String[] suffix = {"", "back"};
		
		for (int i = 0; i < 3; i++)
		{
			ts = getTimeSignal(sourceFilePath, i);
//			if (i == 0) GraphDrawerCart.plotGraph(ts.frequencies.clone(), ts.fftmag.clone())
			for (int j = 0; j < filterRanges[i].length; j++)
			{
//				Robo.wait(10000);
				ts.suppressFrequencies(filterRanges[i][j][0], filterRanges[i][j][1]);
				note += Printer.arrayLnHorizontal(filterRanges[i][j]) + "\t" + TimeSignal.codes[i];
//				if (i == 0) GraphDrawerCart.plotGraph(ts.frequencies.clone(), ts.fftmag.clone());
			}
			pair = ts.getLayersIFFT(blank);
			
			for (int j = 0; j < 2; j++)
			{
				Layer.writeBIN(pair[j], destdir + TimeSignal.codes[i] + suffix[j] + ".bin");
				if (writeImages)
					RHKFileOps.write2Images(destdir + TimeSignal.codes[i] + suffix[j], pair[j].data);
			}
		}
		ColumnIO.writeString(note, destdir + "Filtering note.txt");
		ColumnIO.writeString(sourceFilePath, dest + "_source file.txt");
	}
	
	public static Layer getLayerZFrom3ds(File f)
	{
		Topomap[] t = loadTopomapFrom3ds(f.toString());
		return t[t.length-1].getLayer(0);
	}
	public static Layer getLayerZFromScan(File f)
	{
		Layer[] t = loadLayerFromScan(f.toString(), false);
		String[] names = loadLayerNamesFromScan(f.toString());
		for (int i = 0; i < names.length; i++)
		{
			if (names[i].equalsIgnoreCase("Z"))
			{
				return t[i];
			}
		}
		return null;
	}
	
	public static Topomap getMapFromScanSeries(File[] f, int layerIndex)
	{
		Layer[] layers = new Layer[f.length];
		for (int i = 0; i < f.length; i++)
			layers[i] = loadSpecificLayerFromScan(f[i].toString(), layerIndex);
		return Topomap.newTopomap(layers);
	}
	
	public static int searchForName(String[] names, String name)
	{
		for (int i = 0; i < names.length; i++)
			if (names[i].contains(name))
				return i;
		return -1;
	}
	/**
	 * This does it in the simplest way possible, by looking at the forward and backward lines.
	 * @param f
	 * @param layerCode
	 * @return
	 */
	public static TimeSignal getTimeSignal(String f, int layerCode)
	{
		String[] names = loadLayerNamesFromScan(f);
		double[] times = getLineTimeFromScan(f);
		Layer forward = null, backward = null;
		int forwardI = -1, backwardI = -1;
		if (layerCode == TimeSignal.Z_CODE)
		{
			forwardI = searchForName(names, "Z");
//			backwardI = searchForName(names, "zback");
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.CURRENT_CODE)
		{
			forwardI = searchForName(names, "Current")*2;
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.LOCKIN_CODE)
		{
			forwardI = searchForName(names, "Lockin")*2;
			backwardI = forwardI+1;
		}
		System.out.println(forwardI + "\t" + backwardI);
		//I assume that y is the slow direction and x is the fast one.
//		Layer[] t = loadLayerFromScan(f);
		forward = loadSpecificLayerFromScan(f, forwardI);
		backward = loadSpecificLayerFromScan(f, backwardI);
		int nx = forward.nx, ny = forward.ny;
		double totalTime = (times[0] + times[1])*ny;
		double[] values = new double [2*nx*ny];
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				values[(i+2*nx*j)] = forward.data[i][j];
				values[(i+(2*j+1)*nx)] = backward.data[nx-i-1][j];
			}
		
		return new TimeSignal(values, totalTime);
	}
	public static TimeSignal getTimeSignal(String f, int layerCode, boolean scanWasDown)
	{
		String[] names = loadLayerNamesFromScan(f);
		double[] times = getLineTimeFromScan(f);
		Layer forward = null, backward = null;
		int forwardI = -1, backwardI = -1;
		if (layerCode == TimeSignal.Z_CODE)
		{
			forwardI = searchForName(names, "Z");
//			backwardI = searchForName(names, "zback");
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.CURRENT_CODE)
		{
			forwardI = searchForName(names, "Current")*2;
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.LOCKIN_CODE)
		{
			forwardI = searchForName(names, "Lockin")*2;
			backwardI = forwardI+1;
		}
		System.out.println(forwardI + "\t" + backwardI);
		//I assume that y is the slow direction and x is the fast one.
//		Layer[] t = loadLayerFromScan(f);
		forward = loadSpecificLayerFromScan(f, forwardI);
		backward = loadSpecificLayerFromScan(f, backwardI);
		int nx = forward.nx, ny = forward.ny;
		double totalTime = (times[0] + times[1])*ny;
		double[] values = new double [2*nx*ny];
		if (!scanWasDown)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
				{
					values[(i+2*nx*j)] = forward.data[i][j];
					values[(i+(2*j+1)*nx)] = backward.data[nx-i-1][j];
				}
		else
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
				{
					values[(i+2*nx*j)] = backward.data[i][j];
					values[(i+(2*j+1)*nx)] = forward.data[nx-i-1][j];
				}
		return new TimeSignal(values, totalTime);
	}
	public static double[] getTimeSignalValues(String f, int layerCode)
	{
		String[] names = loadLayerNamesFromScan(f);
		double[] times = getLineTimeFromScan(f);
		Layer forward = null, backward = null;
		int forwardI = -1, backwardI = -1;
		if (layerCode == TimeSignal.Z_CODE)
		{
			forwardI = searchForName(names, "Z");
//			backwardI = searchForName(names, "zback");
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.CURRENT_CODE)
		{
			forwardI = searchForName(names, "Current")*2;
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.LOCKIN_CODE)
		{
			forwardI = searchForName(names, "Lockin_X/R")*2;
			backwardI = forwardI+1;
		}
		System.out.println(forwardI + "\t" + backwardI);
		//I assume that y is the slow direction and x is the fast one.
//		Layer[] t = loadLayerFromScan(f);
		forward = loadSpecificLayerFromScan(f, forwardI);
		backward = loadSpecificLayerFromScan(f, backwardI);
		int nx = forward.nx, ny = forward.ny;
		double totalTime = (times[0] + times[1])*ny;
		double[] values = new double [2*nx*ny];
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				values[(i+2*nx*j)] = forward.data[i][j];
				values[(i+(2*j+1)*nx)] = backward.data[nx-i-1][j];
			}
		
		return values;
	}
	public static TimeSignal getTimeSignalFromValues(String f, int layerCode, double[] data)
	{
		String[] names = loadLayerNamesFromScan(f);
		double[] times = getLineTimeFromScan(f);
		Layer forward = null, backward = null;
		int forwardI = -1, backwardI = -1;
		if (layerCode == TimeSignal.Z_CODE)
		{
			forwardI = searchForName(names, "Z");
//			backwardI = searchForName(names, "zback");
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.CURRENT_CODE)
		{
			forwardI = searchForName(names, "Current")*2;
			backwardI = forwardI+1;
		}
		else if (layerCode == TimeSignal.LOCKIN_CODE)
		{
			forwardI = searchForName(names, "Lockin_X/R")*2;
			backwardI = forwardI+1;
		}
		System.out.println(forwardI + "\t" + backwardI);
		//I assume that y is the slow direction and x is the fast one.
//		Layer[] t = loadLayerFromScan(f);
		forward = loadSpecificLayerFromScan(f, forwardI);
		backward = loadSpecificLayerFromScan(f, backwardI);
		int ny = forward.ny;
		double totalTime = (times[0] + times[1])*ny;
		return new TimeSignal(data, totalTime);
	}
	public static class TimeSignal
	{
		double[] values;
		double totalTime;
		double pixelTime, maxFreq; //maxfreq is the absolute maximum frequency. the real unique maximum is maxfreq/2.
		double[] time;
		double[] fftmag;
		double[] fftmagLog;
		double[] fftmagHalf;
		double[] frequencies;
		double[] freqHalf;
		
		double[] fftComplex;
		double[] fftClone;
		double[] ifftComplex;
		double[] ifftReal;
		static final int Z_CODE = 0;
		static final int CURRENT_CODE = 1;
		static final int LOCKIN_CODE = 2;
		static String[] codes = {"Z", "Current", "LockinX"};
		
		double magmin = 0;
		
		public TimeSignal(double[] values, double totalTime) {
			super();
			this.values = values;
			this.totalTime = totalTime;
			
			time = ArrayOps.generateArrayInclBoth(0, totalTime, values.length);
			pixelTime = totalTime/values.length;
			
			maxFreq = values.length/(totalTime);
			frequencies = ArrayOps.generateArrayInclBoth(0, maxFreq, values.length);
			fftComplex = FFTOps.get1DFFTComplex(values);
			fftClone = fftComplex.clone();
			ifftReal = new double [values.length];
			fftmag = new double [values.length];
			fftmagLog = new double [values.length];
			for (int i = 0; i < values.length; i++)
			{
				fftmag[i] = Complex.mag(fftClone[2*i], fftClone[2*i+1]);
				fftmagLog[i] = Math.log(fftmag[i]);	
			}
			fftmagHalf = new double [values.length/2];
			freqHalf = new double [values.length/2];
			for (int i = 0; i < fftmagHalf.length; i++)
			{
				fftmagHalf[i] = fftmag[i];
				freqHalf[i] = frequencies[i];
			}
			
			magmin = ArrayOps.min(fftmag);
//			fftmag = null;
		}
		public double frequency(int i)
		{
			return i*(maxFreq/(values.length));
		}
		public double getNatrualFrequency(double realFrequency)
		{
			return (realFrequency/maxFreq)*2*Math.PI;
		}
		public void suppressFrequencies(double fmin, double fmax)
		{
			if (fmax >= maxFreq) {
				System.out.println("Error. " + fmax + " was outside of the range. Returning original values.");
				FFTOps.putIFFTReal(fftClone.clone(), ifftReal);
//				return;
			}

			
			double f2max = -fmax + maxFreq;
			double f2min = -fmin + maxFreq;
			double ratio;
			for (int i = 0; i < values.length; i++){
				if ((frequency(i) >= fmin && frequency(i) < fmax) || (frequency(i) >= f2max && frequency(i) < f2min))
				{
					ratio = Complex.mag(fftClone[2*i], fftClone[2*i+1])/magmin;
					fftClone[2*i] /= ratio;
					fftClone[2*i+1] /= ratio;
//					System.out.println(i + "\t" + ratio + "\t" + fmin + "\t" + fmax + "\t" + magmin + "\t" + Complex.mag(fftClone[2*i], fftClone[2*i+1]));
				}
			}
			FFTOps.putIFFTReal(fftClone.clone(), ifftReal);
			for (int i = 0; i < values.length; i++)
			{
				fftmag[i] = Complex.mag(fftClone[2*i], fftClone[2*i+1]);
				fftmagLog[i] = Math.log(fftmag[i]);	
			}
		}
		
		public Layer[] getLayersIFFT(Layer t)
		{
			int nx = t.nx, ny = t.ny;
			double[][][] layerData = new double [2][t.nx][t.ny]; //the [2] are {forward, backward};
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
				{
					 layerData[0][i][j] = ifftReal[(i+2*nx*j)];
					 layerData[1][nx-i-1][j] = ifftReal[(i+(2*j+1)*nx)];
				}
			Layer forward = Layer.newLayer(t, layerData[0]);
			Layer backward = Layer.newLayer(t, layerData[1]);
			return new Layer[] {forward, backward};
		}
		public static TimeSignal getSineWave(double totalTime, double freqHertz, int length)
		{
			double maxFreq = length/(totalTime);
			double k = (freqHertz/maxFreq)*2*Math.PI;
			double[] values = new double [length];
			for (int i = 0; i < values.length; i++)
				values[i] = Math.cos(k*i);
			return new TimeSignal(values, totalTime);
		}
	}
	//This assumes that the number of spectra and sweeps are saved!
	public static PointSpectra getFromDat(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Topomap[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (!s.contains("[DATA]"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String dataColumnLine = null;
		try {
			dataColumnLine = ind.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String[] dataColumns = dataColumnLine.split("\t");
		int nspec; 
		if (!dataColumns[1].contains("[AVG]"))
			nspec = 1;
		else{
			String p = dataColumns[dataColumns.length-1];
			nspec = Integer.parseInt(p.substring(p.indexOf('[')+1, p.indexOf(']')));
		}
		
		double[] x = new double [nspec];
		double[] y = new double [nspec];
		int npts;
		String thing;
		String[] tokens;
		double tokenData;
		
		//The first line is the experiment type. In the header we must find x and y.
		char exp = 0;
		if (headerLines.get(0).contains("bias spectroscopy")) exp = 'b';
		else if (headerLines.get(0).contains("Z spectroscopy")) exp = 'z';
		for (int i = 0; i < headerLines.size(); i++)
		{
			if (headerLines.get(i).startsWith("X (m)"))
			{
				thing = headerLines.get(i).trim();
				tokens = thing.split("\t");
				tokenData = Double.parseDouble(tokens[1]);
				for (int j = 0; j < nspec; j++)
					x[j] = tokenData;
			}
			if (headerLines.get(i).startsWith("Y (m)"))
			{
				thing = headerLines.get(i).trim();
				tokens = thing.split("\t");
				tokenData = Double.parseDouble(tokens[1]);
				for (int j = 0; j < nspec; j++)
					y[j] = tokenData;
			}
		}
		//Now we load the spectral data. We don't know how many voltages there are.
		//We will find out presently.
		ArrayList <String> dataLines = new ArrayList<String>();
		s = "test";
		try {
			while (! (s == null || s.trim().equals("")))
			{
				dataLines.add(ind.readLine());
				s = dataLines.get(dataLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if (dataLines.get(dataLines.size()-1) == null)
			dataLines.remove(dataLines.size()-1);
		if (dataLines.get(dataLines.size()-1).trim().equals(""))
			dataLines.remove(dataLines.size()-1);
		
		npts = dataLines.size();
		String[][] dataMatrix = new String[dataLines.size()][dataColumns.length];
		for (int i = 0; i < dataLines.size(); i++)
			dataMatrix[i] = dataLines.get(i).split("\t");
		
		//Now:
		double[] v = new double [npts];
		//To assign the voltage we must pick the "Bias from Adder column"
		//NOTE this is only true on the Nanonis system with external adder
		int biasI = -1;
		for (int i = 0; i < dataColumns.length; i++)
			if (exp == 'b' && dataColumns[i].contains("from") && (nspec == 1 || dataColumns[i].contains("[AVG]")))
					biasI = i;
		
		if (exp == 'z') biasI = 0;
			
		ArrayList<Integer> dataIndices = new ArrayList<Integer>();
		if (exp == 'b')
			for (int i = 0; i < dataColumns.length; i++)
				if (dataColumns[i].contains("Lockin X") && (nspec == 1 || (dataColumns[i].contains("[") && dataColumns[i].contains("]") && !dataColumns[i].contains("AVG"))))
					dataIndices.add(i);
		
		if (exp == 'z')
			for (int i = 0; i < dataColumns.length; i++)
				if (dataColumns[i].contains("Current") && (nspec == 1 || (dataColumns[i].contains("[") && dataColumns[i].contains("]") && !dataColumns[i].contains("AVG"))))
					dataIndices.add(i);
		System.out.println(filepath);
		double[][] data = new double [nspec][npts];
		for (int i = 0; i < nspec; i++)
			for (int j = 0; j < npts; j++)
			{
				data[i][j] = Double.parseDouble(dataMatrix[j][dataIndices.get(i)]);
			}
		
		for (int i = 0; i < npts; i++)
			v[i] = Double.parseDouble(dataMatrix[i][biasI]);
		
		return new PointSpectra(data, v, x, y);
	}
	/**
	 * This is intended to get FFT spectra exported from the LT spectrum viewer. X and Y are zero, v is frequency.
	 * @param filepath
	 * @return
	 */
	public static PointSpectra getFromASC(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);

		String s = "";
		ArrayList<String> headerLines = new ArrayList<String> ();
		try {
			while (!s.contains(":HEADER_END:"))
			{
				headerLines.add(ind.readLine());
				s = headerLines.get(headerLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String dataColumnLine = null;
		try {
			dataColumnLine = ind.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String[] firstLine = dataColumnLine.split("\t");
		int npts = firstLine.length;
		double[] v = new double [npts];
		for (int i = 0; i < npts; i++)
			v[i] = Double.parseDouble(firstLine[i]);
		//Now we load the spectral data. We don't know how many voltages there are.
		//We will find out presently.
		ArrayList <String> dataLines = new ArrayList<String>();
		s = "test";
		try {
			while (! (s == null || s.trim().equals("")))
			{
				dataLines.add(ind.readLine());
				s = dataLines.get(dataLines.size()-1);
//				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if (dataLines.get(dataLines.size()-1) == null)
			dataLines.remove(dataLines.size()-1);
		if (dataLines.get(dataLines.size()-1).trim().equals(""))
			dataLines.remove(dataLines.size()-1);
		
		int nspec = dataLines.size();
		
		//Now:
		double[][] data = new double [nspec][npts];
		for (int i = 0; i < nspec; i++){
			String[] tok = dataLines.get(i).split("\t");
			for (int j = 0; j < npts; j++)
			{
				data[i][j] = Double.parseDouble(tok[j]);
			}
		}
		
		return new PointSpectra(data, v, new double [nspec], new double [nspec]);
	}
}
