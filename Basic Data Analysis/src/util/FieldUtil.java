package util;

import java.io.File;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import drawing.GaussSquareImpurityAngleSuite;

import main.SRAW;
import schrodinger.MovieMaker;
import util.calc.StripCut;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.Map_4D;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFTOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.geom.Mask;
import util.matrix.Matrix;

public class FieldUtil {

	//converts 2d data in bin file to bmp
	public static void makeImage(String in, String out)
	{
		double[][] data = ColumnIO.readSquareTable(in);
		SRAW.writeImage(out, data);
	}
	public static void makeImage(String in, String out, boolean bin)
	{
		double[][] data;
		if (bin) data = ColumnIO.readSquareTable(in);
		else data = SRAW.getData(in);
		SRAW.writeImage(out, data);
	}
	
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser (Topomap.stddir);
		
		Layer p = Layer.openFree(fc);
		Layer q = Layer.openFree(fc);
		FileOps.writeLines(fc, printXYForTwoFields(p.data, q.data));
		
//		Layer big = Layer.open(fc);
//		for (int i = 0; i < big.nx; i++)
//			for (int j = 0; j < big.ny; j++)
//				big.data[i][j] *= Math.PI*2*2*2;
//		Layer.writeBIN(big, fc);
//		Layer big =
		//		Layer small = Layer.open();
//		double[][] smaller = FieldOps.reduce(2, small.data);
//		double[][] map = FieldOps.getCorrelationMap(smaller, big.data);
//		Layer mapL = Layer.getFreeLayer(map);
//		LayerViewer.show(mapL, mapL.nx);
//		Layer.writeBIN(mapL);
		
//		double[] values = ArrayOps.generateArrayInclBoth(0, 10, 500);
//		for (int i = 0; i < values.length; i++)
//			System.out.println(((values[i]) % (2*Math.PI)) - Math.PI);
//		
		//		int n = 768;
//		double[][] f = new double[n][n];
//		Layer l = Layer.getFreeLayer(f);
//		LayerViewer.show(l, n);
		
//		String dir = "C:\\Users\\madhavanlab2011\\Documents\\PPTs\\2\\";
//		dir = "C:\\data\\analysis\\MEM\\MEM_TEST_2\\";
//		
//		Map_4D map = loadMEMOutput_DAT(0.181, 0.001, 0.019, 0.001, dir, false);
//		Map_4D.writeBIN(map);
//		double[][] z = new double [1024][1024];
//		for (int i = 0; i < z.length; i++)
//			for (int j = 0; j < z[0].length; j++)
//			{
//				z[i][j] = i + j + Math.random() + 5;
//			}
//		FieldOps.fitToPlaneSquare(z);
		
//		Layer t = loadLayer(0.180, 0.000156, 0.018, 0.000375, dir);
//		Layer.writeBIN(t);
//		JFileChooser fc = FileOps.getFC();
//		Topomap t = Topomap.open();
//		ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps = GaussSquareImpurityAngleSuite.getImpuritiesFull(FileOps.selectOpen(null));
//		ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps;
//		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText());
//		for (int i = 1; i <= 10; i++)
//		{
//			dir = "C:\\data\\analysis\\SrIrO 327\\241\\images\\resolution dependence\\" + i + "\\";
//			imps = GaussSquareImpurityAngleSuite.getImpuritiesFull(new File(dir + "implist.txt"));
//			doStuffWithImps(null, imps, latt);
//		}
//		double[][][] u = FileOps.open2Comp(true);
//		double[][][] uShift = shiftUOnePixel(u);
//		
//		Layer t = Layer.openFree();
//		t.data = applyUFieldBigger(u, t.data);
//		Layer.writeBIN(t);
//		Topomap t = Topomap.open();
//		Topomap ts = TopomapUtil.applyUField(uShift, t);
//		Topomap ts = TopomapUtil.applyUField(u, t);
//		Topomap.writeBIN(ts);
//		Topomap tss = TopomapUtil.shiftLeftOnePixel(ts);
//		Layer t = Layer.openFree();
//		FileOps.writeTableBIN(applyUFieldLayer(u, t));
//		createTopomapSpectralFunction();
	}

	/**
	 * If the topography is double the size of the u-field. (i.e. if the u-field was derived from a reduced version of the topography).
	 * @param u
	 * @param biggerTopo
	 */
	public static double[][] applyUFieldBigger(double[][][] u, double[][] biggerTopo)
	{
		double[][][] ubig = FieldOps.expand(u);
		//the length scale has also expanded
		FieldOps.multiply(ubig, 2);
		return applyUField(ubig, biggerTopo);
	}
	
	public static void splitOneFileByHorizontalLine()
	{
		double[][] data = FileOps.openBin("");
		double[][][] split = FieldOps.splitByLineHorizontal(data);
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		RHKFileOps.write3Files(fc, split[0]);
		RHKFileOps.write3Files(fc, split[1]);

	}
	public static void doImpurityCorrelations(String dir)
	{
		JFileChooser fc = new JFileChooser(dir);
		File f = null;
		JOptionPane.showMessageDialog(null, "Select the topograph.");
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
			f = fc.getSelectedFile();
		double[][] topo = ColumnIO.readSquareTable(f.toString());
		
		dir = fc.getSelectedFile().toString().substring(0, f.toString().lastIndexOf("\\")+1);
		String output = dir + "corrmaps\\";
		f = new File(output);
		if (!f.exists()) f.mkdir();
		JOptionPane.showMessageDialog(null, "Select the impurity.");
		while(fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			f = fc.getSelectedFile();
			double[][] imp = ColumnIO.readSquareTable(f.toString());
			
			String impname = f.toString().substring(f.toString().lastIndexOf("\\")+1, f.toString().length()-4);

			writeCorrelationMap(imp, topo, output+impname);
			JOptionPane.showMessageDialog(null, "Select the impurity.");
		}
		printPath2Slash(dir);
	}
	public static void printTrueRegions(boolean[][] f)
	{
		printTrueRegions(FieldOps.isolatedTrueRegions(f));
	}
	public static void printTrueRegions(double[][] listRegions)
	{
		for (int i = 0; i < listRegions[0].length; i++)
		{
			System.out.println("Region "+  i +" centered at (" + listRegions[0][i] + ", " + listRegions[1][i] +")   with area " + listRegions[2][i]);
		}
	}
	public static void writeCorrelationMap(double[][] a, double[][] b, String f)
	{
		double[][] answer = FieldOps.getCorrelationMapCent(a, b);
		writeDatAndBMP(answer, f);
	}
	public static boolean[][] writeSplit(double[][] data, double cutoff, String dir, String name)
	{
		boolean[][] x = FieldOps.isGreaterThan(data, cutoff);
		ColumnIO.writeBin(data, dir + name + ".dat");
		SRAW.writeImage(dir + name, x);
		return x;
	}
	public static void makeImagesOfDirBin(String dir)
	{
		File[] files = new File(dir).listFiles();
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat"))
			{
				SRAW.writeImage(files[i].toString().substring(0, files[i].toString().length() - 4) + ".bmp", ColumnIO.readSquareTable(files[i].toString()));
			}
		}
	}
	public static void makeImagesOfDirBin(String dir, String keyword)
	{
		if (keyword.equals("")) {makeImagesOfDirBin(dir); return;}
		File[] files = new File(dir).listFiles();
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
			{
				SRAW.writeImage(files[i].toString().substring(0, files[i].toString().length() - 4) + ".bmp", ColumnIO.readSquareTable(files[i].toString()));
			}
		}
	}
	public static void makeImagesOfDirBin(String dir, String keyword, double min, double max)
	{
		if (keyword.equals("")) {makeImagesOfDirBin(dir); return;}
		File[] files = new File(dir).listFiles();
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
			{
				SRAW.writeImage(files[i].toString().substring(0, files[i].toString().length() - 4) + ".bmp", ColumnIO.readSquareTable(files[i].toString()), min, max);
			}
		}
	}
	public static void makeImagesOfDirBinMiddle(String dir, String keyword, double min, double max, int factor)
	{
		if (keyword.equals("")) {makeImagesOfDirBin(dir); return;}
		File[] files = new File(dir).listFiles();
		int minx, dx, miny, dy;
		double[][] data;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword))
			{	
				data = ColumnIO.readSquareTable(files[i].toString());
				dx = data.length/factor;
				dy = data[0].length/factor;
				minx = (data.length-dx)/2;
				miny = (data[0].length-dy)/2;
				SRAW.writeImage(files[i].toString().substring(0, files[i].toString().length() - 4) + ".bmp", data, min, max, minx, miny, dx, dy);
			}
		}
	}
	
	public static String[] printXYForTwoFields(double[][] x, double[][] y)
	{
		ArrayList<String> ansT = new ArrayList<String>();
		
		ansT.add("" + FieldOps.correlation(x, y));
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				ansT.add("" + x[i][j] + "\t" + y[i][j]);
		
		String[] ans = new String[ansT.size()];
		for (int i = 0 ; i < ans.length; i++)
			ans[i] = ansT.get(i);
		return ans;
	}
	public static void expandTopo(double kmag, double topolength, double bragglength, String dir, String name)
	{
		double klength = topolength/kmag;
		double factor = bragglength/klength;
		System.out.println(klength + "\t" + factor);
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat");
		double[][] result = FieldOps.rescaleXY(data, factor, 16);
		ColumnIO.writeBin(result, dir + name + "_resc.dat");
		SRAW.writeImage(dir + name + "_resc", result);
	}
	
	
	public static void phaseNPic(String in, String out)
	{
		double[][] phase = ColumnIO.readSquareTable(in);
		int[][] n = FieldOps.phaseSteps(phase, phase.length/2, phase.length/2);
		SRAW.writeImage(out, n);
	}
	
	public static void addPhaseWrite(String dir, String base, double[] L)
	{
		String[] names = new String[L.length];
		String[] out = new String [L.length];
		String[] mark = {"x", "y", "z"};
		for (int i = 0; i < L.length; i++)
		{
			if (L[i] != 0){
				names[i] = "devoutbinft\\" + base + "dev" + mark[i] + "L" + L[i] + "phase.dat";
				out[i] = "devoutbinft\\" + base +  mark[i] + "L" + L[i] + "cont";
				addPhaseWrite(dir + names[i], dir + out[i]);
			}
		}
	}
	public static void addPhaseWrite(String in, String out)
	{
		double[][] phase = ColumnIO.readSquareTable(in);
		int[][] n = FieldOps.phaseSteps(phase, phase.length/2, phase.length/2);
		for (int i = 0; i < phase.length; i++)
			for (int j = 0; j < phase.length; j++)
				phase[i][j] += n[i][j]*Math.PI*2;
		ColumnIO.writeBin(phase, out + ".dat");
		SRAW.writeImage(out + ".bmp", phase);
		SRAW.writeImage(out + "_n.bmp", n);
	}
	
	public static void makeFlatImage(String in, String out)
	{
		double[][] data = ColumnIO.readSquareTable(in);
		FieldOps.subtractFracAvg(data, 0.99);
		FieldOps.subtractSlopes1(data);
		SRAW.writeImage(out, data);
	}
	public static void makeFlatDir(String dir, String keyword, double frac)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth;
		String name;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				data = ColumnIO.readSquareTable(files[i].toString());
				FieldOps.subtractFracAvg(data, frac);
				FieldOps.subtractSlopes1(data);
				ColumnIO.writeBin(data, name + "flat.dat");
				SRAW.writeImage(name + "flat", data);
				FFTOps.writeFFTBMPCent(name + "flatfft", FFTOps.obtainFFT(data), true, true);
			}			
		}
	}
	public static void writeLocalAvgDir(String dir, String keyword, int L, int freq)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth;
		String name, subname;
		String newdir = dir + "locavg\\";
		String minusdir = dir + "subloc\\";
//		String divdir = dir + "divloc\\";
		if (!new File(newdir).exists())
			new File(newdir).mkdir();
		if (!new File(minusdir).exists())
			new File(minusdir).mkdir();
//		if (!new File(divdir).exists())
//			new File(divdir).mkdir();
		
		System.out.println("Started. L = " + L);
		double[][] temp;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				subname = name.substring(name.lastIndexOf("\\") + 1);
				data = ColumnIO.readSquareTable(files[i].toString());
				smooth = FieldOps.gaussSmooth(data, L, freq);
				ColumnIO.writeBin(smooth, newdir + subname + ".dat");
				SRAW.writeImage(newdir + subname, smooth);
				//subtract the local avg, just to see what happens:
				temp = FieldOps.minus(data, smooth);
				ColumnIO.writeBin(temp, minusdir + subname +".dat");
				SRAW.writeImage(minusdir + subname, temp);
				FFTOps.writeFFTBMPCent(minusdir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
//				FieldOps.overEquals(data, smooth);
//				ColumnIO.writeBin(data, divdir + subname + "divloc.dat");
//				SRAW.writeImage(divdir + subname + "divloc", data);
//				FFTOps.writeFFTBMPCent(divdir + subname + "divlocfft", FFTOps.obtainFFT(data), true, true);
			}			
		}
		System.out.println("Done. L = " + L);
	}
	public static void writeLocalAvgDirFT(String dir, String keyword, double L, boolean subtractplane)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth;
		double[][][] fftz, nfftz, ift;
		double[][] gauss;
		String name, subname;
		String newdir = dir + "locavg\\";
		String minusdir = dir + "subloc\\";
//		String divdir = dir + "divloc\\";
		if (!new File(newdir).exists())
			new File(newdir).mkdir();
		if (!new File(minusdir).exists())
			new File(minusdir).mkdir();
//		if (!new File(divdir).exists())
//			new File(divdir).mkdir();
		
		int N;
		double[][] temp;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				subname = name.substring(name.lastIndexOf("\\") + 1);
				data = ColumnIO.readSquareTable(files[i].toString());
				if (subtractplane) FieldOps.subtractSlopes1(data);
				N = data.length;
				fftz = new double [data.length][data[0].length][2];
				nfftz = new double [data.length][data[0].length][2];
				ift = new double [data.length][data[0].length][2];
				gauss = new double [data.length][data[0].length];
				FFTOps.obtainFFTCent(data, fftz);
				for (int k = 0; k < data.length; k++)
					for (int l = 0; l < data[0].length; l++)
						gauss[k][l] = Math.exp((-((double)(k*k) + (l*l))/(2*L*L)));
				gauss[0][0] = 0.9999;

				for (int k = 0; k < N; k++)
					for (int l = 0; l < N; l++)
					{
						nfftz[k][l][0] = fftz[k][l][0]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
						nfftz[k][l][1] = fftz[k][l][1]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
					}
				FFTOps.obtainIFFTCent(nfftz, ift);
				ColumnIO.writeBin(ift, 0, newdir + subname + "locavg.dat");
				SRAW.writeImage(newdir + subname + "locavg", ift, 0);
				//subtract the local avg, just to see what happens:
				temp = FieldOps.minus(data, ift, 0);
				ColumnIO.writeBin(temp, minusdir + subname + ".dat");
				SRAW.writeImage(minusdir + subname + "", temp);
				FFTOps.writeFFTBMPCent(minusdir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
			}			
		}
		System.out.println("Done. L = " + L);
	}
	public static void splitRegionsDir(String dir, String keyword, String regname)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		int[][] reg = ColumnIO.readSquareTableInt(regname + ".dat");
		
		String name, subname;
		String m1dir = dir + "rm1\\";
		String odir = dir + "r0\\";
		String p1dir = dir + "rp1\\";
		//		String divdir = dir + "divloc\\";
		if (!new File(m1dir).exists())
			new File(m1dir).mkdir();
		if (!new File(odir).exists())
			new File(odir).mkdir();
		if (!new File(p1dir).exists())
			new File(p1dir).mkdir();
		
		int N;
		double[][] temp;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				subname = name.substring(name.lastIndexOf("\\") + 1);
				data = ColumnIO.readSquareTable(files[i].toString());
				N = data.length;
				temp = new double [N][N];
				double mean = 0;
				mean = FieldOps.getMeanRegion(data, reg, -1);
				FieldOps.mask(data, temp, reg, -1, mean);
				ColumnIO.writeBin(temp, m1dir + subname + ".dat");
				SRAW.writeImage(m1dir + subname, temp);
				FFTOps.writeFFTBMPCent(m1dir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
				
				mean = FieldOps.getMeanRegion(data, reg, 0);
				FieldOps.mask(data, temp, reg, 0, mean);
				ColumnIO.writeBin(temp, odir + subname + ".dat");
				SRAW.writeImage(odir + subname, temp);
				FFTOps.writeFFTBMPCent(odir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
				
				mean = FieldOps.getMeanRegion(data, reg, 1);
				FieldOps.mask(data, temp, reg, 1, mean);
				ColumnIO.writeBin(temp, p1dir + subname + ".dat");
				SRAW.writeImage(p1dir + subname, temp);
				FFTOps.writeFFTBMPCent(p1dir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
			}			
		}
	}
	public static void splitRegionsDirPerc(String dir, String keyword, String topodir, double lowerperc, double upperperc)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		int[][] reg;
		
		if (topodir.equals("")) topodir = dir;
		double[][] topo;
		
		String name, subname;
		String m1dir = dir + "rm1\\";
		String odir = dir + "r0\\";
		String p1dir = dir + "rp1\\";
		//		String divdir = dir + "divloc\\";
		if (!new File(m1dir).exists())
			new File(m1dir).mkdir();
		if (!new File(odir).exists())
			new File(odir).mkdir();
		if (!new File(p1dir).exists())
			new File(p1dir).mkdir();
		
		int N;
		double[][] temp;
		double mean;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				subname = name.substring(name.lastIndexOf("\\") + 1);
				topo = ColumnIO.readSquareTable(topodir + subname + ".dat");
				data = ColumnIO.readSquareTable(files[i].toString());
				N = data.length;
				reg = new int [N][N];
				temp = new double [N][N];
				FieldOps.cutOffExtremes(topo, lowerperc, upperperc, FieldOps.getSortedDump(topo), temp, reg);
				mean = FieldOps.getMeanRegion(data, reg, -1);
				FieldOps.mask(data, temp, reg, -1, mean);
				ColumnIO.writeBin(temp, m1dir + subname + ".dat");
				SRAW.writeImage(m1dir + subname, temp);
				FFTOps.writeFFTBMPCent(m1dir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
				
				mean = FieldOps.getMeanRegion(data, reg, 0);
				FieldOps.mask(data, temp, reg, 0, mean);
				ColumnIO.writeBin(temp, odir + subname + ".dat");
				SRAW.writeImage(odir + subname, temp);
				FFTOps.writeFFTBMPCent(odir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
				
				mean = FieldOps.getMeanRegion(data, reg, 1);
				FieldOps.mask(data, temp, reg, 1, mean);
				ColumnIO.writeBin(temp, p1dir + subname + ".dat");
				SRAW.writeImage(p1dir + subname, temp);
				FFTOps.writeFFTBMPCent(p1dir + subname + "fft", FFTOps.obtainFFT(temp), true, true);
			}			
		}
	}
	
	public static void splitRegionsPercTopomaps(String dir, Topomap topo, Topomap didv, double lowerperc, double upperperc)
	{
		String name, subname;
		String dirname = dir.substring(0, dir.length() - 1);
		String m1dir = dirname + "_rm1\\";
		String odir = dirname + "_r0\\";
		String p1dir = dirname + "_rp1\\";
		//		String divdir = dir + "divloc\\";
		if (!new File(m1dir).exists())
			new File(m1dir).mkdir();
		if (!new File(odir).exists())
			new File(odir).mkdir();
		if (!new File(p1dir).exists())
			new File(p1dir).mkdir();
		
		int N;
		double[][] temp;
		double[][] topotemp;
		double[][] didvtemp;
		double mean;
		int[][] reg;
		for (int i = 0; i < topo.v.length; i++)
		{
			N = topo.x.length;
			reg = new int [N][N];
			temp = new double [N][N];
			topotemp = FieldOps.copy(topo.data[i]);
			didvtemp = FieldOps.copy(didv.data[i]);
			FieldOps.cutOffExtremes(topotemp, lowerperc, upperperc, FieldOps.getSortedDump(topotemp), temp, reg);
			mean = FieldOps.getMeanRegion(didvtemp, reg, -1);
			FieldOps.mask(didvtemp, temp, reg, -1, mean);
			ColumnIO.writeBin(temp, m1dir + didv.names[i] + ".dat");
			SRAW.writeImage(m1dir + didv.names[i], temp);
			FFTOps.writeFFTBMPCent(m1dir + didv.names[i] + "fft", FFTOps.obtainFFT(temp), true, true);
			
			mean = FieldOps.getMeanRegion(didvtemp, reg, 0);
			FieldOps.mask(didvtemp, temp, reg, 0, mean);
			ColumnIO.writeBin(temp, odir + didv.names[i] + ".dat");
			SRAW.writeImage(odir + didv.names[i], temp);
			FFTOps.writeFFTBMPCent(odir + didv.names[i] + "fft", FFTOps.obtainFFT(temp), true, true);
			
			mean = FieldOps.getMeanRegion(didvtemp, reg, 1);
			FieldOps.mask(didvtemp, temp, reg, 1, mean);
			ColumnIO.writeBin(temp, p1dir + didv.names[i] + ".dat");
			SRAW.writeImage(p1dir + didv.names[i], temp);
			FFTOps.writeFFTBMPCent(p1dir + didv.names[i] + "fft", FFTOps.obtainFFT(temp), true, true);
		}
	}
	public static void doCutoffDir(String dir, String keyword, double gamma, String suboutdir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				smooth = ColumnIO.readSquareTable(avgdir + shortname + "locavg.dat");
				answer = FieldOps.cutOffExtremes(data, smooth, gamma);
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);

			}			
		}
	}
	public static void doCutoffDir(String dir, String keyword, double gamma, String suboutdir, String suboutdirmin)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				smooth = ColumnIO.readSquareTable(avgdir + shortname + "locavg.dat");
				answer = FieldOps.cutOffExtremes(data, smooth, gamma);
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);
				
				FieldOps.minusEquals(answer, smooth);
				if (!new File(dir + suboutdirmin).exists())
					new File(dir + suboutdirmin).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdirmin + shortname + ".dat");
				SRAW.writeImage(dir + suboutdirmin + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdirmin + shortname + "fft", FFTOps.obtainFFT(answer), true, true);

			}			
		}
	}
	public static void doCutoffDir(String dir, String keyword, double gammaup, double gammadown, String suboutdir, String suboutdirmin)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		
		String description = "Gamma_up = " + gammaup + "\tGamma_down = " + gammadown + "\r\n";
		description += "Name\tSigma\r\n";
		double sigma;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				sigma = FieldOps.sigma(data);
				description += shortname + "\t" + sigma;
				smooth = ColumnIO.readSquareTable(avgdir + shortname + "locavg.dat");
				answer = FieldOps.cutOffExtremes(data, smooth, gammaup, gammadown, sigma);
				description += "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);
				
				FieldOps.minusEquals(answer, smooth);
				if (!new File(dir + suboutdirmin).exists())
					new File(dir + suboutdirmin).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdirmin + shortname + ".dat");
				SRAW.writeImage(dir + suboutdirmin + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdirmin + shortname + "fft", FFTOps.obtainFFT(answer), true, true);

			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	
	//This one assumes that the local average is 0.
	public static void doCutoffDir(String dir, String keyword, double gammaup, double gammadown, String suboutdir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, answer;
		String name, shortname;
		String description = "Gamma_up = " + gammaup + "\tGamma_down = " + gammadown + "\r\n";
		description += "Name\tSigma\r\n";
		double sigma;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				sigma = FieldOps.sigma(data);
				description += shortname + "\t" + sigma;
				answer = FieldOps.cutOffExtremes(data, gammaup, gammadown, sigma);
				description += "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);
			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	//This one only cuts values from given percentiles.
	
	public static void doCutoffDirPerc(String dir, String keyword, double lowerperc, double upperperc, String suboutdir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name, shortname;
		String description = "Lower Pecentile = " + (lowerperc*100) + "\tUpper Percentile = " + (upperperc*100) + "\r\n";
		description += "Name\tSigma\r\n";
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				FieldOps.cutOffExtremes(data, lowerperc, upperperc);
				description += "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(data, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, data);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(data), true, true);
			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	public static void doReplaceDirPerc(String dir, String keyword, double lowerperc, double upperperc, String suboutdir, double def)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name, shortname;
		String description = "Lower Pecentile = " + (lowerperc*100) + "\tUpper Percentile = " + (upperperc*100) + "\r\n";
		description += "Name\tSigma\r\n";
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				FieldOps.replaceExtremes(data, lowerperc, upperperc, def);
				description += "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(data, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, data);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(data), true, true);
			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	public static void doCutoffDirPerc(String dir, String keyword, double[] lowerperc, double[] upperperc, String suboutdir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name, shortname;
		String description = "";
		description += "Name\tUpper Perc.\tLower Perc.\r\n";
		int n = 0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				FieldOps.cutOffExtremes(data, lowerperc[n], upperperc[n]);
				description += shortname + "\t" + upperperc[n] + "\t" + lowerperc[n] + "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(data, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, data);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(data), true, true);
				n++;
			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	public static void doReplaceDirPerc(String dir, String keyword, double[] lowerperc, double[] upperperc, String suboutdir, double def)
	{
		File[] files = new File(dir).listFiles();
		double[][] data;
		String name, shortname;
		String description = "";
		description += "Name\tSigma\tUpper Perc.\tLower Perc.\r\n";
		description += "Name\tSigma\r\n";
		int n =0;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				FieldOps.replaceExtremes(data, lowerperc[n], upperperc[n], def);
				description += shortname + "\t" + upperperc[n] + "\t" + lowerperc[n] + "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(data, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, data);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(data), true, true);
				n++;
			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	public static void doReplaceDir(String dir, String keyword, double gammaup, double gammadown, String suboutdir, String suboutdirmin)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		
		String description = "Gamma_up = " + gammaup + "\tGamma_down = " + gammadown + "\r\n";
		description += "Name\tSigma\r\n";
		double sigma;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				sigma = FieldOps.sigma(data);
				description += shortname + "\t" + sigma;
				smooth = ColumnIO.readSquareTable(avgdir + shortname + "locavg.dat");
				answer = FieldOps.replaceExtremes(data, smooth, gammaup, gammadown, sigma);
				description += "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);
				
				FieldOps.minusEquals(answer, smooth);
				if (!new File(dir + suboutdirmin).exists())
					new File(dir + suboutdirmin).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdirmin + shortname + ".dat");
				SRAW.writeImage(dir + suboutdirmin + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdirmin + shortname + "fft", FFTOps.obtainFFT(answer), true, true);

			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	//The following replaces extreme values by values equal to locavg[i][j]-gammarep*sigma;
	public static void doReplaceDir(String dir, String keyword, double gammaup, double gammadown, String suboutdir, double gammarep)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		
		String description = "Gamma_up = " + gammaup + "\tGamma_down = " + gammadown + "\r\n";
		description += "Name\tSigma\r\n";
		double sigma;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				sigma = FieldOps.sigma(data);
				description += shortname + "\t" + sigma;
				answer = FieldOps.replaceExtremes(data, gammaup, gammadown, sigma, gammarep);
				description += "\r\n";
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);
			}			
		}
		ColumnIO.writeString(description, dir + suboutdir + "cutoff_log.txt");
	}
	public static void doubleDir(String dir, String keyword, String suboutdir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				answer = FieldOps.expand(data);
//				answer = FieldOps.expandNoSmooth(data, 2);
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + ".dat");
				SRAW.writeImage(dir + suboutdir + shortname, answer);
//				SRAW.writeImage(dir + suboutdir + shortname + "fft", FieldOps.getLog(FFTOps.obtainFFTmagCent(answer)), -8, 8);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "fft", FFTOps.obtainFFT(answer), true, true);
			}			
		}
	}
	public static void double4Dir(String dir, String keyword, String suboutdir)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				answer = FieldOps.expand4(data);
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + "4.dat");
				SRAW.writeImage(dir + suboutdir + shortname + "4", answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "4fft", FFTOps.obtainFFT(answer), true, true);
			}			
		}
	}
	public static void applyLinTransDir(String dir, String keyword, String suboutdir, String transformlogfile)
	{
		File[] files = new File(dir).listFiles();
		double[][] data, smooth, answer;
		String name, shortname;
		String avgdir = dir + "locavg\\";
		
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && !files[i].toString().contains("locavg"))
			{
				name = files[i].toString().substring(0, files[i].toString().length() - 4);
				shortname = name.substring(name.lastIndexOf("\\")+1);
				data = ColumnIO.readSquareTable(files[i].toString());
				answer = new double [data.length][data[0].length];
				applyLinearTransformationTo(data, transformlogfile, 16, answer);
				if (!new File(dir + suboutdir).exists())
					new File(dir + suboutdir).mkdir();
				ColumnIO.writeBin(answer, dir + suboutdir + shortname + "trans.dat");
				SRAW.writeImage(dir + suboutdir + shortname + "trans", answer);
				FFTOps.writeFFTBMPCent(dir + suboutdir + shortname + "transfft", FFTOps.obtainFFT(answer), true, true);
			}			
		}
	}
	
	public static void applyLinearTransformationTo(double[][] data, String transformLogFile, int detail, double[][] target)
	{
		String trans = ColumnIO.readString(transformLogFile);
		String[] lines = trans.split("\r\n");
		double a = 0, b = 0, c = 0;
		for (int i = 0; i < lines.length; i++)
		{
			if (lines[i].startsWith("a = ")) a = Double.parseDouble(lines[i].substring(4));
			if (lines[i].startsWith("b = ")) b = Double.parseDouble(lines[i].substring(4));
			if (lines[i].startsWith("c = ")) c = Double.parseDouble(lines[i].substring(4));
		}
		
		double[][] mat = {{a, c}, {c, b}};
//		double[][] data = ColumnIO.readSquareTable(dir + filename + ".dat");
//		double[][] target = new double [data.length][data[0].length];
		FieldOps.applyLinTransNumSec(data, target, mat, new int [] {data.length/2, data[0].length/2}, detail, detail);
//		ColumnIO.writeBin(target, dir + name + "trans.dat");
//		SRAW.writeImage(dir + name + "trans", target);
//		FFTOps.writeFFTBMPCent(dir + name + "transfft", FFTOps.obtainFFT(target), true, true);
	}
	public static class RadiallyAveragedStripCut
	{
		StripCut s;
		int nphipts;
		double dphi;
		double[] average = null;
		double[][] allresults = null;
		public RadiallyAveragedStripCut(StripCut s, int nphipts)
		{
			this.s = s;
			this.nphipts = nphipts;
			dphi = 2*Math.PI/(nphipts+1);
			allresults = new double [nphipts][s.result.length]; 
			average = new double [s.result.length];
		}
		
		public void calculate()
		{
			for (int i = 0; i < nphipts; i++)
			{
				s.setPhi(i*dphi);
				s.makeCut();
				for (int j = 0; j < s.result.length; j++)
				{
					allresults[i][j] = s.result[j];
					average[j] += allresults[i][j];
				}
			}
			for (int j = 0; j < s.result.length; j++)
			{
				average[j] /= nphipts;
//				average[j] *= j;
			}
		}
		
		public static class ReturnClass
		{
			public double[] average;
			public double[][] allresults;
		}
	}
	
	public static void writeDatAndBMP(double[][] data, String fullname)
	{
		ColumnIO.writeBin(data, fullname + ".dat");
		SRAW.writeImage(fullname, data);
	}
	
	//This will have be garaunteed to have values between 1 and 0:
	//centers is of length [2][n].
	public static double[][] getSumOfGaussianRings(double[][] centers, double radius, double thickness, int Nx, int Ny)
	{
		double[][] ans = new double [Nx][Ny];
		double r;
		double g;
		System.out.print("Doing (out of " + Nx + "): ");
		for (int i = 0; i < Nx; i++){ System.out.println(" " + (i+1));
			for (int j = 0; j < Ny; j++){
				for (int k = 0; k < centers[0].length; k++)
				{
					r = Math.sqrt(Math.pow(i-centers[0][k], 2) + Math.pow(j-centers[1][k], 2));
					g = Math.exp(-(r-radius)*(r-radius)/(thickness*thickness));
					ans[i][j] += g;
				}
				ans[i][j] = Math.min(ans[i][j], 1);
			}
		}
		System.out.println();
		return ans;
	}
	public static void getSumOfGaussianRingsQuick(double[][] centers, double radius, double thickness, double[][] ans)
	{
		int Nx = ans.length, Ny = ans[0].length;
		double r;
		double g;
		double rangemax = radius + 4*thickness;
		System.out.print("Doing (out of " + Nx + "): ");
		for (int i = 0; i < Nx; i++){ if (i % 100 == 0) System.out.println(" " + (i+1));
			for (int j = 0; j < Ny; j++){ ans[i][j] = 0;
				for (int k = 0; k < centers[0].length; k++)
				{
					if (Math.abs(i - centers[0][k]) < rangemax && Math.abs(j - centers[1][k]) < rangemax)
					{
						r = Math.sqrt(Math.pow(i-centers[0][k], 2) + Math.pow(j-centers[1][k], 2));
						g = Math.exp(-(r-radius)*(r-radius)/(thickness*thickness));
						ans[i][j] += g;
					}
				}
				ans[i][j] = Math.min(ans[i][j], 1);
			}
		}
		System.out.println();
	}
	public static double[][] getSumOfGaussianRingsQuick(double[][] centers, double radius, double thickness, int Nx, int Ny)
	{
		double[][] ans = new double [Nx][Ny];
		getSumOfGaussianRingsQuick(centers, radius, thickness, ans);
		return ans;
	}
	public static void printPath2Slash(String inpath)
	{
		String out = "";
		String[] chars = inpath.split("");
		for (int i = 0; i < chars.length; i++)
		{
			if (chars[i].equals("\\"))
				chars[i] += "\\";
			out += chars[i];
		}
		System.out.println(out);
	}
	
	public static void writeRotations(String dir, int npts)
	{
		double theta = 0, dth = 2*Math.PI/(npts-1);
//		double[][] source = FileOps.openBin(dir);
		double[][] source = Mask.rectMaskCent(16, 16, 8, 8);
		int nx = source.length, ny = source[0].length;
		double[] origin = new double []{nx/2, ny/2};
		double[][] target = new double [nx][ny];
		double[][] matrix = new double [2][2];
		for (int i = 0; i < npts; i++)
		{
			theta = 0 + dth*i;
			Matrix.putRotationMatrix(theta, matrix);
			FieldOps.applyLinearTransformation(source, matrix, origin, target);
//			ColumnIO.writeBin(target, dir + "rot" + i + ".dat");
			SRAW.writeImage(MovieMaker.avidir + "rot" + MovieMaker.fromInt(i), target);
		}
	}
	public static void doRotationCorrelations(String dir, int npts)
	{
		double theta = 0, dth = (Math.PI/2)/(npts-1);
		double[][] data = FileOps.openBin(dir);
		double[][] source = Mask.rectMaskCent(16, 16, 8, 8);
		int nx = source.length, ny = source[0].length;
		double[] origin = new double []{nx/2, ny/2};
		double[][] target = new double [nx][ny];
		double[][] matrix = new double [2][2];
		double[] c = new double [npts];
		for (int i = 0; i < npts; i++)
		{
			theta = 0 + dth*i;
			Matrix.putRotationMatrix(theta, matrix);
			FieldOps.applyLinearTransformation(source, matrix, origin, target);
			c[i] = FieldOps.correlation(target, data);
			System.out.println(theta + "\t" + c[i]);
			//			ColumnIO.writeBin(target, dir + "rot" + i + ".dat");
//			SRAW.writeImage(MovieMaker.avidir + "rot" + MovieMaker.fromInt(i), target);
		}
	}
	
	public static void readWrite(String dir)
	{
		double[][] data = FileOps.openBin(dir);
		SRAW.writeImage(dir + "pic", data);
		
		ColumnIO.writeBin(data, dir + "pic.dat");
	}
	public static void flipX(String dir)
	{
		double[][] data = FileOps.openBin(dir);
		ArrayOps.flipX(data);
		FileOps.writeTableBIN(data);
	}
	public static void doLineSubtract(String dir, boolean truei)
	{
		JFileChooser fc = new JFileChooser(dir);
		double[][] data = FileOps.openTable(fc);
//		FieldOps.subtactLineAvg(data, truei);
//		double[][] ans = FieldOps.normalizeEachRow(data, false);
//		RHKFileOps.write3Files(fc, data);
//		FileOps.writeTableASCII(fc, ans);
	}
	
	public static void doStuffWithImps(Topomap t, ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps, AtomicCoordinatesSet latt)
	{
		
//		ArrayList<GaussSquareImpurityAngleSuite.Impurity> posImps = GaussSquareImpurityAngleSuite.getImpurities(FileOps.selectOpen(null));
//	
//		for (int i = 0; i < imps.size(); i++)
//		{
//			imps.get(i).position = posImps.get(i).position;
//			imps.get(i).latticePos = latt.getAtomicCoords(imps.get(i).position);
//		}
//		FileOps.writeString(GaussSquareImpurityAngleSuite.getImpText(imps));
		double x, y;
		int xm, ym;
		double[][] pixelAt = new double [4][2];
		double xp, yp;
//		int[] badI = new int[] {33, 36, 49, 53, 54, 55};
		int[] badI = new int[] {36, 37, 49, 51, 53, 55};
		ArrayList<GaussSquareImpurityAngleSuite.Impurity> goodImps = ArrayOps.copyExcept(imps, badI);
		
		for (int i = 0; i < goodImps.size(); i++)
			if (goodImps.get(i).latticePos[0] == 0) goodImps.get(i).latticePos = latt.getAtomicCoords(goodImps.get(i).position);
		
		
		int[] isL, isR, isC;
		isR = FieldUtil.indexesBetweenCutoffs(imps, badI, 31, Double.MAX_VALUE);
		isL = FieldUtil.indexesBetweenCutoffs(imps, badI, -Double.MAX_VALUE, 16);
		isC = FieldUtil.indexesBetweenCutoffs(imps, badI, 16, 31);
		ArrayList<GaussSquareImpurityAngleSuite.Impurity> rightImps, leftImps, centImps;
		rightImps = ArrayOps.copyOnly(imps, isR);
		leftImps = ArrayOps.copyOnly(imps, isL);
		centImps = ArrayOps.copyOnly(imps, isC);
//		
		double[] angR = new double [rightImps.size()];
		double[] angL = new double [leftImps.size()];
		double[] angC = new double [centImps.size()];
		for (int i = 0; i < leftImps.size(); i++)
			angL[i] = leftImps.get(i).angle;
		for (int i = 0; i < rightImps.size(); i++)
			angR[i] = rightImps.get(i).angle;
		for (int i = 0; i < centImps.size(); i++)
			angC[i] = centImps.get(i).angle;
//		
		double meanR = ArrayOps.mean(angR), meanL = ArrayOps.mean(angL), meanC = ArrayOps.mean(angC);
		double sigR = ArrayOps.sigma(angR), sigL = ArrayOps.sigma(angL), sigC = ArrayOps.sigma(angC);
		double braggA = Math.toDegrees(latt.getAngleBetween(new double[] {1, 0}, 0));
		double diffR = meanR - braggA;
		double diffL = braggA - meanL;
		
		double meanD = (diffR + diffL)/2;
		System.out.println(braggA);
		System.out.println(meanR + "\t" + sigR);
		System.out.println(meanL + "\t" + sigL);
		System.out.println(meanC + "\t" + sigC);
		System.out.println(diffR);
		System.out.println(diffL);
		System.out.println(meanD);
		System.out.println(braggA + "\t" + angL.length + "\t" + meanL + "\t" + angR.length + "\t" + meanR + "\t" + angC.length + "\t" + diffR + "\t" + diffL);
		
		System.out.println(Printer.arrayVertical(angL));
		System.out.println(Printer.arrayVertical(angR));
		
				double[] lattO = latt.getOrigin();
//		double[] trans = latt.getAtomicCoords(new double[] {-1 + lattO[0], lattO[1]});
//		goodImps = ArrayOps.copyOnly(imps, is);
		
//		goodImps = ArrayOps.copyExcept(imps, badI);
//		for (int i = 0; i < goodImps.size(); i++)
//			System.out.println(i + "\t" + goodImps.get(i).angle);
		
		int sub;
//		boolean[] subA = new boolean [goodImps.size()];
//		double[][] emergenceIs = FileOps.openTable("");
//		double[] index = emergenceIs[1];
//		double[] diffGood = ArrayOps.copyExcept(index, badI);
//		double[][] ds = new double[imps.size()][imps.size()]; //the inter-impurity distances
//		double avg = 0;
		//The index of the positive Fermi energy is 198.
		int indexMin = 198;
		double[] sum = new double [imps.size()];
//		double[][][] specOut = new double [imps.size()][5][t.nlayers];
//		double[][] averageSpec = new double [imps.size()][t.nlayers];
////		for (int n = 1; n <= 5; n++)
//		for (int i = 0; i < leftImps.size(); i++)
//			System.out.println(i + "\t" + isL[i] + "\t" + leftImps.get(i).angle + "\t" + leftImps.get(i).maxCorr + "\t" + (angL[i] - meanL) + "\t" + sum[isL[i]] + "\t" + (Math.abs(angL[i] - braggA) - meanD));
//		for (int i = 0; i < rightImps.size(); i++)
//			System.out.println(i + "\t" + isR[i] + "\t" + rightImps.get(i).angle + "\t" + rightImps.get(i).maxCorr + "\t" + (angR[i] - meanR) + "\t" + sum[isR[i]] + "\t" + (Math.abs(angR[i] - braggA) - meanD));
//		for (int i = 0; i < centImps.size(); i++)
//			System.out.println(i + "\t" + isC[i] + "\t" + centImps.get(i).angle + "\t" + centImps.get(i).maxCorr + "\t" + (angC[i] - meanC) + "\t" + sum[isC[i]] + "\t" + (Math.abs(angC[i] - braggA) - meanD));
		
//		double[] sumGoodOnly = ArrayOps.copyExcept(sum, badI);
//		double[][] averageSpecGoodOnly = ArrayOps.copyExcept(averageSpec, badI);
//		double[][] averageBin = new double [8][t.nlayers];
//		
//		ArrayList<Integer> binIs = new ArrayList<Integer>();
//		double[] binMins = ArrayOps.generateArrayNotIncl(200, 1000, 8);
//		for (int i = 0; i < 8; i++)
//		{
//			binIs = new ArrayList<Integer>();
//			for (int j = 0; j < goodImps.size(); j++)
//				if (sumGoodOnly[j] > binMins[i] && sumGoodOnly[j] < binMins[i] + 100)
//					binIs.add(j);
//			for (int k = 0; k < t.nlayers; k++){
//				for (int j = 0; j < binIs.size(); j++)
//				{
//					averageBin[i][k] += averageSpecGoodOnly[binIs.get(j)][k];
//				}
//				averageBin[i][k] /= binIs.size();
//			}
//			//Finally, print them
//			for (int j = 0; j < t.nlayers; j++)
//			{
//				System.out.print(i + "\t" + t.v[j] + "\t");
//				for (int k = 0; k < binIs.size(); k++)
//					System.out.print(averageSpecGoodOnly[binIs.get(k)][j] + "\t");
//				System.out.println();
//			}
//		}
//		GraphDrawerCart g;
//		double[][][] plotData = new double [4][2][t.nlayers];
//		for (int i = 0; i < t.nlayers; i++)
//		{
//			plotData[0][0][i] = t.v[i];
//			plotData[1][0][i] = t.v[i];
//			plotData[2][0][i] = t.v[i];
//			plotData[3][0][i] = t.v[i];
//		}
//		for (int i = 0; i < imps.size(); i++)
//		{
////			g = new GraphDrawerCart("", t.v, averageSpec[i]);
//////			g.repaint();
////			g.showWindow();
//			for (int k = 0; k < 4; k++)
//				plotData[k][1] = specOut[i][k+1];
////			SRAW.writeImage(Topomap.stddir + "impcornerplots\\imp" + MovieMaker.fromInt(i), GraphDrawerCart.drawPlot(t.v, averageSpec[i]));
//			SRAW.writeImage(Topomap.stddir + "impcornerplots\\imp" + MovieMaker.fromInt(i), GraphDrawerCart.drawPlot(plotData));
////			g.setVisible(false);
//		}
//		System.exit(0);
//		SRAW.writeImage(fc.getCurrentDirectory().toString() + "\\graph" + MovieMaker.fromInt(currentImp), (BufferedImage)spec.dbimage);
//		for (int i = 0; i < imps.size(); i++)
//		{
//			sum[i] = 0;
//			x = imps.get(i).latticePos[0];
//			y = imps.get(i).latticePos[1];
//			xm = FieldOps.roundDown(x);
//			ym = FieldOps.roundDown(y);
//			pixelAt[0] = latt.getPixelCoords(new double [] {xm, ym});
//			pixelAt[1] = latt.getPixelCoords(new double [] {xm+1, ym});
//			pixelAt[2] = latt.getPixelCoords(new double [] {xm, ym+1});
//			pixelAt[3] = latt.getPixelCoords(new double [] {xm+1, ym+1});
//			
//			for (int j = 0; j < t.nlayers; j++)
//			{
//				specOut[i][0][j] = t.v[j];
//				for (int k = 0; k < 4; k++)
//					specOut[i][k+1][j] = FieldOps.getValueAt(t.data[j], pixelAt[k][0]+1, pixelAt[k][1], 0);  //1 pixel horizontal shift
//			}
//			for (int j = 0; j < t.nlayers; j++)
//				averageSpec[i][j] = (specOut[i][1][j] + specOut[i][2][j] + specOut[i][3][j] + specOut[i][4][j])/4;
//			for (int j = indexMin; j < t.nlayers; j++)
//			sum[i] += averageSpec[i][j];
//			
//			ColumnIO.writeTable(specOut, Topomap.stddir + "lists\\imp corner spectra\\imp" + MovieMaker.fromInt(i) + ".txt");
//			
//			for (int j = 0; j < imps.size(); j++)
//			{
//				ds[i][j] = Distance.distance(imps.get(i).position[0] - imps.get(j).position[0], imps.get(i).position[1] - imps.get(j).position[1]);
//			}
//			ArrayOps.quicksort(ds[i]);
//			avg = 0;
//			for (int j = 1; j <= n; j++)
//			{
//				avg += ds[i][j];
//			}
//			avg /= n;
//			
//				for (int j = 0; j < imps.size(); j++)
//				{
//					System.out.println(i + "\t" + ds[i][j]);
//				}
//				x = Ops.modulo(goodImps.get(i).latticePos[0], 1) + trans[0];
//				y = Ops.modulo(goodImps.get(i).latticePos[1], 1) + trans[1];
//				x = goodImps.get(i).latticePos[0];
//				y = goodImps.get(i).latticePos[1];
//				xm = FieldOps.roundDown(x);
//				ym = FieldOps.roundDown(y);
//				subA[i] = (xm+ym)%2 == 0;
//				if (!ArrayOps.contains(badI, i))
//			System.out.println(t.v[(int)index[i]] + "\t" + sum[i]);
//						System.out.println(goodImps.get(i).angle + "\t" + (subA[i] ? 1 : 0));
//	
//						System.out.println(x + "\t" + y);
//		}

		
			
		
//		FileOps.writeString(GaussSquareImpurityAngleSuite.getImpText(imps));
		
		double[] max = new double [imps.size()];
		for (int i = 0; i < imps.size(); i++)
			max[i] = imps.get(i).angle - braggA;
////		
////		int[] badI = new int[] {33, 36, 49, 53, 54, 55};
//		
//		
//		
		double[] maxCopy = ArrayOps.copyExcept(max, badI);
//		double[][] anglePlot;
////		int[] is = FieldUtil.indexesBetweenCutoffs(imps, badI, 31, Double.MAX_VALUE);
////		Printer.printlnHorizontal(is);
////		anglePlot = getAnglePlot(imps, is);
////		FileOps.writeBINandImage(Topomap.stddir + "anglePlot\\" + "right", anglePlot);
////		
////		is = FieldUtil.indexesBetweenCutoffs(imps, badI, -Double.MAX_VALUE, 15.5);
////		anglePlot = getAnglePlot(imps, is);
////		FileOps.writeBINandImage(Topomap.stddir + "anglePlot\\" + "left", anglePlot);
////		Printer.printlnHorizontal(is);
////		is = FieldUtil.indexesBetweenCutoffs(imps, badI, 15.5, 31);
////		anglePlot = getSpecialPlot(goodImps, null, diffGood);
////		FileOps.writeBINandImage(Topomap.stddir + "anglePlot\\" + "difference", anglePlot);
////		Printer.printlnHorizontal(is);
//				
		double[] binBottoms = ArrayOps.generateArrayNotInclLower(-20, 20, 50);
		int[] hist = ArrayOps.getHistogram(maxCopy, binBottoms);
//		
		for (int i = 0; i < binBottoms.length; i++)
			System.out.println(binBottoms[i] + "\t" + hist[i]);
//		Printer.printlnVertical(hist);
	}
	public static double[][] applyUField(double[][][] u, double[][] topo)
	{
		double[][] after = new double[topo.length][topo[0].length];
		FieldOps.applyUFieldSmooth(topo, u, 2, 2, after);
		FieldOps.changeZeroToAverage(after);
		return after;
	}
	
	public static boolean[] splitImpsWRTRoot2(Layer t, ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps, AtomicCoordinatesSet latt)
	{
		for (int i = 0; i < imps.size(); i++)
			if (imps.get(i).latticePos[0] == 0) imps.get(i).latticePos = latt.getAtomicCoords(imps.get(i).position);

		double x, y;
		int xm, ym;
		double[][] pixelAt = new double [4][2];
		boolean[] subA = new boolean[imps.size()];
		String[] output = new String[] {"A", "B"};
		for (int i = 0; i < imps.size(); i++)
		{
			x = imps.get(i).latticePos[0];
			y = imps.get(i).latticePos[1];
			xm = FieldOps.roundDown(x);
			ym = FieldOps.roundDown(y);
			subA[i] = (xm+ym)%2 == 0;
		}
		
		return subA;

	}

	public static double[][] getAnglePlot(ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps, int[] is)
	{
		int N = 256;
		double[] angle = new double [imps.size()];
		for (int i = 0; i < imps.size(); i++)
			angle[i] = imps.get(i).angle;
		
		if (is == null)
		{
			is = new int [imps.size()];
			for (int i = 0; i < imps.size(); i++)
				is[i] = i;
		}
	
		double mean = ArrayOps.mean(angle, is);
		double[][] plot = new double [N][N];
		double[] ds = new double [is.length];
		double dAngle;
		double gauss;
		int nearestImpIndex;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < is.length; k++)
					ds[k] = Distance.distance(i-imps.get(is[k]).position[0], j-imps.get(is[k]).position[1]);
				nearestImpIndex = ArrayOps.minIndex(ds);
				dAngle = imps.get(is[nearestImpIndex]).angle - mean;
				gauss = ds[nearestImpIndex] < 6 ? 1 : 0;
				plot[i][j] = mean + dAngle*gauss;
			}
		return plot;
	}
	public static double[][] getSpecialPlot(ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps, int[] is, double[] special)
	{
		int N = 256;
		
		if (is == null)
		{
			is = new int [imps.size()];
			for (int i = 0; i < imps.size(); i++)
				is[i] = i;
		}
	
		double mean = ArrayOps.mean(special, is);
		System.out.println(mean);
		double[][] plot = new double [N][N];
		double[] ds = new double [is.length];
		double dAngle;
		double gauss;
		int nearestImpIndex;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < is.length; k++)
					ds[k] = Distance.distance(i-imps.get(is[k]).position[0], j-imps.get(is[k]).position[1]);
				nearestImpIndex = ArrayOps.minIndex(ds);
				dAngle = special[is[nearestImpIndex]] - mean;
				gauss = ds[nearestImpIndex] < 6 ? 1 : 0;
				plot[i][j] = mean + dAngle*gauss;
			}
		return plot;
	}
	public static int[] indexesBetweenCutoffs(ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps, int[] excludedI, double lowCut, double highCut)
	{
		ArrayList<Integer> indexes = new ArrayList<Integer>();
		for (int i = 0; i < imps.size(); i++)
		{
			if (imps.get(i).angle > lowCut && imps.get(i).angle < highCut && !ArrayOps.contains(excludedI, i))
				indexes.add(i);
		}
		int[] ans = new int [indexes.size()];
		for (int i = 0; i < indexes.size(); i++)
			ans[i] = indexes.get(i);
		return ans;
	}
	public static double[][][] shiftUOnePixel(double[][][] u)
	{
		//this shifts the u field by +1, 0;
		double[][][] us = new double [u.length][u[0].length][2];
		for (int i = 1; i < u.length; i++)
			for (int j = 0; j < u[0].length; j++)
			{
				us[i][j][0] = u[i-1][j][0];
				us[i][j][1] = u[i-1][j][1];
			}
		
		//the edge of the shifted field is blank by default, we will just copy the first row in u
		for (int j = 0; j < u[0].length; j++){
			us[0][j][0] = u[0][j][0];
			us[0][j][1] = u[0][j][1];
		}
		return us;
	}
	
	
	//This method is for analyzing the Iridates topograph in run 234. It loads a list of impurities, etc....
	public static void method234(double[][] topo)
	{
		ArrayList<GaussSquareImpurityAngleSuite.Impurity> imps = GaussSquareImpurityAngleSuite.getImpurities(FileOps.selectOpen(null));
		
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText());
		double x, y;
		int xm, ym;
		double sum0, sum1;
		double[][] zeroCoord = new double [4][2];
		double[][] oneCoord = new double [4][2];
		int[] chirality = new int[imps.size()];
		int[] sublatt = new int [imps.size()];
		for (int i = 0; i < imps.size(); i++){
			imps.get(i).latticePos = latt.getAtomicCoords(imps.get(i).position);
			x = imps.get(i).position[0];
			y = imps.get(i).position[1];
//			System.out.println(Ops.modulo(imps.get(i).latticePos[0], 1) + "\t" + Ops.modulo(imps.get(i).latticePos[1], 1));
			xm = FieldOps.roundDown(imps.get(i).latticePos[0]);
			ym = FieldOps.roundDown(imps.get(i).latticePos[1]);			
			zeroCoord[0] = latt.getPixelCoords(xm+2, ym);
			zeroCoord[1] = latt.getPixelCoords(xm+1, ym+2);
			zeroCoord[2] = latt.getPixelCoords(xm-1, ym+1);
			zeroCoord[3] = latt.getPixelCoords(xm, ym-1);
			oneCoord[0] = latt.getPixelCoords(xm+2, ym+1);
			oneCoord[1] = latt.getPixelCoords(xm, ym+2);
			oneCoord[2] = latt.getPixelCoords(xm-1, ym);
			oneCoord[3] = latt.getPixelCoords(xm+1, ym-1);
		
			sum0 = 0; sum1 = 0;
			for (int j = 0; j < 4; j++)
			{
				sum0 += FieldOps.getValueAt(topo, zeroCoord[j][0], zeroCoord[j][1], 0);
				sum1 += FieldOps.getValueAt(topo, oneCoord[j][0], oneCoord[j][1], 0);
			}
			chirality[i] = sum0 > sum1 ? 0 : 1;
			sublatt[i] = (xm+ym)%2 == 0 ? 0 : 1;
			if (chirality[i] == 0)
				SRAW.writeImage("C:\\data\\analysis\\SrIrO 327\\234\\imps\\left\\imp" + MovieMaker.fromInt(i), FieldOps.expandBi(topo, 4, FieldOps.round(x-16), 32, FieldOps.round(y-16), 32));
			else
				SRAW.writeImage("C:\\data\\analysis\\SrIrO 327\\234\\imps\\right\\imp" + MovieMaker.fromInt(i), FieldOps.expandBi(topo, 4, FieldOps.round(x-16), 32, FieldOps.round(y-16), 32));
		}
		for (int i = 0; i < imps.size(); i++)
			System.out.println(chirality[i] + "\t" + sublatt[i]);
	}
	
	public static Layer evaluateLayerAtTopomapPoints(Topomap t, Layer y)
	{
		double[][] data = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				data[i][j] = y.evaluateAtMetric(t.x[i], t.y[j]);
			}
		return new Layer (data, t.x, t.y, y.v, y.current);
	}
	public static double[][] applyUFieldLayer(double[][][] u, Layer t)
	{
		double[][] after = new double[t.nx][t.ny];
			FieldOps.applyUFieldSmooth(t.data, u, 2, 2, after);
			FieldOps.changeZeroToAverage(after);
		return after;
	}
	
	public static Map_4D loadMEMOutput_DAT(double a1, double da1, double a2, double da2, String dir, boolean DAT_true_else_SPT)
	{
		String[] pieceNames;
		
		
		if (DAT_true_else_SPT) pieceNames = new String[] {"Re_Sig_Initial","Error Bars","Re_Sig_Calc","Im_Sig_Calc"};
		else pieceNames = new String[] {"Alpha^2 F", "Model A^2F"}; 
		
		
		int npieces = pieceNames.length;
		
		String keyword = DAT_true_else_SPT ? "DAT" : "SPT";
		
		File[] f = new File(dir).listFiles();
		String fstr;
		ArrayList<File> selected = new ArrayList<File>();
		String name;
		String[] elements;
		int imax = 0, jmax = 0;
		int ith, jth;
		double[][] filedata;
		int nlayers = 0;
		for (int n = 0; n < f.length; n++)
		{
			fstr = f[n].toString();
			name = fstr.substring(fstr.lastIndexOf("\\")+1);
			if (name.endsWith(".dat") && name.contains(keyword) && !name.contains("MIN_CHI"))
			{
				selected.add(f[n]);
				elements = name.split("_");
				ith = Integer.parseInt(elements[2]);
				jth = Integer.parseInt(elements[3].substring(0, elements[3].indexOf(".")));
				imax = Math.max(ith, imax);
				jmax = Math.max(jth, jmax);
				
				filedata = ColumnIO.readAllColumns_Special_Alert(f[n], null);
				nlayers = filedata[0].length;
				System.out.print(" " + n);
			}
		}
		int nx = imax+1, ny = jmax+1;
		double[][][][] dataset = new double [npieces][nlayers][nx][ny];
		
		double[] v = null;
		System.out.println();
		for (int n = 0; n < selected.size(); n++)
		{
			fstr = selected.get(n).toString();
			name = fstr.substring(fstr.lastIndexOf("\\")+1);
			elements = name.split("_");
			ith = Integer.parseInt(elements[2]);
			jth = Integer.parseInt(elements[3].substring(0, elements[3].indexOf(".")));

			filedata = ColumnIO.readAllColumns_Special_Alert(selected.get(n), null);
			for (int m = 0; m < npieces; m++)
				for (int p = 0; p < nlayers; p++)
					dataset[m][p][ith][jth] = filedata[m+1][p];
			v = filedata[0];
			System.out.print(" " + n);
		}
		double[] x, y;
		x = ArrayOps.generateArray(a1, da1, nx);
		y = ArrayOps.generateArray(a2, da2, ny);
		return new Map_4D(dataset, v, x, y, null, pieceNames);
	}
	
	public static void createTopomapSpectralFunction()
	{
		int nx = 512, ny = 512, nv = 128;
		
		double[] kx = ArrayOps.generateArrayInclBoth(-1, 1, nx);
		double[] ky = ArrayOps.generateArrayInclBoth(-1, 1, ny);
		double[] v = ArrayOps.generateArrayInclBoth(-0.5, 0, nv);
		
		double[] re = ArrayOps.generateArrayInclBoth(0, 0, nv);
		double[] im = ArrayOps.generateArrayInclBoth(0.003, 0.003, nv);
		
		Topomap t = TopomapUtil.getSpectralFunction(v, kx, ky, re, im);
		Topomap.writeBIN(t);
	}
	
	public static Layer loadLayer(double a1, double da1, double a2, double da2, String dir)
	{
		double[][] data = FileOps.openTable(dir);
		double[] x, y;
		int nx = data.length, ny = data[0].length;
		x = ArrayOps.generateArray(a1, da1, nx);
		y = ArrayOps.generateArray(a2, da2, ny);
		return new Layer(data, x, y, 1, 1);
	}

}
