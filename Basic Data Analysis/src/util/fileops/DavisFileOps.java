package util.fileops;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Date;

import javax.swing.JFileChooser;

import drawing.LayerViewer;
import drawing.SpectraDrawer;
import drawing.TopomapViewer;
import util.ArrayOps;
import util.Printer;
import util.TopomapUtil;

public class DavisFileOps {
	
	/**
	 * Lockin array is 
	 * [0] - freq
	 * 1 - amp
	 * 2 - phase
	 * 3 - sens
	 * 4 - tconst
	 * 5 - rollof
	 * 6 - reserve
	 * 7 - filters
	 * 8 - harmonic
	 * 9 - expand
	 */
	
	public static final int DIDV_MAP = 0;
	public static final int CURRENT_MAP = 1;
	public static final int TOPO_MAP = 2;
	public static final int FDBK_MAP = 3;
	
	public static String davisDir = "E:\\Nb_Bi2Se3\\";
	
	/**
	 * 
	 * @param filepath
	 * @return
	 */
	public static DavisLayerInfo loadFromFile(String filepath)
	{
		ByteBuffer bb = getFromFile(filepath);
		//Get the date
		int rOffset = 256;
		int start = 24;
		int len = 20;
		char[] desc = new char [len];
		for (int i = 0; i < len; i++)
			desc[i] = bb.getChar(rOffset+start-1+i);
		String descs = new String(desc);
		System.out.println(descs);
		
		start = 179;
		double wFact = bb.getFloat(rOffset+start-1);
		System.out.println("WFact: " + wFact);
		start = 183;
		double wZero = bb.getFloat(rOffset+start-1);
		System.out.println("WZero: " + wZero);
		
		start = 151;
		int nx = bb.getInt(rOffset+start-1);
		System.out.println(nx);
		
		start = 155;
		int ny = bb.getInt(rOffset+start-1);
		System.out.println(ny);
		
		start = 225;
		int nz = bb.getInt(rOffset+start-1);
		System.out.println(nz);
		//nz should equal one for a topography
		
		String unit = "";
		start = 195;
		for (int i = 0; i < 10; i++)
			unit += bb.getChar(rOffset+start-1+i);
		System.out.println(unit);

		String xyunit = "";
		start = 205;
		for (int i = 0; i < 10; i++)
			xyunit += bb.getChar(rOffset+start-1+i);
		System.out.println("xyunit: " + xyunit);
		boolean useNM = xyunit.contains("m");

		double xmin;
		start = 163;
		xmin = bb.getFloat(rOffset+start-1);
		System.out.println(xmin);
		double xmax;
		start = 167;
		xmax = bb.getFloat(rOffset+start-1);
		System.out.println(xmax);
		
		double[] x = ArrayOps.generateArrayInclBoth(xmin, xmax, nx);
		ArrayOps.multiply(x, useNM ? 1e-9 : 1e-10); //assuming units are angstroms.
		
		double ymin;
		start = 171;
		ymin = bb.getFloat(rOffset+start-1);
		System.out.println(ymin);
		double ymax;
		start = 175;
		ymax = bb.getFloat(rOffset+start-1);
		System.out.println(ymax);
		
		double[] y = ArrayOps.generateArrayInclBoth(xmin, xmax, nx);
		ArrayOps.multiply(y, useNM ? 1e-9 : 1e-10); //assuming units are angstroms.
		
		//topomap parameters;
		rOffset = 1280;
		start = 1;
		double vstart = bb.getFloat(rOffset+start-1);
		System.out.println("vstart = " + vstart);
		start = 5;
		double vstop = bb.getFloat(rOffset+start-1);
		System.out.println("vstop = " + vstop);
		double[] v = (nz == 1 ? new double [] {0} : ArrayOps.generateArrayInclBoth(vstart, vstop, nz));
		
		rOffset = 1356;
		double[] li = new double [10];
		start = 1;
		for (int i = 0; i < 10; i++)
			li[i] = bb.getFloat(rOffset+start+i*4-1);
		
		System.out.println(Printer.arrayLnHorizontal(li));
		
		rOffset = 1384; start = 23;
		double vset = bb.getFloat(rOffset+start-1);
		System.out.println(vset); // in mV;
		vset = vset/1000.0;
		
		start = 27;
		double iset = bb.getFloat(rOffset+start-1);
		System.out.println(iset); // in nanoamps I think
		iset = iset/1e9;
		
		return new DavisLayerInfo(nx, ny, nz, x, y, v, vset, iset, li, wFact, wZero);
	}
	
	public static PointSpectra loadPointSpectrum(String filepath)
	{
		ByteBuffer bb = getFromFile(filepath);
		double[] allData = new double [bb.array().length/4];
		for (int i = 0; i < allData.length; i++)
			allData[i] = bb.getFloat();
		
		double[] vTemp = new double [allData.length/2];
		for (int i = 0; i < vTemp.length; i++)
			vTemp[i] = allData[2*i];

		double[] vProper = new double [vTemp.length-264];
		for (int i = 0; i < vProper.length; i++)
			vProper[i] = vTemp[i+264];
		
		//Now we read the forward and backward.
		double[] rawFB = new double [vTemp.length];
		for (int i = 0; i < vTemp.length; i++)
			rawFB[i] = allData[2*i + 1];
		double[] actualFB = new double [vProper.length];
		for (int i = 0; i < actualFB.length; i++)
			actualFB[i] = rawFB[i+264];
		
		double[] vForward = new double [vProper.length/2];
		double[] dataForward = new double [vProper.length/2];
		double[] dataBackward = new double [vProper.length/2];
		for (int i = 0; i < vForward.length; i++)
		{
			vForward[i] = vProper[i];
			dataForward[i] = actualFB[i];
			dataBackward[i] = actualFB[i+vForward.length];
		}
		
		return new PointSpectra(new double[][] {dataForward, dataBackward},vForward, new double[] {0, 0}, new double[] {0, 0});
			
	}
	
	public static ByteBuffer getFromFile(String filepath)
	{
		File file = new File(filepath);
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer[] t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);
		
		byte[] contents = new byte[(int)file.length()];
		try {
			ind.read(contents);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ByteBuffer bb = ByteBuffer.wrap(contents);
		bb.order(ByteOrder.LITTLE_ENDIAN);
		try {
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return bb;
	}
	
	public static double[][][] loadFromFile3D(String filepath, DavisLayerInfo info, int type)
	{
		int nx = info.nx, ny = info.ny, nz = info.nz;
		double[][][] ans = new double [nz][nx][ny];
		
		ByteBuffer bb = getFromFile(filepath);
		
		int[] a = new int [nx*ny*nz];
		int start = 2112;
		for (int i = 0; i < nx*ny*nz; i++)
		{
			a[i] = bb.getShort(start+i*2);
			a[i] += 65536;
			a[i] %= 65536;
		}
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int k = 0; k < nz; k++)
				{
					int n = i + j*nx + k*nx*ny;
					ans[k][i][j] = a[n]*info.wFact + info.wZero;
				}

//		if (type == DIDV_MAP){
//			for (int i = 0; i < nx; i++)
//				for (int j = 0; j < ny; j++)
//					for (int k = 0; k < nz; k++)
//					{
//						int n = i + j*nx + k*nx*ny;
//						ans[k][i][j] = a[n]*0.000305180 - 10;
//					}
//		}
//		else if (type == CURRENT_MAP){
//			for (int i = 0; i < nx; i++)
//				for (int j = 0; j < ny; j++)
//					for (int k = 0; k < nz; k++)
//					{
//						int n = i + j*nx + k*nx*ny;
//						ans[k][i][j] = -(a[n]*0.000305180 - 10);
//					}
//		}
//		else if (type == TOPO_MAP){
//			for (int i = 0; i < nx; i++)
//				for (int j = 0; j < ny; j++)
//					for (int k = 0; k < nz; k++)
//					{
//						int n = i + j*nx + k*nx*ny;
//						ans[k][i][j] = a[n]*0.000305180 - 10;
//					}
//		}
//		else
//			System.out.println("Type unknown.");
		return ans;
	}
	public static double[][] loadFromFile2D(String filepath, DavisLayerInfo info)
	{
		return loadFromFile3D(filepath, info, TOPO_MAP)[0];
	}
	public static class DavisLayerInfo
	{

		int nx; int ny; int nz;
		double[] x; double[] y; double[] v;
		double vset; double iset;
		double[] li;
		double wFact, wZero;
	
		public DavisLayerInfo(int nx, int ny, int nz, double[] x, double[] y, double[] v, double vset,
				double iset, double[] li, double wFact, double wZero) {
			super();
			this.nx = nx;
			this.ny = ny;
			this.nz = nz;
			this.x = x;
			this.y = y;
			this.v = v;
			this.vset = vset;
			this.iset = iset;
			this.wFact = wFact;
			this.wZero = wZero;
		}
	}
//	private static double readFloat(byte[] content, int start)
//	{
//		ByteBuffer bb = ByteBuffer.allocate(8);
//		bb.order(ByteOrder.BIG_ENDIAN);
//		for (int i = 0; i < 4; i++)
//			bb.put(content[start+i]);
//		return (double) bb.getFloat();
//	}
//	private static double readFloat(byte[] content, int start)
//	{
//		ByteBuffer bb = ByteBuffer.allocate(8);
//		bb.
//		bb.order(ByteOrder.BIG_ENDIAN);
//		for (int i = 0; i < 4; i++)
//			bb.put(content[start+i]);
//		return (double) bb.getFloat();
//	}
	
	public static Topomap getTopomapFromFile(String filepath)
	{
		DavisLayerInfo info = loadFromFile(filepath);
		int type = filepath.endsWith(".1FL") ? DIDV_MAP : CURRENT_MAP;
		double[][][] data = loadFromFile3D(filepath, info, type);
		Topomap t = new Topomap(data, info.v, info.x, info.y, null);
		if (t.v[0] > t.v[t.nlayers-1]) TopomapUtil.flipEnergy(t);
//		for (int i = 0; i < t.nlayers; i++)
//		{
//			ArrayOps.flipX(t.data[i]);
//		}
		return t;
	}
	public static void writeTopomap(String source, String dest)
	{
		Topomap.writeBIN(getTopomapFromFile(source), dest);
	}
	public static void writeLayer(String source, String dest)
	{
		Layer.writeBIN(getLayerFromFile(source), dest);
	}
	public static Layer getLayerFromFile(String filepath)
	{
		DavisLayerInfo info = loadFromFile(filepath);
		double[][] data = loadFromFile2D(filepath, info);
		Layer t = new Layer(data, info.x, info.y, info.vset, info.iset);
		RHKFileOps.doFitting(t, 10);
		return t;
	}

	public static void openSomeFile()
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(davisDir);
		String s = FileOps.selectOpen(fc).toString();
		int type = getType(s);
		if (type == DIDV_MAP || type == CURRENT_MAP){
			Topomap t = getTopomapFromFile(s);
			new TopomapViewer(t, fc.getCurrentDirectory().toString(), 512);
		}
		else if (type == TOPO_MAP || type == FDBK_MAP){
			Layer t = getLayerFromFile(s);
			LayerViewer.show(t, 1024, true);
		}
		else{
			PointSpectra ps = loadPointSpectrum(s);
			new SpectraDrawer(ps, null);
		}
	}
	public static void main(String[] args)
	{
//		Topomap.setStdDir();
//		JFileChooser fc = new JFileChooser(Topomap.stddir);
//		Topomap t = getTopomapFromFile(FileOps.selectSave(fc).toString());
//		new TopomapViewer(t, fc.getCurrentDirectory().toString(), 512);
		openSomeFile();
	}
	
	public static int getType(String filepath)
	{
		if (filepath.endsWith(".1FL"))
			return DIDV_MAP;
		if (filepath.endsWith(".FFL"))
			return CURRENT_MAP;
		if (filepath.endsWith(".TFL"))
			return TOPO_MAP;
		if (filepath.endsWith(".FFR"))
			return FDBK_MAP;
		return -1;
	}
	
	public static void doDavisInfo()
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		loadFromFile(FileOps.selectOpen(fc).toString());
		
	}
	
}
