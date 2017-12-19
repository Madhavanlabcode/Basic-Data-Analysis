package util.fileops;

import image.ImageEditing;


//Added vectors to store previous versions of zooms
import java.util.Vector;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolator;

//import flanagan.interpolation.BiCubicInterpolation;



import util.ArrayOps;
import util.FieldOps;
import util.geom.MapRectCoordSystem;

public class Layer  extends MapRectCoordSystem{

	public double v, current;
	public double[][] data;
	private double mean;
	
	//Added a vector for storing previous zooms
	public Vector<double[][]> prevZooms;
	
	public double lineTime = 0; //Nanonis only; not saved in the file.
	
	public boolean fileIsComplete = true;
	
	//This is necessary for bicubic interpolation
//	private BiCubicInterpolation interp = null;
	private BicubicSplineInterpolatingFunction interp = null;
	
	public Layer(double[][] data, double[] x, double[] y, double v, double current)
	{
		super(x, y);
		this.data = data;
		this.x = x;
		this.y = y;
		this.v = v;
		this.current = current;
		this.mean = FieldOps.mean(data);
		
		//Added
		this.prevZooms= new Vector();
	}
	public int getByteLength()
	{
		int l = 0;
		l += 24;
		l += 8*(nx + ny);
		l += 8*(nx*ny);
		l += 16; //for the origin
		return l;
	}
	////////////////////////////////////////////////////////
	public void zoomIn(){
		String inputCoord=JOptionPane.showInputDialog("Enter the Coordinates of Upper Left Corner of Zoom separated by a comma");
		String[] pair=inputCoord.split(",");
		int xStart=Integer.parseInt(pair[0].trim());
		int yStart=Integer.parseInt(pair[1].trim());
		String inputZoom=JOptionPane.showInputDialog("Enter the final size of the Zoom square in pixels ");
		int zoomSize=Integer.parseInt(inputZoom.trim());
		
		if((xStart+zoomSize)>(this.nx-1) || (yStart+zoomSize)>(this.ny - 1)){
			System.out.println("Your zoom is out of the bounds of the image");
		}
		else{
			prevZooms.addElement(this.data);
			double[][] tempdata= new double[zoomSize][zoomSize];
			for(int i=xStart; i<(xStart+zoomSize); i++){
				for(int j=yStart; j<(yStart+zoomSize); j++){
					System.out.println(xStart + yStart + zoomSize);
					tempdata[i-xStart][j-yStart]=data[i][j];
				}
			}
			this.data=tempdata;
		}
		
	}
	
	
	////////////////////////////////////////////////////////
	//This method accomplishes the 0.5 pixel shift which is in FieldOps.getValueAt().
	public double evaluateAt(double[] pixPt)
	{
		return FieldOps.getValueAt(data, pixPt[0]+0.5, pixPt[1]+0.5, mean);
	}
	public double evaluateAt(double px, double py)
	{
		return FieldOps.getValueAt(data, px+0.5, py+0.5, mean);
	}
	public double evaluateAtMetric(double[] metricPt)
	{
		double[] temp = this.getPixelCoords(metricPt);
		return FieldOps.getValueAt(data, temp[0]+0.5, temp[1]+0.5, mean);
	}
	public double evaluateAtMetric(double mx, double my)
	{
		double[] temp = this.getPixelCoords(mx, my);
		return FieldOps.getValueAt(data, temp[0]+0.5, temp[1]+0.5, mean);
	}
	public double evaluateAtMetricInt(double mx, double my)
	{
		double[] temp = this.getPixelCoords(mx, my);
		if (FieldOps.withinBounds(FieldOps.round(temp[0]), FieldOps.round(temp[1]), data))
			return data[FieldOps.round(temp[0])][FieldOps.round(temp[1])];
		else return mean;
	}
	
	public double evaluateBiCubic(double[] pixPt)
	{
		if (interp == null)
			initializeInterpolator();
		double[] metric = this.getMetricCoords(pixPt);
		restrictPt(metric);
		return interp.value(metric[0], metric[1]);
	}
	public double evaluateBiCubic(double px, double py)
	{
		if (interp == null)
			initializeInterpolator();
		double[] metric = this.getMetricCoords(px, py);
		restrictPt(metric);
		return interp.value(metric[0], metric[1]);
	}
	public double evaluateBiCubicMetric(double[] metricPt)
	{
		if (interp == null)
			initializeInterpolator();
		double[] pt = restrictPtClone(metricPt);
		return interp.value(pt[0], pt[1]);
	}
	public double evaluateBiCubicMetric(double mx, double my)
	{
		if (interp == null)
			initializeInterpolator();
		double[] pt = restrictPt(mx, my);
		return interp.value(pt[0], pt[1]);
	}
	
	private void initializeInterpolator() {
		// TODO Auto-generated method stub
		BicubicSplineInterpolator erp = new BicubicSplineInterpolator();
		interp = erp.interpolate(x, y, data);
	}
	public void flipX()
	{
		ArrayOps.flipX(data);
		ArrayOps.flip(x);
		resetCoords();
	}
	public void flipY()
	{
		ArrayOps.flipY(data);
		y = ArrayOps.getFlip(y);
//		ArrayOps.flip(x);
		resetCoords();
	}
	public PointSpectra toSpectra()
	{
		return new PointSpectra(FieldOps.transpose(data), x, y, new double [y.length]);
	}
	
	public static void writeBIN(Layer t, String filepath)
	{
		File file = new File(filepath);
		
		FileOutputStream outf = null;
		BufferedOutputStream outbuff = null;
		DataOutputStream outd = null;
		Topomap.putStdDir(file.getParent() + "\\");

		try {
			outf = new FileOutputStream(file);
			outbuff = new BufferedOutputStream(outf, t.getByteLength());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		outd = new DataOutputStream(outbuff);
		try {
			outd.writeInt(t.nx);
			outd.writeInt(t.ny);
			outd.writeDouble(t.v);
			outd.writeDouble(t.current);
			for (int i = 0; i < t.nx; i++)
				outd.writeDouble(t.x[i]);
			for (int i = 0; i < t.ny; i++)
				outd.writeDouble(t.y[i]);
			for (int j = 0; j < t.nx; j++)
				for (int k = 0; k < t.ny; k++)
					outd.writeDouble(t.data[j][k]);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
//			outbuff.flush();
//			outf.flush();
			outd.close();
			outbuff.close();
			outf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static Layer readBIN(String filepath)
	{
		int nx, ny;
		double[] x = null, y = null;
		double v = 0, current = 0;
		double[][] data = null;
		double[] origin = new double[2];
		File file = new File(filepath);
		Topomap.putStdDir(file.getParent() + "\\");
		FileInputStream inf = null;
		BufferedInputStream inbuff = null;
		DataInputStream ind = null;
		Layer t = null;
		try {
			inf = new FileInputStream(file);
			inbuff = new BufferedInputStream(inf, (int)file.length());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		ind = new DataInputStream(inbuff);
		try {
			nx = ind.readInt();
			ny = ind.readInt();

			//read 8 bytes: if the first 4 bytes are an int with length 1:
			//open only layer of topomap as a Layer
			//else continue as if it's a Layer
			byte[] testBytes = new byte[8];
			for(int i=0;i<8;i++){
				testBytes[i]=ind.readByte();
			}
			
			if(ByteBuffer.wrap(testBytes).getInt()==1){
				System.out.println("Opening single-layer Topomap as Layer");
				ind.close();
				inbuff.close();
				inf.close();
				Topomap t1=Topomap.readBIN(filepath);
				return t1.getLayer(0);
			}
			v = ByteBuffer.wrap(testBytes).getDouble();
			
			current = ind.readDouble();
			x = new double[nx]; y = new double [ny];
			data = new double [nx][ny];
			
			for (int i = 0; i < nx; i++)
				x[i] = ind.readDouble();
			for (int i = 0; i < ny; i++)
				y[i] = ind.readDouble();
			for (int j = 0; j < nx; j++)
				for (int k = 0; k < ny; k++)
					data[j][k] = ind.readDouble();
			t = new Layer(data, x, y, v, current);
			
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return t;
	}
	public static void writeBIN(Layer t)
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			Layer.writeBIN(t, s);
		}
	}
	public static void writeBIN(Layer t, JFileChooser fc)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
		{
			String s = fc.getSelectedFile().toString();
			if (!s.endsWith(".bin")) s += ".bin";
			Layer.writeBIN(t, s);
		}
	}
	public static Layer open()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
			Topomap.currentDir = fc.getCurrentDirectory().toString() + "\\";
			if (fc.getSelectedFile().toString().endsWith(".bin"))
				return Layer.readBIN(fc.getSelectedFile().toString());
			else if (fc.getSelectedFile().toString().endsWith(".txt"))
				return RHKFileOps.getLayer(fc.getSelectedFile());
			else;
		}
		return null;
	}
	public static Layer open(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			if (fc.getSelectedFile().toString().endsWith(".bin"))
				return Layer.readBIN(fc.getSelectedFile().toString());
			else if (fc.getSelectedFile().toString().endsWith(".txt"))
				return RHKFileOps.getLayer(fc.getSelectedFile());
			else;
		return null;
	}
	//This replaces the data block of a topomap while keeping the voltages and horizontal positions the same
	public static Layer openFree()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		return openFree(fc);
	}
	public static Layer openFree(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION){
			Topomap.currentDir = fc.getCurrentDirectory().toString() + "\\";
			return openFree(fc.getSelectedFile());
		}
		return null;
		
	}
	public static Layer openFree(File fullpath)
	{
			if (fullpath.toString().endsWith(".bin"))
				return Layer.readBIN(fullpath.toString());
			else if (fullpath.toString().endsWith(".txt"))
				return RHKFileOps.getLayer(fullpath);
			else if (fullpath.toString().endsWith(".dat")) //then there is only the nxn topographic map with no extra data
				return Layer.getFreeLayer(ColumnIO.readSquareTable(fullpath.toString()));
			else if (fullpath.toString().endsWith(".sxm")){
				String message = "What data do you want to view?\r\n";
				String[] theNames = NanonisFileOps.loadLayerNamesFromScan(fullpath.toString());
				for (int i = 0; i < theNames.length; i++)
					message += "" + i + " - " + theNames[i] + "\r\n";
				int o = Integer.parseInt(JOptionPane.showInputDialog(message).trim());
				Layer[] m = NanonisFileOps.loadLayerFromScan(fullpath.toString(), true);
				return m[o];
			}
			else if (fullpath.toString().endsWith(".bmp"))
				return Layer.getFreeLayer(ImageEditing.getFromImage(ImageEditing.open(fullpath.toString())));
			else;
		
		return null;
	}
	public static Layer getFreeLayer(double[][] data)
	{
		double[] x = new double [data.length];
		double[] y = new double [data[0].length];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int j = 0; j < y.length; j++)
			y[j] = j;
		
		double v = 1;
		double current = 1;
		return new Layer(data, x, y, v, current);
	}
	public static Layer getCropped(Layer t, int xi, int xf, int yi, int yf)
	{
		int nx = xf-xi;
		int ny = yf-yi;
		double[] x = new double[nx];
		double[] y = new double[ny];
		
		double[][] data = new double [nx][ny];
		
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				data[i][j] = t.data[xi+i][yi+j];
		for (int i = 0; i < nx; i++)
			x[i] = t.x[i+xi];
		for (int i = 0; i < ny; i++)
			y[i] = t.y[i+yi];
		
		return new Layer(data,x, y, 1, 1);
	}

	
	public static Layer newLayer(Layer t, double[][] data)
	{
		return new Layer(data, t.x, t.y, t.v, t.current);
	}
	public static Layer newLayer(Topomap t, double[][] data)
	{
		return new Layer(data, t.x, t.y, t.v[0], 1.0);
	}

	/**
	 * returns a new layer, double the size of the old one
	 * @param t
	 * @return
	 */
	public static Layer embiggen(Layer t)
	{
		int nx = 2*t.nx;
		int ny = 2*t.ny;
		double[][] data = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				data[i][j] = t.evaluateAt(i/2.0, j/2.0);
			}
		double[] x = new double [nx];
		double[] y = new double [ny];
		for (int i = 0; i < nx; i++)
			x[i] = t.x[0] + i*t.xLength/nx;
		for (int i = 0; i < ny; i++)
			y[i] = t.y[0] + i*t.yLength/ny;
		return new Layer(data, x, y, t.v, t.current);
	}
	public void setMean() {
		mean = FieldOps.mean(data);
	}


}
