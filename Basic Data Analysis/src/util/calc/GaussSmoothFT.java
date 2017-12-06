package util.calc;

import java.io.File;

import javax.swing.JFileChooser;

import main.SRAW;
import util.ArrayOps;
import util.FieldOps;
import util.fileops.ColumnIO;
import util.fourier.FFTOps;

public class GaussSmoothFT extends util.movie.BasicMovie implements Calculation2D_1P {

	int N;
	int power;
	double[][] topo;
	double[][][] fftz, nfftz, ift;
	double[][] nfftlogmag;

	double[][] smooth;
	double[][] subloc;
	double[][] gauss;
	double[][] newfftlogmag;
	
	int current = 0;
	final int DISP_FIELDS = 4;
	
	public GaussSmoothFT(double[][] topo, int power)
	{
		this.power = power;
		N = topo.length;
		this.topo = topo;
		fftz = new double[N][N][2];
		nfftz = new double [N][N][2];
		nfftlogmag = new double [N][N];
		ift = new double [N][N][2];
		gauss = new double [N][N];
		smooth = new double [N][N];
		subloc = new double [N][N];
		newfftlogmag = new double [N][N];
		FFTOps.obtainFFTCent(topo, fftz);
	}
	
	@Override
	public double[][] realField() {
		switch(current)
		{
		case 0:
			return smooth;
		case 1:
			return nfftlogmag;
		case 2:
			return subloc;
		case 3:
			return newfftlogmag;
		}
		return null;
	}

	@Override
	public double[][][] complexField() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getInfo(int x, int y) {
		// TODO Auto-generated method stub
		return "Placeholder \r\n ";
	}

	@Override
	public String getInfo(double x, double y) {
		// TODO Auto-generated method stub
		return "";
	}

	@Override
	public double[] getMinMax(boolean real) {
		// TODO Auto-generated method stub
		return new double [] {ArrayOps.min(realField()), ArrayOps.max(realField())};
	}

	@Override
	public void redoCalc(double p) {
		// TODO Auto-generated method stub
		double r;
		if ( p != 2)
			for (int k = 0; k < topo.length; k++)
				for (int l = 0; l < topo[0].length; l++){
					r = Math.sqrt(((double)(k*k) + (l*l)));
					gauss[k][l] = Math.exp(-Math.pow(r/p, power));
				}
		else
			for (int k = 0; k < topo.length; k++)
				for (int l = 0; l < topo[0].length; l++){
					r = Math.sqrt(((double)(k*k) + (l*l)));
					gauss[k][l] = Math.exp(-(r*r)/(2*p*p));
				}
			
		
		gauss[0][0] = 0.999;

		for (int k = 0; k < N; k++)
			for (int l = 0; l < N; l++)
			{
				nfftz[k][l][0] = fftz[k][l][0]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
				nfftz[k][l][1] = fftz[k][l][1]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
			}
		FFTOps.obtainIFFTCent(nfftz, ift);
		FieldOps.getIndex(ift, smooth, 0);
		FieldOps.minus(topo, smooth, subloc);
		setSpecialCurrent();
	}
	private void setSpecialCurrent()
	{
		if (current == 3){
			FFTOps.obtainFFTmagCent(subloc, newfftlogmag);
			FieldOps.log(newfftlogmag);
		}
		else if (current == 1)
		{
			FieldOps.magnitude(nfftz, nfftlogmag);
//			FieldOps.log(nfftlogmag);
		}
	}

	@Override
	public void resize(double xmin, double xmax, double ymin, double ymax) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double[] getFieldBounds() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void switchDisplayField() {
		// TODO Auto-generated method stub
		current++;
		current %= DISP_FIELDS;
		setSpecialCurrent();
	}

	@Override
	public void save(String dir) {
		// TODO Auto-generated method stub
		JFileChooser fc = new JFileChooser(dir);
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION){
			FFTOps.obtainFFTmagCent(subloc, newfftlogmag);
//			FieldOps.log(newfftlogmag);
			FieldOps.magnitude(nfftz, nfftlogmag);
//				FieldOps.log(nfftlogmag);

			
			ColumnIO.writeBin(subloc, fc.getSelectedFile().toString() + "subloc.dat");
			SRAW.writeImage(fc.getSelectedFile().toString() + "subloc", subloc);
			SRAW.writeImage(fc.getSelectedFile().toString() + "sublocfft", newfftlogmag);
			ColumnIO.writeBin(smooth, fc.getSelectedFile().toString() + "locavg.dat");
			SRAW.writeImage(fc.getSelectedFile().toString() + "locavg", smooth);
			SRAW.writeImage(fc.getSelectedFile().toString() + "locavgfft", nfftlogmag);
		}
		
	}

	@Override
	public void run(String dir, String name) {
		// TODO Auto-generated method stub
		
	}
	public double[][] getRField(double p){
		redoCalc(p);
		return realField();
	}
	public double[][][] getCField(double p){
		redoCalc(p);
		return complexField();
	}
	
	public static GaussSmoothFT getNew(String dir, String name)
	{
		return new GaussSmoothFT(ColumnIO.readSquareTable(dir + name + ".dat"), 2);
	}
	public static GaussSmoothFT getNew(String dir)
	{
		JFileChooser fc = new JFileChooser(dir);
		File f = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
			f = fc.getSelectedFile();
		return new GaussSmoothFT(ColumnIO.readSquareTable(f.toString()), 8);
	}
	public static GaussSmoothFT getNew(String dir, String name, int p)
	{
		return new GaussSmoothFT(ColumnIO.readSquareTable(dir + name + ".dat"), p);
	}

	@Override
	public String getInfo(double p) {
		// TODO Auto-generated method stub
		return null;
	}

}
