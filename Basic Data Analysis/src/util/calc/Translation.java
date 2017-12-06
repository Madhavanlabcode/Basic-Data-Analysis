package util.calc;

import util.ArrayOps;
import util.FieldOps;
import util.fileops.ColumnIO;
import util.fourier.FFTOps;

public class Translation implements Calculation2D{
	
	double[][] x, y, cx, cy, deltax, deltay, big, small, xcut, ycut;
	int[][] xmain, ymain, xoth, yoth;
	double[][] source, target;
	double[][][] dr;
	int N;
	
	int currentDisplayField = 0;
	int DISP_FIELDS = 11;
	
	int poff = 100000;
	
	public Translation(double[][] source, double[][][] dr)
	{
		N = source.length;
		this.source = source;
		target = new double [N][N];
		this.dr = dr;
		x = new double[N][N];
		y = new double[N][N];
		cx = new double[N][N];
		cy = new double[N][N];
		deltax = new double[N][N];
		deltay = new double[N][N];
		big = new double[N][N];
		small = new double[N][N];
		xcut = new double[N][N];
		ycut = new double[N][N];
		xmain = new int [N][N];
		ymain = new int [N][N];
		xoth = new int [N][N];
		yoth = new int [N][N];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				x[i][j] = i+dr[i][j][0]+0.5;
				y[i][j] = j+dr[i][j][1]+0.5;
				cx[i][j] = (int)(x[i][j] + poff) + 0.5;
				cy[i][j] = (int)(y[i][j] + poff) + 0.5;
				xmain[i][j] = (int)cx[i][j];
				ymain[i][j] = (int)cy[i][j];
				cx[i][j] -= poff; cy[i][j] -= poff;
				xmain[i][j] -= poff; ymain[i][j] -= poff;
				deltax[i][j] = Math.abs(x[i][j] - cx[i][j]);
				deltay[i][j] = Math.abs(y[i][j] - cy[i][j]);
				big[i][j] = (1-deltax[i][j])*(1-deltay[i][j]);
				small[i][j] = deltax[i][j]*deltay[i][j];
				xcut[i][j] = deltax[i][j]*(1-deltay[i][j]);
				ycut[i][j] = deltay[i][j]*(1-deltax[i][j]);
				xoth[i][j] = x[i][j] > cx[i][j] ? xmain[i][j]+1 : xmain[i][j]-1;
				yoth[i][j] = y[i][j] > cy[i][j] ? ymain[i][j]+1 : ymain[i][j]-1;
				try{
					target[xmain[i][j]][ymain[i][j]] += source[i][j]*big[i][j];
				}
				catch(ArrayIndexOutOfBoundsException e){;}
				try{
					target[xoth[i][j]][ymain[i][j]] += source[i][j]*xcut[i][j];
				}
				catch(ArrayIndexOutOfBoundsException e){;}
				try{
					target[xmain[i][j]][yoth[i][j]] += source[i][j]*ycut[i][j];
				}
				catch(ArrayIndexOutOfBoundsException e){;}
				try{
					target[xoth[i][j]][yoth[i][j]] += source[i][j]*small[i][j];
				}
				catch(ArrayIndexOutOfBoundsException e){;}
			}
		
	}

	@Override
	public double[][] realField() {
		switch(currentDisplayField){
			case 0: return target;
			case 1:	return x;
			case 2:	return y;
			case 3:	return cx;
			case 4:	return cy;
			case 5: return deltax;
			case 6: return deltay;
			case 7: return big;
			case 8: return small;
			case 9: return xcut;
			case 10: return ycut;
		}
		return null;
	}
	String displaying()
	{
		switch(currentDisplayField){
		case 0: return "target";
		case 1:	return "x";
		case 2:	return "y";
		case 3:	return "cx";
		case 4:	return "cy";
		case 5: return "deltax";
		case 6: return "deltay";
		case 7: return "big";
		case 8: return "small";
		case 9: return "xcut";
		case 10: return "ycut";
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
		return "pixel = [" + x + ", " + y + "] \r\n" 
		+ "shift = " + vec(dr[x][y]) + "\r\n"
		+ "x = " + this.x[x][y] + "\r\n"
		+ "y = " + this.y[x][y] + "\r\n"
		+ "cx = " + cx[x][y] + "\r\n"
		+ "cy = " + cy[x][y] + "\r\n"
		+ "xmain = " + xmain[x][y] + "\r\n"
		+ "ymain = " + ymain[x][y] + "\r\n"
		+ "xoth = " + xoth[x][y] + "\r\n"
		+ "yoth = " + yoth[x][y] + "\r\n"
		+ "deltax = " + deltax[x][y] + "\r\n"
		+ "deltay = " + deltay[x][y] + "\r\n"
		+ "big = " + big[x][y] + "\r\n"
		+ "small = " + small[x][y] + "\r\n"
		+ "xcut = " + xcut[x][y] + "\r\n"
		+ "ycut = " + ycut[x][y] + "\r\n"
		+ "result = " + target[x][y] + "\r\n"
		+ "currently showing: " + displaying() + "\r\n";
	}

	@Override
	public String getInfo(double x, double y) {
//		return "K1 = " + vec(bragg[0]) + "\r\n"
//		+ "K2 = " + vec(bragg[1]) + "\r\n"
//		+ "K3 = " + vec(bragg[2]) + "\r\n"
//		+ "ratio = " + deform[y] + "\r\n"
//		+ "angle = " + Math.toDegrees(phi[x]) + " deg.\r\n"
//		+ "K1' = " + vec(braggmod[x][y][0]) + "\r\n"
//		+ "K2' = " + vec(braggmod[x][y][1]) + "\r\n"
//		+ "K3' = " + vec(braggmod[x][y][2]) + "\r\n"
//		+ "|K1'| = " + mag1[x][y] + "\r\n"
//		+ "|K2'| = " + mag2[x][y] + "\r\n"
//		+ "|K3'| = " + mag3[x][y] + "\r\n"
//		+ "theta12 = " + Math.toDegrees(th12[x][y]) + " deg.\r\n"
//		+ "theta23 = " + Math.toDegrees(th23[x][y]) + " deg.\r\n"
//		+ "theta13 = " + Math.toDegrees(th13[x][y]) + " deg.\r\n"
//		+ "Error12 = " + braggError[0][x][y] + "\r\n"
//		+ "Error23 = " + braggError[1][x][y] + "\r\n"
//		+ "Error31 = " + braggError[2][x][y] + "\r\n"
//		+ "Sum = " + sumErrors[x][y] + "\r\n";
//	// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getMinMax() {
		return new double[] {ArrayOps.min(realField()), ArrayOps.max(realField())};
	}


	@Override
	public void redoCalc() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void resize(double xmin, double xmax, double ymin, double ymax) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double[] getFieldBounds() {
		return new double[]{0, 512, 0, 512};
		
	}

	@Override
	public void switchDisplayField() {
		// TODO Auto-generated method stub
		currentDisplayField++;
		currentDisplayField %= DISP_FIELDS;
		
	}

	@Override
	public void save(String dir) {
		// TODO Auto-generated method stub
		
	}
	private static String vec(double[] x)
	{
		return "[" + x[0] + ", " + x[1] + "]";
	}
	
	public static Translation makeTranslation(String topo, double phi, double deform)
	{
		double[][] data = ColumnIO.readSquareTable(topo);
		double[][] fftmag = FFTOps.obtainFFTmagCent(data);
		double alpha = Math.cos(phi)*Math.cos(phi) + deform*Math.sin(phi)*Math.sin(phi);
		double beta = deform*Math.cos(phi)*Math.cos(phi) + Math.sin(phi)*Math.sin(phi);
		double gamma = (1 - deform)*Math.sin(phi)*Math.cos(phi);
		int N = data.length;
		FieldOps.log(fftmag);
		double[][][] dr = FieldOps.generateLinearTransformation(fftmag, new double[][] {{alpha, gamma}, {gamma, beta}}, new int []{N/2,N/2});
		return new Translation(fftmag, dr);
	}

	@Override
	public void run(String dir, String name) {
		// TODO Auto-generated method stub
		
	}

}
