package util.calc;

import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.geom.Distance;

public class StripCut {

	public double[][] map;
	int N;
	public int[] origin;
	
	double mean;
	
	public double[] result;
	
	int detail = 2;
	public int npts ;//in the linecut
	
	//the strip is on the radius vector, considering values from r0 to r1, with inclination angle theta.
	public double width, r0, r1, phi;
	public double[] rhat = new double[2], phihat = new double [2];
	double drstrip;
	
	public StripCut(double[][] map, double phi, double r0, double r1, double width, int npts)
	{
		this.map = map;
		N = map.length;
		origin = new int []{N/2, N/2};
		this.npts = npts;
		result = new double [npts];
		this.width = width;
		this.r0 = r0;
		this.r1 = r1;
		this.rhat = new double [2];
		this.phi = phi;
		rhat[0] = Math.cos(phi); rhat[1] = Math.sin(phi);
		phihat[0] = -Math.sin(phi); phihat[1] = Math.cos(phi);
		phihat = new double [2];
		drstrip = (r1-r0)/(npts-1);
		mean = FieldOps.mean(map);
		
	}
	
	//default:
	public StripCut(double[][] map, double width)
	{
		this.map = map;
		N = map.length;
		origin = new int []{N/2, N/2};
		this.npts = N/2;
		result = new double [npts];
		this.width = width;
		this.r0 = 0;
		this.r1 = N/2;
		this.rhat = new double [2];
		phihat = new double [2];
		this.phi = 0;
		rhat[0] = Math.cos(phi); rhat[1] = Math.sin(phi);
		phihat[0] = -Math.sin(phi); phihat[1] = Math.cos(phi);
		drstrip = (r1-r0)/(npts-1);
		mean = FieldOps.mean(map);
	}
	public StripCut(double[][] map, double width, int[] origin, int[] destination)
	{
		this.map = map;
		N = map.length;
		this.origin = origin;
		this.npts = N/2;
		result = new double [npts];
		this.width = width;
		this.r0 = 0;
		this.r1 = Distance.distance(origin, destination);
		this.rhat = new double [2];
		phihat = new double [2];
		this.phi = Complex.phase(new double [] {destination[0] - origin[0], destination[1] - origin[1]});
		rhat[0] = Math.cos(phi); rhat[1] = Math.sin(phi);
		phihat[0] = -Math.sin(phi); phihat[1] = Math.cos(phi);
		drstrip = (r1-r0)/(npts-1);
		mean = FieldOps.mean(map);
	}
	public void changeMap(double[][] map)
	{
		this.map = map;
		mean = FieldOps.mean(map);
	}
	
	
	public void makeCut()
	{
		drstrip = (r1-r0)/(npts-1);
		for (int i = 0; i < npts; i++)
			result[i] = 0;
		
		double[] rs = ArrayOps.generateArrayInclBoth(r0, r1, npts);
		double[] ws = ArrayOps.generateArrayInclBoth(-width, width, (int)(detail*width));
		
		double x, y, xw, yw;
		
		for (int i = 0; i < npts; i++){
			x = rs[i]*Math.cos(phi) + (origin[0]);
			y = rs[i]*Math.sin(phi) + (origin[1]);
			for (int j = 0; j < ws.length; j++)
			{
				xw = x - ws[j]*Math.sin(phi);
				yw = y + ws[j]*Math.cos(phi);
				result[i] += FieldOps.getValueAt(map, xw, yw, mean);
			}
			result[i] /= ws.length;
		}	
		
	}
	public void setWidth(double width) {
		this.width = width;
		makeCut();
	}
	public void setR0(double r0) {
		this.r0 = r0;
		makeCut();
	}
	public void setR1(double r1) {
		this.r1 = r1;
//		makeCut();
	}
	public void setPhi(double phi) {
		this.phi = phi;
		rhat[0] = Math.cos(phi); rhat[1] = Math.sin(phi);
		phihat[0] = -Math.sin(phi); phihat[1] = Math.cos(phi);
//		makeCut();
	}
	public void setAll(double[] all) {
		phi = all[0];
		r0 = all[1];
		r1 = all[2];
		width = all[3];
		makeCut();
	}
	public int[][] getFourCorners()
	{
		int[][] fc = new int [2][4];
		fc[0][0] = (int)(r0*rhat[0] - (width/2)*phihat[0]);
		fc[1][0] = (int)(r0*rhat[1] - (width/2)*phihat[1]);
		fc[0][1] = (int)(r0*rhat[0] + (width/2)*phihat[0]);
		fc[1][1] = (int)(r0*rhat[1] + (width/2)*phihat[1]);
		fc[0][2] = (int)(r1*rhat[0] + (width/2)*phihat[0]);
		fc[1][2] = (int)(r1*rhat[1] + (width/2)*phihat[1]);
		fc[0][3] = (int)(r1*rhat[0] - (width/2)*phihat[0]);
		fc[1][3] = (int)(r1*rhat[1] - (width/2)*phihat[1]);
		return fc;
	}
	public String getOutput()
	{
		String it = "r\tIntensity\r\n";
		for (int i = 0; i < npts; i++)
		{
			it += (r0 + i*drstrip) + "\t" + (result[i]) + "\r\n";
		}
		return it;
	}
	public static double dot(double[] a, double[] b)
	{
		return a[0]*b[0] + a[1]*b[1];
	}

	public void setNPTS(int npts) {
		// TODO Auto-generated method stub
		this.npts = npts;
		result = new double [npts];
		
	}	
	
}
