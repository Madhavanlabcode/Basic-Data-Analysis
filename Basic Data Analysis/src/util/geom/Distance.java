package util.geom;

import java.util.ArrayList;

import util.ArrayOps;
import util.Complex;

public class Distance {

	
	//the square extends from xc-sqSide/2 to xc+sqSide/2

	public static double distanceToSquareEdge(double sqSide, double xc, double yc, double x, double y)
	{
		double xm = xc - sqSide/2;
		double xp = xc + sqSide/2;
		double ym = yc - sqSide/2;
		double yp = yc + sqSide/2;
		
		//The square splits space into twelve regions.
		//The first four are the regions "outside" the four corners of the square.
		if (x >= xp && y >= yp)
			return distance(x-xp, y-yp);
		else if (x >= xp && y <= ym)
			return distance(x-xp, y-ym);
		else if (x <= xm && y >= yp)
			return distance(x-xm, y-yp);
		else if (x <= xm && y <= ym)
			return distance(x-xm, y-ym);

		//the next four regions are the horizontal stripes outward from the square's edges.
		else if (y <= ym && within(x, xm, xp))
			return distance(y-ym, 0);
		else if (y >= yp && within(x, xm, xp))
			return distance(y-yp, 0);
		else if (x <= xm && within(y, ym, yp))
			return distance(x-xm, 0);
		else if (x >= xp && within(y, ym, yp))
			return distance(x-xp, 0);
		
		//There are four differnent possible answers inside the square, and each such region is a right triangle with a side for the hypotaneuse.
		//but we won't consider the geometry
		double dxp = Math.abs(x-xp);
		double dxm = Math.abs(x-xm);
		double dyp = Math.abs(y-yp);
		double dym = Math.abs(y-ym);
		return Math.min(Math.min(dxp, dxm),Math.min(dyp, dym));
	}
	public static double distanceToSquareBody(double sqSide, double xc, double yc, double x, double y)
	{
		double xm = xc - sqSide/2;
		double xp = xc + sqSide/2;
		double ym = yc - sqSide/2;
		double yp = yc + sqSide/2;
		
		//The square splits space into twelve regions.
		//The first four are the regions "outside" the four corners of the square.
		if (x >= xp && y >= yp)
			return distance(x-xp, y-yp);
		else if (x >= xp && y <= ym)
			return distance(x-xp, y-ym);
		else if (x <= xm && y >= yp)
			return distance(x-xm, y-yp);
		else if (x <= xm && y <= ym)
			return distance(x-xm, y-ym);

		//the next four regions are the horizontal stripes outward from the square's edges.
		else if (y <= ym && within(x, xm, xp))
			return distance(y-ym, 0);
		else if (y >= yp && within(x, xm, xp))
			return distance(y-yp, 0);
		else if (x <= xm && within(y, ym, yp))
			return distance(x-xm, 0);
		else if (x >= xp && within(y, ym, yp))
			return distance(x-xp, 0);
		
		//inside the square...
		return 0;
	}
	public static double dot(double ax, double ay, double bx, double by)
	{
		return ax*bx + ay*by;
	}
	public static double[] unitVector(double[] a)
	{
		if (a[0] == 0 && a[1] == 0) return new double [2];
		return new double [] {a[0]/Complex.mag(a), a[1]/Complex.mag(a)};
	}
	
	public static double distance(double x, double y)
	{
		return Math.sqrt(x*x + y*y);
	}
	public static double distance(double[] b, double[] a)
	{
		return Math.sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
	}
	public static double distance(int[] i, double[] d)
	{
		return Math.sqrt((i[0]-d[0])*(i[0]-d[0]) + (i[1]-d[1])*(i[1]-d[1]));
	}
	public static double distance(int[] b, int[] a)
	{
		return Math.sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
	}
	public static double getAngleBetween(int[] b, int[] a)
	{
		return Math.acos(dot(a[0], a[1], b[0], b[1])/(Distance.distance(a[0], a[1])*Distance.distance(b[0], b[1])));
	}
	public static boolean within(double x, double xlower, double xupper)
	{
		return x >= xlower && x <= xupper;
	}
	/**
	 * x and y are assumed less than nx, ny
	 * @param x
	 * @param y
	 * @param nx
	 * @param ny
	 * @return
	 */
	public static double distancePeriodic(double x, double y, int nx, int ny)
	{
		if (x < nx/2 && y < ny/2)
			return Math.sqrt(x*x + y*y);
		else if (y < ny/2)
			return Math.sqrt(y*y + (nx-x)*(nx-x));
		else if (x < nx/2)
			return Math.sqrt(x*x + (ny-y)*(ny-y));
		else
			return Math.sqrt((nx-x)*(nx-x) + (ny-y)*(ny-y));

	}
	/**
	 * In order to avoid a long calculation we add all true pixels to an ArrayList. (Needs to be sparse therefore.)
	 * @param bs
	 * @param ds
	 */
	public static void putDistanceFromTrueSparse(boolean[][] bs, double[][] ds) {
		ArrayList<double[]> trues = new ArrayList<double[]>();
		for (int i = 0; i < bs.length; i++)
			for (int j = 0; j < bs[0].length; j++)
				if (bs[i][j]) trues.add(new double[] {i, j});
		
		double[] r = new double [2];
		for (int i = 0; i < bs.length; i++)
			for (int j = 0; j < bs[0].length; j++)
			{
				if (bs[i][j]) ds[i][j] = 0;
				else{
					ds[i][j] = Double.MAX_VALUE;
					r[0] = i;
					r[1] = j;
					for (int k = 0; k < trues.size(); k++)
						ds[i][j] = Math.min(ds[i][j], Distance.distance(r, trues.get(k)));
					
				}
			}
	}
	
	
	public static void putDistanceFromTrueSparseCenter(boolean[][] bs, double[][] ds, double cutoff) {
		ArrayList<double[]> trues = new ArrayList<double[]>();
		for (int i = 0; i < bs.length; i++)
			for (int j = 0; j < bs[0].length; j++)
				if (bs[i][j]) trues.add(new double[] {i, j});
		
		double dCent;
		double[] cent = new double[] {bs.length/2, bs[0].length/2};
		double[] r = new double [2];
		for (int i = 0; i < bs.length; i++)
			for (int j = 0; j < bs[0].length; j++)
			{
				if (bs[i][j]) ds[i][j] = 0;
				else{
					r[0] = i;
					r[1] = j;
					dCent = Distance.distance(r, cent);
					if (dCent > cutoff)
						ds[i][j] = cutoff;
					else{
						ds[i][j] = Double.MAX_VALUE;
						for (int k = 0; k < trues.size(); k++)
							ds[i][j] = Math.min(ds[i][j], Distance.distance(r, trues.get(k)));
					}
				}
			}
	}
	public static void putDistanceFromPointsInListI(ArrayList<int[]> pts, double[][] ds) {
		double[] r = new double [2];
		for (int i = 0; i < ds.length; i++){System.out.print(" " + i);
			for (int j = 0; j < ds[0].length; j++)
			{
				ds[i][j] = Double.MAX_VALUE;
				r[0] = i;
				r[1] = j;
				for (int k = 0; k < pts.size(); k++)
					ds[i][j] = Math.min(ds[i][j], Distance.distance(pts.get(k), r));
					
			}
		}
		System.out.println();
	}
	
	/**
	 * This returns an array of 4 numbers. The first number is the actual distance. The 2nd number is the abscissa of the closest point on the line segment.
	 * The 3rd number is the ordinate. The fourth number is 1 if the distance is within the line segment, 0 if is on the line beyond either end.
	 *
	 * @param xi
	 * @param xf
	 * @param yi
	 * @param yf
	 * @return
	 */
	public static double[] distanceToLineSegment(double xi, double xf, double yi, double yf, double x0, double y0, boolean forceToBeOnSegment)
	{
		double[] ans = new double [4];
		if (xi-xf == 0) //vertical line.
		{
			ans[0] = Math.abs(x0-xi);
			ans[1] = xi;
			ans[2] = y0;
			ans[3] = isBetween(y0, yi, yf) ? 1 : 0;
		}
		else
		{
			double m = (yf-yi)/(xf-xi);
			double b = yi - m*xi;
			ans[1] = (x0 + m*y0 - m*b)/(1 + m*m);
			ans[2] = m*ans[1] + b;
			ans[0] = Distance.distance(ans[1]-x0, ans[2]-y0);
			ans[3] = isBetween(ans[1], xi, xf) ? 1 : 0;
		}
		if (forceToBeOnSegment && ans[3] == 0)
		{
			ans[3] = 1;
			double di = Distance.distance(x0-xi, y0-yi);
			double df = Distance.distance(x0-xf, y0-yf);
			ans[0] = Math.min(di, df);
			ans[1] = df > di ? xf : xi;
			ans[2] = df > di ? yf : yi;
		}
		return ans;
	}
	/**
	 * This is brute force method costing an amount of time proportional to the number of points in the function
	 * @param x
	 * @param y
	 * @param x0
	 * @param y0
	 * @return
	 */
	public static double minimumDistanceToFunctionLinear(double[] x, double [] y, double x0, double y0)
	{
		double[] distance = new double [x.length-1];
		double[] distanceSegment;
		for (int i = 0; i < x.length-1; i++)
		{
			distanceSegment = distanceToLineSegment(x[i], x[i+1], y[i], y[i+1], x0, y0, true);
			distance[i] = distanceSegment[0];
		}
		return ArrayOps.min(distance);
	}

	/**
	 * This is the four-pointed shape, the inverse of the circle inscribed in a square.
	 * @param ds
	 * @param r
	 */
	public static void putDistanceFromPointyShapeCenter(double[][] ds, double r)
	{
		double[] c = {ds.length/2, ds[0].length/2};
		double[] rvec = new double [2];
		double[] right = {c[0]+r, c[1]};
		double[] left = {c[0]-r, c[1]};
		double[] top = {c[0], c[1]+r};
		double[] bottom = {c[0], c[1]-r};
		double[] upright = {c[0]+r, c[1]+r};
		double[] downright = {c[0]+r, c[0]-r};
		double[] upleft = {c[0]-r, c[1]+r};
		double[] downleft = {c[0]-r, c[1]-r};
		
		double[] d = new double [4];
		double dmin;
		
		double dx, dy;
		for (int i = 0; i < ds.length; i++)
			for (int j = 0; j < ds[0].length; j++)
			{
				rvec[0] = i; rvec[1] = j;
				dx = i-c[0];
				dy = j-c[1];
				
				//first four choices: outside the box:
				if (Math.abs(dx) > r || Math.abs(dy) > r){
					if (Math.abs(dy) < Math.abs(dx)){
						if (dx > 0) ds[i][j] = Distance.distance(rvec, right);
						else ds[i][j] = Distance.distance(rvec, left);
					}
					else {
						if (dy > 0) ds[i][j] = Distance.distance(rvec, top);
						else ds[i][j] = Distance.distance(rvec, bottom);
					}
				}
				else //if we are inside the box. 
				{
					d[0] = Distance.distance(rvec, upright);
					d[1] = Distance.distance(rvec, upleft);
					d[2] = Distance.distance(rvec, downright);
					d[3] = Distance.distance(rvec, downleft);
					dmin = ArrayOps.min(d);
					if (dmin >= r) ds[i][j] = 0;
					else ds[i][j] = r-dmin;
				}
			}
	}
	private static boolean isBetween(double value, double a, double b)
	{
		return (value >= a && value <= b) || (value <= a && value >= b);
	}
}
