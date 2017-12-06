package util.geom;

//This includes the conversion between meters and pixels for the Topograph and Layer classes.
public class MapRectCoordSystem {

	public double[] x;
	public double[] y;
	public int nx, ny;
	public double xLength, yLength;
	private double dx, dy;
	
	private double dkx, dky;
	
	//the origin of the coordinate system is automatically (x[0], y[0]).
	
	public MapRectCoordSystem(double[] x, double[] y) {
		super();
		this.x = x;
		this.y = y;
		nx = x.length;
		ny = y.length;
		xLength = ((x[x.length-1]-x[0])/(x.length-1))*nx;
		yLength = ((y[y.length-1]-y[0])/(y.length-1))*ny;
		dx = xLength/nx;
		dy = yLength/ny;
		dkx = 2*Math.PI/xLength;
		dky = 2*Math.PI/yLength;
	}
	
	
	public void resetCoords()
	{
		nx = x.length;
		ny = y.length;
		xLength = ((x[x.length-1]-x[0])/(x.length-1))*nx;
		yLength = ((y[y.length-1]-y[0])/(y.length-1))*ny;
		dx = xLength/nx;
		dy = yLength/ny;
		dkx = 2*Math.PI/xLength;
		dky = 2*Math.PI/yLength;
	}

	public void flipYSign()
	{
		for (int i = 0; i < nx; i++)
			y[i] = -y[i];
		resetCoords();
	}
	public void flipXSign()
	{
		for (int i = 0; i < nx; i++)
			y[i] = -y[i];
		resetCoords();
	}
	/**
	 * This is if a stupid mistake has been made in the size of things
	 * @param factor
	 */
	public void resizeSimple(double factor)
	{
		for (int i = 0; i < nx; i++)
			x[i] *= factor;
		for (int i = 0; i < ny; i++)
			y[i] *= factor;
		resetCoords();
	}
	//The center here enters in nanometers.
	public void makeOrigin(double[] center)
	{
		x[0] = center[0]*Math.pow(10, -9) - xLength/2;
		y[0] = center[1]*Math.pow(10, -9) - yLength/2;
		for (int i = 0; i < nx; i++)
			x[i] = x[0] + i*dx;
		for (int i = 0; i < ny; i++)
			y[i] = y[0] + i*dy;
	}

	public void putMetricCoords(double[] pixPt, double[] atomPt)
	{
		putMetricCoords(pixPt[0], pixPt[1], atomPt);
	}
	public void putMetricCoords(double pixX, double pixY, double[] metricPt)
	{
		metricPt[0] = pixX*dx + x[0];
		metricPt[1] = pixY*dy + y[0];
	}
	public double[] getMetricCoords(double[] pixPt){
		double[] atomPt = new double [2];
		putMetricCoords(pixPt, atomPt);
		return atomPt;
	}
	public double[] getMetricCoords(double pixX, double pixY)
	{
		double[] atomPt = new double [2];
		putMetricCoords(pixX, pixY, atomPt);
		return atomPt;
	}
	
	public void putPixelCoords(double[] at, double[] ans){
		putPixelCoords(at[0], at[1], ans);
	}
	public void putPixelCoords(double metx, double mety, double[] ans)
	{
		ans[0] = (metx-x[0])/dx;
		ans[1] = (mety-y[0])/dy;
	}
	public double[] getPixelCoords(double[] at)
	{
		double[] ans = new double [2];
		putPixelCoords(at[0], at[1], ans);
		return ans;
	}
	public void putPixelCoordsCent(double metx, double mety, double[] ans)
	{
		ans[0] = 0.5+(metx-x[0])/dx;
		ans[1] = 0.5+(mety-y[0])/dy;
	}
	public double[] getPixelCoordsCent(double[] at)
	{
		double[] ans = new double [2];
		putPixelCoordsCent(at[0], at[1], ans);
		return ans;
	}
	public double[] getPixelCoords(double atx, double aty)
	{
		double[] ans = new double [2];
		putPixelCoords(atx, aty, ans);
		return ans;
	}
	public double[] getPixelCoordsCent(double atx, double aty)
	{
		double[] ans = new double [2];
		putPixelCoordsCent(atx, aty, ans);
		return ans;
	}
	
	public boolean containsMetric(double mx, double my)
	{
		double[] p = getPixelCoords(mx, my);
		return contains(p[0], p[1]);
	}
	public boolean contains(double px, double py)
	{
		return px >= 0 && px <= nx-1 && py >= 0 && py <= ny-1;
	}
	
	/**
	 * These attempt to locate a point in k-space. The assumption is that (0, 0) is at the center of the place.
	 * @param pixPt
	 * @param atomPt
	 */
	public void putFourierMetricCoords(double[] pixPt, double[] atomPt)
	{
		putFourierMetricCoords(pixPt[0], pixPt[1], atomPt);
	}
	public void putFourierMetricCoords(double pixX, double pixY, double[] metricPt)
	{
		metricPt[0] = (pixX-nx/2)*dkx;
		metricPt[1] = (pixY-ny/2)*dky;;
	}
	public double[] getFourierMetricCoords(double[] pixPt){
		double[] atomPt = new double [2];
		putFourierMetricCoords(pixPt, atomPt);
		return atomPt;
	}
	public double[] getFourierMetricCoords(double pixX, double pixY)
	{
		double[] atomPt = new double [2];
		putFourierMetricCoords(pixX, pixY, atomPt);
		return atomPt;
	}

	//This takes a point that is over the edge and moves it back to the edge.
	//This assumes that dx and dy are positive!
	public void restrictPt(double[] r)
	{
		if (r[0] < x[0]) r[0] = x[0];
		if (r[0] > x[x.length-1]) r[0] = x[x.length-1];
		if (r[1] < y[0]) r[1] = y[0];
		if (r[1] > y[y.length-1]) r[1] = y[y.length-1];
	}
	public double[] restrictPt(double xp, double yp)
	{
		double[] pt = new double[] {xp, yp};
		if (pt[0] < x[0]) pt[0] = x[0];
		if (pt[0] > x[x.length-1]) pt[0] = x[x.length-1];
		if (pt[1] < y[0]) pt[1] = y[0];
		if (pt[1] > y[y.length-1]) pt[1] = y[y.length-1];
		return pt;
	}
	public double[] restrictPtClone(double[] r)
	{
		double[] pt = r.clone();
		if (pt[0] < x[0]) pt[0] = x[0];
		if (pt[0] > x[x.length-1]) pt[0] = x[x.length-1];
		if (pt[1] < y[0]) pt[1] = y[0];
		if (pt[1] > y[y.length-1]) pt[1] = y[y.length-1];
		return pt;
	}
	
	public void translate(double dx, double dy)
	{
		for (int i = 0; i < nx; i++)
			x[i] += dx;
		for (int i = 0; i < ny; i++)
			y[i] += dy;
	}
}
