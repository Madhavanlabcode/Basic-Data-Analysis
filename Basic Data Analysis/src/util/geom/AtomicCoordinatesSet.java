package util.geom;

import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.matrix.Matrix;
import util.Complex;
import util.FieldOps;
import util.Printer;

public class AtomicCoordinatesSet {

	//a and b are the two lattice vectors
	 double[] a;
	 double[] b;
	 double am, bm;
	 double dotaHatbHat;
	 
	 double alpha, sa, ca; //The angle between them; its sin and cosine.
	 double beta; //the angle supplementary to alpha
	 double[] aHat, bHat;
	 
	 //We will use the law of sines to find the atomic coordinates.
	 double[] origin; //This is the origin of the atomic coordinate system, in pixels.
	 
	 //for the code which determines lattice vectors from the Bragg peaks, see AtomicCoordinatesGenerator.java keyTyped() ' '.
	 
	 
	 private double[] temp = new double [2];
	 private double[] tempProja = new double [2]; //the projection of temp onto a.
	 private double[] tempProjb = new double [2]; //analagously for tempb.
	 private double[] tempReja = new double [2]; 
	 private double[] tempRejb = new double [2];
	 private double trda, trdb;
//	 private double[] tempProjab = new double [2]; //the projection of the projection, onto the other vector.
//	 private double[] tempProjba = new double [2];
	 //	 private double da, db, tam, tbm; //the dot products and temp magnitudes
	 private double tL, tTha, tThb, tda, tdb, tR; //temperoray variables start with t.
	 private double tempRa, tempRb;
	 
	public AtomicCoordinatesSet(double[] a, double[] b, double[] origin) {
		this.a = a;
		this.b = b;
		am = Complex.mag(a);
		bm = Complex.mag(b);
		
		ca = dot(a, b)/(am*bm);
		alpha = Math.acos(ca);
		beta = Math.PI - alpha;
		sa = Math.sin(alpha);
		aHat = new double[] {a[0]/am, a[1]/am};
		bHat = new double[] {b[0]/bm, b[1]/bm};
		

		dotaHatbHat = dot(aHat, bHat);
		
		this.origin = origin;
	}
	public AtomicCoordinatesSet(String s)
	{
		a = new double [2];
		b = new double [2];
		origin = new double [2];
		String[] lines = s.split("\r\n");
		String[] words;
		words = lines[0].split("\t");
		a[0] = Double.parseDouble(words[0]);
		a[1] = Double.parseDouble(words[1]);
		words = lines[1].split("\t");
		b[0] = Double.parseDouble(words[0]);
		b[1] = Double.parseDouble(words[1]);
		words = lines[2].split("\t");
		origin[0] = Double.parseDouble(words[0]);
		origin[1] = Double.parseDouble(words[1]);
		am = Complex.mag(a);
		bm = Complex.mag(b);

		ca = dot(a, b)/(am*bm);
		alpha = Math.acos(ca);
		beta = Math.PI - alpha;
		sa = Math.sin(alpha);
		aHat = new double[] {a[0]/am, a[1]/am};
		bHat = new double[] {b[0]/bm, b[1]/bm};
		dotaHatbHat = dot(aHat, bHat);
	}
	public static AtomicCoordinatesSet open(javax.swing.JFileChooser fc)
	{
		return new AtomicCoordinatesSet (ColumnIO.readString(FileOps.selectOpen(fc).toString()));
	}
	public void setOrigin(double x, double y)
	{
		origin = new double[] {x, y};
	}
	public void moveOrigin(double dx, double dy)
	{
		origin = new double[] {origin[0] + dx, origin[1] + dy};
	}
	public double[] getOrigin()
	{
		return new double[] {origin[0], origin[1]};
	}
	public void putAtomicCoords(double[] pixPt, double[] atomPt)
	{
		putAtomicCoords(pixPt[0], pixPt[1], atomPt);
	}
	public void putAtomicCoordsVector(double[] pixPt, double[] atomPt)
	{
		putAtomicCoordsVector(pixPt[0], pixPt[1], atomPt);
	}
	public void putAtomicCoords3(double pixX, double pixY, double[] atomPt)
	{
		temp[0] = pixX - origin[0];
		temp[1] = pixY - origin[1];
		tR = Complex.mag(temp);
		if (tR == 0)
		{
			atomPt[0] = 0;
			atomPt[1] = 0;
			return;
		}
		tda = dot(temp, aHat);
		tdb = dot(temp, bHat);
//		if (tda > tR) tda = tR;
//		else if (tda < -tR) tda = -tR;
		tempRa = tda/tR;
		tempRb = tdb/tR;
		if (tempRa > 1) tempRa = 1;
		if (tempRa < -1) tempRa = -1;
		if (tempRb > 1) tempRb = 1;
		if (tempRb < -1) tempRb = -1;
		tTha = Math.acos(tempRa);
		tThb = Math.acos(tempRb);
		tL = tR/sa;
		atomPt[0] = tL*Math.sin(tThb)/am;
		atomPt[1] = tL*Math.sin(tTha)/bm;
		//We try to account for the sign by watching the angles.
		if (tTha < alpha && tThb <= alpha); //do nothing; this is the default case.
		else if (tTha >= alpha && tThb < beta) atomPt[0] *= -1;
		else if (tTha > beta && tThb >= beta) {atomPt[0] *= -1; atomPt[1] *= -1;} //"Third quadrant"
		else if (tTha <= beta && tThb > alpha) atomPt[1] *= -1; //"Fourth quadrant."
	}

	/**
	 * Use orthonogonalization to put it properly.
	 * @param pixX
	 * @param pixY
	 * @param atomPt
	 */
	public void putAtomicCoordsVector(double pixX, double pixY, double[] atomPt)
	{
		temp[0] = pixX;
		temp[1] = pixY;
		
		tda = dot(temp, aHat);
		tdb = dot(temp, bHat);
		tempProja[0] = tda*aHat[0];
		tempProja[1] = tda*aHat[1];
		tempProjb[0] = tdb*bHat[0];
		tempProjb[1] = tdb*bHat[1];
		tempReja[0] = temp[0] - tempProja[0];
		tempReja[1] = temp[1] - tempProja[1];
		tempRejb[0] = temp[0] - tempProjb[0];
		tempRejb[1] = temp[1] - tempProjb[1];
		trdb = dot(tempReja, bHat);
		trda = dot(tempRejb, aHat);
		
		atomPt[0] = trda/(am*sa*sa);
		atomPt[1] = trdb/(bm*sa*sa);
	}
	public void putAtomicCoords(double pixX, double pixY, double[] atomPt)
	{
		putAtomicCoordsVector(pixX - origin[0], pixY-origin[1], atomPt);
	}
	public void putAtomicCoordsVector3(double pixX, double pixY, double[] atomPt)
	{
		temp[0] = pixX;
		temp[1] = pixY;
		tR = Complex.mag(temp);
		if (tR == 0)
		{
			atomPt[0] = 0;
			atomPt[1] = 0;
			return;
		}
		tda = dot(temp, aHat);
		tdb = dot(temp, bHat);
//		if (tda > tR) tda = tR;
//		else if (tda < -tR) tda = -tR;
		tempRa = tda/tR;
		tempRb = tdb/tR;
		if (tempRa > 1) tempRa = 1;
		if (tempRa < -1) tempRa = -1;
		if (tempRb > 1) tempRb = 1;
		if (tempRb < -1) tempRb = -1;
		tTha = Math.acos(tempRa);
		tThb = Math.acos(tempRb);
		tL = tR/sa;
		atomPt[0] = tL*Math.sin(tThb)/am;
		atomPt[1] = tL*Math.sin(tTha)/bm;
		//We try to account for the sign by watching the angles.
		if (tTha < alpha && tThb <= alpha); //do nothing; this is the default case.
		else if (tTha >= alpha && tThb < beta) atomPt[0] *= -1;
		else if (tTha > beta && tThb >= beta) {atomPt[0] *= -1; atomPt[1] *= -1;} //"Third quadrant"
		else if (tTha <= beta && tThb > alpha) atomPt[1] *= -1; //"Fourth quadrant."
	}
	public double[] getAtomicCoords(double[] pixPt){
		double[] atomPt = new double [2];
		putAtomicCoords(pixPt, atomPt);
		return atomPt;
	}
	public double[] getAtomicCoords(double pixX, double pixY)
	{
		double[] atomPt = new double [2];
		putAtomicCoords(pixX, pixY, atomPt);
		return atomPt;
	}
	
	public void putPixelCoords(double[] at, double[] ans){
		putPixelCoords(at[0], at[1], ans);
	}
	public void putPixelCoords(double atx, double aty, double[] ans)
	{
		ans[0] = atx*a[0] + aty*b[0] + origin[0];
		ans[1] = atx*a[1] + aty*b[1] + origin[1];
	}
	public double[] getPixelCoords(double[] at)
	{
		double[] ans = new double [2];
		putPixelCoords(at[0], at[1], ans);
		return ans;
	}
	public double[] getPixelCoords(double atx, double aty)
	{
		double[] ans = new double [2];
		putPixelCoords(atx, aty, ans);
		return ans;
	}
	public static double dot(double[] a, double[] b)
	{
		return a[0]*b[0] + a[1]*b[1];
	}
	public static double dot(double ax, double ay, double[] b)
	{
		return ax*b[0] + ay*b[1];
	}
	public double getAngleBetween(double[] vector, int latticeVectorIndex)
	{
		double[] latticeVector = latticeVectorIndex == 0 ? a : b;
		return Math.acos(dot(vector, latticeVector)/(Complex.mag(vector)*Complex.mag(latticeVector)));
	}
	public double getAngleBetweenVectors()
	{
		return Math.acos(dot(a, b)/(am*bm));
	}
	public String forDisplay()
	{
		String s = "";
		s += "The two lattice vectors are " + Printer.vectorP(a) + " and " + Printer.vectorP(b) + ".\r\n";
		s += "They have magnitudes " + am + " and " + bm + " respectively.\r\n";
		s += "The corresponding unit vectors are " + Printer.vectorP(aHat) + " and \r\n";
		s += "" + Printer.vectorP(bHat) + ".";
		s += "Their mutual angle is " + Math.toDegrees(alpha) + " degrees.\r\n";
		s += "The lattice origin is at " + Printer.vectorP(origin) + ".\r\n";
		return s;
	}
	
	public AtomicCoordinatesSet getRotatedCopy(double[][] mat, double[] rotOrigin)
	{
		double[] anew = Matrix.getProductWith(mat, a);
		double[] bnew = Matrix.getProductWith(mat, b);
		temp[0] = origin[0] - rotOrigin[0];
		temp[1] = origin[1] - rotOrigin[1];
		double[] newO = Matrix.getProductWith(mat, temp);
		newO[0] += rotOrigin[0];
		newO[1] += rotOrigin[1];
		return new AtomicCoordinatesSet(anew, bnew, newO);
	}
	/**
	 * Generates a rt2 lattice with lattice vectors a+b and a-b.
	 * @param latt
	 * @return
	 */
	public AtomicCoordinatesSet getRt2Lattice()
	{
		return new AtomicCoordinatesSet(new double[] {a[0] + b[0], a[1] + b[1]}, new double[] {a[0] - b[0], a[1] - b[1]}, origin);
	}
	public AtomicCoordinatesSet getOneOverRt2Lattice()
	{
		return new AtomicCoordinatesSet(new double[] {(a[0] + b[0])/2, (a[1] + b[1])/2}, new double[] {(a[0] - b[0])/2, (a[1] - b[1])/2}, origin);
	}
	public AtomicCoordinatesSet[] getRt2Lattices() //The rt2 and the rt2 shifted by one lattice vector
	{
		return new AtomicCoordinatesSet[] {
				new AtomicCoordinatesSet(new double[] {a[0] + b[0], a[1] + b[1]}, new double[] {a[0] - b[0], a[1] - b[1]}, origin),
				new AtomicCoordinatesSet(new double[] {a[0] + b[0], a[1] + b[1]}, new double[] {a[0] - b[0], a[1] - b[1]}, this.getPixelCoords(1, 0))};
	}
	/**
	 * If the pixel image has been doubled in size this will be the appropriate new lattice
	 * @return
	 */
	public AtomicCoordinatesSet getDoubledLattice()
	{
		return new AtomicCoordinatesSet(new double[] {2*a[0], 2*a[1]}, new double[] {2*b[0], 2*b[1]}, new double[] {2*origin[0], 2*origin[1]});
	}
	
	public String forDisplay2()
	{
		String s = "";
		s += Printer.vectorP(tda/am, tdb/bm) + "\t";
		s += Printer.vectorP(tTha, tThb) + "\t";
		s += Printer.vectorP(Complex.mag(aHat), Complex.mag(bHat)) + "\t";
		s += "tda = " + tda + ",\t" + "tdb = " + tdb + ",\t" + "tR = " + tR + ",\t";
		s += "tda/tR = " + (tempRa) + "\t tdb/tR" + (tempRb) + "\t";
		return s;
	}
	public String toString()
	{
		return Printer.arrayLnHorizontal(a) + Printer.arrayLnHorizontal(b) + Printer.arrayLnHorizontal(origin);
	}
	
	public static AtomicCoordinatesSet generateCentered(double[][] bragg, int N)
	{	
		double m, dot;
		double[][] latticeVectors = new double [2][2];
		for (int i = 0; i < 2; i++)
		{	//This is a clockwise rotation by 90 degrees with respect to the other bragg peak, given that the coordinate system is wrong-handed.
			latticeVectors[i][0] = -bragg[(i+1)%2][1];
			latticeVectors[i][1] = bragg[(i+1)%2][0];
			m = Complex.mag(latticeVectors[i]);
			for (int j = 0; j < 2; j++)
				latticeVectors[i][j] /= m; //which is the same as the magnitude of the Bragg peak of course.
			//Now we have unit vectors.
			dot = AtomicCoordinatesSet.dot(latticeVectors[i], bragg[i]);
			//the dot product should be N. Therefore, enlarge it by N/dot.
			for (int j = 0; j < 2; j++)
				latticeVectors[i][j] *= N/dot; 						
		}
		return new AtomicCoordinatesSet(latticeVectors[0], latticeVectors[1], new double[] {N/2, N/2});
	}
	public static AtomicCoordinatesSet generateCentered(int[][] bragg, int N)
	{	
		double m, dot;
		double[][] braggTrue = new double [2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		double[][] latticeVectors = new double [2][2];
		for (int i = 0; i < 2; i++)
		{	//This is a clockwise rotation by 90 degrees with respect to the other bragg peak, given that the coordinate system is wrong-handed.
			latticeVectors[i][0] = -braggTrue[(i+1)%2][1];
			latticeVectors[i][1] = braggTrue[(i+1)%2][0];
			m = Complex.mag(latticeVectors[i]);
			for (int j = 0; j < 2; j++)
				latticeVectors[i][j] /= m; //which is the same as the magnitude of the Bragg peak of course.
			//Now we have unit vectors.
			dot = AtomicCoordinatesSet.dot(latticeVectors[i], braggTrue[i]);
			//the dot product should be N. Therefore, enlarge it by N/dot.
			for (int j = 0; j < 2; j++)
				latticeVectors[i][j] *= N/dot; 						
		}
		return new AtomicCoordinatesSet(latticeVectors[0], latticeVectors[1], new double[] {N/2, N/2});
	}
	public static int[][] generateBragg(AtomicCoordinatesSet latt, int N)
	{	
		double m, dot;
		double[][] braggVectors = new double [2][2];
		double[][] sourceVectors = new double [2][2];
		sourceVectors[0] = latt.getA();
		sourceVectors[1] = latt.getB();
		for (int i = 0; i < 2; i++)
		{	//This is a clockwise rotation by 90 degrees with respect to the other bragg peak, given that the coordinate system is wrong-handed.
			braggVectors[i][0] = -sourceVectors[(i+1)%2][1];
			braggVectors[i][1] = sourceVectors[(i+1)%2][0];
			m = Complex.mag(braggVectors[i]);
			for (int j = 0; j < 2; j++)
				braggVectors[i][j] /= m; //the bragg vectors are now unit vectors
			//Now we have unit vectors.
			dot = AtomicCoordinatesSet.dot(braggVectors[i], sourceVectors[i]);
			//the dot product should be N. Therefore, enlarge it by N/dot.
			for (int j = 0; j < 2; j++)
				braggVectors[i][j] *= N/dot; 						
		}
		int[][] ans = new int [2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				ans[i][j] = FieldOps.round(braggVectors[i][j]);
		return ans;
	}

	//returns new double [] {aInv, bInv}.
	public double[][] getReciprocal()
	{	
		double dota, dotb;
		double[][] rVec = new double [2][2];
		rVec[0] = new double[] {b[1], -b[0]};
		rVec[1] = new double[] {a[1], -a[0]};
		dota = AtomicCoordinatesSet.dot(rVec[0], a);
		dotb = AtomicCoordinatesSet.dot(rVec[1], b);
		//The reciprocal lattice vectors satisfy rVec[0] * a = 2*PI so
		for (int j = 0; j < 2; j++)
		{
			rVec[0][j] *= 2*Math.PI/(dota); 						
			rVec[1][j] *= 2*Math.PI/(dotb); 						
		}
		return rVec;
	}
	public double[] getA() {
		return a;
	}
	public double[] getB() {
		return b;
	}
	public AtomicCoordinatesSet getReciprocalLattice()
	{	
		double dota, dotb;
		double[][] rVec = new double [2][2];
		rVec[0] = new double[] {b[1], -b[0]};
		rVec[1] = new double[] {a[1], -a[0]};
		dota = AtomicCoordinatesSet.dot(rVec[0], a);
		dotb = AtomicCoordinatesSet.dot(rVec[1], b);
		//The reciprocal lattice vectors satisfy rVec[0] * a = 2*PI so
		for (int j = 0; j < 2; j++)
		{
			rVec[0][j] *= 2*Math.PI/(dota); 						
			rVec[1][j] *= 2*Math.PI/(dotb); 						
		}
		return new AtomicCoordinatesSet(rVec[0], rVec[1], origin.clone());
	}

}
