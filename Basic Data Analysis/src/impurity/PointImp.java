package impurity;

import java.util.ArrayList;

import util.FieldOps;
import util.fileops.ColumnIO;
import util.fileops.Layer;

public class PointImp {



	public double[] pixelPos;
	
	public PointImp(double[] inputPos) {
		super();
		pixelPos = inputPos;
	}
	
	public boolean equals(PointImp b){
		if(b.pixelPos[0]==this.pixelPos[0]&&b.pixelPos[1]==this.pixelPos[1])
			return true;
		else
			return false;
	}
	
	public static PointImp[] readFromGaussSquareFile(java.io.File f)
	{
		double[][] table = ColumnIO.readNColumns(f, 2, 3);
		PointImp[] imps = new PointImp[table[0].length]; 
		for (int i = 0; i < imps.length; i++)
			imps[i] = new PointImp(new double [] {table[0][i], table[1][i]});
		return imps;
	}
	public static int[] readEnergyIndices(java.io.File f)
	{
		double[][] table = ColumnIO.readNColumns(f, 3, 3);
		int[] is = new int[table[0].length]; 
		for (int i = 0; i < is.length; i++)
			is[i] = (int)table[2][i];
		return is;
	}
	public static void writeToFile(PointImp[] imps, java.io.File f)
	{
		String[] lines = new String [imps.length+3];
		lines[0] = "Impurity List";
		lines[1] = "()";
		lines[2] = "X\tY";
		for (int i = 0; i < imps.length; i++)
			lines[i+3] = "" + imps[i].pixelPos[0] + "\t" + imps[i].pixelPos[1];
		
		ColumnIO.writeLines(lines, f.toString());
	}
	
	public static PointImp[] getAboveCutoffContiguous(Layer t, double cutoff)
	{
		boolean[][] possibleNotVisitedImps = FieldOps.isGreaterThan(t.data, cutoff);
		FieldOps.negate(possibleNotVisitedImps);
		ArrayList<ArrayList<int[]>> blobs = FieldOps.getAllContiguousBlobs(possibleNotVisitedImps);
		
		PointImp[] imps = new PointImp[blobs.size()];
		for (int i = 0; i < imps.length; i++)
		{
			imps[i] = new PointImp(FieldOps.getCentroid(blobs.get(i)));
		}
		return imps;
		
	}
}
