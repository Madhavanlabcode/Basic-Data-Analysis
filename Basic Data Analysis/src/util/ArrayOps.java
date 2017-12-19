package util;

import java.util.ArrayList;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import util.fourier.FFTOps;
import flanagan.analysis.Regression;


public class ArrayOps {

	public static double mean(double[] x)
	{
		double mean = 0;
		for (int i = 0; i < x.length; i++)
			mean += x[i];
		return mean/x.length;
	}
	public static double mean(double[] x, int[] includedIs)
	{
		double mean = 0;
		for (int i = 0; i < includedIs.length; i++)
			mean += x[includedIs[i]];
		return mean/includedIs.length;
	}
	//this returns two numbers: the mean below and above the cutoff;
	public static double[] meanSplit(double[] x, double cutoff)
	{
		double[] ans = new double [2];
		int[] n = new int [2];
		for (int i = 0; i < x.length; i++)
			if (x[i] < cutoff) {ans[0] += x[i]; n[0]++;}
			else {ans[1] += x[i]; n[1]++;}
		for (int i = 0; i < 2; i++)
			ans[i]/=n[i];
		return ans;
	}
	public static double mean(double[][] x)
	{
		double mean = 0;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				mean += x[i][j];
		return mean/(x.length*x[0].length);
	}
	public static boolean contains(int[] nums, int value)
	{
		boolean contains = false;
		for (int i = 0; i < nums.length && !contains; i++)
			contains = value == nums[i];
		return contains;
	}
	public static boolean contains(double[] nums, double value)
	{
		boolean contains = false;
		for (int i = 0; i < nums.length && !contains; i++)
			contains = value == nums[i];
		return contains;
	}
	
	//returns the mean of all elements in the array except the ones which are zero.
	public static double meanExceptZero(double[][] x)
	{
		double mean = 0;
		int n = 0;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				if (x[i][j] != 0){
					mean += x[i][j];
					n++;
				}
		return mean/n;
	}
	public static double meanExceptValue(double[][] x, double value)
	{
		double mean = 0;
		int n = 0;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				if (x[i][j] != value){
					mean += x[i][j];
					n++;
				}
		return mean/n;
	}
	public static double meanInsideBox(double[][] x, int xi, int xf, int yi, int yf)
	{
		double mean = 0;
		int n = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
					mean += x[i][j];
					n++;
				}
		return mean/n;
	}
	
	//mask ranges from zero to one
	public static double meanWithoutMask(double[][] data, double[][] mask) {
		double sum = 0;
		double denominator = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++){
				sum += (1-mask[i][j])*data[i][j];
				denominator += mask[i][j];
			}
		return sum / ((data.length*data[0].length) - denominator);
	}

	public static double meanExceptValueString(double[][] x, String value)
	{
		double mean = 0;
		int n = 0;
		String s;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			{
				s = "" + x[i][j];
				if (!s.equals(value)){
					mean += x[i][j];
					n++;
				}
			}
		return mean/n;
	}
	public static double meanBeforeIndex(double[][] x, int ind)
	{
		double mean = 0;
		int n = 0;
		String s;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			{
				if (n++ < ind)
					mean += x[i][j];
				
			}
		return mean/n;
	}
	public static double sigma(double[] x)
	{
		return sigma(x, mean(x));
	}
	public static double sigma(double[] x, double mean)
	{
		return Math.sqrt(variance(x, mean));
	}
	public static double variance(double[] x)
	{
		return variance(x, mean(x));
	}
	public static double variance(double[] x, double mean)
	{
		double sigma = 0;
		for (int i = 0; i < x.length; i++)
			sigma += Math.pow(x[i] - mean, 2);
		return sigma/x.length;
	}
	
	public static double max(double[] x)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			max = Math.max(x[i], max);
		return max;
	}
	public static int maxIndex(double[] x)
	{
		double max = -Double.MAX_VALUE;
		int index = -1;
		for (int i = 0; i < x.length; i++){
			max = Math.max(x[i], max);
			if (max == x[i])
				index = i;
		}
		return index;
	}
	public static double centIndex(double[] x)
	{
		double cent = 0;
		double sum = 0;
		for (int i = 0; i < x.length; i++){
			cent += Math.abs(i*x[i]);
			sum += Math.abs(x[i]);
		}
		return cent/sum;
	}
	public static double sigmaIndex(double[] x)
	{
		double sig = 0;
		double sum = 0;
		for (int i = 0; i < x.length; i++){
			sig += i*i*x[i];
			sum += x[i];
		}
		
		double cent = centIndex(x);
		
		return (Math.sqrt(sig/sum - cent*cent));

	}
	
	public static double indexToValue(double[] x, double i)
	{
		return ((x[x.length-1]-x[0])/(x.length-1))*i + x[0];
	}
	public static int maxIndex(double[] x, int bmin, int bmax)
	{
		double max = -Double.MAX_VALUE;
		int index = -1;
		for (int i = bmin; i < bmax; i++){
			max = Math.max(x[i], max);
			if (max == x[i])
				index = i;
		}
		return index;
	}
	public static double min(double[] x)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			min = Math.min(x[i], min);
		return min;
	}
	public static int minIndex(double[] x)
	{
		double min = Double.MAX_VALUE;
		int index = -1;
		for (int i = 0; i < x.length; i++){
			min = Math.min(x[i], min);
			if (min == x[i])
				index = i;
		}
		return index;
	}

	//The quicksort code below (5 methods) was taken from http://www.cs.princeton.edu/introcs/42sort/QuickSort.java.html
    public static void quicksort(double[] a) {
        shuffle(a);                        // to guard against worst-case
        quicksort(a, 0, a.length - 1);
    }

    // quicksort a[left] to a[right]
    public static void quicksort(double[] a, int left, int right) {
        if (right <= left) return;
        int i = partition(a, left, right);
        quicksort(a, left, i-1);
        quicksort(a, i+1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partition(double[] a, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (a[++i] < a[right])      // find item on left to swap
                ;                               // a[right] acts as sentinel
            while (a[right] < a[--j])     // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, i, j);                      // swap two elements into place
        }
        exch(a, i, right);                      // swap with partition element
        return i;
    }

    // is x < y ?

    // exchange a[i] and a[j]
    private static void exch(double[] a, int i, int j) {
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }

    // shuffle the array a[]
    private static void shuffle(Comparable[] a) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N-i));   // between i and N-1
            exch(a, i, r);
        }
	}
	//The quicksort code below (5 methods) was taken from http://www.cs.princeton.edu/introcs/42sort/QuickSort.java.html
    public static void quicksort(Comparable[] a) {
        shuffle(a);                        // to guard against worst-case
        quicksort(a, 0, a.length - 1);
    }

    // quicksort a[left] to a[right]
    public static void quicksort(Comparable[] a, int left, int right) {
        if (right <= left) return;
        int i = partition(a, left, right);
        quicksort(a, left, i-1);
        quicksort(a, i+1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partition(Comparable[] a, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (a[++i].compareTo(a[right]) < 0)      // find item on left to swap
                ;                               // a[right] acts as sentinel
            while (a[right].compareTo(a[--j]) < 0)     // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, i, j);                      // swap two elements into place
        }
        exch(a, i, right);                      // swap with partition element
        return i;
    }
    
    // is x < y ?

    // exchange a[i] and a[j]
    private static void exch(Comparable[] a, int i, int j) {
        Comparable swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }

    // shuffle the array a[]
    private static void shuffle(double[] a) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N-i));   // between i and N-1
            exch(a, i, r);
        }
	}

    public static double max(double[][] x)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				if (x[i][j] != Double.NaN) max = Math.max(x[i][j], max);
		return max;
	}
    public static double max(double[][][] x)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				for (int k = 0; k < x[i][j].length; k++)
					max = Math.max(x[i][j][k], max);
		return max;
	}
    public static double max(double[][] x, int index, boolean firstIndexSummedOver)
	{
		double max = -Double.MAX_VALUE;
		if (!firstIndexSummedOver)
			for (int i = 0; i < x[0].length; i++)
				max = Math.max(x[index][i], max);
		else
			for (int i = 0; i < x.length; i++)
				max = Math.max(x[i][index], max);
		return max;
	}
    public static double max(double[][][] x, int thirdindex)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			max = Math.max(x[i][j][thirdindex], max);
		return max;
	}
	public static double min(double[][] x)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				if (x[i][j] != Double.NaN) 
					min = Math.min(x[i][j], min);
		return min;
	}
	public static double min(double[][][] x, int thirdindex)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			min = Math.min(x[i][j][thirdindex], min);
		return min;
	}
	public static double min(double[][][] x)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				for (int k = 0; k < x[i][j].length; k++)
					min = Math.min(x[i][j][k], min);
		return min;
	}
    public static int max(int[][] x)
	{
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			max = Math.max(x[i][j], max);
		return max;
	}
	public static int min(int[][] x)
	{
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			min = Math.min(x[i][j], min);
		return min;
	}
	public static double sum(double[] x)
	{
		double sum = 0;
		for (int i = 0; i < x.length; i++)
			sum += x[i];
		return sum;
	}
	public static int sum(int[] x)
	{
		int sum = 0;
		for (int i = 0; i < x.length; i++)
			sum += x[i];
		return sum;
	}
	public static double sum(double[][] x)
	{
		double sum = 0;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			sum += x[i][j];
		return sum;
	}
	public static double sum(double[][][] x)
	{
		double sum = 0;
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				for (int k = 0; k < x[i][j].length; k++)
				sum += x[i][j][k];
		return sum;
	}
	
	public static void star(double[][] source, double[][] target)
	{
		for (int i = 0; i < source.length; i++)
		{
			target[i][0] = source[i][0];
			target[i][1] = -source[i][1];
		}
	}
	
	public static void sqrt(double[] f)
	{
		int N = f.length;
		for (int i = 0; i < N; i++)
			f[i] = Math.sqrt(f[i]);
	}
	public static void sqrt(double[] source, double[] target)
	{
		int N = source.length;
		for (int i = 0; i < N; i++)
			target[i] = Math.sqrt(source[i]);
	}
	
	public static void magnitude(double[][] f, double[] target)
	{
		int N = f.length;
		for (int i = 0; i < N; i++)
			target[i] = Math.sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1]);
	}
	public static void magnitudeSq(double[][] f, double[] target)
	{
		int N = f.length;
		for (int i = 0; i < N; i++)
			target[i] = f[i][0]*f[i][0] + f[i][1]*f[i][1];
	}
	public static void phase(double[][] f, double[] target)
	{
		int N = f.length;
		for (int i = 0; i < N; i++)
			target[i] = FieldOps.atan(f[i][0], f[i][1]);
	}
	public static void flipY(double[][] x)
	{
		double[][] y = new double[x.length][x[0].length];
		for (int i = 0; i < x.length; i++)
		
			for (int j = 0; j < x[i].length; j++)
				y[i][j] = x[i][x[0].length-1-j];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				x[i][j] = y[i][j];
	}
	public static void flipX(double[][] x)
	{
		double[][] y = new double[x.length][x[0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				y[i][j] = x[x.length-1-i][j];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				x[i][j] = y[i][j];
	}
	public static void flip(double[] x)
	{
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = x[x.length-1-i];
		for (int i = 0; i < x.length; i++)
			x[i] = y[i];
	}
	public static double[] getFlip(double[] x)
	{
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = x[x.length-1-i];
		return y;
	}
	
	public static double[] generateArray(double min, double d, int npts)
	{
		double [] a = new double [npts];
		for (int i = 0; i < a.length; i++)
			a[i] = min + i*d;
		return a;
	}

	//This creates an array from min almost up to max according to the most obvious procedure.
	public static double[] generateArrayNotInclUpper(double min, double max, int npts)
	{
		double dm = (max-min)/npts;
		double [] a = new double [npts];
		for (int i = 0; i < a.length; i++)
			a[i] = min + i*dm;
		return a;
	}
	public static double[] generateArrayNotInclLower(double min, double max, int npts)
	{
		double dm = (max-min)/npts;
		double [] a = new double [npts];
		for (int i = 0; i < a.length; i++)
			a[i] = min + (i+1)*dm;
		return a;
	}
	public static double[] generateArrayInclBoth(double min, double max, int npts)
	{
		if (npts == 1) return new double[] {min};
		double dm = (max-min)/(npts-1);
		double [] a = new double [npts];
		for (int i = 0; i < a.length; i++)
			a[i] = min + i*dm;
		return a;
	}
	public static double[] getArrayFrom2ndIndex(int j, double[][] arg)
	{
		double[] ans = new double [arg.length];
		for (int i = 0; i < arg.length; i++)
			ans[i] = arg[i][j];
		return ans;
	}
	public static void putArrayNotInclUpper(double min, double max, int npts, double[] a)
	{
		double dm = (max-min)/npts;
		for (int i = 0; i < a.length; i++)
			a[i] = min + i*dm;
	}
	public static void putArrayNotInclLower(double min, double max, int npts, double[] a)
	{
		double dm = (max-min)/npts;
		for (int i = 0; i < a.length; i++)
			a[i] = min + (i+1)*dm;
	}
	public static void putArrayInclBoth(double min, double max, int npts, double[] a)
	{
		double dm = (max-min)/(npts-1);
		for (int i = 0; i < a.length; i++)
			a[i] = min + i*dm;
	}
	public static int[] getHistogram(double[] values, double[] binmins)
	{
		int[] hist = new int[binmins.length];
		for (int i = 0; i < values.length; i++){
			for (int j = 0; j < binmins.length-1; j++)
			{
				if (values[i] >= binmins[j] && values[i] < binmins[j+1])
					hist[j]++;
			}
			if (values[i] > binmins[binmins.length-1])
				hist[binmins.length-1]++;
		}
		return hist;
	}
	/**
	 * This one bins an int array into bins of size one. [0] is the bin values, [1] is their populations.
	 * @param values
	 * @return
	 */
	public static int[][] getHistogram(int[] values)
	{
		int max = ArrayOps.max(new int[][] {values});
		int min = ArrayOps.min(new int[][] {values});
		int delta = max - min;
		int[] hist = new int[delta+1];
		int[] bins = new int[delta+1];
		for (int i = min; i <= max; i++)
			bins[i-min] = i;
		
		for (int i = 0; i < values.length; i++)
			for (int j = 0; j < bins.length; j++)
				if (values[i] == bins[j]) hist[j]++;

		return new int[][] {bins, hist};
	}
	/**
	 * Values is [n][2].
	 * @param values
	 * @param binminsX
	 * @param binminsY
	 * @return
	 */
	public static double[][] get2DHistogram(double[][] values, double[] binminsX, double[] binminsY)
	{
		double[][] hist = new double[binminsX.length][binminsY.length];
		for (int i = 0; i < values.length; i++){
			for (int k = 0; k < binminsX.length-1; k++)
				for (int j = 0; j < binminsY.length-1; j++)
				{
					if (values[i][0] >= binminsX[k] && values[i][0] < binminsX[k+1] &&
							values[i][1] >= binminsY[j] && values[i][1] < binminsY[j+1])
						hist[k][j]++;
				}
		
			if (values[i][1] > binminsY[binminsY.length-1])
				for (int j = 0; j < binminsX.length-1; j++)
					if (values[i][0] >= binminsX[j] && values[i][0] < binminsX[j+1])
						hist[j][binminsY.length-1]++;
			if (values[i][0] > binminsX[binminsX.length-1])
				for (int j = 0; j < binminsY.length-1; j++)
					if (values[i][0] >= binminsY[j] && values[i][0] < binminsY[j+1])
						hist[binminsX.length-1][j]++;
			
			if (values[i][0] > binminsY[binminsY.length-1] && values[i][0] > binminsX[binminsX.length-1])
				hist[binminsX.length-1][binminsY.length-1]++;
		}

		
		return hist;
	}
	public static double[] copyExcept(double[] data, int[] indicesNotToCopy)
	{
		double[] ans = new double [data.length - indicesNotToCopy.length];
		int nskipped = 0;
		for (int i = 0; i < data.length; i++)
		{
			if (isInSet(i, indicesNotToCopy))
				nskipped++;
			else
				ans[i-nskipped] = data[i];
		}
		return ans;
	}
	public static double[] copyWithin(double[] data, int imin, int imax)
	{
		double[] ans = new double [imax - imin];
		for (int i = imin; i < imax; i++)
		{
			ans[i-imin] = data[i];
		}
		return ans;
	}
	public static double[][] copyExcept(double[][] data, int[] indicesNotToCopy)
	{
		double[][] ans = new double [data.length - indicesNotToCopy.length][data[0].length];
		int nskipped = 0;
		for (int i = 0; i < data.length; i++)
		{
			if (isInSet(i, indicesNotToCopy))
				nskipped++;
			else
				ans[i-nskipped] = data[i];
		}
		return ans;
	}
	public static Object[] copyExcept(Object[] data, int[] indicesNotToCopy)
	{
		Object[] ans = new Object [data.length - indicesNotToCopy.length];
		int nskipped = 0;
		for (int i = 0; i < data.length; i++)
		{
			if (isInSet(i, indicesNotToCopy))
				nskipped++;
			else
				ans[i-nskipped] = data[i];
		}
		return ans;
	}
	
	public static boolean isInSet(int i, int[] set)
	{
		boolean is = false;
		for (int j = 0; j < set.length && !is; j++)
			if (i == set[j])
				is = true;
		return is;
	}
	
	public static int indexOf(double[] sorted, double value, boolean arrayIncreasing)
	{
		if (arrayIncreasing){
			int i = 0;
			for (i = 0; i < sorted.length; i++)
				if (sorted[i] >= value) return i;
			return i;
		}
		else
		{
			int i = 0;
			for (i = 0; i < sorted.length; i++)
				if (sorted[i] <= value) return i;
			return i;
		}
	}
	
	public static ArrayList copyExcept(ArrayList list, int[] excludedI)
	{
		ArrayList<Object> copy = new ArrayList();
		int nskipped = 0;
		for (int i = 0; i < list.size(); i++)
		{
			if (isInSet(i, excludedI))
				nskipped++;
			else
				copy.add(list.get(i));
		}
		return copy;
		
	}
	public static ArrayList copyOnly(ArrayList list, int[] includedI)
	{
		ArrayList<Object> copy = new ArrayList();
		int nskipped = 0;
		for (int i = 0; i < list.size(); i++)
		{
			if (!isInSet(i, includedI))
				nskipped++;
			else
				copy.add(list.get(i));
		}
		return copy;
	}
	
	public static double[] getPercentiles(double[] unsorted)
	{
		double[] perc = new double[101];
		double[] sorted = unsorted.clone();
		quicksort(sorted);
		
		int n = 0;
		for (int i = 0; i < 100; i++)
		{
			n = FieldOps.round(sorted.length*i/100);
			perc[i] = sorted[n];
		}
		perc[100] = sorted[sorted.length-1];
		return perc;
	}
	public static double[] toDouble(int[] x)
	{
		double[] y = new double [x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = x[i];
		return y;
	}
	public static double[][] toDouble(int[][] x)
	{
		double[][] y = new double [x.length][x[0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[0].length; j++)
			y[i][j] = x[i][j];
		return y;
	}
	public static double[][] toDouble(boolean[][] x)
	{
		double[][] y = new double [x.length][x[0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[0].length; j++)
			y[i][j] = x[i][j] ? 1 : 0;
		return y;
	}
	public static void putDouble(boolean[][] x, double[][] y)
	{
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[0].length; j++)
				y[i][j] = x[i][j] ? 1 : 0;
	}
	public static double[][][] toDouble(boolean[][][] x)
	{
		double[][][] y = new double [x.length][x[0].length][x[0][0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				for (int k = 0; k < x[i][j].length; k++)
					y[i][j][k] = x[i][j][k] ? 1 : 0;
				
		return y;
	}
	
	//This
	public static int[] getHistogramExcludeValue(double[] values, double[] binmins, double exc) {
		int[] hist = new int[binmins.length];
		for (int i = 0; i < values.length; i++){
			for (int j = 0; j < binmins.length-1; j++)
			{
				if (values[i] >= binmins[j] && values[i] < binmins[j+1] && values[i] != exc)
					hist[j]++;
			}
			if (values[i] > binmins[binmins.length-1])
				hist[binmins.length-1]++;
		}
		return hist;
	}
	
	/**
	 * Returns the array {a, b, c, error} where y = ax^2 +bx +c
	 * @param x
	 * @param y
	 * @return
	 */
	public static double[] fitToParabola(double[] x, double[] y)
	{
		double sy = ArrayOps.sum(y);
		double sx = sum(x);
		double sxx = sumPower(x, 2);
		double sxxx = sumPower(x, 3);
		double sxxxx = sumPower(x, 4);
		
		int n = x.length;
		
		double sxy = sumPowerTwo(x, y, 1, 1);
		double sxxy = sumPowerTwo(x, y, 2, 1);

		double[][] A = {
				{n, sx, sxx},
				{sx, sxx, sxxx},
				{sxx, sxxx, sxxxx}};
		double[] B = {sy,sxy,sxxy};		
		//use Apahce linear solver:
		RealMatrix ac = new Array2DRowRealMatrix(A, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(B, false);

		RealVector solution = null;
		solution = solver.solve(bs);
		double[] answer = solution.toArray();
		double error = 0;
		for (int i = 0; i < n; i++)
			error += Math.pow(answer[0]*x[i]*x[i] + answer[1]*x[i] + answer[2] - y[i], 2);
		
//		System.out.print(Printer.arrayLnHorizontal(answer));
		
		return new double[] {answer[2], answer[1], answer[0], error/n};
	}
	
	public static double sumPower(double[] data, int power)
	{
		double sum = 0;
		for (int i = 0; i < data.length; i++)
			sum += Math.pow(data[i], power);
		return sum;
	}
	public static double sumPowerTwo(double[] one, double[] two, int powOne, int powTwo)
	{
		double sum = 0;
		for (int i = 0; i < one.length; i++)
			sum += Math.pow(one[i], powOne)*Math.pow(two[i], powTwo);
		return sum;
	}
	
	/**
	 * This assumes that no more than one period needs to be added between any two points.
	 * 
	 * @param phase
	 * @param period
	 */
	public static void renderPhasesContinuous(double[] phase, double period)
	{
		for (int i = 1; i < phase.length; i++)
		{
			if(Math.abs(phase[i]-phase[i-1]) > Math.abs(phase[i]+period-phase[i-1])) phase[i]+= period;
			if(Math.abs(phase[i]-phase[i-1]) > Math.abs(phase[i]-period-phase[i-1])) phase[i]-= period;
		}
	}
	public static double[] add(double[] y, double dh) {
		double[] ans = new double [y.length];
		for (int i = 0; i < ans.length; i++)
		{
			ans[i] = y[i]+dh;
		}
		return ans;
	}
	public static void addEquals(double[] y, double dh) {
		for (int i = 0; i < y.length; i++)
			y[i]+=dh;
	}
	public static double[] clone(double d, int nx) {
		double[] a = new double [nx];
		for (int i = 0; i < nx; i++)
			a[i] = d;
		return a;
	}
	
	public static double[] getGaussDecayingSine(double k, double L, double nL)
	{
		int length = FieldOps.round(L*nL);
		double[] ans = new double [length];
		for (int i = 0; i < length; i++)
			ans[i] = Math.cos(k*i)*Math.exp(-((double)i*i)/(L*L));
		return ans;
	}
	public static double[] getCosine(double k, int length)
	{
		double[] ans = new double [length];
		for (int i = 0; i < length; i++)
			ans[i] = Math.cos(k*i);
		return ans;
	}
	public static double getProductWithEvenArray(double[] longArray, double[] shortArray, int offset)
	{
		double sum = longArray[offset]*shortArray[0];
		double sumShort = 1;
		for (int i = 1; i < shortArray.length; i++)
		{
			if (offset+i < longArray.length){
				sum += longArray[offset+i]*shortArray[i]; sumShort += shortArray[i];}
			if (offset-i > 0){
				sum += longArray[offset-i]*shortArray[i]; sumShort += shortArray[i];}
		}
		return sum/sumShort;
	}
	public static double[] getFilteredWithEvenArray(double[] longArray, double[] shortArray, int offset)
	{
		double[] ans = new double[longArray.length];
		ans[offset] = longArray[offset];
		for (int i = 1; i < shortArray.length; i++)
		{
			if (offset+i < longArray.length)
				ans[offset+i] = longArray[offset+i]*shortArray[i];
			if (offset-i > 0)
				ans[offset-i] = longArray[offset-i]*shortArray[i];
		}
		return ans;
	}
	/**
	 * 
	 * @param z
	 * @param include
	 * @param power
	 * @return
	 */
	public static double[] fitToNPolynomial(double[] x, double[] y, int power)
	{
		

		double[] paramsNoZ = new double [2*power+1];
		double[] paramsZ = new double [power+1];
		
		for (int i = 0; i < y.length; i++)
		{	
			for (int p = 0; p < 2*power+1; p++)
				paramsNoZ[p] += Math.pow(x[i], p);
				
			for (int p = 0; p < power+1; p++)
				paramsZ[p] += y[i]*Math.pow(x[i], p);
		}
				
		double[][] matrix = new double [power+1][power+1];
		for (int i = 0; i < power+1; i++)
			for (int j = 0; j < power+1; j++)
				matrix[i][j] = paramsNoZ[i+j];

		//use Apahce linear solver:
		RealMatrix ac = new Array2DRowRealMatrix(matrix, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(paramsZ, false);

		RealVector solution = null;
		try{solution = solver.solve(bs);
		}
		catch (Exception e)
		{
//			System.out.print("error! "); //return Flanagan's solver.
			return performPolynomialFit(y, x, power);
		}
		double[] answer = solution.toArray();
//		double[] otherAnswer =  performPolynomialFit(y, x, power);
		return answer;
		//		Printer.printlnHorizontal(answer);
//		return answer;
	}
	public static double[][] subtractPolynomialFit(double[] x, double[] y, int power)
	{
		double[] a = fitToNPolynomial(x, y, power);
		double[] fit = new double[y.length];
		double[] ans = new double [y.length];
		for (int i = 0; i < ans.length; i++)
		{
			fit[i] = 0.999*a[0];
			for (int j = 1; j < power+1; j++)
					fit[i] += a[j]*Math.pow(x[i], j);
			ans[i] = y[i] - fit[i];
		}		
		return new double[][] {ans, fit};	
	}
	public static double[][] subtractLineForPeriodicness(double[] y)
	{
		double[] fit = new double[y.length];
		double[] ans = new double [y.length];
		double dy = y[y.length-1] - y[0];
		double dx = y.length-1;
		double m = dy/dx;
		for (int i = 0; i < ans.length; i++)
		{
			fit[i] = m*i + y[0];
			ans[i] = y[i] - fit[i];
		}		
		return new double[][] {ans, fit};	
	}
	public static void multiply(double[][] x, double d) {
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
			x[i][j] *= d;
	}
	public static void multiply(double[] x, double d) {
		for (int i = 0; i < x.length; i++)
			x[i] *= d;
	}
	public static double[] toArrayLineByLine(double[][] data)
	{
		double [] answer = new double [data.length*data[0].length];
		int n = 0;
		for (int j = 0; j < data[0].length; j++)
			for (int i = 0; i < data.length; i++)
				answer[n++] = data[i][j];
		return answer;
	}
	//New x is supposed given in pixel units.
	public static double[] getInterpolatedFourier(double[] y, double[] newX){
		double[] q = new double [y.length];
		double[][] expqr = new double [y.length][2];
		for (int i = 0; i < y.length; i++){
			q[i] = 2*Math.PI*i/y.length;
			if (q[i] >= Math.PI) q[i] -= 2*Math.PI;
		}
		double[] tempZ = new double [2];
		double[] tempFFTZ = new double[2];
		double[] result = new double [newX.length];
//		double[][][] tempNonsense = new double [2*t.nlayers][t.nx][t.ny];
		
		double[] fftz = FFTOps.get1DFFTComplex(y);
		for (int i = 0; i < result.length; i++){
			for (int k = 0; k < q.length; k++)
				{
					expqr[k][0] = Math.cos(q[k]*newX[i]);
					expqr[k][1] = Math.sin(q[k]*newX[i]);
					tempFFTZ[0] = fftz[2*k];
					tempFFTZ[1] = fftz[2*k+1];
					Complex.product(expqr[k], tempFFTZ, tempZ);
//					fftz[2*k] = tempZ[0];
//					fftz[2*k+1] = tempZ[1];
					//Now I have multiplied the FFT by exp(iqr). I must sum upon q:
					result[i] += tempZ[0];
//					tempNonsense[2*k][i][j] = fftz[2*k];
//					tempNonsense[2*k+1][i][j] = fftz[2*k+1];
				}
			result[i] /= q.length;
			}
		
		return result;

	}
	public static double[] gaussSmooth(double[] data, double L)
	{
		int Lp = (int)(L+1);
		int nlayers = data.length;
		double[] gauss = new double [nlayers];
		double[] smooth = new double [nlayers];
			for (int l = 0; l < nlayers; l++)
				gauss[l] = Math.exp(-(l*l)/(2*L*L));
		
		double gsum;
		int x, y;
		int xprime, yprime, xpmin, xpmax, ypmin, ypmax;
		for (int j = 0; j < nlayers; j++)
		{
			gsum = 0;
			smooth[j] = 0;
			x = j;
			//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
			xpmin = Math.max(0, x - 6*Lp);
			xpmax = Math.min(nlayers, x + 6*Lp);
			for (xprime = xpmin; xprime < xpmax; xprime++)
			{
				gsum += gauss[Math.abs(x-xprime)];
				smooth[j] += gauss[Math.abs(x-xprime)]*data[xprime];
			}
			smooth[j] /= gsum;
//			System.out.println("" + j + "\t" + gsum);
		}
		return smooth;

	}
	/**
	 * Returns Flanagan's coefficients a.
	 * @param data
	 * @param v
	 * @param npoly
	 * @return
	 */
	public static double[] performPolynomialFit(double[] data, double[] v, int npoly)
	{
		Regression reg;
		double[] g = new double[data.length];
		double[] a;
		double[] poly = new double [data.length];
		double[] subt = new double [data.length];
		double[][] powers = new double[data.length][npoly+1];
		for (int k = 0; k < data.length; k++)
			for (int i = 0; i < npoly+1; i++)
				powers[k][i] = Math.pow(v[k], i);
			
		
		reg = new Regression(v,data);
		reg.polynomial(npoly);
		a = reg.getBestEstimates();
//					System.out.print(i + ",\t" + j + ",\t");
//					for (int k = 0; k  < a.length; k++)
//						System.out.print(a[k] + ",\t");
//					System.out.println();
		for (int k = 0; k < data.length; k++)
		{
			for (int m = 0; m < npoly+1; m++)
				poly[k] += a[m]*powers[k][m];
			subt[k] = data[k]-poly[k];
		}
		return a;
//		return new double [][] {poly, subt};
	}
	public static double[] getDerivative(double[] f) {
		double[] ans = new double [f.length];
		ans[0] = f[1]-f[0];
		for (int i = 1; i < ans.length - 1; i++)
			ans[i] = (f[i+1]-f[i-1])/2;
		ans[ans.length-1] = f[f.length-1] - f[f.length-2];
		return ans;
	}
	public static void normalizeToRange(double[] bufferd, double newmin, double newmax) {
		// TODO Auto-generated method stub
		double oldmin = min(bufferd);
		double oldmax = max(bufferd);
		
		for (int i = 0; i < bufferd.length; i++)
		{
			double position = (bufferd[i] - oldmin)/(oldmax-oldmin);
			bufferd[i] = (1 - position)*newmin + position*newmax;
		}
		
		
	}
	public static double[] concatenate(double[][] ds) {
		int length = 0;
		for (int i = 0; i < ds.length; i++)
			length += ds[i].length;
		double[] ans = new double [length];
		int n = 0;
		for (int i = 0; i < ds.length; i++)
			for (int j = 0; j < ds[i].length; j++)
				ans[n++] = ds[i][j];
		return ans;
	}

	public static double[][] subtract(double y, double[][] z) {
		double[][] ans = new double[z.length][z[0].length];
		for (int i=0; i < ans.length; i++)
			for(int j=0;j<ans[0].length;j++)
				ans[i][j] = y-z[i][j];
		return ans;
	}
}
