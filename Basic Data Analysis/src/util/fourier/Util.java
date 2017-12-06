package util.fourier;

public class Util {

	public static int fHat2tofHat(int i, int N)
	{
		return ((i - N/2) + N) % N;
	}
	public static int reltofHat(int i, int N)
	{
		return (i + N) % N;
	}
	public static int fHat2toRel(int i, int N)
	{
		return i - N/2;
	}
}
