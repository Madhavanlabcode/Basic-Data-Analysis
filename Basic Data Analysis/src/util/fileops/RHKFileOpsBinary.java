package util.fileops;

import java.io.File;
import java.nio.ByteBuffer;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import drawing.LayerViewer;
import util.ArrayOps;
import util.FieldOps;
import util.TopomapUtil;
import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.MultiFileViewer.TimeStampedLayer;
/**
 * This code is intended to reproduce as literally as possible the matlab code
 * in read_sm4.m.
 * @author Dan
 *
 */
public class RHKFileOpsBinary {
	public static String rhkdir = "C:\\Users\\Charlie\\Desktop\\spm data\\Practice\\processed data\\Topomap - Run #488 SnTe on Bi2Te3 90 mins Area 2 0T\\topomap\\";
	public static String[] defaultSufficesSix = new String[] {"didv_r","didv_f","curr_r","curr_f","topo_r","topo_f"};

	/**
	 * Note
	 * @param args
	 */
	public static void main(String[] args)
	{
		Topomap.setStdDir();
//		exportAnEntireDirectoryTopomapStyle();
		JFileChooser fc = new JFileChooser(rhkdir);
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		String fp = FileOps.selectOpen(fc).toString();
		String fo = FileOps.selectSave(fco).toString();
		exportOneFileQuick(fp,fo);
//		EntireRHKFile rhk = loadFile(fp);
//		Layer[] layers = getLayers(fp);
//		for (int i = 0; i < layers.length; i++)
//			LayerViewer.show(layers[i], 1024, true);
	}
	
	public static void exportAnEntireDirectoryTopomapStyle()
	{
		JFileChooser fco = new JFileChooser(Topomap.stddir);
		JFileChooser fci = new JFileChooser(rhkdir);
		
		
		String dir = FileOps.selectDir(fci);
		File[] f = new File(dir).listFiles();
		ArrayOps.quicksort(f);
		//The files should now be in alphabetical order. if they are indexed properly, this is the same as 
		//chronological order.
		String useri = JOptionPane.showInputDialog("As you know, the RHK bias voltages are probably wrong.\r\n"
				+ "If you want, enter the limits of the topomap below, in CHRONOLOGICAL order, in millivolts.");
		String out = FileOps.selectSave(fco).toString();

		EntireRHKFile[] themAll = new EntireRHKFile[f.length];
		for (int i = 0 ; i < f.length; i++)
			themAll[i] = loadFile(f[i].toString());
		
		double[] v = new double [f.length];
		if (useri != "") v = ArrayOps.generateArrayInclBoth(Double.parseDouble(useri.split(",")[0]), Double.parseDouble(useri.split(",")[1]), f.length);
		else
			for (int i = 0; i < f.length; i++)
				v[i] = themAll[i].headers[0].bias;
		
		
		int npages = themAll[0].getNPages();
		String[] suffices = null;
		if (npages == 6) suffices = defaultSufficesSix;
		else
		{
			JOptionPane.showMessageDialog(null, "There are " + ((npages/2)*2) + " pages and we don't recognize their names.\r\nIn the next box, enter them separated by commas.");
			String useri2 = JOptionPane.showInputDialog("Enter them now.");
			suffices = useri2.split(",");
		}
		
		//We go through the pages one by one...
		int nx = themAll[0].theData.x[0].length;
		int ny = themAll[0].theData.y[0].length;
		
		for (int i = 0; i < npages; i++)
		{
			double[][][] data = new double [f.length][nx][ny];
			for (int k = 0; k < f.length; k++)
				data[k] = themAll[k].theData.z[i];
			double[] x = themAll[0].theData.x[i];
			double[] y = themAll[0].theData.y[i];
			Topomap t = new Topomap(data, v, x, y, null);
			if (t.v[0] > t.v[t.nlayers-1])
				TopomapUtil.flipEnergy(t);
			Topomap.writeBIN(t, out + suffices[i] + ".bin");
		}
		
	}
	
	/**
	 * This method operates on the following assumptions:
	 * If the file has an even number of layers, it is a scan. If an odd number of layers,
	 * it contains spectral data.
	 * If the odd number is greater than one, it is a dI/dV map, otherwise it is spectra alone.
	 * If the non-spectral layers are six in number, they are in this order didv_r, didv_f, curr_r, curr_f,topo_r,topo_f.
	 * (For now, if they are not six in number the used will be asked for the suffices.)
	 * 
	 *  If the file is a dI/dV map, the grid of spectra is square.
	 *  
	 *  In spectral data X represents the voltage, and y represents something else. I am not sure how to accurately record the coordinates of each spectrum and it may be necesssary to use ASCII export fr that.
	 * 
	 *
	 */
	public static void exportOneFileQuick(String input, String output)
	{
		EntireRHKFile f = loadFile(input);
		int npages = f.getNPages();
		String[] suffices;
		if (npages == 1)
		{
			PointSpectra ps = new PointSpectra(FieldOps.transpose(f.theData.z[0]), f.theData.x[0], f.theData.y[0], f.theData.y[0]);
			PointSpectra.writeBIN(ps, output);
			return;
		}
		if (npages/2 == 3)
			suffices = new String[] {"didv_r","didv_f","curr_r","curr_f","topo_r","topo_f"};
		else
		{
			JOptionPane.showMessageDialog(null, "There are " + ((npages/2)*2) + " non-spectral pages and we don't recognize their names.\r\nIn the next box, enter them separated by commas.");
			String useri = JOptionPane.showInputDialog("Enter them now.");
			suffices = useri.split(",");
		}
		
		Layer[] layers = getLayersOnly(f);
		for (int i = 0; i < layers.length; i++) Layer.writeBIN(layers[i], output + suffices[i] + ".bin");
		
//		JOptionPane.showMessageDialog(null, "Just export the spectral data yourself, via ASCII."); 
//		if (npages %2 != 0)
//		{
//			int nx = (int)Math.sqrt(f.theData.y[0].length);
//			int nlayers = f.theData.x[0].length;
//			double[][][] mapData = new double [nlayers][nx][nx];
//			double[] v = new double [nlayers][]
//			for (int k = 0; k < nlayers; k++){
//				for (int i = 0; i < nx; i++)
//					for (int j = 0; j < nx; j++)
//						mapData
//			}
//		}
	}
	
	public static Layer[] getLayers(String fp)
	{
		EntireRHKFile f = loadFile(fp);
		int nlayers = f.headers.length;
		Layer[] ans = new Layer[nlayers];
		for (int i = 0; i < nlayers; i++)
		{
			ans[i] = new Layer(f.theData.z[i], f.theData.x[i], f.theData.y[i], f.headers[i].bias, f.headers[i].current);
		}
		return ans;	
	}
	public static Layer[] getLayers(EntireRHKFile f)
	{
		int nlayers = f.headers.length;
		Layer[] ans = new Layer[nlayers];
		for (int i = 0; i < nlayers; i++)
		{
			ans[i] = new Layer(f.theData.z[i], f.theData.x[i], f.theData.y[i], f.headers[i].bias, f.headers[i].current);
		}
		return ans;	
	}
	/**
	 * This only returns the layers.
	 * @param f
	 * @return
	 */
	public static Layer[] getLayersOnly(EntireRHKFile f)
	{
		int nlayers = f.headers.length;
		if (nlayers %2 == 0) return getLayers(f);
		Layer[] ans = new Layer[nlayers-1];
		for (int i = 1; i < nlayers; i++)
		{
			ans[i-1] = new Layer(f.theData.z[i], f.theData.x[i], f.theData.y[i], f.headers[i].bias, f.headers[i].current);
		}
		return ans;	
	}
	
	public static EntireRHKFile loadFile(String filepath)
	{
		ByteBuffer bb = DavisFileOps.getFromFile(filepath);
		RHKFileHeader header = readFileHeader(bb);
		RHKObject[] objectList = new RHKObject[header.objectListCount];
		for (int i = 0; i < objectList.length; i++)
			objectList[i] = readObjects(bb);
		PageIndexHeader pih = readPageIndexHeader(bb);
		RHKObject object = readObjects(bb);
		PageIndex[] pis = new PageIndex[pih.pageCount];
		RHKObject[][] pageIndexArray = new RHKObject[pih.pageCount][4];
		for (int i = 0; i < pih.pageCount; i++){
			pis[i] = readPageIndex(bb);
			for (int j = 0; j < pis[i].objectListCount; j++)
				pageIndexArray[i][j] = readObjects(bb);
		}
		PageHeader[] headers = new PageHeader[pih.pageCount];
		StringData[] sdata = new StringData[pih.pageCount];
		for (int i = 0; i < pih.pageCount; i++)
		{
			int start = pageIndexArray[i][0].offset-1;
			bb.position(start);
			headers[i] = readPageHeader(bb);
			bb.getInt(); //read 4 bytes as in matlab: 207
			for (int j = 0; j < 9; j++)
				readObjects(bb); //the "object list string" is never seen agan in the matlab file.
			sdata[i] = readStringData(bb);
		}
		RHKData theData = readData(bb, pih, headers, pageIndexArray);
//		System.out.println("This data has: " + theData.z.length + " layers.");
		return new EntireRHKFile(header, objectList, pih, object, pis, pageIndexArray,headers, sdata, theData);
	}
	
	public static class EntireRHKFile{
		RHKFileHeader header;
		RHKObject[] objectList;
		PageIndexHeader pih;
		RHKObject ignorethis;
		PageIndex[] pis;
		RHKObject[][] pageIndexArray;
		PageHeader[] headers;
		StringData[] sdata;
		RHKData theData;
		public EntireRHKFile(RHKFileHeader header, RHKObject[] objectList,
				PageIndexHeader pih, RHKObject ignorethis, PageIndex[] pis,
				RHKObject[][] pageIndexArray, PageHeader[] headers,
				StringData[] sdata, RHKData theData) {
			super();
			this.header = header;
			this.objectList = objectList;
			this.pih = pih;
			this.ignorethis = ignorethis;
			this.pis = pis;
			this.pageIndexArray = pageIndexArray;
			this.headers = headers;
			this.sdata = sdata;
			this.theData = theData;
		}
		
		public int getNPages(){
			return headers.length;
		}
	}
	
	public static RHKData readData(ByteBuffer bb, PageIndexHeader pih, PageHeader[] headers, RHKObject[][] pageIndexArray)
	{
		int nlayers = pih.pageCount;
		double[][][] z = new double[nlayers][][];
		double[][] x = new double [nlayers][];
		double[][] y = new double [nlayers][];
		int[] aux; double[] aux2;
		int nx, ny;
		for (int i = 0; i < nlayers; i++)
		{
			bb.position(pageIndexArray[i][1].offset-1);
			aux = new int [pageIndexArray[i][1].size/4];
			aux2 = new double [pageIndexArray[i][1].size/4];
			for (int j = 0; j < aux.length; j++){
				aux[j] = bb.getInt();
				aux2[j] = headers[i].zOffset + aux[j]*headers[i].zScale/256.0;
			}
			nx = headers[i].xy[2];
			ny = headers[i].xy[3];
			z[i] = new double[nx][ny];
			for (int p = 0; p < nx; p++)
				for (int q = 0; q < ny; q++)
					z[i][p][q] = aux2[p + q*nx];
			
			
			x[i] = ArrayOps.generateArray(headers[i].xOffset, headers[i].xScale, nx);
			y[i] = ArrayOps.generateArray(headers[i].yOffset, headers[i].yScale, ny);
			ArrayOps.flipX(z[i]);
		}
		
		return new RHKData(z, x, y);
	}
	public static class RHKData
	{
		//The first index of each array is the unknown PageCount of the PageIndexHeader.
		double[][][] z;
		double[][] x;
		double[][] y;
		public RHKData(double[][][] z, double[][] x, double[][] y) {
			super();
			this.z = z;
			this.x = x;
			this.y = y;
		}
	}
	
	private static StringData readStringData(ByteBuffer bb){
		String[] s = new String [17];
		for (int i = 0; i < 17; i++)
		{
			int length = bb.getShort();
			s[i] = "";
			for (int j = 0; j < length; j++)
				s[i] += bb.getChar();
		}
		return new StringData (s[0],
				s[1], 
				s[2], 
				s[3], 
				s[4], 
				s[5], 
				s[6], 
				s[7], 
				s[8], 
				s[9], 
				s[10], 
				s[11], 
				s[12], 
				s[13], 
				s[14],
				s[15], 
				s[16] );
	}
	
	public static class StringData{

		String label,
			systemText,
			sessionText,
			userText,
			path,
			date,
			time,
			xUnits,
			yUnits,
			zUnits,
			xLabel,
			yLabel,
			statusChannelText,
			completedLineCount,
			overSampleingCount,
			slicedVoltage,
			pLLProStatus;
			
		public StringData(String label, String systemText, String sessionText, String userText,
				String path, String date, String time, String xUnits,
				String yUnits, String zUnits, String xLabel, String yLabel,
				String statusChannelText, String completedLineCount,
				String overSampleingCount, String slicedVoltage,
				String pLLProStatus) {
			super();
			this.label = label;
			this.systemText = systemText;
			this.sessionText = sessionText;
			this.userText = userText;
			this.path = path;
			this.date = date;
			this.time = time;
			this.xUnits = xUnits;
			this.yUnits = yUnits;
			this.zUnits = zUnits;
			this.xLabel = xLabel;
			this.yLabel = yLabel;
			this.statusChannelText = statusChannelText;
			this.completedLineCount = completedLineCount;
			this.overSampleingCount = overSampleingCount;
			this.slicedVoltage = slicedVoltage;
			this.pLLProStatus = pLLProStatus;
		}
	}
	
	private static PageHeader readPageHeader(ByteBuffer bb)
	{
		bb.get(new byte[3]); //skip 3 bytes
		short stringCount = bb.getShort();
		int pageType = bb.getInt();
		
		int dataSubSource = bb.getInt();
		int lineType = bb.getInt();
		int[] xy = new int [4];
		for (int i = 0; i < 4; i++)
			xy[i] = bb.getInt(); //These are x_corner, y_corner, width, height
		int imageType = bb.getInt();
		int scanDir = bb.getInt();
		int groupID = bb.getInt();
		int pageDataSize = bb.getInt();
		int minZValue = bb.getInt();
		int maxZValue = bb.getInt();
		float xScale = bb.getFloat();
		float yScale = bb.getFloat();
		float zScale = bb.getFloat();
		float xyScale = bb.getFloat();
		float xOffset = bb.getFloat();
		float yOffset = bb.getFloat();
		float zOffset = bb.getFloat();
		float period = bb.getFloat();
		float bias = bb.getFloat();
		float current = bb.getFloat();
		float angle = bb.getFloat();
		int colorInfoCount = bb.getInt();
		int gridXSize = bb.getInt();
		int gridYSize = bb.getInt();
		int[] reserved = new int [16];
		for (int i = 0; i < 16; i++)
			reserved[i] = bb.getInt();
		return new PageHeader(stringCount, pageType, dataSubSource,
			 lineType, xy, imageType, scanDir,
			 groupID, pageDataSize, minZValue, maxZValue,
			 xScale, yScale, zScale, xyScale,
			 xOffset, yOffset, zOffset, period,
			 bias, current, angle, colorInfoCount,
			 gridXSize, gridYSize, reserved);
	}
	public static class PageHeader
	{
		short stringCount;
		int pageType;
		
		int dataSubSource;
		int lineType;
		int[] xy;
		int imageType;
		int scanDir;
		int groupID;
		int pageDataSize;
		int minZValue;
		int maxZValue;
		float xScale;
		float yScale;
		float zScale;
		float xyScale;
		float xOffset;
		float yOffset;
		float zOffset;
		float period;
		float bias;
		float current;
		float angle;
		int colorInfoCount;
		int gridXSize;
		int gridYSize;
		int[] reserved;
		public PageHeader(short stringCount, int pageType, int dataSubSource,
				int lineType, int[] xy, int imageType, int scanDir,
				int groupID, int pageDataSize, int minZValue, int maxZValue,
				float xScale, float yScale, float zScale, float xyScale,
				float xOffset, float yOffset, float zOffset, float period,
				float bias, float current, float angle, int colorInfoCount,
				int gridXSize, int gridYSize, int[] reserved) {
			super();
			this.stringCount = stringCount;
			this.pageType = pageType;
			this.dataSubSource = dataSubSource;
			this.lineType = lineType;
			this.xy = xy;
			this.imageType = imageType;
			this.scanDir = scanDir;
			this.groupID = groupID;
			this.pageDataSize = pageDataSize;
			this.minZValue = minZValue;
			this.maxZValue = maxZValue;
			this.xScale = xScale;
			this.yScale = yScale;
			this.zScale = zScale;
			this.xyScale = xyScale;
			this.xOffset = xOffset;
			this.yOffset = yOffset;
			this.zOffset = zOffset;
			this.period = period;
			this.bias = bias;
			this.current = current;
			this.angle = angle;
			this.colorInfoCount = colorInfoCount;
			this.gridXSize = gridXSize;
			this.gridYSize = gridYSize;
			this.reserved = reserved;
		}
	}
	
	private static PageIndex readPageIndex(ByteBuffer bb)
	{
		short[] pageID = new short[8];
		for (int i = 0; i < 8; i++)
			pageID[i] = bb.getShort();
		int pageDataType = bb.getInt();
		int pageSourceType = bb.getInt();
		int objectListCount = bb.getInt();
		int minorv = bb.getInt();
		return new PageIndex(pageID, pageDataType, pageSourceType,
				objectListCount, minorv);
	}
	
	public static class PageIndex
	{
		short[] pageID;
		int pageDataType;
		int pageSourceType;
		int objectListCount;
		int minorv;
		public PageIndex(short[] pageID, int pageDataType, int pageSourceType,
				int objectListCount, int minorv) {
			super();
			this.pageID = pageID;
			this.pageDataType = pageDataType;
			this.pageSourceType = pageSourceType;
			this.objectListCount = objectListCount;
			this.minorv = minorv;
		}
		
		
	}
	private static PageIndexHeader readPageIndexHeader(ByteBuffer bb)
	{
		int pageCount = bb.getInt();
		int objectListCount = bb.getInt();
		int[] reserved = new int [2];
		reserved[0] = bb.getInt();
		reserved[1] = bb.getInt();
		return new PageIndexHeader(pageCount, objectListCount, reserved);
	}
	
	public static class PageIndexHeader{
		int pageCount;
		int objectListCount;
		int[] reserved;
		
		public PageIndexHeader(int pageCount, int objectListCount,
				int[] reserved) {
			super();
			this.pageCount = pageCount;
			this.objectListCount = objectListCount;
			this.reserved = reserved;
		}
	}
	
	private static RHKObject readObjects(ByteBuffer bb)
	{
		int nameCode = bb.getInt();
		int offset = bb.getInt();
		int size = bb.getInt();
		return new RHKObject (nameCode, offset, size);
	}
	
	public static class RHKObject{
		int nameCode;
		int offset;
		int size;

		public RHKObject(int nameCode, int offset, int size) {
			super();
			this.nameCode = nameCode;
			this.offset = offset;
			this.size = size;
		}
	}
	
	public static RHKFileHeader readFileHeader(ByteBuffer bb)
	{
		bb.rewind();
		int headerSize = bb.getShort();
		String signature = "";
		for (int i = 0; i < 18; i++)
			signature += bb.getChar();

		int nPages = bb.getInt();
		int objectListCount = bb.getInt();
		int objectFieldSize = bb.getInt();
		int[] reserved = new int [2];
		reserved[0] = bb.getInt();
		reserved[1] = bb.getInt();
		
		return new RHKFileHeader(headerSize, signature, nPages, objectListCount, objectFieldSize, reserved);
	}
	
	public static class RHKFileHeader
	{
		int headerSize; String signature;
		int nPages;
		int objectListCount;
		int objectFieldSize;
		int[] reserved;

		public RHKFileHeader(int headerSize, String signature, int nPages,
				int objectListCount, int objectFieldSize, int[] reserved) {
			super();
			this.headerSize = headerSize;
			this.signature = signature;
			this.nPages = nPages;
			this.objectListCount = objectListCount;
			this.objectFieldSize = objectFieldSize;
			this.reserved = reserved;
		}
	}
	

}
