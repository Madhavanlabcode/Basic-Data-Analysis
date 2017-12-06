package util;

import java.util.ArrayList;

public class Sort {
    private static long comparisons = 0;
    private static long exchanges   = 0;

    /***********************************************************************
     *  Bubblesort code, written by VE.
     ***********************************************************************/    
    public static void bubbleSort(double[] A) {
    	 // Sort the array A into increasing order.       
    	int out, in;

        for (out = A.length - 1; out > 1; out--)
          // outer loop (backward)
          for (in = 0; in < out; in++)
            // inner loop (forward)
            if (A[in] > A[in + 1]) // out of order?
              swap(in, in + 1, A); // swap them
      }
    public static void bubbleSort(ArrayList<Double> key, ArrayList<Object> data) {
   	 // Sort the array A into increasing order.       
   	int out, in;

       for (out = key.size() - 1; out > 0; out--)
         // outer loop (backward)
         for (in = 0; in < out; in++)
           // inner loop (forward)
           if (key.get(in).compareTo(key.get(in+1)) > 0) // out of order?
             {
        	   Double temp = key.get(in);
        	   key.set(in, key.get(in+1));
        	   key.set(in+1, temp);
        	   Object tempo = data.get(in);
        	   data.set(in, data.get(in+1));
        	   data.set(in+1, tempo);
             }
     }
//private function needed by bubble sort, performs the swap operation
    
//      private static void swap(int one, int two, ArrayList<Material> m) {
//        Material temp = m.get(one);
//        m.set(one, m.get(two));
//        m.set(two, temp);
//      }
//      public static void bubbleSort(ArrayList<Material> m) {
//     	 // Sort the array A into increasing order.       
//     	int out, in;
//
//         for (out = m.size() - 1; out > 1; out--)
//           // outer loop (backward)
//           for (in = 0; in < out; in++)
//             // inner loop (forward)
//             if (m.get(in).name.compareToIgnoreCase(m.get(in + 1).name) > 0) // out of order?
//               swap(in, in + 1, m); // swap them
//       }
 //private function needed by bubble sort, performs the swap operation
       private static void swap(int one, int two, double[] A) {
         double temp = A[one];
         A[one] = A[two];
         A[two] = temp;
       }

      /***********************************************************************
       *  Insertion Sort code written by VE.
       ***********************************************************************/
    public static void insertionSort(double[] A) {
        // Sort the array A into increasing order.
        
     int itemsSorted; // Number of items that have been sorted so far.

     for (itemsSorted = 1; itemsSorted < A.length; itemsSorted++) {
           // Assume that items A[0], A[1], ... A[itemsSorted-1] 
           // have already been sorted.  Insert A[itemsSorted]
           // into the sorted list.
           
        double temp = A[itemsSorted];  // The item to be inserted.
        int loc = itemsSorted - 1;  // Start at end of list.
        
        while (loc >= 0 && A[loc] > temp) {
           A[loc + 1] = A[loc]; // Bump item from A[loc] up to loc+1.
           loc = loc - 1;       // Go on to next location.
        }
        
        A[loc + 1] = temp; // Put temp in last vacated space.
     }
    }
    public static void insertionSort(String[] A) {
        // Sort the array A into increasing order.
        
     int itemsSorted; // Number of items that have been sorted so far.

     for (itemsSorted = 1; itemsSorted < A.length; itemsSorted++) {
           // Assume that items A[0], A[1], ... A[itemsSorted-1] 
           // have already been sorted.  Insert A[itemsSorted]
           // into the sorted list.
           
         double temp = (double)Integer.parseInt(A[itemsSorted].substring(0, A[itemsSorted].indexOf(',')));  // The item to be inserted.
        String tempString = A[itemsSorted];  // The item to be inserted.
        int loc = itemsSorted - 1;  // Start at end of list.
        
        while (loc >= 0 && (double)Integer.parseInt(A[loc].substring(0, A[loc].indexOf(','))) > temp) {
           A[loc + 1] = A[loc]; // Bump item from A[loc] up to loc+1.
           loc = loc - 1;       // Go on to next location.
        }
        
        A[loc + 1] = tempString; // Put temp in last vacated space.
     }
    }
     
    
    /***********************************************************************
     *  Quicksort code from Sedgewick 7.1, 7.2. - Thanks!
     ***********************************************************************/
    public static void quicksort(double[] a) {
        shuffle(a);                        // to guard against worst-case
        quicksort(a, 0, a.length - 1);
    }
    public static void quicksort(double[] a, int left, int right) {
        if (right <= left) return;
        int i = partition(a, left, right);
        quicksort(a, left, i-1);
        quicksort(a, i+1, right);
    }

    private static int partition(double[] a, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))      // find item on left to swap
                ;                               // a[right] acts as sentinel
            while (less(a[right], a[--j]))      // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, i, j);                      // swap two elements into place
        }
        exch(a, i, right);                      // swap with partition element
        return i;
    }

    // is x < y ?
    private static boolean less(double x, double y) {
        comparisons++;
        return (x < y);
    }

    // exchange a[i] and a[j]
    private static void exch(double[] a, int i, int j) {
        exchanges++;
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }

    // shuffle the array a
    private static void shuffle(double[] a) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N-i));   // between i and N-1
            exch(a, i, r);
        }
    }



    // test client
    public static void main(String[] args) {
    	int N = 200000;

        // generate N random real numbers between 0 and 1
        long start = System.currentTimeMillis();
        double[] a = new double[N];
        for (int i = 0; i < N; i++)
            a[i] = Math.random();
        long stop = System.currentTimeMillis();
        double elapsed = (stop - start) / 1000.0;
        System.out.println("Generating input:  " + elapsed + " seconds");

        // sort them
        start = System.currentTimeMillis();
        bubbleSort(a);
        //insertionSort(a);
        //quicksort(a);
        stop = System.currentTimeMillis();
        elapsed = (stop - start) / 1000.0;
        System.out.println("Bubblesort:   " + elapsed + " seconds");
        //System.out.println("Insertionsort:   " + elapsed + " seconds");
        //System.out.println("Quicksort:   " + elapsed + " seconds");
       
        
    }
}

