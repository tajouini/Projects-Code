package guimodule;

public class main {
	// linear search algorithm
	public static String findStringAirpoirtCode(String toFind, Airport[] airports) {
		boolean found = false;
		int i=0;
		while (!found && i<=airports.length){
			found = airports[i].getCity().equals(toFind);
			i ++;
		}
		if (found){
			return airports[i-1].getCode();
		}
		else {
			return null;
		}
		
	}
// sorting algorithm
	public static void selectionSort(int[] vals)  {
		int i = 0;
		int smallest = vals[i];
		int z = i;
	for (i=0; i<=vals.length-2; i++) {	
		// find the smallest element in i to length-1
		for (int j=i+1; j<=vals.length-1; j++) {
			smallest = vals[i];
			if (smallest>vals[j]) {
				smallest = vals[j]; // value
				z = j; // position
			}
			
		}
		// swap positions
		vals[z]= vals[i];
		vals[i] = smallest;
	}
		
	}

	
}
