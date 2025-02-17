package duke_coursera;

public class WordPlay {

	public static boolean isVowel (char ch) {
		switch(ch) {
		  case 'a':
		    // code block
		    break;
		  case 'A':
		    // code block
			return true;
		  case 'O':  return true;
		  case 'o': return true;
		  case 'U': return true;
		  case 'u': return true;
		  case 'I': return true;
		  case 'i': return true;
		  case 'e': return true;
		  case 'E': return true;

		  default:
		    // code block
			return false;

		}
		return false;
	}
	
	public static String replaceVowels(String phrase, char ch) {
		StringBuilder z = new StringBuilder(phrase);
		for (int i=0; i<phrase.length(); i++) {
			if (isVowel(phrase.charAt(i))) {
				z.setCharAt(i,ch);
			}
		}
		return z.toString();
		
	}
	
	public static String emphasize(String phrase, char ch) {
		String str = Character.toString(ch);
		str =  str.toLowerCase();
		char ch_lower = str.charAt(0);
		StringBuilder z = new StringBuilder(phrase);
		while(phrase.indexOf(ch_lower)>=0) {
			int loc = phrase.toLowerCase().indexOf(ch_lower)+1;
			if (loc % 2 ==0) {
				z.setCharAt(loc-1,'+');
			}
			else {
				z.setCharAt(loc-1,'*');
			}
			phrase = z.toString().toLowerCase();
		}
		return z.toString();
		
	}
	
	
	
	public static void main (String[] args) {
	//	System.out.println(emphasize("dna ctgaaactga",'a'));
	//	System.out.println(emphasize("Marry Bella Abracadabra",'a'));
		String alphabet = "abcdefghijklmnopqrstuvwxyz";
		System.out.print(alphabet.indexOf("f"));

	}
}
