
/**
 * Write a description of CaesarCipher here.
 * 
 * @author (your name) 
 * @version (a version number or a date)
 */
public class CaesarCipher {
	
    public static String encrypt(String input, int key) {
		StringBuilder encrypt = new StringBuilder(input);
		String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		String shifted_alphabet = alphabet.substring(key)+alphabet.substring(0);

		
		for (int i = 0; i<input.length(); i++){
			// get the character to encrypt
			char c  =  input.charAt(i);

			char cUP  =  Character.toUpperCase(c);
			int loc = alphabet.indexOf(cUP);
			if (loc!=-1) {
				// get its location
				
				// choose the corresponding character from the shifted alphabet
				if (Character.isLowerCase(c)) {
					encrypt.setCharAt(i,shifted_alphabet.toLowerCase().charAt(loc));
				}
				else {
					encrypt.setCharAt(i,shifted_alphabet.charAt(loc));

				}
				
			}
		
		
		
		}
		
		
		return encrypt.toString();
	}
	
//	public static void testCeasar() {
//		FileResource fr = new FileResource();
//		String message = fr.asString();
//		String encrypted = encrypt(message,key);
//		System.out.println("key is "+ key +"\n"+ encrypted);
//	}
	
	public static String encryptTwoKeys(String input, int key1, int key2) {
		StringBuilder encrypt = new StringBuilder(input);
		String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		String shifted_alphabet1 = alphabet.substring(key1)+alphabet.substring(0);
		String shifted_alphabet2 = alphabet.substring(key2)+alphabet.substring(0);
		int loc = 0;
		
		for (int i = 0; i<input.length(); i++){
			// get the character to encrypt
			char c  =  input.charAt(i);

			char cUP  =  input.toUpperCase().charAt(i);
			loc = alphabet.indexOf(cUP);
			if (loc!=-1 && (i % 2 ==0)) {
				// get its location
				
				// choose the corresponding character from the shifted alphabet
				if (Character.isLowerCase(c)) {
					encrypt.setCharAt(i,shifted_alphabet1.toLowerCase().charAt(loc));
				}
				else {
					encrypt.setCharAt(i,shifted_alphabet1.charAt(loc));

				}
				
			}
			else {
				loc = alphabet.indexOf(cUP);
				if (loc!=-1 && (i % 2 !=0)) {
					// get its location
					
					// choose the corresponding character from the shifted alphabet
					if (Character.isLowerCase(c)) {
						encrypt.setCharAt(i,shifted_alphabet2.toLowerCase().charAt(loc));
					}
					else {
						encrypt.setCharAt(i,shifted_alphabet2.charAt(loc));

					}
					
			}
		
		
		
		}
		
		}
		
		return encrypt.toString();
	}
	
	
	
	
	
	public static void main (String[] args) {
		System.out.println(encrypt("At noon be in the conference room with your hat on for a surprise party. YELL LOUD!", 8));

	}

}
