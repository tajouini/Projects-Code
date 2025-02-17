
/**
 * Write a description of Decrypt here.
 * 
 * @author (your name) 
 * @version (a version number or a date)
 */
public class Decrypt {
    
    
    public String decrypt(String encrypted){
        CaesarCipher cc = new CaesarCipher();

        int[] freqs = countLetters(encrypted);
        int maxDex = maxIndex(freqs);
        int dkey = maxDex-4;
        if (maxDex<4){
            dkey = 26 -(4-maxDex);
        }
        System.out.println("Key is "+dkey);
        return cc.encrypt(encrypted, 26-dkey);
    }

    // public int[] countLetters(String str){
        // int[] freqs = new int[26];
        // String alphabet = "abcdefghijklmnopqrstuvwxyz";
        // for (int k=0; k<str.length(); k++){
            // char car = str.toLowerCase().charAt(k);
            // if (alphabet.indexOf(car)!=-1){
                // freqs[alphabet.indexOf(car)]+=1;
                // System.out.println("alpha index "+alphabet.indexOf(car));
            // }
        // for (int j = 0; j<freqs.length; j++){
           // // System.out.println("length "+freqs.length);

            // System.out.println(alphabet.charAt(j)+" with freq "+freqs[j]);

        // }
    
        
      //  }
     //   return freqs;
     // }
     
      public int[] countLetters(String message){
        String alph = "abcdefghijklmnoqprstuvwxyz";
        int[] counts = new int[26];
        for (int k = 0; k < message.length(); k++){
            char ch = message.toLowerCase().charAt(k);
            int dex = alph.indexOf(ch);
            if (dex != -1){
                counts[dex] += 1;
                
            }
        }
            
            for (int j = 0; j<counts.length; j++){

            System.out.println(alph.charAt(j)+" with freq "+counts[j]);

             }
        return counts;
    }
     
     
     
      public int maxIndex(int[] vals){
        int maxDex = 0;
        for (int k = 0; k < vals.length; k++){
            if (vals[k] > vals[maxDex]){
                maxDex = k;
            }
            
        }
       // System.out.print(maxDex);
        return maxDex;
    }
    
    public String halfOfString(String message,int start){
        String str = "";
        for (int k = start; k< message.length() ; k+= 2) { 
            str = str + message.charAt(k);    
        }
         
        return str;
    }
    
    
    
    public String decryptTwoKeys(String encrypted){
        CaesarCipher cc = new CaesarCipher();

        String str1 = halfOfString(encrypted, 0);
        String str2 = halfOfString(encrypted, 1);
        int key1 = getKey(str1);
        int key2 = getKey(str2);
        System.out.println("Key1 "+Integer.toString(key1)+" Key2 "
        +Integer.toString(key2));
        //key1 = 14;
        //key2 = 24;
        return cc.encryptTwoKeys(encrypted,26-key1,26-key2); 
    }
    
    public int getKey(String str){
        int key = 0;
        int[] freqs = countLetters(str);
        int maxDex = maxIndex(freqs); 
        int dKey = maxDex-4;
        if (maxDex<4){
            dKey = 26 -(4-maxDex);
        }
        return dKey;
    }
    
    
    public void testDecrypt(){
        String str = "Aal uttx hm aal Qtct Fhljha pl Wbdl. Pvxvxlx!";
       // CaesarCipher cc = new CaesarCipher();
       // String encry = cc.encryptTwoKeys(str,21,8);
       // System.out.println(encry);
       System.out.println(decryptTwoKeys(str));
    }
}
