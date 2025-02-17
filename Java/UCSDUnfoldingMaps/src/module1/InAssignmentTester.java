package module1;

public class InAssignmentTester
{
	public static void main(String[] args)
	{
		int var1 = 17;
		int var2 = var1 + 1;
		var1 = var2 + 1;
		System.out.println("var1: "+var1+" var2: "+var2);
		
		StringBuilder sb = new StringBuilder("start");
		sb.insert(4,"le");
		System.out.println(sb);
		
		String s = "startnow"; 
		s = s.substring(0,6);//+ "le" + s.substring(4);
		System.out.println(s); 
		
	}
	
	
}