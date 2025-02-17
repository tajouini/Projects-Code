 

public class _2TernaryExpression {

    public static void main(String[] args) {

        int x = 1;
        int y = 2;
        int z = x > y ? x++ : y++;
        System.out.println(x);
        System.out.println(y);
    }
}
