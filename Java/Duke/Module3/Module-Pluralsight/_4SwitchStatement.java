 

public class _4SwitchStatement {

    public static void main(String[] args) {
        // correct switch statement
        Month month = Month.OCTOBER;
        switch (month) {
            case DECEMBER:
                System.out.println("Winter");
                break;
            case APRIL:
                System.out.println("Spring");
                break;
            case JUNE:
                System.out.println("Summer");
                break;
            case OCTOBER:
                System.out.println("Fall");
                break;
            default:
                System.out.println("oops!");
                break;
        }

        System.out.println("No break statements");
        Month month2 = Month.JUNE;
        switch (month2) {
            case DECEMBER:
                System.out.println("Winter");
            case APRIL:
                System.out.println("Spring");
            case JUNE:
                System.out.println("Summer");
            case OCTOBER:
                System.out.println("Fall");
            default:
                System.out.println("oops!");
        }
        System.out.println("Multiple cases to a single block");
        Month month3 = Month.JANUARY;
        switch (month3) {
            case DECEMBER:
            case JANUARY:
            case FEBRUARY:
                System.out.println("Winter");
                break;
            case OCTOBER:
                System.out.println("Fall");
                break;
        }
    }

    enum Month {
        DECEMBER, JANUARY, FEBRUARY,
        MARCH, APRIL, MAY,
        JUNE, JULY, AUGUST,
        SEPTEMBER, OCTOBER, NOVEMBER
    }
}
