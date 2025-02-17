// See https://aka.ms/new-console-template for more information
using System;
using System.Collections.Generic;

// increment operation
public class Increment
{
    public double x;
    public double AddOne()
    {
        return x + 1;
    }
}
// decrement operation
public class Decrement
{
    public double x;
    public double MinusOne()
    {
        return x - 1;
    }
}
// double operation
public class Double
{
    public double x;
    public double MultiplyByTwo()
    {
        return x * 2;
    }
}
// add random operation
public class Randadd
{
    public double x;
    public double ChangeNum()
    {
        Random rand = new Random();
        int min = 1;
        int max = 10;
        return x + rand.Next(min, max + 1);
    }

}
// undo operation
public class Undo
{
    public string last_operation = string.Empty;
    // we use the stack according to LIFO (Last In, First Out) 
    public Stack<double> myStack = new Stack<double>();
    
    public void addToMyStack(double x)  // add new value into the stack
    {
        myStack.Push(x);
    }
    // if we did not reach the end return a value upon an undo operation
    public double getNum() 
    {
        // update the array size
        if (myStack.Count>=2)
        {
           // arlist.RemoveAt(arlist.Count-1);
            myStack.Pop();
           
            return myStack.Peek();
        }
        else throw new ArgumentOutOfRangeException();
                            
                        
       
    }
}

// Main Program
public class Program
{ // Main method
    public static void Main()
    {
        double num;
        string operation;
        Undo myUndo = new Undo();

        Console.Write("Enter an initial value: ");
        // Make sure that the entered number is a digit
        if (double.TryParse(Console.ReadLine(), out num))
        {   
            myUndo.addToMyStack(num);
            // prompt the user for commands in an infinite loop
            while (true)
            {   
                Console.Write("Enter a command: ");
                operation = Console.ReadLine();
                // Operations written in UPPER CASE are allowed
                operation = operation.ToLower();

                if (operation == "increment") // increment the value by one
                {
                    Increment myIncrement = new Increment { x = num };
                    num = myIncrement.AddOne();
                    myUndo.last_operation = "increment";
                    myUndo.addToMyStack(num);
                    Console.WriteLine("current value: "+num);
                }
                else if (operation == "decrement") // decrement the value by one
                {
                    Decrement myDecrement = new Decrement { x = num };
                    myUndo.last_operation = "decrement";
                    num = myDecrement.MinusOne();
                    myUndo.addToMyStack(num);
                    Console.WriteLine("current value: "+num);
                }
                else if (operation == "double") // double the value
                {
                    Double myDouble = new Double { x = num };
                    myUndo.last_operation = "double";
                    num = myDouble.MultiplyByTwo();
                    myUndo.addToMyStack(num);
                    Console.WriteLine("current value: "+num);
                }
                else if (operation == "randadd") // add a random number
                {
                    Randadd myRandadd = new Randadd { x = num };
                    myUndo.last_operation = "randadd";
                    num = myRandadd.ChangeNum();
                    myUndo.addToMyStack(num);
                    Console.WriteLine("current value: "+num);
                }
                else if (operation == "undo") // perform undo operation
                {
                    if (string.IsNullOrEmpty(myUndo.last_operation))
                    {
                        Console.WriteLine("Error: cannot start with undo");
                        break;
                    }
                    else
                    {   
                        try{
                        num = myUndo.getNum();
                        Console.WriteLine("current value: "+num);
                       }
                        catch(System.ArgumentOutOfRangeException){
                        // We reached the end of the undo operations       
                        Console.WriteLine("Error: cannot undo more");
                        }
                        
                    }
                }
                else
                {
                    Console.WriteLine("invalid command");
                }
            }
        }
        else
        {
            Console.WriteLine("invalid value");
        }
    }
}


