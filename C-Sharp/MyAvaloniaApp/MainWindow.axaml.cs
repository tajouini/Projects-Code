using System;
using System.Collections.Generic;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;

namespace MyAvaloniaApp
{
    
    public class Increment  // increment operation
    {
        public double x;
        public double AddOne()
        {
            return x + 1;
        }
    }
// decrement operation
    public class Decrement // decrement operation
    {
        public double x;
        public double MinusOne()
        {
            return x - 1;
        }
    }
// double operation
    public class Double // double operation
    {
        public double x;
        public double MultiplyByTwo()
        {
            return x * 2;
        }
    }
// add random operation
    public class Randadd   // add a random number operation
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
    
    public class Undo  // undo the last operation
    {
        public string last_operation = string.Empty;
        //public ArrayList arlist = new ArrayList(); // recommended
        public Stack<double> myStack = new Stack<double>();
    
        public void addToMyStack(double x)
        {
            myStack.Push(x);
        }
        public double getNum()
        {
            // update the array size
            if (myStack.Count>=2)
            {
                myStack.Pop();
           
                return myStack.Peek();
            }
            else throw new ArgumentOutOfRangeException();
                            
                        
       
        }
    }
    
    
    
    
    public partial class MainWindow : Window   // main class
    {
        bool validInput = false;
        private double num;
        Undo myUndo = new Undo();
        public MainWindow()
        {
            InitializeComponent();
           
        }
        
        
        // Ensure this method's signature matches what Avalonia expects
        public void NumberTextBox_TextInput(object sender, RoutedEventArgs e)
        {
            // Only allow numeric input
            if (double.TryParse(NumberTextBox.Text, out double num3))
            {
                num = num3;
                System.Console.WriteLine("Reading number: " + num); // Use Debug.WriteLine for output in a GUI context
                myUndo.myStack = new Stack<double>();
                myUndo.addToMyStack(num);
                validInput = true;
                ResultTextBlock.Text = "Click a button";
            }
            else
            {
                ResultTextBlock.Text = "This is not a number!";
                e.Handled = true; // Ignore the input if it's not a number
            }
            
                
        }
        
        

        private void IncrementClick(object sender, RoutedEventArgs e)
        {
            System.Console.WriteLine("Button clicked!");
            // Implement your button click logic here
            var button = sender as Button;
            if (button != null && validInput)
            {
                // Perform operations based on button content or any other logic
                Increment myIncrement = new Increment { x = num };
                num = myIncrement.AddOne();
                myUndo.last_operation = "increment";
                myUndo.addToMyStack(num);
                ResultTextBlock.Text = num.ToString();
                Console.WriteLine("current value: "+num);
            }
        }
        
        private void DecrementClick(object sender, RoutedEventArgs e)
        {   
            System.Console.WriteLine("Button clicked!");
            // Implement your button click logic here
            var button = sender as Button;
            if (button != null && validInput)
            {
                // Perform operations based on button content or any other logic
                Decrement myDecrement = new Decrement { x = num };
                myUndo.last_operation = "decrement";
                num = myDecrement.MinusOne();
                myUndo.addToMyStack(num);
                ResultTextBlock.Text = num.ToString();
                Console.WriteLine("current value: "+num);
            }
        }
        
        private void DoubleClick(object sender, RoutedEventArgs e)
        {   
           
            System.Console.WriteLine("Button clicked!");
            // Implement your button click logic here
            var button = sender as Button;
            if (button != null && validInput)
            {
                // Perform operations based on button content or any other logic
                Double myDouble = new Double { x = num };
                myUndo.last_operation = "double";
                num = myDouble.MultiplyByTwo();
                myUndo.addToMyStack(num);
                ResultTextBlock.Text = num.ToString();
                Console.WriteLine("current value: "+num);
            }
        }
        private void RandAddClick(object sender, RoutedEventArgs e)
        {   
            System.Console.WriteLine("Button clicked!");
            // Implement your button click logic here
            var button = sender as Button;
            if (button != null && validInput)
            {
                // Perform operations based on button content or any other logic
                Randadd myRandadd = new Randadd { x = num };
                myUndo.last_operation = "randadd";
                num = myRandadd.ChangeNum();
                myUndo.addToMyStack(num);
                ResultTextBlock.Text = num.ToString();
                Console.WriteLine("current value: "+num);
            }
        }
        
        private void UndoClick(object sender, RoutedEventArgs e)
        {   

            System.Console.WriteLine("Button clicked!");
            // Implement your button click logic here
            var button = sender as Button;
            if (button != null && validInput)
            {
                // Perform operations based on button content or any other logic
                if (string.IsNullOrEmpty(myUndo.last_operation))
                {
                    ResultTextBlock.Text = "Error: cannot start with undo";
                    Console.WriteLine("Error: cannot start with undo");
                }
                else
                {   
                    try{
                        num = myUndo.getNum();
                        ResultTextBlock.Text = num.ToString();
                        Console.WriteLine("current value: "+num);
                    }
                    catch(System.ArgumentOutOfRangeException)
                    {
                        ResultTextBlock.Text = "Error: cannot undo more";    
                        Console.WriteLine("Error: cannot undo more");
                    }
                        
                }
            }
        }
        
  
       
     
      
        
      
    }
}