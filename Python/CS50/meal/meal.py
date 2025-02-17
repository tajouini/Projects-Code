#In meal.py, implement a program that prompts the user for a time and outputs whether it’s breakfast time, lunch time, or dinner time.
# If it’s not time for a meal, don’t output anything at all. Assume that the user’s input will be formatted in 24-hour time as #:## or ##:##.
# And assume that each meal’s time range is inclusive. For instance, whether it’s 7:00, 7:01, 7:59, or 8:00, or anytime in between, it’s time for breakfast.

#Structure your program per the below, wherein convert is a function (that can be called by main) that converts time, a str in 24-hour format, to the corresponding number of hours as a float.
# For instance, given a time like "7:30" (i.e., 7 hours and 30 minutes), convert should return 7.5 (i.e., 7.5 hours).

def main():
    time = input("Enter time in the format ##:## ")
    f_time = convert(time)
    if 7.00<=f_time <= 8.00:
        print("breakfast time")
    elif 12.00<= f_time <= 13.00:
        print("lunch time")
    elif 18.00<= f_time <= 19.00:
        print("dinner time")
    else:
        print("")

def convert(time):  # converts time, a str in 24-hour format, to the corresponding number of hours as a float
    time = time.strip()
    for i in range(len(time)):
        if time[i] == ":":
         hours = time[0:i]
         minutes_f = time[i+1:i+3]
         if minutes_f == "30":
            minutes_f = "5"
         elif minutes_f == "15":
            minutes_f = "25"
         elif minutes_f == "45":
            minutes_f = "75"

    return float(str(hours)+"."+str(minutes_f))

if __name__ == "__main__":
    main()
