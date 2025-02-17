#In a file called outdated.py,
# implement a program that prompts the user for a date, anno Domini,
# in month-day-year order, formatted like 9/8/1636 or September 8, 1636,
# wherein the month in the latter might be any of the values in the list below:
#Then output that same date in YYYY-MM-DD format. If the user’s input is not a valid date in either format, prompt the user again.
# Assume that every month has no more than 31 days; no need to validate whether a month has 28, 29, 30, or 31 days.


Months = ["January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December"]
Numbers_months = [ "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
Numbers_days = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", str(range(13,32))]

dict1 = dict(zip(Months, Numbers_months))



while True:
    count1 = 0;
    count2 = 0;
    try:
        user_date = input("Date: ")
        user_date = user_date.strip()
        for d in range(len(user_date)):   # September 8, 1999
            if  user_date[d] == " " and count1 == 0:
                count1 += 1
                month = user_date[0:d]
                if user_date[d+2] == ",":
                    day = int(user_date[d+1])
                    year = user_date[d+4:]
                    if month in Months and Numbers_days[day-1] !="":
                        to_print = year+"-"+dict1[month]+"-"+Numbers_days[day-1]
                        print(to_print)
                    else:
                        break

                elif user_date[d+3] == ",":
                    day = int(user_date[d+1:d+3])
                    year = user_date[d+6:]
                    if Numbers_days[day-1] !="":
                        to_print = year+"-"+dict1[month]+"-"+Numbers_days[day-1]
                        print(to_print)
                    else:
                        break

            elif user_date[d] == "/" and count2 == 0:   # 8/9/1999
                count2 += 1
                month = int(user_date[0:d])
                if user_date[d+2] == "/":
                    day = int(user_date[d+1])
                    year = user_date[d+3:]
                    to_print = year+"-"+Numbers_months[month-1]+"-"+Numbers_days[day-1]
                    print(to_print)


                elif user_date[d+3] == "/":
                    day = int(user_date[d+1:d+3])
                    year = user_date[d+4:]
                    to_print = year+"-"+Numbers_months[month-1]+"-"+Numbers_days[day-1]
                    print(to_print)


    except (NameError, IndexError,ValueError):
        print("REJECTED")
        break
    else:
        break
