def main():
    dollars = dollars_to_float(input("How much was the meal? "))
    percent = percent_to_float(input("What percentage would you like to tip? "))
    tip = dollars * percent
    print(f"Leave ${tip:.2f}")


def dollars_to_float(string):
    new_str = ""
    for i in string:
        if i == "$":
            new_str = new_str +""
        else:
            new_str = new_str +i
    return float(new_str)

def percent_to_float(string):  # format ##%
    arg = ""
    for i in string:
        if i == "%":
            arg = arg +""
        else:
            arg = arg +i
    return float(arg)/100


main()
