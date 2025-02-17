#In a file called coke.py,
# implement a program that prompts the user to insert a coin,
# one at a time, each time informing the user of the amount due.
# Once the user has inputted at least 50 cents, output how many cents in change the user is owed.
# Assume that the user will only input integers, and ignore any integer that isn’t an accepted denomination.




sum = 0

while True:
    coin = input("insert a coin one at a time ")
    if int(coin) in [25, 10, 5]:
        sum = sum + int(coin)
        if sum >= 50:
            print(f"Change Owed: {sum-50}")
            break
        else:
            print(f"Amount Due: {50-sum}")
    else:
            sum = sum
            print(f"Amount Due: {50-sum}")




