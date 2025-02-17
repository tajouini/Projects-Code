import random

level = -1

while level<=0:
    try:
        level = int(input("Level: "))
    except ValueError:
        pass

num = random.randint(1,level)

var = True
while var:
    try:
        guess = int(input("Guess: "))
    except ValueError:
        pass
    else:
        if guess<num:
            print("too small!")
        elif guess>num:
            print("too large!")
        else:
            print("Just right!")
            var = False



