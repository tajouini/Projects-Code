# prompts the user for a level,
#. If the user does not input 1, 2, or 3, the program should prompt again.
import random


def main():
    level = get_level()
    score = str(generate_integer(level))
    print("Score: "+score)

def get_level():
    while True:
        try:
            level = int(input("Level: "))
            if level not in [1,2,3]:
                raise ValueError
        except ValueError:
            pass
        else:
            break
    return level


def generate_integer(level):
    score = 0
    if level ==1:
        start = 0
        end = 10
    elif level ==2:
        start = 10
        end = 99
    else:
        start = 100
        end = 999

    for i in range(10):
        X = random.randint(start, end)
        Y = random.randint(start,end)
        bool = True
        trials = 1
        while bool and trials<=3:
            sum = int(input(str(X)+  " + " +str(Y)+" = "))
            if sum != X+Y:
                print("EEE")
                trials +=1
            else:
                bool = False
        if trials<=3:
            score +=1
        else:
            print(str(X) + " + "+ str(Y)+" = "+str(X+Y))
    return score

if __name__ == "__main__":
    main()
