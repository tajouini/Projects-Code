
def reward(word):
    if word.lower().strip()[0:5] == "hello":
        reward = "$0"
    elif word.lower().strip()[0] == "h":
        reward = "$20"
    else:
        reward = "$100"
    return reward

word = str(input("Greet me! "))
print(reward(word))
