ans = input("What is the Answer to the Great Question of Life, the Universe, and Everything? ")
if str(ans).strip() == "42" or (str(ans).lower().strip() == "forty-two") or (str(ans).lower().strip() == "forty two"):
    print("Yes")
else:
    print("No")
