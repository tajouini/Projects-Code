import requests
import sys
import json
while True:
    try:
        n = float(sys.argv[1])
    except ValueError:
        sys.exit("Command-line argument is not a number")
    except IndexError:
        sys.exit("Missing command-line argument")
    else:
        try:
            response = requests.get("https://api.coindesk.com/v1/bpi/currentprice.json" + str(n))
            o = response.json()
            print(o)
          #  for result in o["bpi"]:
          #      print(result["USD"])
        except requests.RequestException:
            break
        else:
            print(f"${amount:,.4f}")

