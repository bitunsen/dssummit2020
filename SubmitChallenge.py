import requests
url='https://lf8q0kx152.execute-api.us-east-2.amazonaws.com/default/computeFitnessScore'
x=requests.post(url,json={"qconfig":"4 0 7 5 2 6 1 3","userID":311234,"githubLink":"https://github.com/bitunsen/dssummit2020.git"})
print(x.text)
