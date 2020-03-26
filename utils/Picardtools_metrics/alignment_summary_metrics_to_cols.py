import sys

f = open(sys.argv[1], "r")
l = []

i = 0
for x in f:
    if x[0] != "#":
        x = x.strip()
        if len(x) != 0:
            l.append(x.split('\t'))
            i+=1

for t,v in zip(l[0], l[1]):
    print('{}\t{}'.format(t,v))
