od = open("isolate2.py", "w")
fd = open("isolate.py", "r")
lines = fd.readlines()
fd.close()
writen = []
for l in lines:
  if l not in writen:
    od.write(l)
    writen.append(l)
  else:
    od.write(l.rstrip() + "  # duplicate\n")
od.close() 
