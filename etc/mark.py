# Mark any duplicate lines in a plaintext file

od = open("times2.txt", "w")  # Copy of input file with duplicate lines marked
fd = open("times.txt", "r")  # Input file
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
