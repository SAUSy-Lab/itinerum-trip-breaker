in_fd = "../inputs/SECcoordinates.csv"
user = "73fcaa05-20af-434b-a8d2-7d079ba35a51"
out = "../inputs/SECiso-coordinates.csv"
fd = open(in_fd, "r")
out_fd = open(out, "w")
line = fd.readline()
out_fd.write(line)
for l in fd:
  if user in l:
    out_fd.write(l)
