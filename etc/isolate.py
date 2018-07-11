in_fd = "/inputs/SECcoordinates.csv"
users = ["73fcaa05-20af-434b-a8d2-7d079ba35a51", "1EDF2943-E2C1-4608-8CD2-FA53A54C340E"]
out = "/inputs/SECiso-coordinates.csv"
fd = open(in_fd, "r")
out_fd = open(out, "w")
line = fd.readline()
out_fd.write(line)
for l in fd:
  for user in users:
    if user in l:
      out_fd.write(l)
