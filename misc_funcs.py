#
# This file defines functions not associated with object classes
# 


def ts_str(date_time):
	"""
	Return a timestamp string from a timezone aware datetime object.
	e.g. '2017-09-08T16:54:16-04:00'
	"""
	return date_time.isoformat()


def read_headers(fname):
	""" (str) -> dict
	Return a dictionary mapping header names to column indices.

	Removes the need to hard coding column numbers when reading files.
	"""
	fd = open(fname)
	d = {}
	header = fd.readline()
	titles = header.split(",")
	for i in range(len(titles)):
		d[titles[i].strip()] = i
	fd.close()
	return d
