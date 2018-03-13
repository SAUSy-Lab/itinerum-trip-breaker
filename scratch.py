def project(longitude,latitude,projection_string='epsg:3347'):
	"""Project lat-lon values. Default of 3347 is StatsCan Lambert.
		Units in meters."""
	from pyproj import Proj, transform
	inProj = Proj( init = 'epsg:4326' )
	outProj = Proj( init = projection_string )
	x,y = transform( inProj, outProj, longitude, latitude )
	return x,y

def 
