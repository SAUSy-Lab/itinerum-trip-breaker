/*
create an index giving some estimate of how much this point looks like a wild spatial error
*/

﻿--ALTER TABLE it_coords ADD COLUMN error_index real;
WITH sub AS (
	SELECT 
		uid,
		user_id,
		ST_X(ST_Transform(geom,32617)) AS x,
		ST_Y(ST_Transform(geom,32617)) AS y,
		geom::geography AS geog,
		geom,
		row_number() OVER ( PARTITION BY user_id ORDER BY report_time) AS num
	FROM it_coords 
	WHERE 
		(geom::geography <-> ST_MakePoint(-79.3793,43.6567)::geography) / 1000 < 1000
		--AND user_id = '00e70d32-92e0-4f38-9065-ca8c7ed1664c'
), sub2 AS (
	SELECT 
		p2.uid,
		degrees(acos(round(
			((
				((p2.x-p1.x)^2 + (p2.y-p1.y)^2) + -- a
				((p2.x-p3.x)^2 + (p2.y-p3.y)^2) - -- b
				((p3.x-p1.x)^2 + (p3.y-p1.y)^2) -- c
			) /
			sqrt(
				4 * 
				((p2.x-p1.x)^2 + (p2.y-p1.y)^2) * -- a
				((p2.x-p3.x)^2 + (p2.y-p3.y)^2) -- b
			))::numeric
		,8))) AS a,
		least( (p1.geog<->p2.geog) , (p2.geog<->p3.geog)) AS d
	FROM sub AS p1 
	JOIN sub AS p2 ON p1.num +1 = p2.num AND p1.user_id = p2.user_id
	JOIN sub AS p3 ON p2.num +1 = p3.num AND p2.user_id = p3.user_id
	WHERE 
		NOT p1.geom~p2.geom AND NOT p2.geom~p3.geom
)
UPDATE it_coords SET 
	error_index = 1/((a+1)/d)
FROM sub2 
WHERE sub2.uid = it_coords.uid;
