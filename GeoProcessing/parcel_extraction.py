#!/usr/bin/env python
# PURPOSE       : To create urban parcel from road networks
# Ref: Liu, Xingjian, and Ying Long. “Automated Identification and Characterization of Parcels with OpenStreetMap and Points of Interest.” 
#           Environment and Planning B-Planning & Design, vol. 43, no. 2, 2016, pp. 341–360.


import grass.script as gs
import re


# Display the mapset content
vectors = gs.read_command("g.list", type='vect')
print vectors

roads = raw_input("Please enter the name of the road layer : ")
river = raw_input("Please enter the name of the river layer : ")
boundary = raw_input("Please enter the name of the urban boundary layer : ") 
roadspace_name = raw_input("Please enter the name of the output road space layer : ") 
parcel_name = raw_input("Please enter the name of the output urban parcel layer : ") 


# Set the buffer zone width for different type of roads and rivers
roads_width = {"motorway": 15, "trunk": 15, "primary": 15, "secondary": 12, "tertiary": 10,
            "unclassified": 4, "residential": 6, "living_street": 6, "pedestrian": 6,
            "service": 6, "track": 4, "track_grade1": 4, "track_grade2": 4, "track_grade3": 4, "track_grade4": 4, "track_grade5": 4,
            "bridleway": 4, "cycleway": 4, "footway": 4, "path": 2, "steps": 2, "unknown": 2}

river_width = {"river": 12, "stream": 8, "canal": 10, "drain": 6}

area_min = 300.0
ff1_min = 0.10
ff2_min = 0.20
area_1 = 25000
area_2 = 2500

# --1. clean the topology of road and river networks
#####################################################################################################
# m.snaplp only support DBF as database driver 
# Before running this command, use: "db.connect -p" to check the database connection
# If not use DBF, use command db.connect, db.copy, v.db.connect to modify the connection
#####################################################################################################
gs.run_command("db.connect", driver="dbf", database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/")

gs.run_command("v.clean", flags="c", input=roads, output="temp_roads_clean", type="line", tool="snap,rmdangle,prune,rmline,rmsa",thresh="5, 250",overwrite=True)
gs.run_command("v.clean", flags="c", input=river, output="temp_river_clean", type="line", tool="snap,rmdangle,prune,rmline,rmsa",thresh="5, 250",overwrite=True)

'''
gs.run_command("db.copy", from_driver="sqlite", from_database="$GISDBASE/$LOCATION_NAME/$MAPSET/sqlite/sqlite.db", from_table="temp_roads_clean", to_driver="dbf", to_database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/", to_table="temp_roads_clean")
gs.run_command("db.copy", from_driver="sqlite", from_database="$GISDBASE/$LOCATION_NAME/$MAPSET/sqlite/sqlite.db", from_table="temp_river_clean", to_driver="dbf", to_database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/", to_table="temp_river_clean")
gs.run_command("db.copy", from_driver="sqlite", from_database="$GISDBASE/$LOCATION_NAME/$MAPSET/sqlite/sqlite.db", from_table=boundary, to_driver="dbf", to_database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/", to_table=boundary)
gs.run_command("v.db.connect", map="temp_roads_clean", table="temp_roads_clean", driver=dbf, database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/", overwrite=True)
gs.run_command("v.db.connect", map="temp_river_clean", table="temp_river_clean", driver=dbf, database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/", overwrite=True)
gs.run_command("v.db.connect", map=boundary, table="temp_river_clean", driver=dbf, database="$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/", overwrite=True)
'''

# ###ref to Geo: m.snaplp
gs.run_command("m.snaplp_tmp", input="temp_roads_clean", output="temp_roads_clean_2", polygon=boundary, snap=2.0, overwrite=True)
gs.run_command("m.snaplp_tmp", input="temp_river_clean", output="temp_river_clean_2", polygon=boundary, snap=2.0, overwrite=True)

# --2. create the buffer zone for each road
num_road_type = {}
for road_type in roads_width.keys():
    num_road_type[road_type]=len(re.findall(r"\S+", gs.read_command("v.db.select", map="temp_roads_clean", columns="fclass", where="fclass = '%s'"%road_type, flags="c")))

num_river_type = {}
for river_type in river_width.keys():
    num_river_type[river_type]=len(re.findall(r"\S+", gs.read_command("v.db.select", map="temp_river_clean", columns="fclass", where="fclass = '%s'"%river_type, flags="c")))

print "----------Create the buffer zone for each road and river----------"
bf_road_list = []
for road_type, width in roads_width.items():
    if num_road_type[road_type] > 0:
        gs.run_command("v.buffer", input="temp_roads_clean", layer=1, where="fclass = '%s'"%road_type, output="bf_"+road_type, type="line", distance=width)
        bf_road_list.append("bf_"+road_type)

bf_river_list = []
for river_type, width in river_width.items():
    if num_river_type[river_type] > 0:
        gs.run_command("v.buffer", input="temp_river_clean", layer=1, where="fclass = '%s'"%river_type, output="bf_"+river_type, type="line", distance=width)
        bf_river_list.append("bf_"+river_type)

# --3. create the road and river buffer zone by overlay operation i.e., union
l_ro = len(bf_road_list)
gs.run_command("v.overlay", ainput=bf_road_list[0], binput=bf_road_list[1], output="temp_road"+str(0), operator="or", overwrite=True)
for i in range(1, l_ro-1):
    gs.run_command("v.overlay", ainput="temp_road"+str(i-1), binput=bf_road_list[i+1], output="temp_road"+str(i), operator="or", overwrite=True)

l_ri = len(bf_river_list)
gs.run_command("v.overlay", ainput=bf_river_list[0], binput=bf_river_list[1], output="temp_river"+str(0), operator="or", overwrite=True)
for i in range(1, l_ri-1):
    gs.run_command("v.overlay", ainput="temp_river"+str(i-1), binput=bf_river_list[i+1], output="temp_river"+str(i), operator="or", overwrite=True)

# combine road buffer zone with river buffer zone 
gs.run_command("v.overlay", ainput="temp_road"+str(l_ro-2), binput="temp_river"+str(l_ri-2), output=roadspace_name, operator="or", overwrite=True)

# --4. create the parcel by reversing the road zone 
gs.run_command("v.overlay", ainput=boundary, binput=roadspace_name, output=parcel_name, operator="not", overwrite=True)
# clean the temporary layer
gs.run_command("g.remove", flags="f", type="vector", pattern="bf*")
gs.run_command("g.remove", flags="f", type="vector", pattern="temp*")

# --5. clean the topology of urban parcel
# ------------adapted from GeoPUMMA:p.A1.clean_topology.py
print "----------Clean the topology of urban parcel----------"
gs.run_command("v.category", input=parcel_name, output='clean2', type='boundary', option='add', overwrite=True)
gs.run_command("v.extract", flags='t', input='clean2', output='clean3', type='boundary', overwrite=True)
gs.run_command("v.type", input='clean3', output='clean4', from_type='boundary', to_type='line', overwrite=True)

gs.run_command("m.snaplp", input="clean4", output="clean5", polygon=boundary, snap=2, overwrite=True)

gs.run_command("v.clean", flags="c", input='clean4', output='clean5', type='line', tool="snap,rmdangle,prune,rmline,rmsa", thresh="3, 150", overwrite=True)
gs.run_command("v.type", input='clean5', output='clean6', from_type='line', to_type='boundary', overwrite=True)
gs.run_command("v.centroids", input='clean6', output='clean7', overwrite=True)
gs.run_command("v.category", input='clean7', output='clean8', type='boundary', option='del', overwrite=True)
gs.run_command("v.clean", input='clean8', output='clean9', tool='rmarea', thresh='0', overwrite=True)
gs.run_command("v.category", input='clean9', output='clean10', option='del', type='boundary', overwrite=True)
gs.run_command("v.build.polylines", input='clean10', output='clean11', cats='multi', overwrite=True)
    
gs.run_command("v.db.addtable", map='clean11', col='b_cat INT', layer='1', overwrite=True)
gs.run_command("v.distance", from_='clean11', from_type='centroid', from_layer='1', to=parcel_name, upload='cat',column='b_cat', overwrite=True)
    
gs.run_command('v.db.addcolumn',map='clean11',columns='c_cat INT')
gs.run_command('v.db.update', map='clean11',column='c_cat',value='b_cat')
gs.run_command("v.reclass", input='clean11', output='clean12', column='c_cat', overwrite=True)
gs.run_command("db.copy", from_table=parcel_name, to_table='clean12', overwrite=True)
gs.run_command("v.db.connect", map='clean12', table='clean12', layer='1', overwrite=True)
    
#Extracting only features with category > 0
gs.run_command("v.extract",input='clean12',output=parcel_name+"_temp", where="cat>0", overwrite=True)
    
gs.run_command("g.remove", flags="f", type="vector", pattern="clean*")

gs.run_command("v.db.droptable", flags="f", map=parcel_name+"_temp")
gs.run_command("v.db.addtable", map=parcel_name+"_temp")

gs.run_command("v.db.addcolumn", map=parcel_name+"_temp", col="cat int,area double precision,perimeter double precision,ff double precision")

gs.run_command("v.to.db", map=parcel_name+"_temp", option="cat", col="cat")
gs.run_command("v.to.db", map=parcel_name+"_temp", option="area", col="area")
gs.run_command("v.to.db", map=parcel_name+"_temp", option="perimeter", col="perimeter")
gs.run_command("v.db.update", map=parcel_name+"_temp", col="ff", query_column="16*area/(perimeter*perimeter)")

gs.run_command("v.extract", flags="d", input=parcel_name+"_temp", output=parcel_name+"_cleaned", type="area", where="(area>=%s) AND ((ff>=%s) OR (area>=%s)) AND ((ff>=%s) OR (area>=%s))"%(area_min, ff1_min, area_1, ff2_min, area_2), overwrite=True)

gs.run_command("v.overlay", ainput=boundary, binput=parcel_name+"_cleaned", output=roadspace_name, operator="not", overwrite=True)

gs.run_command("g.remove", flags="f", type="vector", name=parcel_name+"_temp")
