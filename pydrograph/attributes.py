#  site_no         -- Site identification number
#  station_nm      -- Site name
#  site_tp_cd      -- Site type
#  lat_va          -- DMS latitude
#  long_va         -- DMS longitude
#  dec_lat_va      -- Decimal latitude
#  dec_long_va     -- Decimal longitude
#  coord_meth_cd   -- Latitude-longitude method
#  coord_acy_cd    -- Latitude-longitude accuracy
#  coord_datum_cd  -- Latitude-longitude datum
#  dec_coord_datum_cd -- Decimal Latitude-longitude datum
#  district_cd     -- District code
#  state_cd        -- State code
#  county_cd       -- County code
#  country_cd      -- Country code
#  land_net_ds     -- Land net location description
#  map_nm          -- Name of location map
#  map_scale_fc    -- Scale of location map
#  alt_va          -- Altitude of Gage/land surface
#  alt_meth_cd     -- Method altitude determined
#  alt_acy_va      -- Altitude accuracy
#  alt_datum_cd    -- Altitude datum
#  huc_cd          -- Hydrologic unit code
#  basin_cd        -- Drainage basin code
#  topo_cd         -- Topographic setting code
#  data_types_cd   -- Flags for the type of data collected
#  instruments_cd  -- Flags for instruments at site
#  construction_dt -- Date of first construction
#  inventory_dt    -- Date site established or inventoried
#  drain_area_va   -- Drainage area
#  contrib_drain_area_va -- Contributing drainage area
#  tz_cd           -- Mean Greenwich time offset
#  local_time_fg   -- Local standard time flag
#  reliability_cd  -- Data reliability code
#  gw_file_cd      -- Data-other GW files
#  nat_aqfr_cd     -- National aquifer code
#  aqfr_cd         -- Local aquifer code
#  aqfr_type_cd    -- Local aquifer type code
#  well_depth_va   -- Well depth
#  hole_depth_va   -- Hole depth
#  depth_src_cd    -- Source of depth data
#  project_no      -- Project number
#  rt_bol          -- Real-time data flag
#  peak_begin_date -- Peak-streamflow data begin date
#  peak_end_date   -- Peak-streamflow data end date
#  peak_count_nu   -- Peak-streamflow data count
#  qw_begin_date   -- Water-quality data begin date
#  qw_end_date     -- Water-quality data end date
#  qw_count_nu     -- Water-quality data count
#  gw_begin_date   -- Field water-level measurements begin date
#  gw_end_date     -- Field water-level measurements end date
#  gw_count_nu     -- Field water-level measurements count
#  sv_begin_date   -- Site-visit data begin date
#  sv_end_date     -- Site-visit data end date
#  sv_count_nu     -- Site-visit data count

streamflow_attributes = [ \
'site_no',         
'station_nm',      
'site_tp_cd',             
'dec_lat_va',      
'dec_long_va',     
'coord_meth_cd',   
'coord_acy_cd',    
'coord_datum_cd',  
'dec_coord_datum_cd',
'district_cd',     
'state_cd',        
'county_cd',       
'country_cd',      
'land_net_ds',     
'map_nm',          
'map_scale_fc',    
'alt_va',          
'alt_meth_cd',     
'alt_acy_va',      
'alt_datum_cd',    
'huc_cd',          
'basin_cd',        
'topo_cd',         
'inventory_dt',    
'drain_area_va',   
'contrib_drain_area_va',
'tz_cd',           
'local_time_fg',   
'reliability_cd',    
'project_no',     
'rt_bol',         
'peak_begin_date', 
'peak_end_date',   
'peak_count_nu',   
'qw_begin_date',   
'qw_end_date',     
'qw_count_nu',        
'sv_begin_date',   
'sv_end_date',     
'sv_count_nu']

iv_attributes = [ \
'agency_cd',    
'site_no',  
'station_nm',   
'site_tp_cd',   
'dec_lat_va',   
'dec_long_va',  
'coord_acy_cd', 
'dec_coord_datum_cd',   
'alt_va',   
'alt_acy_va',   
'alt_datum_cd', 
'huc_cd',   
'data_type_cd', 
'parm_cd',  
'stat_cd',   
'ts_id',    
'loc_web_ds',   
'medium_grp_cd',    
'parm_grp_cd',  
'srs_id',   
'access_cd',    
'begin_date',   
'end_date', 
'count_nu']

gw_attributes = [
'site_no',
'station_nm',
'site_tp_cd',
'dec_lat_va',
'dec_long_va',
'coord_meth_cd',
'coord_acy_cd',
'coord_datum_cd',
'dec_coord_datum_cd',
'district_cd',
'state_cd',
'county_cd',
'country_cd',
'land_net_ds',
'map_nm',
'map_scale_fc',
'alt_va',
'alt_meth_cd',
'alt_acy_va',
'alt_datum_cd',
'huc_cd',
'basin_cd',
'topo_cd',
'data_types_cd',
'instruments_cd',
'construction_dt',
'inventory_dt',
'tz_cd',
'local_time_fg',
'reliability_cd',
'gw_file_cd',
'nat_aqfr_cd',
'aqfr_cd',
'aqfr_type_cd',
'well_depth_va',
'hole_depth_va',
'depth_src_cd',
'project_no',
'rt_bol',
'peak_begin_date',
'peak_end_date',
'peak_count_nu',
'qw_begin_date',
'qw_end_date',
'qw_count_nu',
'gw_begin_date',
'gw_end_date',
'gw_count_nu',
'sv_begin_date',
'sv_end_date',
'sv_count_nu'
]
