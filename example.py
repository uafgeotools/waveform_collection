from obspy import UTCDateTime
from waveform_collection import gather_waveforms, gather_waveforms_bulk

# Time window of data to gather for all examples
STARTTIME = UTCDateTime(2019, 9, 22, 6)
ENDTIME   = UTCDateTime(2019, 9, 22, 6, 10)

#%% Example 1 - Gather all infrasound records from Dillingham infrasound array

st = gather_waveforms(source='IRIS', network='AV', station='DLL',
                      location='*', channel='*', starttime=STARTTIME,
                      endtime=ENDTIME)

#%% Example 2 - Gather all BHN channel seismic records from within a 200 km radius of Iliamna volcano, Alaska

st_bulk = gather_waveforms_bulk(lon_0=-153.0918, lat_0=60.0319, max_radius=200,
                                starttime=STARTTIME, endtime=ENDTIME,
                                channel='BHN')
