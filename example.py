from obspy import UTCDateTime
from waveform_collection import gather_waveforms, gather_waveforms_bulk

# Time window of data to gather for first two examples
STARTTIME = UTCDateTime(2019, 9, 22, 6)
ENDTIME   = UTCDateTime(2019, 9, 22, 6, 10)
#%% Example 1 - Gather all infrasound records from Dillingham infrasound array

st = gather_waveforms(source='IRIS', network='AV', station='DLL',
                      location='*', channel='*', starttime=STARTTIME,
                      endtime=ENDTIME)

#%% Example 2 - Gather all BHN channel seismic records from within a 200 km radius of Iliamna volcano, Alaska

st_bulk = gather_waveforms_bulk(lon_0=-153.0918, lat_0=60.0319, max_radius=200,
                                starttime=STARTTIME, endtime=ENDTIME,
                                channel='BHN', parallel=False, cores=7)
# gather_waveforms_bulk can also be run in parallel, faster for ~>1 day of data

#%% Example 3 - Gather 10 days of 3-component seismic data from SSLS in parallel using 6 cores.
STARTTIME = UTCDateTime(2023, 9, 1)
ENDTIME   = UTCDateTime(2023, 9, 11)  # 10 days

# The parallel feature is best used to improve the speed when gathering especially long periods of data (days to months).
st_long = gather_waveforms(source='IRIS', network='AV', station='SSLS',
                      location='*', channel='BH*', starttime=STARTTIME,
                      endtime=ENDTIME, parallel=True, cores=6)