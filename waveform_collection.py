from obspy.clients.fdsn import Client as FDSN_Client
from obspy.clients.earthworm import Client as EW_Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.geodetics import gps2dist_azimuth
from obspy import Stream, read, UTCDateTime
import numpy as np
from xarray import DataArray
import os
import glob
import json
import fnmatch
import warnings
from .grid import calculate_time_buffer


# Always show warnings
warnings.filterwarnings('always')

# Get location of folder containing this script
dirname = os.path.dirname(__file__)

# Load AVO infrasound station calibration values (units are Pa/ct)
with open(os.path.join(dirname, 'avo_json', 'avo_infra_calibs.json')) as f:
    AVO_INFRA_CALIBS = json.load(f)

# Load AVO infrasound station coordinates (elevation units are meters)
with open(os.path.join(dirname, 'avo_json', 'avo_infra_coords.json')) as f:
    AVO_INFRA_COORDS = json.load(f)

# Define IRIS and AVO clients (define WATC client within function)
iris_client = FDSN_Client('IRIS')
avo_client = EW_Client('pubavo1.wr.usgs.gov', port=16023)  # 16023 is long-term

# Channels to use in data requests - covering all the bases here!
CHANNELS = 'BDF,BDG,BDH,BDI,BDJ,BDK,HDF,DDF'

# Define some conversion factors
KM2M = 1000     # [m/km]
SEC2MIN = 1/60  # [min/s]
HR2SEC = 3600   # [s/hr]


def gather_waveforms(source, network, station, starttime, endtime,
                     time_buffer=0, remove_response=False,
                     return_failed_stations=False, watc_username=None,
                     watc_password=None):
    """
    Gather infrasound waveforms from IRIS or WATC FDSN, or AVO Winston, and
    output a Stream object with station/element coordinates attached.
    Optionally remove the sensitivity.

    NOTE:
        Usual RTM usage is to specify a starttime/endtime that brackets the
        estimated source origin time. Then time_buffer is used to download
        enough extra data to account for the time required for an infrasound
        signal to propagate to the farthest station. Because this buffer is so
        critical, this function issues a warning if it remains set to its
        default of 0 s.

    Args:
        source: Which source to gather waveforms from - options are:
                'IRIS' <-- IRIS FDSN
                'WATC' <-- WATC FDSN
                'AVO'  <-- AVO Winston
        network: SEED network code
        station: SEED station code
        starttime: Start time for data request (UTCDateTime)
        endtime: End time for data request (UTCDateTime)
        time_buffer: [s] Extra amount of data to download after endtime
                     (default: 0)
        remove_response: Toggle conversion to Pa via remove_sensitivity() if
                         available, else just do a simple scalar multiplication
                         (default: False)
        return_failed_stations: If True, returns a list of station codes that
                                were requested but not downloaded. This
                                disables the standard failed station warning
                                message (default: False)
        watc_username: Username for WATC FDSN server (default: None)
        watc_password: Password for WATC FDSN server (default: None)
    Returns:
        st_out: Stream containing gathered waveforms
        failed_stations: (Optional) List containing station codes that were
                         requested but not downloaded
    """

    print('--------------')
    print('GATHERING DATA')
    print('--------------')

    # Warn if buffer is set to 0 s
    if time_buffer == 0:
        warnings.warn('Time buffer is set to 0 seconds. Are you sure you\'ve '
                      'downloaded enough data for RTM?')

    # IRIS FDSN
    if source == 'IRIS':

        print('Reading data from IRIS FDSN...')
        try:
            st_out = iris_client.get_waveforms(network, station, '*', CHANNELS,
                                               starttime, endtime + time_buffer,
                                               attach_response=remove_response)
        except FDSNNoDataException:
            st_out = Stream()  # Just create an empty Stream object

    # WATC FDSN
    elif source == 'WATC':

        print('Connecting to WATC FDSN...')
        watc_client = FDSN_Client('http://10.30.5.10:8080',
                                  user=watc_username,
                                  password=watc_password)

        print('Successfully connected. Reading data from WATC FDSN...')
        try:
            st_out = watc_client.get_waveforms(network, station, '*', CHANNELS,
                                               starttime, endtime + time_buffer,
                                               attach_response=remove_response)
        except FDSNNoDataException:
            st_out = Stream()  # Just create an empty Stream object

    # AVO Winston
    elif source == 'AVO':

        print('Reading data from AVO Winston...')
        try:
            # Array case
            if station in ['ADKI', 'AKS', 'DLL', 'OKIF', 'SDPI']:

                # Select the correct channel
                if station in ['DLL', 'OKIF']:
                    channel = 'HDF'
                else:
                    channel = 'BDF'

                st_out = Stream()  # Make an empty Stream object to populate

                # Deal with funky channel naming convention for AKS (for all
                # other arrays, six numbered elements are assumed)
                if station == 'AKS':
                    for channel in ['BDF', 'BDG', 'BDH', 'BDI', 'BDJ', 'BDK']:
                        st_out += avo_client.get_waveforms(network, station,
                                                           '--', channel,
                                                           starttime,
                                                           endtime + time_buffer)
                else:
                    for location in ['01', '02', '03', '04', '05', '06']:
                        st_out += avo_client.get_waveforms(network, station,
                                                           location, channel,
                                                           starttime,
                                                           endtime + time_buffer)

            # Single station case
            else:
                st_out = avo_client.get_waveforms(network, station, '--',
                                                  'BDF', starttime,
                                                  endtime + time_buffer)

                # Special case for single stations with HDF channel
                if station in ['CLES1', 'CLES2', 'HAG', 'PVV', 'SSLN']:
                    st_out += avo_client.get_waveforms(network, station, '--',
                                                       'HDF', starttime,
                                                       endtime + time_buffer)

        # KeyError means that the station is not on AVO Winston for ANY time
        # period, OR that the user didn't format the request (e.g., station
        # string) appropriately
        except KeyError:
            st_out = Stream()  # Just create an empty Stream object

    else:
        raise ValueError('Unrecognized source. Valid options are \'IRIS\', '
                         '\'WATC\', or \'AVO\'.')

    st_out.sort()

    # Check that all requested stations are present in Stream
    requested_stations = station.split(',')
    downloaded_stations = [tr.stats.station for tr in st_out]
    failed_stations = []
    for sta in requested_stations:
        # The below check works with wildcards, but obviously cannot detect if
        # ALL stations corresponding to a given wildcard (e.g., O??K) were
        # downloaded. Thus, if careful station selection is desired, specify
        # each station explicitly and the below check will then be effective.
        if not fnmatch.filter(downloaded_stations, sta):
            if not return_failed_stations:
                # If we're not returning the failed stations, then show this
                # warning message to alert the user
                warnings.warn(f'Station {sta} not downloaded from {source} '
                              'server for this time period.')
            failed_stations.append(sta)

    # If the Stream is empty, then we can stop here
    if st_out.count() == 0:
        print('No data downloaded.')
        if return_failed_stations:
            return st_out, failed_stations
        else:
            return st_out

    # Otherwise, show what the Stream contains
    print(st_out.__str__(extended=True))  # This syntax prints the WHOLE Stream

    # Add zeros to ensure all Traces have same length
    st_out.trim(starttime, endtime + time_buffer, pad=True, fill_value=0)

    print('Assigning coordinates...')

    # Assign coordinates using IRIS FDSN regardless of data source
    try:
        inv = iris_client.get_stations(network=network, station=station,
                                       starttime=starttime,
                                       endtime=endtime + time_buffer,
                                       level='channel')
    except FDSNNoDataException:
        inv = []

    for tr in st_out:
        for nw in inv:
            for sta in nw:
                for cha in sta:
                    # Being very thorough to check if everything matches!
                    if (tr.stats.network == nw.code and
                            tr.stats.station == sta.code and
                            tr.stats.location == cha.location_code and
                            tr.stats.channel == cha.code):

                        tr.stats.longitude = cha.longitude
                        tr.stats.latitude = cha.latitude
                        tr.stats.elevation = cha.elevation

    # Check if any Trace did NOT get coordinates assigned, and try to use JSON
    # coordinates if available
    for tr in st_out:
        try:
            tr.stats.longitude, tr.stats.latitude, tr.stats.elevation
        except AttributeError:
            try:
                tr.stats.latitude, tr.stats.longitude,\
                    tr.stats.elevation = AVO_INFRA_COORDS[tr.stats.station]
                warnings.warn(f'Using coordinates from JSON file for {tr.id}.')
            except KeyError:
                print(f'No coordinates available for {tr.id}. Stopping.')
                raise

    # Remove sensitivity
    if remove_response:

        print('Removing sensitivity...')

        for tr in st_out:
            try:
                # Just removing sensitivity for now. remove_response() can lead
                # to errors. This should be sufficient for now. Plus some
                # IRIS-AVO responses are wonky.
                tr.remove_sensitivity()
            except ValueError:
                try:
                    calib = AVO_INFRA_CALIBS[tr.stats.station]
                    tr.data = tr.data * calib
                    tr.stats.processing.append('RTM: Data multiplied by '
                                               f'calibration value of {calib} '
                                               'Pa/ct')
                    warnings.warn('Using calibration value from JSON file for '
                                  f'{tr.id}.')
                except KeyError:
                    print('No calibration value available for {tr.id}. '
                          'Stopping.')
                    raise

    print('Done')

    # Return the Stream with coordinates attached (and responses removed if
    # specified)
    if return_failed_stations:
        return st_out, failed_stations
    else:
        return st_out


def gather_waveforms_bulk(lon_0, lat_0, max_radius, starttime, endtime,
                          time_buffer=0, remove_response=False,
                          watc_username=None, watc_password=None):
    """
    Bulk gather infrasound waveforms within a specified maximum radius of a
    specified location. Waveforms are gathered from IRIS (and optionally WATC)
    FDSN, and AVO Winston. Outputs a Stream object with station/element
    coordinates attached. Optionally removes the sensitivity. [Output Stream
    has the same properties as output Stream from gather_waveforms().]

    NOTE 1:
        WATC database will NOT be used for station search NOR data download
        unless BOTH watc_username and watc_password are set.

    NOTE 2:
        Usual RTM usage is to specify a starttime/endtime that brackets the
        estimated source origin time. Then time_buffer is used to download
        enough extra data to account for the time required for an infrasound
        signal to propagate to the farthest station. This function can
        automatically calculate an appropriate buffer time (it assumes that the
        station search center and source grid center are identical, which in
        practice should be the case since the grid center should be used as the
        station search center).

    Args:
        lon_0: [deg] Longitude of search center
        lat_0: [deg] Latitude of search center
        max_radius: [km] Maximum radius to search for stations within
        starttime: Start time for data request (UTCDateTime)
        endtime: End time for data request (UTCDateTime)
        time_buffer: Either a buffer time in s or an RTM grid (i.e., an
                     xarray.DataArray output from define_grid() for this
                     event). If a grid is specified, the buffer time in s is
                     automatically calculated based upon the grid params and
                     this function's station locations. This is the extra
                     amount of data to download after endtime, and is simply
                     passed on to the calls to gather_waveforms() (default: 0)
        remove_response: Toggle conversion to Pa via remove_sensitivity() if
                         available, else just do a simple scalar multiplication
                         (default: False)
        watc_username: Username for WATC FDSN server (default: None)
        watc_password: Password for WATC FDSN server (default: None)
    Returns:
        st_out: Stream containing bulk gathered waveforms
    """

    print('-------------------')
    print('BULK GATHERING DATA')
    print('-------------------')

    print('Creating station list...')

    # Grab IRIS inventory - not accounting for buffer here
    iris_inv = iris_client.get_stations(starttime=starttime, endtime=endtime,
                                        channel=CHANNELS, level='channel')

    inventories = [iris_inv]  # Add IRIS inventory to list

    # If the user supplied both a WATC password and WATC username, then search
    # through WATC database
    if watc_username and watc_password:

        print('Connecting to WATC FDSN...')
        watc_client = FDSN_Client('http://10.30.5.10:8080', user=watc_username,
                                  password=watc_password)
        print('Successfully connected.')

        # Grab WATC inventory - not accounting for buffer here
        watc_inv = watc_client.get_stations(starttime=starttime,
                                            endtime=endtime, channel=CHANNELS,
                                            level='channel')

        inventories.append(watc_inv)  # Add WATC inventory to list

    requested_station_list = []  # Initialize list of stations to request

    max_station_dist = 0  # [m] Keep track of the most distant station

    # Big loop through all channels in all inventories!
    for inv in inventories:
        for nw in inv:
            for stn in nw:
                for cha in stn:
                    dist, _, _ = gps2dist_azimuth(lat_0, lon_0, cha.latitude,
                                                  cha.longitude)  # [m]
                    if dist <= max_radius * KM2M:
                        requested_station_list.append(stn.code)
                        # Keep track of most distant station (within radius)
                        if dist > max_station_dist:
                            max_station_dist = dist

    # Loop through each entry in AVO infrasound station coordinates JSON file
    for sta, coord in AVO_INFRA_COORDS.items():
        dist, _, _ = gps2dist_azimuth(lat_0, lon_0, *coord[0:2])  # [m]
        if dist <= max_radius * KM2M:
            requested_station_list.append(sta)
            # Keep track of most distant station (within radius)
            if dist > max_station_dist:
                max_station_dist = dist

    if not requested_station_list:
        raise ValueError('Station list is empty. Expand the station search '
                         'and try again.')

    # Put into the correct format for ObsPy (e.g., 'HOM,O22K,DLL')
    requested_stations = ','.join(np.unique(requested_station_list))

    print('Done')

    # Check if time_buffer is an xarray.DataArray - if so, the user wants a
    # buffer time to be automatically calculated from this grid
    if type(time_buffer) == DataArray:
        time_buffer = calculate_time_buffer(grid=time_buffer,
                                            max_station_dist=max_station_dist)  # [s]

    if time_buffer != 0:
        print(f'Using time buffer of {time_buffer:.1f} s '
              f'(~{time_buffer * SEC2MIN:.0f} min)')

    print('Making calls to gather_waveforms()...')

    st_out = Stream()  # Initialize empty Stream to populate

    # Gather waveforms from IRIS
    iris_st, iris_failed = gather_waveforms(source='IRIS', network='*',
                                            station=requested_stations,
                                            starttime=starttime,
                                            endtime=endtime,
                                            time_buffer=time_buffer,
                                            remove_response=remove_response,
                                            return_failed_stations=True)
    st_out += iris_st

    # If IRIS couldn't grab all stations in requested station list, try WATC
    # (if the user set username and password)
    if iris_failed:

        if watc_username and watc_password:
            # Gather waveforms from WATC
            watc_st, watc_failed = gather_waveforms(source='WATC', network='*',
                                                    station=','.join(iris_failed),
                                                    starttime=starttime,
                                                    endtime=endtime,
                                                    time_buffer=time_buffer,
                                                    remove_response=remove_response,
                                                    return_failed_stations=True,
                                                    watc_username=watc_username,
                                                    watc_password=watc_password)
        else:
            # Return an empty Stream and same failed stations
            watc_st, watc_failed = Stream(), iris_failed

        st_out += watc_st

        # If WATC couldn't grab all stations missed by IRIS, try AVO
        if watc_failed:

            # Gather waveforms from AVO
            remaining_failed = []
            for sta in watc_failed:
                avo_st, avo_failed = gather_waveforms(source='AVO',
                                                      network='AV',
                                                      station=sta,
                                                      starttime=starttime,
                                                      endtime=endtime,
                                                      time_buffer=time_buffer,
                                                      remove_response=remove_response,
                                                      return_failed_stations=True)

                st_out += avo_st
                remaining_failed += avo_failed

            if remaining_failed:
                print('--------------')
                for sta in remaining_failed:
                    warnings.warn(f'Station {sta} found in radius search but '
                                  'no data found.')

    print('--------------')
    print('Finishing gathering waveforms from station list. Check warnings '
          'for any missed stations.')

    return st_out


def read_local(data_dir, coord_file, network, station, starttime, endtime):
    """
    Read in waveforms from "local" 1-hour, IRIS-compliant miniSEED files, and
    output a Stream object with station/element coordinates attached.

    NOTE 1:
        The expected naming convention for the miniSEED files is:
        <network>.<station>.<location>.<channel>.<year>.<julian_day>.<hour>

    NOTE 2:
        This function assumes that the response has been removed from the
        waveforms in the input miniSEED files. This is usually the case.

    Args:
        data_dir: Directory containing miniSEED files
        coord_file: JSON file containing coordinates for local stations (full
                    path required)
        network: SEED network code
        station: SEED station code
        starttime: Start time for data request (UTCDateTime)
        endtime: End time for data request (UTCDateTime)

    Returns:
        st_out: Stream containing gathered waveforms
    """

    print('-----------------------------')
    print('GATHERING LOCAL MINISEED DATA')
    print('-----------------------------')

    # Take (hour) floor of starttime
    starttime_hr = UTCDateTime(starttime.year, starttime.month, starttime.day,
                               starttime.hour)

    # Take (hour) floor of endtime - this ensures we check this miniSEED file
    endtime_hr = UTCDateTime(endtime.year, endtime.month, endtime.day,
                             endtime.hour)

    # Define filename template (don't check location or channel!)
    template = f'{network}.{station}.*.*.{{}}.{{}}.{{}}'

    # Initialize Stream object
    st_out = Stream()

    # Initialize the starting hour
    tmp_time = starttime_hr

    # Cycle forward in time, advancing hour by hour through miniSEED files
    while tmp_time <= endtime_hr:

        pattern = template.format(tmp_time.strftime('%Y'),
                                  tmp_time.strftime('%j'),
                                  tmp_time.strftime('%H'))

        files = glob.glob(os.path.join(data_dir, pattern))

        for file in files:
            st_out += read(file)

        tmp_time += HR2SEC  # Add an hour!

    st_out.merge()  # Merge traces with the same ID
    st_out.sort()

    # If the Stream is empty, then we can stop here
    if st_out.count() == 0:
        print('No data downloaded.')
        return st_out

    # Otherwise, show what the Stream contains
    print(st_out.__str__(extended=True))  # This syntax prints the WHOLE Stream

    # Add zeros to ensure all Traces have same length
    st_out.trim(starttime, endtime, pad=True, fill_value=0)

    print('Assigning coordinates...')

    # Assign coordinates by searching through user-supplied JSON file
    with open(coord_file) as f:
        local_coords = json.load(f)
    for tr in st_out:
        try:
            tr.stats.latitude, tr.stats.longitude,\
                tr.stats.elevation = local_coords[tr.stats.station]
        except KeyError:
            print(f'No coordinates available for {tr.id}. Stopping.')
            raise

    print('Done')

    # Return the Stream with coordinates attached
    return st_out
