from obspy import Stream, read, UTCDateTime
import glob
import os
from .common import load_json_file


# Define some conversion factors
HR2SEC = 3600   # [s/hr]


def read_local(data_dir, coord_file, network, station, location, channel,
               starttime, endtime, merge=True):
    """
    Read in waveforms from "local" 1-hour, IRIS-compliant miniSEED files, and
    output a Stream object with station/element coordinates attached.

    NOTE 1:
        The expected naming convention for the miniSEED files is:
        <network>.<station>.<location>.<channel>.<year>.<julian_day>.<hour>

    NOTE 2:
        This function assumes that the response has been removed from the
        waveforms in the input miniSEED files.

    Args:
        data_dir: Directory containing miniSEED files
        coord_file: JSON file containing coordinates for local stations (full
                    path required)
        network: SEED network code [wildcards (*, ?) accepted]
        station: SEED station code [wildcards (*, ?) accepted]
        location: SEED location code [wildcards (*, ?) accepted]
        channel: SEED channel code [wildcards (*, ?) accepted]
        starttime: Start time for data request (UTCDateTime)
        endtime: End time for data request (UTCDateTime)
        merge: Toggle merging of Traces with identical IDs (default: True)

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

    # Define filename template
    template = f'{network}.{station}.{location}.{channel}.{{}}.{{}}.{{}}'

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

    if merge:
        st_out.merge()  # Merge Traces with the same ID
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
    local_coords = load_json_file(coord_file)
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
