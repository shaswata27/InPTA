#!/usr/bin/env python3

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.dates as mdates
from matplotlib.ticker import tor
from astropy.time import Time
from astropy.coordinates import get_sun, SkyCoord, EarthLocation, AltAz, Angle
from datetime import datetime, timedelta
import pandas as pd


def create_or_clear_directory(directory_path):
    if os.path.exists(directory_path):
        # If the directory exists, remove all its contents
        for filename in os.listdir(directory_path):
            file_path = os.path.join(directory_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                    print(f"Deleted directory: {file_path}")  # Debugging statement
            except Exception as e:
                print(f'Failed to delete {file_path}. Reason: {e}')
    else:
        # If the directory does not exist, create it
        os.makedirs(directory_path)
        
    #freshly opening the summary file
    summaryfile = f'{directory_path}/summary.txt'
    with open(summaryfile, 'w') as file:
        file.write(f"# This file contains the Summary of separtion angle between the sun and the sources during the observation time: \n") 
        file.write("----------------------------------------------------------------------------------------- \n \n")
   
    return summaryfile
               
def convert_ist_to_utc(ist_time):
    """Convert IST to UTC by subtracting 5 hours 30 minutes."""
    return (datetime.strptime(ist_time, '%Y-%m-%d %H:%M:%S') - timedelta(hours=5, minutes=30)).strftime('%Y-%m-%d %H:%M:%S')

def generate_time_range(start_time, end_time, interval_minutes):
    """Generate a list of times from start to end at a given interval."""
    start = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
    end = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    print(f"start = {start}")
    print(f"end = {end}")
    times = [start + timedelta(minutes=i) for i in range(0, int((end - start).total_seconds() / 60) + 1, interval_minutes)]
    print(f"times = {times}")
    return [t.strftime('%Y-%m-%d %H:%M:%S') for t in times]

def sepang_calc(time, psr_coord):
 
    suncoord = get_sun(time)
    print(suncoord)
    sep = suncoord.separation(psr_coord)
    print(sep)
    print(f'sepdeg = {sep.deg}')
    print("\n")

    return sep.deg
    
def get_positions(times, gmrt_location, RA, DEC):
    """Get the positions of the Sun and Pulsar for each timestamp."""
    sun_positions = []
    pulsar_positions = []
    sep_ang = [] #new

    for time in times:
        formatted_time = Time(time, format='iso', scale='utc')
        print(f"Formatted time = {formatted_time}")
        # Get Sun's position
        suncoord = get_sun(formatted_time)
        sun_altaz = suncoord.transform_to(AltAz(obstime=formatted_time, location=gmrt_location))
        
        # Get Pulsar's position
        pulsar_coord = SkyCoord(RA, DEC, frame='icrs')
        pulsar_altaz = pulsar_coord.transform_to(AltAz(obstime=formatted_time, location=gmrt_location))
        
        # Store azimuth and altitude (to plot on sky map)
        sun_positions.append((sun_altaz.az.deg, sun_altaz.alt.deg))
        pulsar_positions.append((pulsar_altaz.az.deg, pulsar_altaz.alt.deg))
        
        #sep_ang.append(sepang_calc(time, RA, DEC))) #new
        sep_ang.append(sepang_calc(formatted_time, pulsar_coord))

    return np.array(sun_positions), np.array(pulsar_positions), np.array(sep_ang) #new
    #return np.array(sun_positions), np.array(pulsar_positions)
     
# Function to plot the separation angle as a function of time
def plot_separation_angle(times, output_folder, separation_angles, targetname, threshold):
    
    # Convert times to datetime objects if they are not already in datetime format
    times = pd.to_datetime(times) #this time is in utc
    
    # Add 5 hours and 30 minutes to convert to IST
    times_ist = times + timedelta(hours=5, minutes=30)

    
    #print(f"times = {times}")
    #print(f"times[0] = {times[0]}")
    
    #print(f"times_ist = {times_ist}")
    #print(f"times_ist[0] = {times_ist[0]}")
    #sys.exit(1)
    
    #print(f"sepangles = {separation_angles}")
    #print(f"sepangles[0] = {separation_angles[0]}")
    #sys.exit(1)
    
    overlapflag = 0
    colorcode = None
    for sep in separation_angles: #looking for if at any point of time the separtion angle is smaller than the threshold
        if sep < threshold:
            overlapflag = 1
            break
    
    if overlapflag == 1:
        colorcode = 'red'
    else:
        colorcode = 'blue'
        
    fig, ax = plt.subplots(figsize=(8, 6))
    #ax.plot(times, separation_angles, label='Separation Angle', color=colorcode)
    ax.plot(times_ist, separation_angles, label='Separation Angle', color=colorcode)
    #print(f"times = {times}")
    #print(f"length of times = {len(times)}")
    #ax.set(xlabel='Time (UTC)', ylabel='Separation Angle (degrees)', title=f'Sun-Pulsar Separation Angle Timeseries - {targetname}')
    ax.set(xlabel='Time (IST)', ylabel='Separation Angle (degrees)', title=f'Sun-Pulsar Separation Angle Timeseries - {targetname}')
    
    # Set x-axis major locator to hourly intervals and formatter
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))  # Ticks every 1 hour
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H-%M'))  # Format as Hour-Minute
    
    # Use MaxNLocator to limit the number of ticks
    ax.xaxis.set_major_locator(MaxNLocator(nbins=9))  # Limit the number of ticks to a reasonable value
    
    
    #plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    
    figname = f"{output_folder}/{targetname}_separation_angle_vs_time.pdf"
    #Save to file (pdf or jpeg)
    plt.savefig(figname, format='pdf')  # You can change the format here
    
    #plt.close(fig)

def endtimecalc(startime, obs_time):
    # Convert the string to a datetime object
    time_obj = datetime.strptime(startime, "%Y-%m-%d %H:%M:%S")

    # Add obs_time minutes
    new_time_obj = time_obj + timedelta(hours=obs_time)

    # Convert back to string in the same format
    new_time_str = new_time_obj.strftime("%Y-%m-%d %H:%M:%S")
    
    return new_time_str

def parse_angle(strngval):
    # Parse the angle using astropy
    try:
        dec_angle = Angle(strngval)
        print(dec_angle)
    except ValueError as e:
        print(f"Error parsing angle: {e}")
        sys.exit(1)
    return dec_angle
    

def observatory_coord(obsrv_coord_file, obsname):
    with open(obsrv_coord_file, 'r') as file:
        for line in file:
            print(line)
            row = re.split(r'\s+', line.strip())  # Split by any amount of whitespace
            print(row)
            print(len(row))
            if len(row) != 3:
                print("Something is wrong in the psrcoord file \n")
                sys.exit(1)
            print(f"0: {row[0]},  1: {row[1]},  2: {row[2]}")
            if isinstance(row[0], str):
                print(f"{row[0]} is a string")
            if isinstance(row[1], str):
                print(f"{row[1]} is a string")
            if isinstance(row[2], str):
                print(f"{row[2]} is a string")
            
            if row[0] == obsname:
                lat = row[1]
                long = row[2]
                print(f"The observatory coordinate found")
                return lat, long
    print(f"The observatory was not elisted in the observatory coord list. Please Check the file!")
    sys.exit(1) #code execution stopped
    return

def is_leap_year(year):
    """Check if a given year is a leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

def validate_date(dd, mm, yyyy):
    """Validate the date considering month ranges and leap years."""
    if mm < 1 or mm > 12:
        raise ValueError("Month must be between 1 and 12.")
    
    if dd < 1:
        raise ValueError("Day must be at least 1.")
    
    days_in_month = [31, 29 if is_leap_year(yyyy) else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if dd > days_in_month[mm - 1]:
        raise ValueError(f"Day exceeds the number of days in month {mm}.")
    
    return True

def prompt_for_date():
    """Prompt the user for the date in DD-MM-YYYY format."""
    while True:
        date_str = input("Enter the date (DD-MM-YYYY): ")
        try:
            dd, mm, yyyy = map(int, date_str.split('-'))
            validate_date(dd, mm, yyyy)
            return f"{yyyy:04d}-{mm:02d}-{dd:02d}"  # Convert to YYYY-MM-DD format
        except ValueError as e:
            print(f"Invalid date: {e}. Please try again.")

def prompt_for_time():
    """Prompt the user for the time in HH:MM:SS format."""
    while True:
        time_str = input("Enter the time (HH:MM:SS): ")
        try:
            datetime.strptime(time_str, '%H:%M:%S')  # Validate time format
            return time_str
        except ValueError:
            print("Invalid time format. Please use HH:MM:SS.")

def prompt_for_duration():
    """Prompt the user for the duration in minutes."""
    while True:
        duration_str = input("Enter the duration (in hours): ")
        try:
            duration = int(duration_str)
            if duration < 0 or duration > 24:
                raise ValueError("Duration in hours cannot be negative or more than 24.")
            return duration
        except ValueError as e:
            print(f"Invalid duration: {e}. Please enter a valid number.")

def prompt_for_sourcelist():
    """Prompt the user for the pulsar name."""
    while True:
        sourcelistfile = input("Enter the source list with full path: ")
        #chekcing whether the above file exist    
        if os.path.exists(sourcelistfile) and sourcelistfile.endswith('.txt'):
            print("source list file in the path exists.")
            return sourcelistfile
        else:
            print("source list file in the path does NOT exists.")
            sys.exit(1)


def prompt_for_sepangthreshold():
    """Prompt the user for the separation angle threshold like 9 degrees below which the solar effects can enter the side beams."""
    while True:
        SAthreshold_str = input("Enter the separation angle threshold (in degrees): ")
        try:
            SAthreshold = int(SAthreshold_str)
            if SAthreshold < 0 or SAthreshold > 180:
                raise ValueError("SAthreshold cannot be negative or more than 180 degrees.")
            return SAthreshold
        except ValueError as e:
            print(f"Invalid threshold: {e}. Please enter a valid threshold (degrees).")
        
        
def main(obsrv_coord_file, outputfolder, summary, src_list_file, start_time_ist, obs_time, threshold, obsname):
    
    end_time_ist = endtimecalc(start_time_ist, obs_time)
    # Convert start and end times from IST to UTC
    start_time_utc = convert_ist_to_utc(start_time_ist)
    end_time_utc = convert_ist_to_utc(end_time_ist)
    
    # Generate times for the given interval (every 10 minutes)
    interval_minutes = 10
    times = generate_time_range(start_time_utc, end_time_utc, interval_minutes)
    print(f"first times generation in main: {times}")
    #print(f"# of times = {len(times)}")
    #sys.exit(1)
    latitude, longitude = observatory_coord(obsrv_coord_file, obsname)   
    # Set GMRT location
    gmrt_location = EarthLocation(lat=latitude, lon=longitude)
    
    
    
    with open(src_list_file, 'r') as file:
        for line in file:
            print(line)
            row = re.split(r'\s+', line.strip())  # Split by any amount of whitespace
            print(row)
            print(len(row))
            if len(row) < 3:
                print("Something is wrong in the Srclist file \n")
                sys.exit(1)
            print(f"0: {row[0]},  1: {row[1]},  2: {row[2]}")
            if isinstance(row[0], str):
                print(f"{row[0]} is a string")
            if isinstance(row[1], str):
                print(f"{row[1]} is a string")
            if isinstance(row[2], str):
                print(f"{row[2]} is a string")
            
            target_name = row[0]
            RA = parse_angle(row[1]) #converts into degrees
            print(f"{RA}\n")
            DEC = parse_angle(row[2]) #converts into degrees

            # Get the positions of the Sun and Pulsar at each time step
            sun_positions, target_positions, sep_ang_series = get_positions(times, gmrt_location, RA, DEC)
            
            if not times or sep_ang_series is None or len(times) != len(sep_ang_series) or not sep_ang_series.any():
                print(f"Error: Invalid data for target {target_name}. Skipping...")
                sys.exit(1)
            
            # plot the separation angle timeseries
            plot_separation_angle(times, outputfolder, sep_ang_series, target_name, threshold)
            
            #writing to a text file
            with open(summary, 'a') as file:
                file.write("----------------------------------------------------------------------------------------- \n")
                file.write(f"{target_name}        {row[1]}        {row[2]} \n")
                file.write("----------------------------------------------------------------------------------------- \n")
                file.write("Source        Obs Time         Separation Angle \n")
                
                for t, s in zip(times, sep_ang_series):
                    if s <= threshold:
                        file.write(f"{target_name}        {t}        {s} \n")
                        
                file.write("############################################################################################# \n \n")
            

if __name__ == "__main__":

    # Ensure correct number of arguments given by the user
    if len(sys.argv) != 1:
        print("Usage: ./<script>.py")
        sys.exit(1)
    
    
    #prompt for source list file
    srclistfile = prompt_for_sourcelist() #give the full path when prompted
    
    # Prompt for date and time separately
    date_part = prompt_for_date()
    time_part = prompt_for_time()
    # Concatenate date and time
    start_time_ist = f"{date_part} {time_part}"
    print(f"start_time_ist: {start_time_ist}")

    #prompting user for obs time
    obstime = int(prompt_for_duration())  #complete observation duration for that epoch  
    
    #files to be kept in sacrosanct locations
    obsrvcoordfile = "ObservatoryCoord.txt"    
    
    #output files to be kept in sacrosanct locations
    outputfolder = "output_final"    

    #prompt for the separation angle threshold
    threshold = prompt_for_sepangthreshold()
      
    if os.path.exists(obsrvcoordfile):
        print("Observatory Coordinate File exists.")
    else:
        print("Observatory Coordinate File does not exist.")
        sys.exit(1)
        
    if os.path.exists(outputfolder):
        print("Output Folder exists.")
    else:
        print("Output Folder does not exist.")
        
    summary_file = create_or_clear_directory(outputfolder)
    
    obs_name = "GMRT" # you might include this as well in the command-line argument in future when we deal with other observatories. Te obsrvcoord file and the script is already prepared to take account of this
    
    with open(summary_file, 'a') as file:
        file.write(f"Observatory Name: {obs_name} \n")
        file.write(f"Start Time: {start_time_ist} \n")
        file.write(f"Observation Duration: {obstime} \n")
    
    
    main(obsrvcoordfile, outputfolder, summary_file, srclistfile, start_time_ist, obstime, threshold, obs_name)
